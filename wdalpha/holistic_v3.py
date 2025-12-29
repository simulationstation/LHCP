"""
Holistic V3 Pipeline - Blend-aware line-group modeling with calibrated uncertainties

Key improvements over V2:
1. GROUP fitting: Fit multiple absorption components per window instead of one line
2. Robust z_guess: Cross-correlation or HLSP-anchored shift estimation
3. Uncertainty calibration: Empirical inflation factor from injection/recovery
4. Identifiability protection: Orthogonalize distortion vs sensitivity regressor
"""
from __future__ import annotations

import json
import logging
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.interpolate import UnivariateSpline
from scipy.signal import correlate


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class V3Config:
    """Configuration for holistic V3 pipeline."""
    seed: int = 42
    target: str = "G191-B2B"
    data_root: Path = field(default_factory=lambda: Path("data"))
    results_root: Path = field(default_factory=lambda: Path("results/holistic_v3"))
    jobs: int = 4

    # Group fitting parameters
    window_half_width_aa: float = 0.4  # Window half-width in Angstroms
    max_components: int = 3  # Maximum blend components to fit
    min_component_depth: float = 0.02  # Minimum depth to keep component
    shared_width: bool = True  # Share width across components in group

    # Analysis thresholds (relaxed to keep more lines)
    analysis_sigma_max: float = 0.10  # Max centroid uncertainty (Angstroms)
    analysis_chi2_max: float = 50.0  # Max local reduced chi2
    analysis_snr_min: float = 2.0  # Minimum SNR proxy
    analysis_asymmetry_max: float = 2.0  # Relaxed asymmetry threshold
    analysis_contam_max: float = 0.9  # Max window contamination fraction
    analysis_excluded_flags: List[str] = field(default_factory=lambda: ['EDGE'])

    # Gold thresholds (strict subset)
    gold_sigma_max: float = 0.02
    gold_chi2_max: float = 5.0
    gold_snr_min: float = 5.0

    # Inference settings
    use_student_t: bool = True
    student_t_df: float = 4.0
    distortion_degree: int = 1
    condition_number_threshold: float = 100.0
    auto_orthogonalize: bool = True
    ridge_lambda: float = 1e-6  # Ridge penalty for ill-conditioning

    # Uncertainty calibration
    calibration_n_trials: int = 50
    target_median_pull: float = 0.6745  # Gaussian reference

    # Edge buffer
    edge_buffer_aa: float = 5.0

    # Continuum settings for stress tests
    continuum_spline_s_values: List[float] = field(default_factory=lambda: [0.0005, 0.001, 0.002])


# =============================================================================
# Logging utilities
# =============================================================================

def setup_logging(output_dir: Path) -> Path:
    """Setup logging to file."""
    log_path = output_dir / "logs.txt"
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler()
        ]
    )
    return log_path


def log_progress(msg: str):
    """Log progress message."""
    logging.info(msg)


# =============================================================================
# Robust z_guess estimation
# =============================================================================

def estimate_z_guess_xcorr(
    wavelength: np.ndarray,
    flux: np.ndarray,
    line_wavelengths: np.ndarray,
    z_range: Tuple[float, float] = (-0.001, 0.003),
    n_steps: int = 500
) -> Tuple[float, str, float]:
    """
    Estimate global redshift via cross-correlation with synthetic line comb.

    Returns: (z_guess, method, uncertainty_proxy)
    """
    # Create derivative of observed spectrum (enhances absorption lines)
    dflux = np.gradient(flux, wavelength)
    dflux = dflux / np.std(dflux)

    # Test different z values
    z_test = np.linspace(z_range[0], z_range[1], n_steps)
    xcorr_values = []

    for z in z_test:
        # Create synthetic comb at shifted wavelengths
        comb = np.zeros_like(wavelength)
        for lam0 in line_wavelengths:
            lam_shifted = lam0 * (1 + z)
            # Add Gaussian kernel at line position
            idx_nearest = np.argmin(np.abs(wavelength - lam_shifted))
            if 0 <= idx_nearest < len(wavelength):
                width = 10  # points
                i_lo = max(0, idx_nearest - width)
                i_hi = min(len(wavelength), idx_nearest + width)
                x = np.arange(i_lo, i_hi) - idx_nearest
                comb[i_lo:i_hi] += np.exp(-0.5 * (x / 3)**2)

        # Cross-correlation value (negative because absorption = flux dip)
        xcorr = -np.sum(dflux * comb)
        xcorr_values.append(xcorr)

    xcorr_values = np.array(xcorr_values)
    best_idx = np.argmax(xcorr_values)
    z_guess = z_test[best_idx]

    # Uncertainty proxy from peak width
    peak_height = xcorr_values[best_idx]
    half_max = 0.5 * (peak_height + np.median(xcorr_values))
    above_half = xcorr_values > half_max
    if np.sum(above_half) > 1:
        z_uncertainty = (z_test[above_half][-1] - z_test[above_half][0]) / 2
    else:
        z_uncertainty = (z_range[1] - z_range[0]) / n_steps

    return z_guess, "xcorr", z_uncertainty


def estimate_z_guess_hlsp(
    hlsp_df: pd.DataFrame,
    atomic_wavelengths: np.ndarray,
    match_tolerance_aa: float = 0.3
) -> Tuple[float, str, float]:
    """
    Estimate z_guess from HLSP linelist by matching observed to rest wavelengths.

    Returns: (z_guess, method, uncertainty_proxy)
    """
    hlsp_obs = hlsp_df['WAVEOBS'].values
    hlsp_rest = hlsp_df['WAVEREST'].values

    # Use HLSP rest wavelengths directly if available and non-zero
    valid_mask = (hlsp_rest > 0) & np.isfinite(hlsp_rest)
    if np.sum(valid_mask) < 5:
        raise ValueError("Not enough valid HLSP lines for z estimation")

    z_per_line = (hlsp_obs[valid_mask] / hlsp_rest[valid_mask]) - 1

    # Robust median
    z_guess = np.median(z_per_line)
    z_mad = np.median(np.abs(z_per_line - z_guess))
    z_uncertainty = 1.4826 * z_mad / np.sqrt(np.sum(valid_mask))

    return z_guess, "hlsp", z_uncertainty


# =============================================================================
# Group fitting model
# =============================================================================

def gaussian_absorption(wavelength: np.ndarray, depth: float, center: float, sigma: float) -> np.ndarray:
    """Single Gaussian absorption component (0 = full absorption, 1 = continuum)."""
    return 1.0 - depth * np.exp(-0.5 * ((wavelength - center) / sigma)**2)


def multi_component_model(
    wavelength: np.ndarray,
    params: np.ndarray,
    n_components: int,
    continuum_degree: int = 1,
    shared_width: bool = True
) -> np.ndarray:
    """
    Multi-component absorption model.

    Parameters layout:
    - If shared_width: [continuum_coeffs..., shared_sigma, depth1, center1, depth2, center2, ...]
    - If not shared_width: [continuum_coeffs..., depth1, center1, sigma1, depth2, ...]
    """
    # Extract continuum polynomial
    n_cont = continuum_degree + 1
    cont_coeffs = params[:n_cont]

    # Normalize wavelength for polynomial stability
    lam_mid = (wavelength.min() + wavelength.max()) / 2
    lam_norm = wavelength - lam_mid

    continuum = np.polyval(cont_coeffs[::-1], lam_norm)

    # Extract component parameters
    remaining = params[n_cont:]

    if shared_width:
        shared_sigma = remaining[0]
        comp_params = remaining[1:]
        model = np.ones_like(wavelength)
        for i in range(n_components):
            depth = comp_params[2*i]
            center = comp_params[2*i + 1]
            if depth > 0:
                model *= gaussian_absorption(wavelength, depth, center, shared_sigma)
    else:
        model = np.ones_like(wavelength)
        for i in range(n_components):
            depth = remaining[3*i]
            center = remaining[3*i + 1]
            sigma = remaining[3*i + 2]
            if depth > 0 and sigma > 0:
                model *= gaussian_absorption(wavelength, depth, center, sigma)

    return continuum * model


def fit_group(
    wavelength: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    candidate_centers: List[float],
    n_components: int,
    continuum_degree: int = 1,
    shared_width: bool = True
) -> Dict[str, Any]:
    """
    Fit a group of absorption components to a spectral window.

    Returns fit results including centroids, depths, uncertainties, and diagnostics.
    """
    n_cont = continuum_degree + 1

    # Initial guesses
    cont_init = [1.0] + [0.0] * continuum_degree  # Start with flat continuum

    # Sort candidates by expected position
    candidate_centers = sorted(candidate_centers)[:n_components]
    while len(candidate_centers) < n_components:
        # Pad with dummy positions
        candidate_centers.append(np.mean(wavelength))

    # Initial width estimate from data resolution
    sigma_init = 0.05  # Angstroms, typical for E140H

    if shared_width:
        # params: [cont_coeffs..., shared_sigma, depth1, center1, depth2, center2, ...]
        p0 = cont_init + [sigma_init]
        bounds_lo = [-np.inf] * n_cont + [0.01]
        bounds_hi = [np.inf] * n_cont + [0.5]

        for center in candidate_centers:
            p0.extend([0.2, center])  # depth, center
            bounds_lo.extend([0.0, wavelength.min()])
            bounds_hi.extend([1.0, wavelength.max()])
    else:
        p0 = cont_init
        bounds_lo = [-np.inf] * n_cont
        bounds_hi = [np.inf] * n_cont

        for center in candidate_centers:
            p0.extend([0.2, center, sigma_init])
            bounds_lo.extend([0.0, wavelength.min(), 0.01])
            bounds_hi.extend([1.0, wavelength.max(), 0.5])

    p0 = np.array(p0)
    bounds = (bounds_lo, bounds_hi)

    # Objective function
    def residuals(params):
        model = multi_component_model(wavelength, params, n_components, continuum_degree, shared_width)
        return (flux - model) / error

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = optimize.least_squares(
                residuals, p0, bounds=bounds,
                method='trf', max_nfev=1000
            )

        if not result.success:
            return {'converged': False, 'message': result.message}

        params_opt = result.x
        residuals_final = result.fun
        chi2 = np.sum(residuals_final**2)
        dof = len(wavelength) - len(params_opt)
        reduced_chi2 = chi2 / max(1, dof)

        # Estimate parameter uncertainties from Jacobian
        try:
            J = result.jac
            cov = np.linalg.inv(J.T @ J) * (chi2 / max(1, dof))
            param_errors = np.sqrt(np.diag(cov))
        except:
            param_errors = np.zeros_like(params_opt)

        # Extract component results
        components = []
        if shared_width:
            shared_sigma = params_opt[n_cont]
            shared_sigma_err = param_errors[n_cont] if len(param_errors) > n_cont else 0
            comp_start = n_cont + 1

            for i in range(n_components):
                depth = params_opt[comp_start + 2*i]
                center = params_opt[comp_start + 2*i + 1]
                depth_err = param_errors[comp_start + 2*i] if len(param_errors) > comp_start + 2*i else 0
                center_err = param_errors[comp_start + 2*i + 1] if len(param_errors) > comp_start + 2*i + 1 else 0

                components.append({
                    'depth': depth,
                    'center': center,
                    'sigma': shared_sigma,
                    'depth_err': depth_err,
                    'center_err': center_err,
                    'sigma_err': shared_sigma_err
                })
        else:
            comp_start = n_cont
            for i in range(n_components):
                depth = params_opt[comp_start + 3*i]
                center = params_opt[comp_start + 3*i + 1]
                sigma = params_opt[comp_start + 3*i + 2]
                depth_err = param_errors[comp_start + 3*i] if len(param_errors) > comp_start + 3*i else 0
                center_err = param_errors[comp_start + 3*i + 1] if len(param_errors) > comp_start + 3*i + 1 else 0
                sigma_err = param_errors[comp_start + 3*i + 2] if len(param_errors) > comp_start + 3*i + 2 else 0

                components.append({
                    'depth': depth,
                    'center': center,
                    'sigma': sigma,
                    'depth_err': depth_err,
                    'center_err': center_err,
                    'sigma_err': sigma_err
                })

        # Model for diagnostics
        model_best = multi_component_model(wavelength, params_opt, n_components, continuum_degree, shared_width)

        return {
            'converged': True,
            'params': params_opt,
            'param_errors': param_errors,
            'chi2': chi2,
            'dof': dof,
            'reduced_chi2': reduced_chi2,
            'components': components,
            'continuum_coeffs': params_opt[:n_cont],
            'model': model_best,
            'residuals': residuals_final
        }

    except Exception as e:
        return {'converged': False, 'message': str(e)}


def select_best_model(
    wavelength: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    candidate_centers: List[float],
    max_components: int = 3,
    continuum_degree: int = 1,
    shared_width: bool = True,
    min_depth: float = 0.02
) -> Dict[str, Any]:
    """
    Fit models with K=1..max_components and select best by BIC.
    """
    best_bic = np.inf
    best_result = None
    best_k = 0
    all_results = {}

    n_data = len(wavelength)

    for k in range(1, max_components + 1):
        result = fit_group(
            wavelength, flux, error,
            candidate_centers[:k] if len(candidate_centers) >= k else candidate_centers,
            k, continuum_degree, shared_width
        )

        if result['converged']:
            # Count effective parameters (components with significant depth)
            n_params = len(result['params'])
            significant = sum(1 for c in result['components'] if c['depth'] > min_depth)

            # BIC = chi2 + k*ln(n)
            bic = result['chi2'] + n_params * np.log(n_data)
            all_results[k] = {'result': result, 'bic': bic, 'n_significant': significant}

            if bic < best_bic:
                best_bic = bic
                best_result = result
                best_k = k

    if best_result is None:
        return {'converged': False, 'message': 'No model converged'}

    # Add model selection diagnostics
    best_result['best_k'] = best_k
    best_result['bic'] = best_bic
    best_result['all_models'] = all_results

    # Delta BIC vs K=1
    if 1 in all_results:
        best_result['delta_bic_vs_k1'] = all_results[1]['bic'] - best_bic
    else:
        best_result['delta_bic_vs_k1'] = 0

    # Prune negligible components
    best_result['components'] = [c for c in best_result['components'] if c['depth'] > min_depth]

    return best_result


def fit_all_groups(
    wavelength: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    atomic_df: pd.DataFrame,
    hlsp_df: Optional[pd.DataFrame],
    z_guess: float,
    config: V3Config
) -> pd.DataFrame:
    """
    Fit all target lines using group fitting approach.

    For each target line, find neighbors and fit a multi-component model.
    """
    results = []

    # Build list of all potential absorption wavelengths (atomic + HLSP)
    all_wavelengths = list(atomic_df['lambda0_ang'].values * (1 + z_guess))
    if hlsp_df is not None and 'WAVEOBS' in hlsp_df.columns:
        all_wavelengths.extend(hlsp_df['WAVEOBS'].values)
    all_wavelengths = np.array(sorted(set(all_wavelengths)))

    n_lines = len(atomic_df)
    log_progress(f"Group fitting {n_lines} target lines...")

    for idx, row in atomic_df.iterrows():
        line_id = f"line_{idx:04d}"
        species = row['species']
        lambda0 = row['lambda0_ang']
        q_cm1 = row.get('q_cm-1', row.get('q_cm1', 0))
        lambda0_unc = row.get('sigma_lambda0_ang', row.get('lambda0_unc_ang', 0.005))

        # Expected observed wavelength
        lambda_expected = lambda0 * (1 + z_guess)

        # Check if in range
        if lambda_expected < wavelength.min() + config.edge_buffer_aa:
            results.append({
                'line_id': line_id, 'species': species, 'lambda0_ang': lambda0,
                'q_cm1': q_cm1, 'lambda0_unc_ang': lambda0_unc,
                'converged': False, 'flags': 'EDGE'
            })
            continue
        if lambda_expected > wavelength.max() - config.edge_buffer_aa:
            results.append({
                'line_id': line_id, 'species': species, 'lambda0_ang': lambda0,
                'q_cm1': q_cm1, 'lambda0_unc_ang': lambda0_unc,
                'converged': False, 'flags': 'EDGE'
            })
            continue

        # Define window
        w_lo = lambda_expected - config.window_half_width_aa
        w_hi = lambda_expected + config.window_half_width_aa

        # Extract window data
        mask = (wavelength >= w_lo) & (wavelength <= w_hi)
        if np.sum(mask) < 10:
            results.append({
                'line_id': line_id, 'species': species, 'lambda0_ang': lambda0,
                'q_cm1': q_cm1, 'lambda0_unc_ang': lambda0_unc,
                'converged': False, 'flags': 'INSUFFICIENT_DATA'
            })
            continue

        w_wave = wavelength[mask]
        w_flux = flux[mask]
        w_error = error[mask]

        # Find candidate centers (including target line)
        candidates = [lambda_expected]
        for lam in all_wavelengths:
            if w_lo < lam < w_hi and abs(lam - lambda_expected) > 0.01:
                candidates.append(lam)
        candidates = sorted(set(candidates))

        # Fit group model
        fit_result = select_best_model(
            w_wave, w_flux, w_error,
            candidates,
            max_components=min(config.max_components, len(candidates)),
            continuum_degree=1,
            shared_width=config.shared_width,
            min_depth=config.min_component_depth
        )

        if not fit_result['converged']:
            results.append({
                'line_id': line_id, 'species': species, 'lambda0_ang': lambda0,
                'q_cm1': q_cm1, 'lambda0_unc_ang': lambda0_unc,
                'converged': False, 'flags': 'FIT_FAILED'
            })
            continue

        # Find the component closest to the target wavelength
        target_comp = None
        min_dist = np.inf
        for comp in fit_result['components']:
            dist = abs(comp['center'] - lambda_expected)
            if dist < min_dist:
                min_dist = dist
                target_comp = comp

        if target_comp is None or min_dist > config.window_half_width_aa:
            results.append({
                'line_id': line_id, 'species': species, 'lambda0_ang': lambda0,
                'q_cm1': q_cm1, 'lambda0_unc_ang': lambda0_unc,
                'converged': False, 'flags': 'NO_COMPONENT_FOUND'
            })
            continue

        # Build flags
        flags = []
        n_comps = len(fit_result['components'])
        if n_comps > 1:
            # Check if blend is resolved
            separations = []
            for i, c1 in enumerate(fit_result['components']):
                for c2 in fit_result['components'][i+1:]:
                    sep = abs(c1['center'] - c2['center'])
                    sigma_avg = (c1['sigma'] + c2['sigma']) / 2
                    separations.append(sep / sigma_avg)

            if separations and min(separations) > 2.0:
                flags.append('BLEND_RESOLVED')
            else:
                flags.append('BLEND_MODELED')

        if fit_result['reduced_chi2'] > config.analysis_chi2_max:
            flags.append('HIGH_CHI2')

        # SNR proxy
        snr = target_comp['depth'] / max(target_comp['depth_err'], 1e-10) if target_comp['depth_err'] > 0 else target_comp['depth'] / np.median(w_error)
        if snr < config.analysis_snr_min:
            flags.append('LOW_SNR')

        # Asymmetry check (compare residuals left vs right of center)
        center = target_comp['center']
        left_mask = w_wave < center
        right_mask = w_wave > center
        if np.sum(left_mask) > 0 and np.sum(right_mask) > 0:
            left_resid = np.mean(fit_result['residuals'][left_mask])
            right_resid = np.mean(fit_result['residuals'][right_mask])
            asymmetry = abs(left_resid - right_resid) / np.std(fit_result['residuals'])
        else:
            asymmetry = 0

        if asymmetry > config.analysis_asymmetry_max:
            flags.append('ASYMMETRIC')

        # Window contamination
        contam = (n_comps - 1) / config.max_components
        if contam > config.analysis_contam_max:
            flags.append('HIGH_CONTAM')

        results.append({
            'line_id': line_id,
            'species': species,
            'lambda0_ang': lambda0,
            'q_cm1': q_cm1,
            'lambda0_unc_ang': lambda0_unc,
            'lambda_obs': target_comp['center'],
            'sigma_lambda_obs': target_comp['center_err'],
            'depth': target_comp['depth'],
            'width': target_comp['sigma'],
            'depth_err': target_comp['depth_err'],
            'chi2_local': fit_result['chi2'],
            'reduced_chi2_local': fit_result['reduced_chi2'],
            'dof_local': fit_result['dof'],
            'n_components': n_comps,
            'best_k': fit_result['best_k'],
            'delta_bic_vs_k1': fit_result['delta_bic_vs_k1'],
            'asymmetry': asymmetry,
            'snr': snr,
            'window_contam': contam,
            'converged': True,
            'flags': ';'.join(flags) if flags else ''
        })

        if (idx + 1) % 50 == 0:
            log_progress(f"  Fitted {idx + 1}/{n_lines} lines")

    log_progress(f"Group fitting complete: {sum(1 for r in results if r['converged'])} converged")
    return pd.DataFrame(results)


# =============================================================================
# Analysis set selection
# =============================================================================

def select_analysis_set(
    measurements: pd.DataFrame,
    config: V3Config
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Select analysis set (relaxed) and gold set (strict) from measurements.
    """
    # Start with converged fits
    df = measurements[measurements['converged'] == True].copy()

    # Analysis set criteria (relaxed)
    analysis_mask = (
        (df['sigma_lambda_obs'] > 0) &  # Must have valid uncertainty
        (df['sigma_lambda_obs'] <= config.analysis_sigma_max) &
        (df['reduced_chi2_local'] <= config.analysis_chi2_max) &
        (df['snr'] >= config.analysis_snr_min) &
        (df['asymmetry'] <= config.analysis_asymmetry_max) &
        (df['window_contam'] <= config.analysis_contam_max)
    )

    # Exclude lines with hard-exclude flags
    for flag in config.analysis_excluded_flags:
        analysis_mask &= ~df['flags'].str.contains(flag, na=False)

    analysis_df = df[analysis_mask].copy()

    # Gold set (stricter)
    gold_mask = (
        analysis_mask &
        (df['sigma_lambda_obs'] <= config.gold_sigma_max) &
        (df['reduced_chi2_local'] <= config.gold_chi2_max) &
        (df['snr'] >= config.gold_snr_min)
    )
    gold_df = df[gold_mask].copy()

    return analysis_df, gold_df


# =============================================================================
# Inference with identifiability protection
# =============================================================================

def build_design_matrix(
    lambda_obs: np.ndarray,
    sigma_lambda: np.ndarray,
    q_cm1: np.ndarray,
    distortion_degree: int = 1,
    condition_threshold: float = 100.0,
    auto_orthogonalize: bool = True,
    ridge_lambda: float = 1e-6
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Build design matrix with identifiability protection.

    Returns: (design_matrix, weights, identifiability_info)
    """
    n = len(lambda_obs)

    # Sensitivity to alpha: S_i = -2 * q_i / omega_0 (omega_0 in cm^-1)
    omega_0 = 1e8 / lambda_obs  # Convert Angstroms to cm^-1
    S = -2.0 * q_cm1 / omega_0

    # Normalize wavelength for polynomial
    lam_mean = np.mean(lambda_obs)
    lam_std = np.std(lambda_obs) if np.std(lambda_obs) > 0 else 1.0
    lam_norm = (lambda_obs - lam_mean) / lam_std

    # Build distortion basis
    distortion_cols = []
    for d in range(1, distortion_degree + 1):
        distortion_cols.append(lam_norm ** d)

    # Full design matrix: [S | distortion]
    if distortion_cols:
        X = np.column_stack([S] + distortion_cols)
    else:
        X = S.reshape(-1, 1)

    # Weights from centroid uncertainties
    w = 1.0 / np.clip(sigma_lambda, 1e-10, None)
    W = np.diag(w)

    # Check condition number
    XtWX = X.T @ W @ X
    try:
        eigenvalues = np.linalg.eigvalsh(XtWX)
        cond_number = eigenvalues.max() / max(eigenvalues.min(), 1e-12)
    except:
        cond_number = np.inf

    ident_info = {
        'condition_number': cond_number,
        'orthogonalized': False,
        'ridge_applied': False,
        'distortion_degree_used': distortion_degree
    }

    # Orthogonalize distortion to S if ill-conditioned
    if auto_orthogonalize and cond_number > condition_threshold and distortion_degree > 0:
        log_progress(f"  Condition number {cond_number:.1f} > {condition_threshold}, orthogonalizing...")

        # Gram-Schmidt: orthogonalize distortion columns to S
        S_norm = S / np.linalg.norm(S)
        new_distortion = []
        for col in distortion_cols:
            # Remove component parallel to S
            col_orth = col - np.dot(col, S_norm) * S_norm
            if np.linalg.norm(col_orth) > 1e-10:
                new_distortion.append(col_orth)

        if new_distortion:
            X = np.column_stack([S] + new_distortion)
        else:
            X = S.reshape(-1, 1)
            ident_info['distortion_degree_used'] = 0

        ident_info['orthogonalized'] = True

        # Recompute condition number
        XtWX = X.T @ W @ X
        try:
            eigenvalues = np.linalg.eigvalsh(XtWX)
            cond_number = eigenvalues.max() / max(eigenvalues.min(), 1e-12)
            ident_info['condition_number'] = cond_number
        except:
            pass

    # Apply ridge regularization if still ill-conditioned
    if cond_number > condition_threshold and ridge_lambda > 0:
        ident_info['ridge_applied'] = True
        ident_info['ridge_lambda'] = ridge_lambda

    return X, w, ident_info


def infer_delta_alpha(
    lambda_obs: np.ndarray,
    sigma_lambda: np.ndarray,
    q_cm1: np.ndarray,
    lambda0: np.ndarray,
    config: V3Config,
    sigma_inflation: float = 1.0
) -> Dict[str, Any]:
    """
    Infer delta_alpha/alpha with robust likelihood and identifiability protection.
    """
    n = len(lambda_obs)

    # Apply sigma inflation and clamp to minimum
    sigma_min = 0.001  # Minimum 1 mA uncertainty
    sigma_eff = np.clip(sigma_lambda * sigma_inflation, sigma_min, None)

    # Build design matrix
    X, w, ident_info = build_design_matrix(
        lambda_obs, sigma_eff, q_cm1,
        distortion_degree=config.distortion_degree,
        condition_threshold=config.condition_number_threshold,
        auto_orthogonalize=config.auto_orthogonalize,
        ridge_lambda=config.ridge_lambda
    )

    # Response: delta_lambda = lambda_obs - lambda0*(1+z0)
    # We'll fit z0 along with delta_alpha

    # Sensitivity regressor
    omega_0 = 1e8 / lambda_obs
    S = -2.0 * q_cm1 / omega_0

    def neg_log_likelihood(params):
        """Negative log-likelihood for Student-t or Gaussian."""
        delta_alpha = params[0]
        z0 = params[1]
        jitter = np.exp(params[2]) if len(params) > 2 else 0

        # Distortion nuisance parameters
        distortion = np.zeros(n)
        if len(params) > 3:
            lam_mean = np.mean(lambda_obs)
            lam_std = np.std(lambda_obs) if np.std(lambda_obs) > 0 else 1.0
            lam_norm = (lambda_obs - lam_mean) / lam_std

            if ident_info.get('orthogonalized', False):
                # Use orthogonalized basis
                S_norm = S / np.linalg.norm(S)
                for d in range(1, ident_info['distortion_degree_used'] + 1):
                    col = lam_norm ** d
                    col_orth = col - np.dot(col, S_norm) * S_norm
                    if np.linalg.norm(col_orth) > 1e-10 and d <= len(params) - 3:
                        distortion += params[2 + d] * col_orth
            else:
                for d in range(1, config.distortion_degree + 1):
                    if d + 2 < len(params):
                        distortion += params[2 + d] * (lam_norm ** d)

        # Predicted wavelength shift
        delta_lambda_pred = lambda0 * z0 + lambda0 * delta_alpha * S + distortion

        # Observed shift
        delta_lambda_obs = lambda_obs - lambda0

        # Residuals
        resid = delta_lambda_obs - delta_lambda_pred
        var = sigma_eff**2 + jitter**2

        if config.use_student_t:
            df = config.student_t_df
            nll = 0.5 * (df + 1) * np.sum(np.log(1 + resid**2 / (df * var)))
            nll += 0.5 * np.sum(np.log(var))
        else:
            nll = 0.5 * np.sum(resid**2 / var + np.log(var))

        # Ridge penalty on distortion
        if ident_info.get('ridge_applied', False) and len(params) > 3:
            ridge_pen = config.ridge_lambda * np.sum(params[3:]**2)
            nll += ridge_pen

        return nll

    # Initial guess
    z0_init = np.median(lambda_obs / lambda0 - 1)
    n_dist_params = ident_info.get('distortion_degree_used', config.distortion_degree)
    p0 = [0.0, z0_init, np.log(1e-4)] + [0.0] * n_dist_params

    # Optimize
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = optimize.minimize(neg_log_likelihood, p0, method='L-BFGS-B')

        if not result.success:
            return {'converged': False, 'message': result.message}

        params_opt = result.x
        delta_alpha = params_opt[0]
        z0 = params_opt[1]
        jitter = np.exp(params_opt[2])

        # Estimate uncertainties via Hessian
        try:
            hess = result.hess_inv.todense() if hasattr(result.hess_inv, 'todense') else result.hess_inv
            if isinstance(hess, np.ndarray):
                param_errors = np.sqrt(np.diag(hess))
            else:
                # Numerical Hessian
                from scipy.optimize import approx_fprime
                eps = 1e-6
                hess_num = np.zeros((len(params_opt), len(params_opt)))
                for i in range(len(params_opt)):
                    def grad_i(x):
                        return approx_fprime(x, neg_log_likelihood, eps)[i]
                    hess_num[i] = approx_fprime(params_opt, grad_i, eps)
                try:
                    cov = np.linalg.inv(hess_num)
                    param_errors = np.sqrt(np.abs(np.diag(cov)))
                except:
                    param_errors = np.zeros(len(params_opt))
        except:
            param_errors = np.zeros(len(params_opt))

        delta_alpha_err = param_errors[0] if len(param_errors) > 0 else 0
        z0_err = param_errors[1] if len(param_errors) > 1 else 0

        # Compute chi2
        delta_lambda_pred = lambda0 * z0 + lambda0 * delta_alpha * S
        if len(params_opt) > 3:
            lam_mean = np.mean(lambda_obs)
            lam_std = np.std(lambda_obs) if np.std(lambda_obs) > 0 else 1.0
            lam_norm = (lambda_obs - lam_mean) / lam_std
            for d in range(1, n_dist_params + 1):
                if d + 2 < len(params_opt):
                    delta_lambda_pred += params_opt[2 + d] * (lam_norm ** d)

        delta_lambda_obs = lambda_obs - lambda0
        resid = delta_lambda_obs - delta_lambda_pred
        chi2 = np.sum(resid**2 / sigma_eff**2)
        dof = n - len(params_opt)

        return {
            'converged': True,
            'delta_alpha': delta_alpha,
            'delta_alpha_err': delta_alpha_err,
            'z0': z0,
            'z0_err': z0_err,
            'jitter': jitter,
            'chi2': chi2,
            'dof': dof,
            'reduced_chi2': chi2 / max(1, dof),
            'n_lines': n,
            'condition_number': ident_info['condition_number'],
            'orthogonalized': ident_info['orthogonalized'],
            'distortion_degree_used': ident_info['distortion_degree_used'],
            'sigma_inflation': sigma_inflation
        }

    except Exception as e:
        return {'converged': False, 'message': str(e)}


# =============================================================================
# Uncertainty calibration
# =============================================================================

def calibrate_uncertainties(
    analysis_df: pd.DataFrame,
    config: V3Config,
    n_trials: int = 50
) -> Tuple[float, pd.DataFrame]:
    """
    Calibrate uncertainty inflation factor via injection/recovery.

    Returns: (inflation_factor, injection_results_df)
    """
    log_progress(f"Calibrating uncertainties via {n_trials} injection trials...")

    lambda_obs = analysis_df['lambda_obs'].values
    sigma_lambda = analysis_df['sigma_lambda_obs'].values
    q_cm1 = analysis_df['q_cm1'].values
    lambda0 = analysis_df['lambda0_ang'].values

    injection_results = []

    # Inject delta_alpha = 0 and measure pulls
    inject_value = 0.0

    np.random.seed(config.seed)

    for trial in range(n_trials):
        # Perturb observations by their uncertainties
        noise = np.random.normal(0, sigma_lambda)
        lambda_obs_perturbed = lambda_obs + noise

        # Fit
        result = infer_delta_alpha(
            lambda_obs_perturbed, sigma_lambda, q_cm1, lambda0,
            config, sigma_inflation=1.0
        )

        if result['converged']:
            recovered = result['delta_alpha']
            sigma = result['delta_alpha_err']
            bias = recovered - inject_value
            pull = bias / sigma if sigma > 0 else 0

            injection_results.append({
                'trial': trial,
                'inject_value': inject_value,
                'recovered': recovered,
                'sigma': sigma,
                'bias': bias,
                'pull': pull
            })

    inj_df = pd.DataFrame(injection_results)

    if len(inj_df) == 0:
        log_progress("  Warning: No successful injection trials")
        return 1.0, pd.DataFrame()

    # Compute inflation factor
    # Target: median(|pull|) = 0.6745 (Gaussian reference)
    pulls = inj_df['pull'].values
    median_abs_pull = np.median(np.abs(pulls))

    if median_abs_pull > 0:
        inflation_factor = median_abs_pull / config.target_median_pull
    else:
        inflation_factor = 1.0

    # Clamp to reasonable range
    inflation_factor = np.clip(inflation_factor, 0.5, 20.0)

    log_progress(f"  Median |pull| = {median_abs_pull:.3f}")
    log_progress(f"  Inflation factor = {inflation_factor:.3f}")

    inj_df['inflation_factor'] = inflation_factor

    return inflation_factor, inj_df


# =============================================================================
# Stress harness
# =============================================================================

def run_stress_harness(
    analysis_df: pd.DataFrame,
    config: V3Config,
    sigma_inflation: float
) -> pd.DataFrame:
    """
    Run controlled stress test grid.
    """
    log_progress("Running stress harness...")

    results = []

    # Species splits
    species_splits = ['Fe', 'Ni', 'both']

    # Wavelength splits
    wave_splits = ['full', 'blue', 'red']

    # Distortion degrees
    dist_degrees = [0, 1, 2]

    # Likelihoods
    likelihoods = ['student_t', 'gaussian']

    lambda_median = analysis_df['lambda_obs'].median()

    for species in species_splits:
        # Filter by species
        if species == 'both':
            species_df = analysis_df
        elif species == 'Fe':
            species_df = analysis_df[analysis_df['species'].str.contains('Fe')]
        elif species == 'Ni':
            species_df = analysis_df[analysis_df['species'].str.contains('Ni')]
        else:
            continue

        if len(species_df) < 3:
            continue

        for wave_split in wave_splits:
            # Filter by wavelength
            if wave_split == 'full':
                split_df = species_df
            elif wave_split == 'blue':
                split_df = species_df[species_df['lambda_obs'] < lambda_median]
            elif wave_split == 'red':
                split_df = species_df[species_df['lambda_obs'] >= lambda_median]
            else:
                continue

            if len(split_df) < 3:
                continue

            for dist_deg in dist_degrees:
                for likelihood in likelihoods:
                    # Create modified config
                    test_config = V3Config(
                        seed=config.seed,
                        distortion_degree=dist_deg,
                        use_student_t=(likelihood == 'student_t'),
                        student_t_df=config.student_t_df,
                        condition_number_threshold=config.condition_number_threshold,
                        auto_orthogonalize=config.auto_orthogonalize,
                        ridge_lambda=config.ridge_lambda
                    )

                    result = infer_delta_alpha(
                        split_df['lambda_obs'].values,
                        split_df['sigma_lambda_obs'].values,
                        split_df['q_cm1'].values,
                        split_df['lambda0_ang'].values,
                        test_config,
                        sigma_inflation=sigma_inflation
                    )

                    if result['converged']:
                        results.append({
                            'species': species,
                            'wavelength_split': wave_split,
                            'distortion_degree': dist_deg,
                            'likelihood': likelihood,
                            'n_lines': result['n_lines'],
                            'delta_alpha': result['delta_alpha'],
                            'delta_alpha_err': result['delta_alpha_err'],
                            'chi2': result['chi2'],
                            'dof': result['dof'],
                            'reduced_chi2': result['reduced_chi2'],
                            'condition_number': result['condition_number'],
                            'orthogonalized': result['orthogonalized'],
                            'sigma_inflation': sigma_inflation
                        })

    stress_df = pd.DataFrame(results)
    log_progress(f"Stress harness complete: {len(stress_df)} configurations")

    return stress_df


# =============================================================================
# Leave-one-out analysis
# =============================================================================

def leave_one_out_analysis(
    analysis_df: pd.DataFrame,
    config: V3Config,
    sigma_inflation: float,
    baseline_result: Dict[str, Any]
) -> pd.DataFrame:
    """
    Identify most influential lines via leave-one-out.
    """
    log_progress("Running leave-one-out analysis...")

    baseline_da = baseline_result['delta_alpha']
    results = []

    for idx in range(len(analysis_df)):
        # Leave out line idx
        loo_df = analysis_df.drop(analysis_df.index[idx])

        result = infer_delta_alpha(
            loo_df['lambda_obs'].values,
            loo_df['sigma_lambda_obs'].values,
            loo_df['q_cm1'].values,
            loo_df['lambda0_ang'].values,
            config,
            sigma_inflation=sigma_inflation
        )

        if result['converged']:
            line_row = analysis_df.iloc[idx]
            delta_from_baseline = result['delta_alpha'] - baseline_da

            results.append({
                'line_id': line_row['line_id'],
                'species': line_row['species'],
                'lambda0_ang': line_row['lambda0_ang'],
                'delta_alpha_loo': result['delta_alpha'],
                'delta_from_baseline': delta_from_baseline,
                'influence': abs(delta_from_baseline)
            })

    loo_df = pd.DataFrame(results)
    if len(loo_df) > 0:
        loo_df = loo_df.sort_values('influence', ascending=False)

    log_progress(f"LOO analysis complete: {len(loo_df)} lines analyzed")

    return loo_df


# =============================================================================
# Attrition audit
# =============================================================================

def generate_attrition_audit(
    atomic_df: pd.DataFrame,
    measurements_df: pd.DataFrame,
    analysis_df: pd.DataFrame,
    gold_df: pd.DataFrame,
    config: V3Config,
    output_dir: Path
):
    """Generate attrition waterfall and exclusion reasons."""

    # Waterfall stages
    stages = [
        ('S0', 'Total atomic lines', len(atomic_df), 0, ''),
    ]

    n_total = len(atomic_df)

    # In E140H range (check which lines were even attempted)
    n_attempted = len(measurements_df)
    n_out_of_range = n_total - n_attempted
    stages.append(('S1', 'Attempted fits', n_attempted, n_out_of_range, 'OUT_OF_RANGE_OR_EDGE'))

    # Converged
    n_converged = len(measurements_df[measurements_df['converged'] == True])
    n_failed = n_attempted - n_converged
    stages.append(('S2', 'Fit converged', n_converged, n_failed, 'FIT_FAILED'))

    # Analysis set
    n_analysis = len(analysis_df)
    n_quality_cut = n_converged - n_analysis
    stages.append(('S3', 'Analysis set', n_analysis, n_quality_cut, 'QUALITY_CUTS'))

    # Gold set
    n_gold = len(gold_df)
    n_strict_cut = n_analysis - n_gold
    stages.append(('S4', 'Gold set', n_gold, n_strict_cut, 'STRICT_CUTS'))

    # Save waterfall
    waterfall_df = pd.DataFrame(stages, columns=['stage', 'description', 'count', 'dropped', 'reason'])
    waterfall_df.to_csv(output_dir / 'audit_attrition.csv', index=False)

    # Exclusion reasons (for converged fits that didn't make analysis)
    converged_df = measurements_df[measurements_df['converged'] == True].copy()

    reasons = {}

    # Check each criterion
    high_sigma = converged_df['sigma_lambda_obs'] > config.analysis_sigma_max
    reasons['sigma_exceeded'] = high_sigma.sum()

    high_chi2 = converged_df['reduced_chi2_local'] > config.analysis_chi2_max
    reasons['chi2_exceeded'] = high_chi2.sum()

    low_snr = converged_df['snr'] < config.analysis_snr_min
    reasons['low_snr'] = low_snr.sum()

    high_asym = converged_df['asymmetry'] > config.analysis_asymmetry_max
    reasons['asymmetry_exceeded'] = high_asym.sum()

    high_contam = converged_df['window_contam'] > config.analysis_contam_max
    reasons['contamination_exceeded'] = high_contam.sum()

    # Flag-based exclusions
    for flag in config.analysis_excluded_flags:
        has_flag = converged_df['flags'].str.contains(flag, na=False)
        reasons[f'has_{flag.lower()}_flag'] = has_flag.sum()

    # Blend flags (informational, not exclusion)
    blend_modeled = converged_df['flags'].str.contains('BLEND_MODELED', na=False)
    reasons['blend_modeled_kept'] = blend_modeled.sum()

    blend_resolved = converged_df['flags'].str.contains('BLEND_RESOLVED', na=False)
    reasons['blend_resolved_kept'] = blend_resolved.sum()

    # Save exclusion reasons
    reasons_df = pd.DataFrame([
        {'reason': k, 'count': v, 'pct': 100 * v / max(1, n_converged)}
        for k, v in sorted(reasons.items(), key=lambda x: -x[1])
    ])
    reasons_df.to_csv(output_dir / 'audit_exclusion_reasons.csv', index=False)

    # Generate markdown summary
    with open(output_dir / 'audit_examples.md', 'w') as f:
        f.write("# V3 Attrition Audit Report\n\n")
        f.write("## Waterfall Summary\n\n")
        f.write("| Stage | Description | Count | Dropped | Reason |\n")
        f.write("|-------|-------------|-------|---------|--------|\n")
        for _, row in waterfall_df.iterrows():
            f.write(f"| {row['stage']} | {row['description']} | {row['count']} | {row['dropped']} | {row['reason']} |\n")

        f.write("\n## Key Improvements in V3\n\n")
        f.write(f"- **Blend-aware fitting**: {reasons.get('blend_modeled_kept', 0)} lines with blends successfully modeled\n")
        f.write(f"- **Analysis set size**: {n_analysis} lines (target >= 50)\n")
        f.write(f"- **Relaxed thresholds**: sigma_max={config.analysis_sigma_max}, chi2_max={config.analysis_chi2_max}\n")

        f.write("\n## Top Exclusion Reasons (from converged fits)\n\n")
        for _, row in reasons_df.head(5).iterrows():
            f.write(f"- **{row['reason']}**: {row['count']} ({row['pct']:.1f}%)\n")

    log_progress(f"Attrition audit saved to {output_dir}")


# =============================================================================
# Gates evaluation
# =============================================================================

def evaluate_gates(
    analysis_df: pd.DataFrame,
    baseline_result: Dict[str, Any],
    stress_df: pd.DataFrame,
    sigma_inflation: float,
    config: V3Config
) -> Dict[str, Dict[str, Any]]:
    """
    Evaluate quality gates.
    """
    gates = {}

    # Gate 1: N_analysis >= 50 (target >= 100)
    n_analysis = len(analysis_df)
    gates['gate_n'] = {
        'status': 'PASS' if n_analysis >= 50 else 'FAIL',
        'value': n_analysis,
        'threshold': 50,
        'note': f'N_analysis={n_analysis}, target>=100'
    }

    # Gate 2: chi2/dof reasonable after calibration
    reduced_chi2 = baseline_result.get('reduced_chi2', np.inf)
    gates['gate_chi2'] = {
        'status': 'PASS' if reduced_chi2 < 5.0 else 'FAIL',
        'value': reduced_chi2,
        'threshold': 5.0,
        'note': f'chi2/dof={reduced_chi2:.2f} (calibrated sigma)'
    }

    # Gate 3: Stability across stress tests
    if len(stress_df) > 0:
        # Focus on "both" species, full wavelength configurations
        baseline_configs = stress_df[
            (stress_df['species'] == 'both') &
            (stress_df['wavelength_split'] == 'full')
        ]
        if len(baseline_configs) > 0:
            da_values = baseline_configs['delta_alpha'].values
            da_range = da_values.max() - da_values.min()
            baseline_sigma = baseline_result.get('delta_alpha_err', 1) * sigma_inflation
            drift_ratio = da_range / (3 * baseline_sigma)
            gates['gate_stability'] = {
                'status': 'PASS' if drift_ratio < 1.0 else 'FAIL',
                'value': da_range,
                'threshold': 3 * baseline_sigma,
                'note': f'drift={da_range:.2e}, 3*sigma={3*baseline_sigma:.2e}'
            }
        else:
            gates['gate_stability'] = {'status': 'SKIP', 'note': 'No baseline configs'}
    else:
        gates['gate_stability'] = {'status': 'SKIP', 'note': 'No stress results'}

    # Gate 4: Species consistency
    fe_results = stress_df[
        (stress_df['species'] == 'Fe') &
        (stress_df['wavelength_split'] == 'full') &
        (stress_df['distortion_degree'] == config.distortion_degree)
    ]
    ni_results = stress_df[
        (stress_df['species'] == 'Ni') &
        (stress_df['wavelength_split'] == 'full') &
        (stress_df['distortion_degree'] == config.distortion_degree)
    ]

    if len(fe_results) > 0 and len(ni_results) > 0:
        fe_row = fe_results.iloc[0]
        ni_row = ni_results.iloc[0]
        fe_da = fe_row['delta_alpha']
        fe_err = fe_row['delta_alpha_err'] * sigma_inflation
        ni_da = ni_row['delta_alpha']
        ni_err = ni_row['delta_alpha_err'] * sigma_inflation

        tension = abs(fe_da - ni_da) / np.sqrt(fe_err**2 + ni_err**2)
        gates['gate_species'] = {
            'status': 'PASS' if tension < 3.0 else 'FAIL',
            'value': tension,
            'threshold': 3.0,
            'note': f'Fe={fe_da:.2e}+/-{fe_err:.2e}, Ni={ni_da:.2e}+/-{ni_err:.2e}, tension={tension:.1f}sigma'
        }
    else:
        gates['gate_species'] = {'status': 'SKIP', 'note': 'Missing species data'}

    return gates


# =============================================================================
# Executive summary generation
# =============================================================================

def generate_exec_summary(
    config: V3Config,
    baseline_result: Dict[str, Any],
    stress_df: pd.DataFrame,
    loo_df: pd.DataFrame,
    gates: Dict[str, Dict[str, Any]],
    sigma_inflation: float,
    z_guess_info: Dict[str, Any],
    n_analysis: int,
    n_gold: int,
    output_dir: Path
):
    """Generate executive summary markdown."""

    with open(output_dir / 'EXEC_SUMMARY.md', 'w') as f:
        f.write("# EXEC_SUMMARY - Holistic V3 Pipeline\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## Configuration\n\n")
        f.write(f"- Target: {config.target}\n")
        f.write(f"- Seed: {config.seed}\n")
        f.write(f"- Group fitting: max_components={config.max_components}, shared_width={config.shared_width}\n")
        f.write(f"- Distortion degree: {config.distortion_degree}\n")
        f.write(f"- Likelihood: {'Student-t (df={})'.format(config.student_t_df) if config.use_student_t else 'Gaussian'}\n")
        f.write(f"- Uncertainty inflation: {sigma_inflation:.3f}\n\n")

        f.write("## Global Shift Estimation\n\n")
        f.write(f"- Method: {z_guess_info.get('method', 'unknown')}\n")
        f.write(f"- z_guess: {z_guess_info.get('z_guess', 0):.6f}\n")
        f.write(f"- Uncertainty: {z_guess_info.get('uncertainty', 0):.6f}\n\n")

        f.write("## Baseline Result (calibrated uncertainties)\n\n")
        f.write(f"- **N_analysis = {n_analysis}**, N_gold = {n_gold}\n")
        if baseline_result.get('converged'):
            da = baseline_result['delta_alpha']
            da_err = baseline_result['delta_alpha_err'] * sigma_inflation
            f.write(f"- **delta_alpha/alpha = {da:.4e} +/- {da_err:.4e}**\n")
            f.write(f"- z0 = {baseline_result['z0']:.6e} +/- {baseline_result['z0_err']:.6e}\n")
            f.write(f"- chi2/dof = {baseline_result['chi2']:.1f} / {baseline_result['dof']} = {baseline_result['reduced_chi2']:.2f}\n")
            f.write(f"- Condition number = {baseline_result['condition_number']:.1f}\n")
            f.write(f"- Orthogonalized: {baseline_result['orthogonalized']}\n\n")
        else:
            f.write("- Baseline fit did not converge!\n\n")

        f.write("## Species Breakdown\n\n")
        f.write("| Species | N | delta_alpha | sigma (calibrated) |\n")
        f.write("|---------|---|-------------|--------------------|\n")
        for species in ['Fe', 'Ni', 'both']:
            subset = stress_df[
                (stress_df['species'] == species) &
                (stress_df['wavelength_split'] == 'full') &
                (stress_df['distortion_degree'] == config.distortion_degree) &
                (stress_df['likelihood'] == ('student_t' if config.use_student_t else 'gaussian'))
            ]
            if len(subset) > 0:
                row = subset.iloc[0]
                f.write(f"| {species} | {row['n_lines']} | {row['delta_alpha']:.4e} | {row['delta_alpha_err']*sigma_inflation:.4e} |\n")
        f.write("\n")

        f.write("## Stress Test Summary\n\n")
        if len(stress_df) > 0:
            # Calculate max drift
            both_full = stress_df[
                (stress_df['species'] == 'both') &
                (stress_df['wavelength_split'] == 'full')
            ]
            if len(both_full) > 0:
                da_min = both_full['delta_alpha'].min()
                da_max = both_full['delta_alpha'].max()
                f.write(f"- Max drift (both, full): {da_max - da_min:.4e}\n")
            f.write(f"- Total configurations tested: {len(stress_df)}\n\n")

        f.write("## Most Influential Lines (LOO)\n\n")
        if len(loo_df) > 0:
            f.write("| Line ID | Species | lambda0 | delta_from_baseline |\n")
            f.write("|---------|---------|---------|---------------------|\n")
            for _, row in loo_df.head(5).iterrows():
                f.write(f"| {row['line_id']} | {row['species']} | {row['lambda0_ang']:.3f} | {row['delta_from_baseline']:.4e} |\n")
            f.write("\n")

        f.write("## Gates\n\n")
        f.write("| Gate | Status | Note |\n")
        f.write("|------|--------|------|\n")
        for gate_name, gate_info in gates.items():
            status = gate_info['status']
            note = gate_info.get('note', '')
            f.write(f"| {gate_name} | **{status}** | {note} |\n")
        f.write("\n")

        f.write("## Artifacts\n\n")
        f.write("- line_measurements_all.csv\n")
        f.write("- line_measurements_analysis.csv\n")
        f.write("- line_measurements_gold.csv\n")
        f.write("- summary_delta_alpha.csv\n")
        f.write("- injection_recovery.csv\n")
        f.write("- loo_analysis.csv\n")
        f.write("- audit_attrition.csv\n")
        f.write("- audit_exclusion_reasons.csv\n")
        f.write("- identifiability_report.md\n")
        f.write("- diagnostic_plots/\n")

    log_progress(f"Executive summary saved to {output_dir / 'EXEC_SUMMARY.md'}")


def generate_identifiability_report(
    baseline_result: Dict[str, Any],
    config: V3Config,
    output_dir: Path
):
    """Generate identifiability report."""

    with open(output_dir / 'identifiability_report.md', 'w') as f:
        f.write("# Identifiability Report\n\n")
        f.write("## Configuration\n\n")
        f.write(f"- Requested distortion degree: {config.distortion_degree}\n")
        f.write(f"- Condition number threshold: {config.condition_number_threshold}\n")
        f.write(f"- Auto-orthogonalize: {config.auto_orthogonalize}\n")
        f.write(f"- Ridge lambda: {config.ridge_lambda}\n\n")

        f.write("## Results\n\n")
        if baseline_result.get('converged'):
            f.write(f"- Final condition number: {baseline_result['condition_number']:.1f}\n")
            f.write(f"- Orthogonalization applied: {baseline_result['orthogonalized']}\n")
            f.write(f"- Distortion degree used: {baseline_result.get('distortion_degree_used', config.distortion_degree)}\n\n")

            if baseline_result['orthogonalized']:
                f.write("### Action Taken\n\n")
                f.write("Distortion polynomial was orthogonalized to the sensitivity regressor S_i\n")
                f.write("to remove correlation and improve numerical stability.\n")
        else:
            f.write("Baseline fit did not converge.\n")

    log_progress(f"Identifiability report saved to {output_dir / 'identifiability_report.md'}")
