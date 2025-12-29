"""
Holistic V2 Pipeline - Integrated redesign addressing attrition issues

Key changes from V1:
1. Relaxed "analysis set" selection (50+ lines) vs strict "gold" (conservative)
2. Robust Student-t likelihood for outlier-resistant inference
3. Identifiability protection for distortion vs sensitivity
4. Signal injection/recovery test
5. Comprehensive progress logging
"""
from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import minimize
from scipy.special import gammaln
from concurrent.futures import ThreadPoolExecutor

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


@dataclass
class V2Config:
    """Configuration for holistic v2 pipeline"""
    seed: int = 42
    target: str = "G191-B2B"
    data_root: Path = field(default_factory=lambda: Path("data"))
    results_root: Path = field(default_factory=lambda: Path("results/holistic_v2"))
    spectrum_path: Optional[Path] = None
    atomic_path: Optional[Path] = None
    jobs: int = 4

    # Line fitting
    window_aa: float = 0.6
    max_components: int = 2
    edge_buffer_aa: float = 0.2

    # Analysis set thresholds (relaxed from gold)
    analysis_sigma_max: float = 0.05  # Relaxed from 0.01
    analysis_chi2_max: float = 20.0   # Relaxed from 5.0
    analysis_snr_min: float = 2.0     # Relaxed from 3.0
    analysis_asymmetry_max: float = 1.0  # Relaxed from 0.5
    analysis_contam_max: float = 0.8  # Relaxed from 0.5
    # Only exclude EDGE - allow BLEND_SUSPECT and POOR_FIT with downweighting
    analysis_excluded_flags: List[str] = field(default_factory=lambda: ['EDGE'])

    # Gold thresholds (strict, for comparison)
    gold_sigma_max: float = 0.01
    gold_chi2_max: float = 5.0
    gold_snr_min: float = 5.0
    gold_excluded_flags: List[str] = field(default_factory=lambda: ['POOR_FIT', 'EDGE', 'BLEND_SUSPECT'])

    # Inference
    include_lab_uncertainty: bool = True
    distortion_model: str = "poly"
    distortion_degree: int = 1
    jitter_init: float = 3e-7
    max_iter: int = 2000

    # Robust likelihood
    use_student_t: bool = True
    student_t_df: float = 4.0  # Degrees of freedom (lower = more robust)

    # Identifiability
    condition_number_threshold: float = 100.0
    auto_orthogonalize: bool = True

    # Minimum lines gate
    min_analysis_lines: int = 50
    min_gold_lines: int = 5


@dataclass
class InferenceResultV2:
    """Inference result with additional diagnostics"""
    delta_alpha: float
    delta_alpha_err: float
    z0: float
    z0_err: float
    jitter: float
    jitter_err: float
    chi2: float
    dof: int
    n_lines: int
    distortion_coeffs: np.ndarray
    condition_number: float
    orthogonalized: bool
    likelihood_type: str
    convergence_info: str


def setup_logging(run_dir: Path) -> Path:
    """Setup logging file"""
    log_path = run_dir / "logs.txt"
    return log_path


def log_progress(log_path: Path, message: str):
    """Log progress with timestamp"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}\n"
    with open(log_path, 'a') as f:
        f.write(line)
    print(line.strip())


def omega0_from_lambda(lambda0_ang: np.ndarray) -> np.ndarray:
    """Convert wavelength to wavenumber"""
    lambda_cm = lambda0_ang * 1e-8
    return 1.0 / lambda_cm


def sensitivity_factor(q_cm1: np.ndarray, omega0_cm1: np.ndarray) -> np.ndarray:
    """Compute sensitivity to alpha variation"""
    return -2.0 * (q_cm1 / omega0_cm1)


def build_design_matrix(
    lambda0: np.ndarray,
    q_cm1: np.ndarray,
    distortion_degree: int,
    orthogonalize: bool = False,
    condition_threshold: float = 100.0
) -> Tuple[np.ndarray, np.ndarray, float, bool]:
    """
    Build design matrix with identifiability protection.

    Returns: (sensitivity, distortion_matrix, condition_number, was_orthogonalized)
    """
    omega0 = omega0_from_lambda(lambda0)
    sensitivity = sensitivity_factor(q_cm1, omega0)

    # Build distortion matrix
    if distortion_degree <= 0:
        distortion = np.zeros((lambda0.size, 0))
    else:
        x = (lambda0 - np.mean(lambda0)) / (np.max(lambda0) - np.min(lambda0) + 1e-10)
        cols = [x ** k for k in range(1, distortion_degree + 1)]
        distortion = np.vstack(cols).T

    # Check condition number
    if distortion.shape[1] > 0:
        design = np.column_stack([sensitivity, distortion])
        cond = np.linalg.cond(design)
    else:
        cond = 1.0

    was_orthogonalized = False

    # Orthogonalize distortion to sensitivity if ill-conditioned
    if orthogonalize and cond > condition_threshold and distortion.shape[1] > 0:
        # Gram-Schmidt: make distortion orthogonal to sensitivity
        sens_norm = sensitivity / (np.linalg.norm(sensitivity) + 1e-10)
        for i in range(distortion.shape[1]):
            proj = np.dot(distortion[:, i], sens_norm) * sens_norm
            distortion[:, i] = distortion[:, i] - proj
        was_orthogonalized = True

        # Recompute condition number
        design = np.column_stack([sensitivity, distortion])
        cond = np.linalg.cond(design)

    return sensitivity, distortion, cond, was_orthogonalized


def student_t_nll(residuals: np.ndarray, sigma: np.ndarray, df: float) -> float:
    """
    Student-t negative log-likelihood for robust inference.

    NLL = sum[ 0.5*(df+1)*log(1 + (r/s)^2/df) + log(s) + const ]
    """
    z = residuals / sigma
    nll = 0.5 * (df + 1) * np.sum(np.log(1 + z**2 / df))
    nll += np.sum(np.log(sigma))
    # Constant terms (don't affect optimization)
    # nll += n * (gammaln(0.5*df) - gammaln(0.5*(df+1)) + 0.5*np.log(df*np.pi))
    return float(nll)


def gaussian_nll(residuals: np.ndarray, sigma: np.ndarray) -> float:
    """
    Gaussian negative log-likelihood with log-variance penalty.
    """
    z = residuals / sigma
    nll = 0.5 * np.sum(z**2 + np.log(sigma**2))
    return float(nll)


def infer_delta_alpha_v2(
    lines: pd.DataFrame,
    atomic: pd.DataFrame,
    config: V2Config,
    use_student_t: bool = True,
) -> InferenceResultV2:
    """
    Robust inference of delta_alpha with Student-t likelihood option.
    """
    # Merge data
    merged = lines.merge(atomic, on="line_id", suffixes=("_obs", "_lab"))

    # Handle column names after merge
    if "lambda0_ang" in merged.columns:
        lambda0 = merged["lambda0_ang"].to_numpy()
    elif "lambda0_ang_obs" in merged.columns:
        lambda0 = merged["lambda0_ang_obs"].to_numpy()
    else:
        lambda0 = merged["lambda0_ang_lab"].to_numpy()

    lambda_obs = merged["lambda_obs"].to_numpy()
    sigma_lambda = merged["sigma_lambda_obs"].to_numpy()

    if "q_cm1" in merged.columns:
        q_cm1 = merged["q_cm1"].to_numpy()
    elif "q_cm1_obs" in merged.columns:
        q_cm1 = merged["q_cm1_obs"].to_numpy()
    else:
        q_cm1 = merged["q_cm1_lab"].to_numpy()

    # Compute observed redshift
    z_obs = lambda_obs / lambda0 - 1.0
    sigma_z = sigma_lambda / lambda0

    # Include lab uncertainty if available
    if config.include_lab_uncertainty:
        for col in ["lambda0_unc_ang", "lambda0_unc_ang_lab"]:
            if col in merged.columns:
                lab_unc = merged[col].to_numpy()
                sigma_z = np.sqrt(sigma_z**2 + (lab_unc / lambda0)**2)
                break

    # Build design matrix with identifiability protection
    sensitivity, distortion, cond_num, was_ortho = build_design_matrix(
        lambda0, q_cm1, config.distortion_degree,
        orthogonalize=config.auto_orthogonalize,
        condition_threshold=config.condition_number_threshold
    )
    n_distortion = distortion.shape[1]

    # Parameter packing
    def pack_params(delta_alpha, z0, dist_coeffs, log_jitter):
        params = [delta_alpha, z0]
        if n_distortion > 0:
            params.extend(dist_coeffs.tolist())
        params.append(log_jitter)
        return np.array(params)

    def unpack_params(params):
        delta_alpha = params[0]
        z0 = params[1]
        dist_coeffs = params[2:2+n_distortion] if n_distortion > 0 else np.array([])
        log_jitter = params[2+n_distortion]
        return delta_alpha, z0, dist_coeffs, log_jitter

    # Initial guess
    params0 = pack_params(
        0.0,
        float(np.median(z_obs)),
        np.zeros(n_distortion),
        np.log(config.jitter_init)
    )

    # Bounds for jitter
    n_params = len(params0)
    lower_bounds = np.full(n_params, -np.inf)
    upper_bounds = np.full(n_params, np.inf)
    jitter_idx = 2 + n_distortion
    lower_bounds[jitter_idx] = np.log(1e-15)
    upper_bounds[jitter_idx] = np.log(10.0 * np.median(sigma_z))

    # Objective function
    df = config.student_t_df

    def objective(params):
        delta_alpha, z0, dist_coeffs, log_jitter = unpack_params(params)
        jitter = np.exp(log_jitter)

        # Model prediction
        dist_term = distortion @ dist_coeffs if n_distortion > 0 else 0.0
        z_model = z0 + dist_term + sensitivity * delta_alpha

        # Effective sigma
        sigma_eff = np.sqrt(sigma_z**2 + jitter**2)

        # Residuals
        residuals = z_obs - z_model

        if use_student_t:
            return student_t_nll(residuals, sigma_eff, df)
        else:
            return gaussian_nll(residuals, sigma_eff)

    # Optimize
    result = minimize(
        objective,
        params0,
        method='L-BFGS-B',
        bounds=list(zip(lower_bounds, upper_bounds)),
        options={'maxiter': config.max_iter, 'disp': False}
    )

    # Extract results
    delta_alpha, z0, dist_coeffs, log_jitter = unpack_params(result.x)
    jitter = float(np.exp(log_jitter))

    # Compute chi2 for diagnostics
    sigma_eff = np.sqrt(sigma_z**2 + jitter**2)
    dist_term = distortion @ dist_coeffs if n_distortion > 0 else 0.0
    z_model = z0 + dist_term + sensitivity * delta_alpha
    chi2 = float(np.sum(((z_obs - z_model) / sigma_eff)**2))
    dof = max(1, len(z_obs) - len(result.x))

    # Error estimation via numerical Hessian
    eps = 1e-8
    def residuals_func(params):
        da, z, dc, lj = unpack_params(params)
        j = np.exp(lj)
        s_eff = np.sqrt(sigma_z**2 + j**2)
        dt = distortion @ dc if n_distortion > 0 else 0.0
        zm = z + dt + sensitivity * da
        return (z_obs - zm) / s_eff

    resid0 = residuals_func(result.x)
    jac_cols = []
    for i in range(len(result.x)):
        params_plus = result.x.copy()
        params_plus[i] += eps
        resid_plus = residuals_func(params_plus)
        jac_cols.append((resid_plus - resid0) / eps)
    jac = np.column_stack(jac_cols)

    try:
        cov = np.linalg.pinv(jac.T @ jac)
        errs = np.sqrt(np.diag(cov))
    except:
        errs = np.full(len(result.x), np.nan)

    delta_alpha_err = float(errs[0]) if len(errs) > 0 else np.nan
    z0_err = float(errs[1]) if len(errs) > 1 else np.nan
    jitter_err = float(jitter * errs[2+n_distortion]) if len(errs) > 2+n_distortion else np.nan

    return InferenceResultV2(
        delta_alpha=float(delta_alpha),
        delta_alpha_err=delta_alpha_err,
        z0=float(z0),
        z0_err=z0_err,
        jitter=jitter,
        jitter_err=jitter_err,
        chi2=chi2,
        dof=dof,
        n_lines=len(z_obs),
        distortion_coeffs=dist_coeffs,
        condition_number=cond_num,
        orthogonalized=was_ortho,
        likelihood_type="student_t" if use_student_t else "gaussian",
        convergence_info=f"success={result.success}, nit={result.nit}"
    )


def run_injection_recovery(
    lines: pd.DataFrame,
    atomic: pd.DataFrame,
    config: V2Config,
    inject_values: List[float] = [0.0, 1e-4, 3e-4],
    n_trials: int = 10,
) -> pd.DataFrame:
    """
    Run injection/recovery test on the actual line distribution.
    """
    np.random.seed(config.seed)

    results = []

    for inject_da in inject_values:
        for trial in range(n_trials):
            # Create simulated data
            merged = lines.merge(atomic, on="line_id", suffixes=("_obs", "_lab"))

            if "lambda0_ang" in merged.columns:
                lambda0 = merged["lambda0_ang"].to_numpy()
            elif "lambda0_ang_obs" in merged.columns:
                lambda0 = merged["lambda0_ang_obs"].to_numpy()
            else:
                lambda0 = merged["lambda0_ang_lab"].to_numpy()

            if "q_cm1" in merged.columns:
                q_cm1 = merged["q_cm1"].to_numpy()
            elif "q_cm1_obs" in merged.columns:
                q_cm1 = merged["q_cm1_obs"].to_numpy()
            else:
                q_cm1 = merged["q_cm1_lab"].to_numpy()

            sigma_lambda = merged["sigma_lambda_obs"].to_numpy()

            # True model
            omega0 = omega0_from_lambda(lambda0)
            sensitivity = sensitivity_factor(q_cm1, omega0)
            z0_true = 5e-5  # Small true redshift
            z_true = z0_true + sensitivity * inject_da

            # Add noise based on measured errors and typical jitter
            jitter = 5e-5
            sigma_z = sigma_lambda / lambda0
            sigma_eff = np.sqrt(sigma_z**2 + jitter**2)
            z_noisy = z_true + np.random.normal(0, sigma_eff)

            # Create simulated lines dataframe
            sim_lines = lines.copy()
            sim_lines['lambda_obs'] = lambda0 * (1 + z_noisy)

            # Run inference
            try:
                res = infer_delta_alpha_v2(sim_lines, atomic, config, use_student_t=True)

                results.append({
                    'inject_delta_alpha': inject_da,
                    'trial': trial,
                    'recovered_delta_alpha': res.delta_alpha,
                    'recovered_err': res.delta_alpha_err,
                    'bias': res.delta_alpha - inject_da,
                    'pull': (res.delta_alpha - inject_da) / res.delta_alpha_err if res.delta_alpha_err > 0 else np.nan,
                    'chi2': res.chi2,
                    'dof': res.dof,
                    'n_lines': res.n_lines
                })
            except Exception as e:
                results.append({
                    'inject_delta_alpha': inject_da,
                    'trial': trial,
                    'recovered_delta_alpha': np.nan,
                    'recovered_err': np.nan,
                    'bias': np.nan,
                    'pull': np.nan,
                    'chi2': np.nan,
                    'dof': np.nan,
                    'n_lines': 0
                })

    return pd.DataFrame(results)


def select_analysis_set(
    all_lines: pd.DataFrame,
    config: V2Config
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Select analysis set (relaxed) and gold set (strict).

    Returns: (analysis_set, gold_set)
    """
    # Filter converged fits
    converged = all_lines[all_lines['lambda_obs'].notna()].copy()

    # Analysis set (relaxed thresholds)
    analysis_mask = (
        (converged['sigma_lambda_obs'] <= config.analysis_sigma_max) &
        (converged['reduced_chi2_local'] <= config.analysis_chi2_max) &
        (converged['snr'] >= config.analysis_snr_min) &
        (converged['asymmetry'].abs() <= config.analysis_asymmetry_max) &
        (converged['window_contam'] <= config.analysis_contam_max)
    )

    # Check excluded flags
    for flag in config.analysis_excluded_flags:
        analysis_mask &= ~converged['flags'].str.contains(flag, na=False)

    analysis_set = converged[analysis_mask].copy()

    # Gold set (strict thresholds)
    gold_mask = (
        (converged['sigma_lambda_obs'] <= config.gold_sigma_max) &
        (converged['reduced_chi2_local'] <= config.gold_chi2_max) &
        (converged['snr'] >= config.gold_snr_min)
    )

    for flag in config.gold_excluded_flags:
        gold_mask &= ~converged['flags'].str.contains(flag, na=False)

    gold_set = converged[gold_mask].copy()

    return analysis_set, gold_set


def run_stress_harness(
    analysis_set: pd.DataFrame,
    atomic: pd.DataFrame,
    config: V2Config,
    log_path: Path
) -> pd.DataFrame:
    """
    Run stress test harness with multiple configurations.
    """
    results = []

    # Configurations to test
    distortion_degrees = [0, 1, 2]
    likelihood_types = ['student_t', 'gaussian']
    species_sets = ['Fe', 'Ni', 'both']
    wavelength_splits = ['full', 'blue', 'red']

    log_progress(log_path, f"Starting stress harness with {len(analysis_set)} analysis lines")

    for dist_deg in distortion_degrees:
        for lik_type in likelihood_types:
            for species in species_sets:
                for wave_split in wavelength_splits:
                    # Select subset
                    subset = analysis_set.copy()

                    # Species filter
                    if species == 'Fe':
                        subset = subset[subset['species'].str.contains('Fe', na=False)]
                    elif species == 'Ni':
                        subset = subset[subset['species'].str.contains('Ni', na=False)]

                    # Wavelength split
                    if wave_split == 'blue':
                        med_wave = subset['lambda0_ang'].median()
                        subset = subset[subset['lambda0_ang'] < med_wave]
                    elif wave_split == 'red':
                        med_wave = subset['lambda0_ang'].median()
                        subset = subset[subset['lambda0_ang'] >= med_wave]

                    # Skip if too few lines
                    if len(subset) < 2:
                        results.append({
                            'distortion_degree': dist_deg,
                            'likelihood': lik_type,
                            'species': species,
                            'wavelength_split': wave_split,
                            'n_lines': len(subset),
                            'delta_alpha': np.nan,
                            'delta_alpha_err': np.nan,
                            'chi2': np.nan,
                            'dof': np.nan,
                            'condition_number': np.nan,
                            'orthogonalized': False,
                            'notes': 'insufficient_lines'
                        })
                        continue

                    # Run inference
                    try:
                        test_config = V2Config(
                            **{k: v for k, v in config.__dict__.items()
                               if not k.startswith('_')}
                        )
                        test_config.distortion_degree = dist_deg

                        use_t = (lik_type == 'student_t')
                        res = infer_delta_alpha_v2(subset, atomic, test_config, use_student_t=use_t)

                        results.append({
                            'distortion_degree': dist_deg,
                            'likelihood': lik_type,
                            'species': species,
                            'wavelength_split': wave_split,
                            'n_lines': res.n_lines,
                            'delta_alpha': res.delta_alpha,
                            'delta_alpha_err': res.delta_alpha_err,
                            'chi2': res.chi2,
                            'dof': res.dof,
                            'condition_number': res.condition_number,
                            'orthogonalized': res.orthogonalized,
                            'notes': ''
                        })
                    except Exception as e:
                        results.append({
                            'distortion_degree': dist_deg,
                            'likelihood': lik_type,
                            'species': species,
                            'wavelength_split': wave_split,
                            'n_lines': len(subset),
                            'delta_alpha': np.nan,
                            'delta_alpha_err': np.nan,
                            'chi2': np.nan,
                            'dof': np.nan,
                            'condition_number': np.nan,
                            'orthogonalized': False,
                            'notes': str(e)[:50]
                        })

    return pd.DataFrame(results)


def generate_exec_summary(
    run_dir: Path,
    config: V2Config,
    baseline_result: InferenceResultV2,
    stress_results: pd.DataFrame,
    n_analysis: int,
    n_gold: int,
    injection_results: pd.DataFrame,
) -> Dict[str, str]:
    """
    Generate executive summary with gates.
    """
    gates = {}

    # Gate N: minimum analysis lines
    if n_analysis >= config.min_analysis_lines:
        gates['gate_n'] = 'PASS'
        gate_n_note = f"N_analysis={n_analysis} >= {config.min_analysis_lines}"
    else:
        gates['gate_n'] = 'FAIL'
        gate_n_note = f"N_analysis={n_analysis} < {config.min_analysis_lines}"

    # Gate chi2: reduced chi2 reasonable
    if baseline_result.dof > 0:
        reduced_chi2 = baseline_result.chi2 / baseline_result.dof
        if reduced_chi2 < 3.0:
            gates['gate_chi2'] = 'PASS'
            gate_chi2_note = f"reduced_chi2={reduced_chi2:.2f} < 3.0"
        else:
            gates['gate_chi2'] = 'FAIL'
            gate_chi2_note = f"reduced_chi2={reduced_chi2:.2f} >= 3.0"
    else:
        gates['gate_chi2'] = 'FAIL'
        gate_chi2_note = "dof=0"

    # Gate stability: drift across stress tests
    valid_stress = stress_results[stress_results['delta_alpha'].notna()]
    if len(valid_stress) > 0:
        max_da = valid_stress['delta_alpha'].max()
        min_da = valid_stress['delta_alpha'].min()
        drift = max_da - min_da
        if drift <= baseline_result.delta_alpha_err * 3:
            gates['gate_stability'] = 'PASS'
            gate_stab_note = f"drift={drift:.2e} <= 3*sigma={3*baseline_result.delta_alpha_err:.2e}"
        else:
            gates['gate_stability'] = 'FAIL'
            gate_stab_note = f"drift={drift:.2e} > 3*sigma={3*baseline_result.delta_alpha_err:.2e}"
    else:
        gates['gate_stability'] = 'FAIL'
        gate_stab_note = "no valid stress results"

    # Gate species: Fe vs Ni consistency
    fe_results = stress_results[(stress_results['species'] == 'Fe') &
                                 (stress_results['wavelength_split'] == 'full') &
                                 stress_results['delta_alpha'].notna()]
    ni_results = stress_results[(stress_results['species'] == 'Ni') &
                                 (stress_results['wavelength_split'] == 'full') &
                                 stress_results['delta_alpha'].notna()]

    if len(fe_results) > 0 and len(ni_results) > 0:
        fe_da = fe_results['delta_alpha'].median()
        ni_da = ni_results['delta_alpha'].median()
        fe_err = fe_results['delta_alpha_err'].median()
        ni_err = ni_results['delta_alpha_err'].median()
        combined_err = np.sqrt(fe_err**2 + ni_err**2)

        # Check if they agree in sign or overlap within errors
        if np.sign(fe_da) == np.sign(ni_da) or abs(fe_da - ni_da) < 2 * combined_err:
            gates['gate_species'] = 'PASS'
            gate_species_note = f"Fe={fe_da:.2e}, Ni={ni_da:.2e}, consistent"
        else:
            gates['gate_species'] = 'FAIL'
            gate_species_note = f"Fe={fe_da:.2e}, Ni={ni_da:.2e}, tension"
    else:
        gates['gate_species'] = 'WARN'
        gate_species_note = "insufficient data for species comparison"

    # Write summary
    with open(run_dir / 'EXEC_SUMMARY.md', 'w') as f:
        f.write("# EXEC_SUMMARY - Holistic V2 Pipeline\n\n")

        f.write("## Baseline Result (Student-t, Both species, Full wavelength)\n\n")
        f.write(f"- **Δα/α** = {baseline_result.delta_alpha:.4e} ± {baseline_result.delta_alpha_err:.4e}\n")
        f.write(f"- z0 = {baseline_result.z0:.4e} ± {baseline_result.z0_err:.4e}\n")
        f.write(f"- jitter = {baseline_result.jitter:.4e}\n")
        f.write(f"- χ² = {baseline_result.chi2:.2f}, dof = {baseline_result.dof}\n")
        f.write(f"- N_analysis = {n_analysis}, N_gold = {n_gold}\n")
        f.write(f"- Condition number = {baseline_result.condition_number:.1f}\n")
        f.write(f"- Orthogonalized = {baseline_result.orthogonalized}\n")
        f.write(f"- Likelihood = {baseline_result.likelihood_type}\n\n")

        f.write("## Gates\n\n")
        f.write(f"| Gate | Status | Note |\n")
        f.write(f"|------|--------|------|\n")
        f.write(f"| N_analysis ≥ 50 | **{gates['gate_n']}** | {gate_n_note} |\n")
        f.write(f"| χ²/dof < 3 | **{gates['gate_chi2']}** | {gate_chi2_note} |\n")
        f.write(f"| Stability | **{gates['gate_stability']}** | {gate_stab_note} |\n")
        f.write(f"| Species | **{gates['gate_species']}** | {gate_species_note} |\n\n")

        f.write("## Injection Recovery Summary\n\n")
        if len(injection_results) > 0:
            for inject_val in injection_results['inject_delta_alpha'].unique():
                subset = injection_results[injection_results['inject_delta_alpha'] == inject_val]
                mean_bias = subset['bias'].mean()
                std_bias = subset['bias'].std()
                mean_pull = subset['pull'].mean()
                f.write(f"- Inject {inject_val:.1e}: bias = {mean_bias:.2e} ± {std_bias:.2e}, mean_pull = {mean_pull:.2f}\n")

        f.write("\n## Stress Test Summary (first 10 rows)\n\n")
        f.write("| Dist | Likelihood | Species | Split | N | Δα/α | σ | Cond |\n")
        f.write("|------|------------|---------|-------|---|------|---|------|\n")
        for _, row in stress_results.head(10).iterrows():
            da_str = f"{row['delta_alpha']:.2e}" if pd.notna(row['delta_alpha']) else "N/A"
            err_str = f"{row['delta_alpha_err']:.2e}" if pd.notna(row['delta_alpha_err']) else "N/A"
            cond_str = f"{row['condition_number']:.1f}" if pd.notna(row['condition_number']) else "N/A"
            f.write(f"| {row['distortion_degree']} | {row['likelihood']} | {row['species']} | {row['wavelength_split']} | {row['n_lines']} | {da_str} | {err_str} | {cond_str} |\n")

        f.write("\n## Artifacts\n\n")
        f.write(f"- line_measurements_all.csv\n")
        f.write(f"- line_measurements_analysis.csv\n")
        f.write(f"- line_measurements_gold.csv\n")
        f.write(f"- summary_delta_alpha.csv\n")
        f.write(f"- inferred_alpha.json\n")
        f.write(f"- injection_recovery.csv\n")
        f.write(f"- diagnostic_plots/\n")

    return gates
