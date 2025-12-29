from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize
from concurrent.futures import ThreadPoolExecutor

from wdalpha.utils.units import omega0_from_lambda, sensitivity_factor


@dataclass
class InferenceResult:
    delta_alpha: float
    delta_alpha_err: float
    z0: float
    z0_err: float
    jitter: float
    jitter_err: float
    chi2: float
    dof: int
    distortion_coeffs: np.ndarray
    loss: str


def _distortion_matrix(lambda0: np.ndarray, degree: int) -> np.ndarray:
    if degree <= 0:
        return np.zeros((lambda0.size, 0))
    x = (lambda0 - np.mean(lambda0)) / (np.max(lambda0) - np.min(lambda0))
    cols = [x ** k for k in range(1, degree + 1)]
    return np.vstack(cols).T


def _distortion_binned(lambda0: np.ndarray, n_bins: int) -> np.ndarray:
    if n_bins <= 1:
        return np.zeros((lambda0.size, 0))
    bins = np.linspace(np.min(lambda0), np.max(lambda0), n_bins + 1)
    mat = np.zeros((lambda0.size, n_bins))
    indices = np.digitize(lambda0, bins) - 1
    for i in range(n_bins):
        mat[:, i] = indices == i
    return mat[:, :-1]


def build_design(
    lambda0: np.ndarray,
    q_cm1: np.ndarray,
    model: str,
    degree: int,
    n_bins: int,
) -> Tuple[np.ndarray, np.ndarray]:
    omega0 = omega0_from_lambda(lambda0)
    sensitivity = sensitivity_factor(q_cm1, omega0)
    if model == "binned":
        distortion = _distortion_binned(lambda0, n_bins)
    else:
        distortion = _distortion_matrix(lambda0, degree)
    return sensitivity, distortion


def _pack_params(
    delta_alpha: float,
    z0: float,
    distortion: np.ndarray,
    log_jitter: float,
    log_jitter_species: np.ndarray | None = None,
) -> np.ndarray:
    parts = [delta_alpha, z0]
    if distortion.size:
        parts.extend(distortion.tolist())
    parts.append(log_jitter)
    if log_jitter_species is not None:
        parts.extend(log_jitter_species.tolist())
    return np.array(parts, dtype=float)


def _unpack_params(
    params: np.ndarray,
    n_distortion: int,
    n_species: int,
) -> Tuple[float, float, np.ndarray, float, np.ndarray | None]:
    delta_alpha = params[0]
    z0 = params[1]
    idx = 2
    distortion = params[idx : idx + n_distortion] if n_distortion else np.array([])
    idx += n_distortion
    log_jitter = params[idx]
    idx += 1
    log_jitter_species = None
    if n_species > 0:
        log_jitter_species = params[idx : idx + n_species]
    return delta_alpha, z0, distortion, log_jitter, log_jitter_species


def infer_delta_alpha(
    lines: pd.DataFrame,
    atomic: pd.DataFrame,
    distortion_model: str,
    distortion_degree: int,
    distortion_bins: int,
    robust_loss: str,
    huber_scale: float,
    include_lab_uncertainty: bool,
    jitter_init: float,
    jitter_per_species: bool,
    max_iter: int,
) -> InferenceResult:
    merged = lines.merge(atomic, on="line_id", suffixes=("_obs", "_lab"))
    # Handle lambda0_ang column which may be suffixed after merge
    if "lambda0_ang" in merged.columns:
        lambda0 = merged["lambda0_ang"].to_numpy()
    elif "lambda0_ang_obs" in merged.columns:
        lambda0 = merged["lambda0_ang_obs"].to_numpy()
    elif "lambda0_ang_lab" in merged.columns:
        lambda0 = merged["lambda0_ang_lab"].to_numpy()
    else:
        raise KeyError("lambda0_ang not found in merged data")
    lambda_obs = merged["lambda_obs"].to_numpy()
    sigma_lambda = merged["sigma_lambda_obs"].to_numpy()
    # Handle q_cm1 column which may be suffixed after merge
    if "q_cm1" in merged.columns:
        q_cm1 = merged["q_cm1"].to_numpy()
    elif "q_cm1_obs" in merged.columns:
        q_cm1 = merged["q_cm1_obs"].to_numpy()
    elif "q_cm1_lab" in merged.columns:
        q_cm1 = merged["q_cm1_lab"].to_numpy()
    else:
        raise KeyError("q_cm1 not found in merged data")
    # Handle species column which may be suffixed after merge
    if "species" in merged.columns:
        species = merged["species"].to_numpy()
    elif "species_obs" in merged.columns:
        species = merged["species_obs"].to_numpy()
    elif "species_lab" in merged.columns:
        species = merged["species_lab"].to_numpy()
    else:
        species = np.array(["unknown"] * len(merged))

    z_obs = lambda_obs / lambda0 - 1.0
    sigma_z = sigma_lambda / lambda0
    if include_lab_uncertainty and "lambda0_unc_ang" in merged.columns:
        sigma_z = np.sqrt(sigma_z**2 + (merged["lambda0_unc_ang"].to_numpy() / lambda0) ** 2)

    sensitivity, distortion = build_design(lambda0, q_cm1, distortion_model, distortion_degree, distortion_bins)
    n_distortion = distortion.shape[1]
    unique_species = np.unique(species) if jitter_per_species else np.array([])
    n_species = unique_species.size if jitter_per_species else 0
    species_index = {sp: idx for idx, sp in enumerate(unique_species)}

    params0 = _pack_params(
        0.0,
        float(np.median(z_obs)),
        np.zeros(n_distortion),
        np.log(jitter_init),
        np.log(np.full(n_species, jitter_init)) if jitter_per_species else None,
    )

    # Build bounds to prevent jitter from exploding to infinity
    # Jitter must stay within reasonable limits relative to sigma_z
    log_jitter_lower = np.log(1e-15)
    log_jitter_upper = np.log(10.0 * np.median(sigma_z))  # upper bound: 10x median error
    n_params = len(params0)
    lower_bounds = np.full(n_params, -np.inf)
    upper_bounds = np.full(n_params, np.inf)
    # Set bounds for log_jitter (at index 2 + n_distortion)
    jitter_idx = 2 + n_distortion
    lower_bounds[jitter_idx] = log_jitter_lower
    upper_bounds[jitter_idx] = log_jitter_upper
    # Set bounds for per-species jitter if applicable
    if n_species > 0:
        for i in range(n_species):
            lower_bounds[jitter_idx + 1 + i] = log_jitter_lower
            upper_bounds[jitter_idx + 1 + i] = log_jitter_upper

    def model(params: np.ndarray) -> np.ndarray:
        delta_alpha, z0, distortion_coeffs, _, _ = _unpack_params(params, n_distortion, n_species)
        distortion_term = distortion @ distortion_coeffs if n_distortion else 0.0
        return z0 + distortion_term + sensitivity * delta_alpha

    def neg_log_likelihood(params: np.ndarray) -> float:
        """Proper Gaussian NLL with log-variance penalty to prevent jitter degeneracy."""
        delta_alpha, z0, distortion_coeffs, log_jitter, log_jitter_species = _unpack_params(
            params, n_distortion, n_species
        )
        base = z0 + (distortion @ distortion_coeffs if n_distortion else 0.0) + sensitivity * delta_alpha
        jitter = np.exp(log_jitter)
        if jitter_per_species and log_jitter_species is not None:
            jitter_species = np.exp(log_jitter_species)
            sigma_eff_sq = sigma_z**2 + np.array([jitter_species[species_index[sp]] for sp in species]) ** 2
        else:
            sigma_eff_sq = sigma_z**2 + jitter**2
        # Gaussian NLL: 0.5 * sum[(z - model)^2 / sigma^2 + log(sigma^2)]
        # The log(sigma^2) term prevents jitter from going to infinity
        resid_sq = (z_obs - base) ** 2
        nll = 0.5 * np.sum(resid_sq / sigma_eff_sq + np.log(sigma_eff_sq))
        return float(nll)

    def residuals(params: np.ndarray) -> np.ndarray:
        """Residuals for chi2 calculation and Jacobian estimation."""
        delta_alpha, z0, distortion_coeffs, log_jitter, log_jitter_species = _unpack_params(
            params, n_distortion, n_species
        )
        base = z0 + (distortion @ distortion_coeffs if n_distortion else 0.0) + sensitivity * delta_alpha
        jitter = np.exp(log_jitter)
        sigma = sigma_z.copy()
        if jitter_per_species and log_jitter_species is not None:
            jitter_species = np.exp(log_jitter_species)
            sigma = np.sqrt(sigma**2 + np.array([jitter_species[species_index[sp]] for sp in species]) ** 2)
        else:
            sigma = np.sqrt(sigma**2 + jitter**2)
        return (z_obs - base) / sigma

    # Use minimize with proper NLL instead of least_squares
    # This prevents the jitter-inflation degeneracy
    opt_bounds = list(zip(lower_bounds, upper_bounds))
    result = minimize(
        neg_log_likelihood,
        params0,
        method='L-BFGS-B',
        bounds=opt_bounds,
        options={'maxiter': max_iter, 'disp': False},
    )

    delta_alpha, z0, distortion_coeffs, log_jitter, _ = _unpack_params(result.x, n_distortion, n_species)
    jitter = float(np.exp(log_jitter))
    chi2 = float(np.sum(residuals(result.x) ** 2))
    dof = max(1, z_obs.size - result.x.size)

    # Compute numerical Jacobian of residuals for error estimation
    eps = 1e-8
    resid0 = residuals(result.x)
    jac_cols = []
    for i in range(len(result.x)):
        params_plus = result.x.copy()
        params_plus[i] += eps
        resid_plus = residuals(params_plus)
        jac_cols.append((resid_plus - resid0) / eps)
    jac = np.column_stack(jac_cols)
    cov = np.linalg.pinv(jac.T @ jac)
    errs = np.sqrt(np.diag(cov))
    delta_alpha_err = float(errs[0]) if errs.size > 0 else np.nan
    z0_err = float(errs[1]) if errs.size > 1 else np.nan
    jitter_err = float(jitter * errs[2 + n_distortion]) if errs.size > (2 + n_distortion) else np.nan

    return InferenceResult(
        delta_alpha=float(delta_alpha),
        delta_alpha_err=delta_alpha_err,
        z0=float(z0),
        z0_err=z0_err,
        jitter=jitter,
        jitter_err=jitter_err,
        chi2=chi2,
        dof=dof,
        distortion_coeffs=distortion_coeffs,
        loss=robust_loss,
    )


def build_influence(
    lines: pd.DataFrame,
    atomic: pd.DataFrame,
    config: Dict[str, object],
    jobs: int,
) -> pd.DataFrame:
    # Reset index to ensure 0-based contiguous indices for drop() to work
    lines = lines.reset_index(drop=True)
    rows = []
    line_ids = lines["line_id"].to_numpy()

    def _run(idx: int) -> Dict[str, object]:
        subset = lines.drop(index=idx)
        res = infer_delta_alpha(subset, atomic, **config)
        return {"line_id": line_ids[idx], "delta_alpha": res.delta_alpha, "sigma": res.delta_alpha_err}

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for row in executor.map(_run, range(len(lines))):
            rows.append(row)
    return pd.DataFrame(rows)
