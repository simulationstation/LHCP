from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

C_KMS = 299792.458


@dataclass
class InferenceResult:
    delta_alpha: float
    delta_alpha_err: float
    z: float
    z_err: float
    chi2: float
    dof: int


def _model(wavenumber_obs, x, z):
    return (wavenumber_obs * (1 + z)) - x


def _omega_obs(lambda_obs):
    return 1.0 / (lambda_obs * 1e-8)


def infer_delta_alpha(lines_df: pd.DataFrame, atomic_df: pd.DataFrame, include_lab_uncertainty: bool = True) -> InferenceResult:
    if "line_id" in lines_df.columns and "line_id" in atomic_df.columns:
        merged = lines_df.merge(atomic_df, on="line_id", suffixes=("_obs", "_lab"))
    else:
        merged = lines_df.merge(atomic_df, on="species", suffixes=("_obs", "_lab"))
    omega_obs = _omega_obs(merged["lambda_obs"].to_numpy())
    omega_0 = merged["wavenumber_cm1"].to_numpy()
    q = merged["q_cm1"].to_numpy()

    x_guess = 0.0
    z_guess = np.median(merged["lambda_obs"].to_numpy() / merged["wavelength_aa"].to_numpy() - 1.0)

    sigma_lambda = merged["sigma_lambda_obs"].to_numpy()
    sigma_omega = np.abs(_omega_obs(merged["lambda_obs"].to_numpy() + sigma_lambda) - omega_obs)
    if include_lab_uncertainty:
        sigma_omega = np.sqrt(sigma_omega**2 + (merged["wavelength_unc_aa"].to_numpy() / merged["wavelength_aa"].to_numpy()) ** 2 * omega_0**2)

    def model(params, omega0, qval):
        x, z = params
        return (omega0 + qval * x) / (1 + z)

    def residuals(params, omega_obs_val, omega0, qval):
        return omega_obs_val - model(params, omega0, qval)

    def wrapped(omega0, qval, x, z):
        return model((x, z), omega0, qval)

    # curve_fit signature: f(xdata, *params) -> ydata
    # We use a dummy xdata (indices) since our model uses omega_0/q directly
    n_lines = len(omega_obs)
    if n_lines < 2:
        raise ValueError(f"Need at least 2 lines for inference, got {n_lines}")

    popt, pcov = curve_fit(
        lambda _, x, z: wrapped(omega_0, q, x, z),
        np.arange(n_lines),  # dummy xdata
        omega_obs,  # ydata we're fitting
        p0=[x_guess, z_guess],
        sigma=sigma_omega,
        absolute_sigma=True,
        maxfev=10000,
    )
    x, z = popt
    perr = np.sqrt(np.diag(pcov))
    model_vals = model((x, z), omega_0, q)
    chi2 = float(np.sum(((omega_obs - model_vals) / sigma_omega) ** 2))
    dof = len(omega_obs) - 2
    delta_alpha = x / 2.0
    delta_alpha_err = perr[0] / 2.0
    return InferenceResult(
        delta_alpha=delta_alpha,
        delta_alpha_err=delta_alpha_err,
        z=z,
        z_err=perr[1],
        chi2=chi2,
        dof=dof,
    )


def save_result(result: InferenceResult, out_path: Path, metadata: Dict | None = None) -> None:
    payload = {
        "delta_alpha": result.delta_alpha,
        "delta_alpha_err": result.delta_alpha_err,
        "z": result.z,
        "z_err": result.z_err,
        "chi2": result.chi2,
        "dof": result.dof,
    }
    if metadata:
        payload["metadata"] = metadata
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2))
