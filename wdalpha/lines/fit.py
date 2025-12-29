from __future__ import annotations

from dataclasses import dataclass
from typing import List

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.optimize import curve_fit


def gaussian(x, amp, cen, sigma, offset):
    return offset + amp * np.exp(-0.5 * ((x - cen) / sigma) ** 2)


def multi_gaussian(x, *params):
    offset = params[-1]
    n = (len(params) - 1) // 3
    model = np.full_like(x, offset, dtype=float)
    for i in range(n):
        amp, cen, sigma = params[3 * i : 3 * i + 3]
        model += amp * np.exp(-0.5 * ((x - cen) / sigma) ** 2)
    return model


@dataclass
class LineFitResult:
    line_id: str
    species: str
    lambda_obs: float
    sigma_lambda_obs: float
    chi2: float
    n_components: int


def fit_line_window(wave, flux, error, center_guess, n_components=1) -> LineFitResult | None:
    """Fit absorption line(s) with Gaussian profile(s).

    Returns None if fitting fails.
    """
    wave_range = wave.max() - wave.min()

    # Initial guesses - ensure they're valid for absorption lines
    amp_guess = np.min(flux) - np.median(flux)
    # Ensure amplitude is negative (absorption)
    amp_guess = min(amp_guess, -0.01)

    sigma_guess = wave_range / (6 * n_components)
    # Ensure sigma is within bounds
    sigma_guess = np.clip(sigma_guess, 1e-3, wave_range * 0.4)

    offset_guess = np.median(flux)
    # Ensure offset is positive
    offset_guess = max(offset_guess, 0.1)

    params = []
    for i in range(n_components):
        cen_offset = (i - (n_components - 1) / 2) * sigma_guess
        cen_guess = np.clip(center_guess + cen_offset, wave.min() + 0.01, wave.max() - 0.01)
        params.extend([amp_guess, cen_guess, sigma_guess])
    params.append(offset_guess)

    bounds_low = []
    bounds_high = []
    for _ in range(n_components):
        bounds_low.extend([-np.inf, wave.min(), 1e-4])
        bounds_high.extend([0.0, wave.max(), wave_range])
    bounds_low.append(0.0)
    bounds_high.append(np.inf)

    try:
        popt, pcov = curve_fit(
            multi_gaussian,
            wave,
            flux,
            p0=params,
            sigma=error,
            absolute_sigma=True,
            bounds=(bounds_low, bounds_high),
            maxfev=10000,
        )
    except (ValueError, RuntimeError) as e:
        # Fitting failed - return None
        return None

    # Check for valid covariance
    if np.any(np.isinf(pcov)) or np.any(np.isnan(pcov)):
        return None

    perr = np.sqrt(np.diag(pcov))
    centers = popt[1::3][:n_components]
    center_errors = perr[1::3][:n_components]
    best_idx = np.argmin(np.abs(centers - center_guess))
    model = multi_gaussian(wave, *popt)
    chi2 = float(np.sum(((flux - model) / error) ** 2))
    return LineFitResult(
        line_id="",
        species="",
        lambda_obs=float(centers[best_idx]),
        sigma_lambda_obs=float(center_errors[best_idx]),
        chi2=chi2,
        n_components=n_components,
    )


def fit_lines(wavelength, flux, error, windows_df, max_components=1) -> List[LineFitResult]:
    rows = list(windows_df.to_dict(orient="records"))

    def _fit_row(row):
        mask = (wavelength >= row["window_min"]) & (wavelength <= row["window_max"])
        if mask.sum() < 5:
            return None
        wave = wavelength[mask]
        fl = flux[mask]
        err = error[mask]
        res = fit_line_window(
            wave,
            fl,
            err,
            row["wavelength_aa"],
            n_components=max_components if row.get("blend_flag", False) else 1,
        )
        if res is None:
            return None
        res.line_id = str(row.get("line_id", row.get("species", "")))
        res.species = str(row["species"])
        return res

    results: List[LineFitResult] = []
    with ThreadPoolExecutor() as executor:
        for res in executor.map(_fit_row, rows):
            if res is not None:
                results.append(res)
    return results
