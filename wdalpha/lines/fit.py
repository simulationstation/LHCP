from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Sequence

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.optimize import curve_fit


def _continuum(x: np.ndarray, c0: float, c1: float, x0: float) -> np.ndarray:
    return c0 + c1 * (x - x0)


def _gaussian(x: np.ndarray, amp: float, cen: float, sigma: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - cen) / sigma) ** 2)


def _model(x: np.ndarray, params: Sequence[float], n_components: int, x0: float) -> np.ndarray:
    c0, c1 = params[0], params[1]
    model = _continuum(x, c0, c1, x0)
    offset = 2
    for i in range(n_components):
        amp, cen, sigma = params[offset + 3 * i : offset + 3 * i + 3]
        model += _gaussian(x, amp, cen, sigma)
    return model


@dataclass
class LineFitResult:
    line_id: str
    species: str
    lambda0: float
    lambda_obs: float
    sigma_lambda_obs: float
    depth: float
    width: float
    continuum_c0: float
    continuum_c1: float
    chi2_local: float
    reduced_chi2_local: float
    asymmetry: float
    snr: float
    window_contam: float
    flags: List[str] = field(default_factory=list)
    n_components: int = 1


def estimate_z_guess(
    wavelength: np.ndarray,
    flux: np.ndarray,
    line_centers: Iterable[float],
    window_aa: float,
) -> float:
    centers = []
    for center in line_centers:
        mask = (wavelength >= center - window_aa) & (wavelength <= center + window_aa)
        if mask.sum() < 5:
            continue
        local_wave = wavelength[mask]
        local_flux = flux[mask]
        idx = np.argmin(local_flux)
        lambda_obs = local_wave[idx]
        centers.append(lambda_obs / center - 1.0)
    if not centers:
        return 0.0
    return float(np.median(centers))


def _initial_center(wave: np.ndarray, flux: np.ndarray, center_guess: float) -> float:
    idx = np.argmin(flux)
    return float(wave[idx]) if wave.size else center_guess


def _count_minima(flux: np.ndarray) -> int:
    if flux.size < 5:
        return 1
    grad = np.diff(flux)
    minima = np.where((grad[:-1] < 0) & (grad[1:] > 0))[0]
    return int(minima.size)


def _asymmetry_proxy(wave: np.ndarray, flux: np.ndarray, center: float) -> float:
    left = flux[wave < center]
    right = flux[wave >= center]
    if left.size == 0 or right.size == 0:
        return 0.0
    return float((np.mean(left) - np.mean(right)) / max(1e-6, np.std(flux)))


def fit_line_window(
    wave: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    center_guess: float,
    max_components: int = 2,
    edge_buffer_aa: float = 0.2,
) -> LineFitResult | None:
    wave_range = wave.max() - wave.min()
    x0 = center_guess
    center_init = _initial_center(wave, flux, center_guess)
    amp_guess = np.min(flux) - np.median(flux)
    amp_guess = min(amp_guess, -0.01)
    sigma_guess = max(1e-3, wave_range / 8.0)
    c0_guess = np.median(flux)
    c1_guess = 0.0

    n_minima = _count_minima(flux)
    n_components = 1 if n_minima <= 1 else min(max_components, 2)
    params = [c0_guess, c1_guess]
    bounds_low = [-np.inf, -np.inf]
    bounds_high = [np.inf, np.inf]
    for i in range(n_components):
        cen_offset = (i - (n_components - 1) / 2) * sigma_guess
        params.extend([amp_guess, center_init + cen_offset, sigma_guess])
        bounds_low.extend([-np.inf, wave.min(), 1e-4])
        bounds_high.extend([0.0, wave.max(), wave_range])

    try:
        popt, pcov = curve_fit(
            lambda x, *p: _model(x, p, n_components, x0),
            wave,
            flux,
            p0=params,
            sigma=error,
            absolute_sigma=True,
            bounds=(bounds_low, bounds_high),
            maxfev=20000,
        )
    except (ValueError, RuntimeError):
        return None

    if np.any(np.isnan(pcov)) or np.any(np.isinf(pcov)):
        return None

    perr = np.sqrt(np.diag(pcov))
    centers = popt[3::3]
    center_errors = perr[3::3]
    best_idx = int(np.argmin(np.abs(centers - center_guess)))
    model = _model(wave, popt, n_components, x0)
    chi2 = float(np.sum(((flux - model) / error) ** 2))
    dof = max(1, wave.size - len(popt))
    reduced_chi2 = chi2 / dof
    depth = float(np.min(model) - np.median(model))
    width = float(np.abs(popt[4 + 3 * best_idx]))
    asymmetry = _asymmetry_proxy(wave, flux, centers[best_idx])
    snr = float(np.abs(depth) / max(1e-6, np.median(error)))
    window_contam = float(np.std(flux - model) / max(1e-6, np.std(flux)))

    flags = []
    if n_minima > 1:
        flags.append("MULTI_MINIMA")
    if n_components > 1:
        flags.append("BLEND_SUSPECT")
    if centers[best_idx] - wave.min() < edge_buffer_aa or wave.max() - centers[best_idx] < edge_buffer_aa:
        flags.append("EDGE")
    return LineFitResult(
        line_id="",
        species="",
        lambda0=center_guess,
        lambda_obs=float(centers[best_idx]),
        sigma_lambda_obs=float(center_errors[best_idx]),
        depth=float(depth),
        width=float(width),
        continuum_c0=float(popt[0]),
        continuum_c1=float(popt[1]),
        chi2_local=chi2,
        reduced_chi2_local=float(reduced_chi2),
        asymmetry=float(asymmetry),
        snr=float(snr),
        window_contam=float(window_contam),
        flags=flags,
        n_components=n_components,
    )


def fit_lines(
    wavelength: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    line_table: List[dict],
    window_aa: float,
    max_components: int,
    edge_buffer_aa: float,
    jobs: int = 4,
) -> List[LineFitResult]:
    rows = list(line_table)

    def _fit_row(row: dict) -> LineFitResult | None:
        center = row["lambda_expected"]
        mask = (wavelength >= center - window_aa) & (wavelength <= center + window_aa)
        if mask.sum() < 8:
            return None
        wave = wavelength[mask]
        fl = flux[mask]
        err = error[mask]
        res = fit_line_window(
            wave,
            fl,
            err,
            center,
            max_components=max_components,
            edge_buffer_aa=edge_buffer_aa,
        )
        if res is None:
            return None
        res.line_id = str(row["line_id"])
        res.species = str(row["species"])
        res.lambda0 = float(row["lambda0"])
        if row.get("edge_flag", False):
            res.flags.append("EDGE")
        return res

    results: List[LineFitResult] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for res in executor.map(_fit_row, rows):
            if res is not None:
                results.append(res)
    return results


def apply_quality_flags(
    result: LineFitResult,
    min_snr: float,
    max_reduced_chi2: float,
    max_asymmetry: float,
    saturation_threshold: float,
) -> LineFitResult:
    if result.snr < min_snr:
        result.flags.append("LOW_SNR")
    if result.reduced_chi2_local > max_reduced_chi2:
        result.flags.append("POOR_FIT")
    if abs(result.asymmetry) > max_asymmetry:
        result.flags.append("BLEND_SUSPECT")
    if result.depth < -saturation_threshold:
        result.flags.append("SATURATED")
    result.flags = sorted(set(result.flags))
    return result
