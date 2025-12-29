from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np
from scipy.interpolate import UnivariateSpline


def build_mask(wavelength: np.ndarray, line_centers: Iterable[float], half_width: float) -> np.ndarray:
    mask = np.ones_like(wavelength, dtype=bool)
    for center in line_centers:
        mask &= np.abs(wavelength - center) > half_width
    return mask


def normalize_spectrum(
    wavelength: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    mask: np.ndarray,
    spline_s: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Spline continuum normalize spectrum.

    Returns normalized flux, error, and continuum.
    """
    weights = 1.0 / np.clip(error, 1e-12, None)
    spline = UnivariateSpline(wavelength[mask], flux[mask], w=weights[mask], s=spline_s)
    continuum = spline(wavelength)
    norm_flux = flux / continuum
    norm_error = error / continuum
    return norm_flux, norm_error, continuum


def build_continuum(
    wavelength: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    line_centers: Iterable[float],
    mask_width_aa: float,
    spline_s: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mask = build_mask(wavelength, line_centers, mask_width_aa)
    norm_flux, norm_error, continuum = normalize_spectrum(
        wavelength,
        flux,
        error,
        mask,
        spline_s,
    )
    return norm_flux, norm_error, continuum, mask
