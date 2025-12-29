from __future__ import annotations

import numpy as np


def linear_rebin(wavelength: np.ndarray, flux: np.ndarray, error: np.ndarray, new_wave: np.ndarray):
    flux_new = np.interp(new_wave, wavelength, flux)
    error_new = np.interp(new_wave, wavelength, error)
    return new_wave, flux_new, error_new
