from __future__ import annotations

import numpy as np


def detect_lines(wavelength: np.ndarray, flux: np.ndarray, threshold: float = 5.0):
    """Simple line detection placeholder."""
    median = np.median(flux)
    sigma = np.std(flux)
    mask = flux < median - threshold * sigma
    return wavelength[mask]
