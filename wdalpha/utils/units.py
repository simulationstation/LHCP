from __future__ import annotations

import numpy as np


def angstrom_to_cm(value: float | np.ndarray) -> float | np.ndarray:
    return np.asarray(value) * 1e-8


def cm_to_angstrom(value: float | np.ndarray) -> float | np.ndarray:
    return np.asarray(value) / 1e-8


def omega0_from_lambda(lambda0_ang: np.ndarray) -> np.ndarray:
    lambda_cm = angstrom_to_cm(lambda0_ang)
    return 1.0 / lambda_cm


def sensitivity_factor(q_cm1: np.ndarray, omega0_cm1: np.ndarray) -> np.ndarray:
    return -2.0 * (q_cm1 / omega0_cm1)
