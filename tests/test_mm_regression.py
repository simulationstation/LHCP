import numpy as np
import pandas as pd

from wdalpha.inference.mm_regression import infer_delta_alpha


def test_mm_regression_recovery():
    delta_alpha = 2e-6
    x = 2 * delta_alpha
    z = 1e-4
    atomic = pd.DataFrame(
        {
            "species": ["FeII", "MgII", "SiII"],
            "wavelength_aa": [2600.0, 2800.0, 1526.0],
            "wavelength_unc_aa": [0.005, 0.005, 0.005],
            "wavenumber_cm1": [1 / (2600e-8), 1 / (2800e-8), 1 / (1526e-8)],
            "q_cm1": [1300.0, 250.0, -1500.0],
        }
    )
    omega_0 = atomic["wavenumber_cm1"].to_numpy()
    q = atomic["q_cm1"].to_numpy()
    omega_obs = (omega_0 + q * x) / (1 + z)
    lambda_obs = 1 / omega_obs / 1e-8
    lines = pd.DataFrame(
        {
            "species": atomic["species"],
            "lambda_obs": lambda_obs,
            "sigma_lambda_obs": np.full(3, 0.002),
        }
    )
    result = infer_delta_alpha(lines, atomic, include_lab_uncertainty=False)
    assert np.isclose(result.delta_alpha, delta_alpha, rtol=1e-2)
