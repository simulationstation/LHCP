import numpy as np
import pandas as pd

from wdalpha.inference.model import infer_delta_alpha


def test_inference_recovery_with_distortion():
    rng = np.random.default_rng(42)
    lambda0 = np.linspace(1400.0, 1500.0, 8)
    # Use realistic q values that don't correlate with lambda0
    # Real q-coefficients vary non-monotonically across species/transitions
    q = np.array([300.0, 1200.0, 500.0, 800.0, 1400.0, 200.0, 1000.0, 600.0])
    # Use delta_alpha that gives adequate SNR for detection
    delta_alpha_true = 1e-4
    z0_true = 2e-4
    distortion_coeff = 3e-6
    x = (lambda0 - np.mean(lambda0)) / (np.max(lambda0) - np.min(lambda0))
    sensitivity = -2.0 * q / (1 / (lambda0 * 1e-8))
    z_obs = z0_true + distortion_coeff * x + sensitivity * delta_alpha_true
    z_obs += rng.normal(scale=3e-7, size=lambda0.size)
    lambda_obs = lambda0 * (1 + z_obs)
    sigma_lambda = np.full_like(lambda0, 0.001)
    atomic = pd.DataFrame(
        {
            "line_id": [f"line_{i}" for i in range(lambda0.size)],
            "species": ["FeII"] * lambda0.size,
            "lambda0_ang": lambda0,
            "q_cm1": q,
        }
    )
    lines = pd.DataFrame(
        {
            "line_id": atomic["line_id"],
            "species": atomic["species"],
            "lambda_obs": lambda_obs,
            "sigma_lambda_obs": sigma_lambda,
        }
    )
    result = infer_delta_alpha(
        lines,
        atomic,
        distortion_model="poly",
        distortion_degree=1,
        distortion_bins=4,
        robust_loss="linear",
        huber_scale=1.5,
        include_lab_uncertainty=False,
        jitter_init=3e-7,
        jitter_per_species=False,
        max_iter=2000,
    )
    assert np.isclose(result.delta_alpha, delta_alpha_true, rtol=0.2), \
        f"delta_alpha={result.delta_alpha:.2e} vs true={delta_alpha_true:.2e}"
