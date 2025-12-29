import numpy as np
from wdalpha.lines.fit import fit_line_window


def test_fit_line_window():
    wave = np.linspace(5000, 5001, 200)
    true_center = 5000.5
    flux = 1.0 - 0.4 * np.exp(-0.5 * ((wave - true_center) / 0.05) ** 2)
    error = np.full_like(wave, 0.02)
    result = fit_line_window(wave, flux, error, true_center, n_components=1)
    assert abs(result.lambda_obs - true_center) < 0.01
    assert result.sigma_lambda_obs > 0
