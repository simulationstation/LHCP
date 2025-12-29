import numpy as np

from wdalpha.utils.units import angstrom_to_cm, cm_to_angstrom, omega0_from_lambda, sensitivity_factor


def test_unit_conversions(tmp_path):
    lam = np.array([2600.0])
    cm = angstrom_to_cm(lam)
    assert np.isclose(cm_to_angstrom(cm), lam)
    omega0 = omega0_from_lambda(lam)
    assert omega0 > 0
    sens = sensitivity_factor(np.array([1300.0]), omega0)
    assert sens < 0
