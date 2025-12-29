import pandas as pd
from wdalpha.io.tables import build_atomic_table


def test_unit_conversions(tmp_path):
    nist = tmp_path / "nist.csv"
    df = pd.DataFrame(
        {
            "species": ["FeII"],
            "wavelength": [2600.0],
            "wavelength_unc": [0.01],
            "q": [1300.0],
            "q_unc": [50.0],
            "wavelength_unit": ["AA"],
            "q_unit": ["cm-1"],
        }
    )
    df.to_csv(nist, index=False)
    out = tmp_path / "combined.csv"
    built = build_atomic_table(nist, [], out)
    assert "wavelength_aa" in built.columns
    assert built["wavenumber_cm1"].iloc[0] > 0
