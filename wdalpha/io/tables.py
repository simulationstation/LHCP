from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd
from astropy import units as u


REQUIRED_COLUMNS = {"species", "wavelength", "wavelength_unc", "q", "q_unc"}


def _read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")
    return df


def _convert_units(df: pd.DataFrame) -> pd.DataFrame:
    wave_unit = df.get("wavelength_unit", "AA")
    q_unit = df.get("q_unit", "cm-1")
    wave = (df["wavelength"].to_numpy() * u.Unit(wave_unit)).to(u.AA).value
    wave_unc = (df["wavelength_unc"].to_numpy() * u.Unit(wave_unit)).to(u.AA).value
    wavenumber = (1 / (wave * u.AA)).to(1 / u.cm).value
    q = (df["q"].to_numpy() * u.Unit(q_unit)).to(1 / u.cm).value
    q_unc = (df["q_unc"].to_numpy() * u.Unit(q_unit)).to(1 / u.cm).value
    out = df.copy()
    out["wavelength_aa"] = wave
    out["wavelength_unc_aa"] = wave_unc
    out["wavenumber_cm1"] = wavenumber
    out["q_cm1"] = q
    out["q_unc_cm1"] = q_unc
    return out


def build_atomic_table(nist_path: Path, supplemental_paths: List[Path], out_path: Path) -> pd.DataFrame:
    frames = [_read_table(nist_path)]
    for path in supplemental_paths:
        frames.append(_read_table(path))
    df = pd.concat(frames, ignore_index=True)
    df = _convert_units(df)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return df
