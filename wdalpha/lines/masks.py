from __future__ import annotations

from typing import Tuple

import pandas as pd


def build_line_windows(atomic_df: pd.DataFrame, half_width: float) -> pd.DataFrame:
    windows = atomic_df.copy()
    windows["window_min"] = windows["wavelength_aa"] - half_width
    windows["window_max"] = windows["wavelength_aa"] + half_width
    windows["blend_flag"] = False
    return windows


def find_window(windows: pd.DataFrame, wavelength: float) -> Tuple[float, float]:
    row = windows.iloc[(windows["wavelength_aa"] - wavelength).abs().argsort()[:1]]
    return float(row["window_min"].values[0]), float(row["window_max"].values[0])
