from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_summary(run_dir: Path) -> pd.DataFrame:
    lines = pd.read_csv(run_dir / "line_measurements.csv")
    inferred = pd.read_json(run_dir / "inferred_alpha.json", typ="series")
    summary = lines.describe(include="all")
    summary["delta_alpha"] = inferred.get("delta_alpha")
    summary["delta_alpha_err"] = inferred.get("delta_alpha_err")
    summary_path = run_dir / "summary.csv"
    summary.to_csv(summary_path)
    return summary
