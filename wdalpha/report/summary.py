from __future__ import annotations

from pathlib import Path
from typing import Dict

import pandas as pd


def build_exec_summary(
    run_dir: Path,
    result_payload: Dict[str, object],
    gating: Dict[str, str],
    summary_table: pd.DataFrame,
) -> None:
    lines = [
        "# EXEC_SUMMARY",
        "",
        "## Δα/α Result",
        f"- Δα/α = {result_payload['delta_alpha']:.3e} ± {result_payload['delta_alpha_err']:.3e}",
        f"- z0 = {result_payload['z0']:.3e} ± {result_payload['z0_err']:.3e}",
        f"- jitter = {result_payload['jitter']:.3e}",
        "",
        "## Gates",
    ]
    for name, status in gating.items():
        lines.append(f"- {name}: {status}")
    lines.extend(
        [
            "",
            "## Artifacts",
            f"- Line measurements (all): {run_dir / 'line_measurements_all.csv'}",
            f"- Line measurements (gold): {run_dir / 'line_measurements_gold.csv'}",
            f"- Summary table: {run_dir / 'summary_delta_alpha.csv'}",
            f"- Plots: {run_dir / 'diagnostic_plots'}",
        ]
    )
    lines.extend(["", "## Stability Table (head)", ""])
    lines.append(summary_table.head(10).to_string(index=False))
    (run_dir / "EXEC_SUMMARY.md").write_text("\n".join(lines))
