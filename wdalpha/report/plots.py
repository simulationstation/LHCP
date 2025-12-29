from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_q_vs_shift(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.errorbar(lines_df["q_cm1"], lines_df["delta_lambda"], yerr=lines_df["sigma_lambda_obs"], fmt="o")
    ax.set_xlabel("q (cm$^{-1}$)")
    ax.set_ylabel("Observed shift (Å)")
    ax.grid(True, alpha=0.3)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_residuals(residuals, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.hist(residuals, bins=20, alpha=0.7)
    ax.set_xlabel("Residual (cm$^{-1}$)")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.3)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
