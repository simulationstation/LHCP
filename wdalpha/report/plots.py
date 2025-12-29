from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def _save(fig: plt.Figure, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_residual_vs_wavelength(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.scatter(lines_df["lambda0_ang"], lines_df["residual_z"], s=18, alpha=0.8)
    ax.axhline(0, color="k", ls="--", alpha=0.4)
    ax.set_xlabel("Lambda0 (Å)")
    ax.set_ylabel("Residual (z)")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)


def plot_residual_vs_q(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.scatter(lines_df["q_cm1"], lines_df["residual_z"], s=18, alpha=0.8)
    ax.axhline(0, color="k", ls="--", alpha=0.4)
    ax.set_xlabel("q (cm$^{-1}$)")
    ax.set_ylabel("Residual (z)")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)


def plot_residual_vs_sensitivity(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.scatter(lines_df["sensitivity"], lines_df["residual_z"], s=18, alpha=0.8)
    ax.axhline(0, color="k", ls="--", alpha=0.4)
    ax.set_xlabel("Sensitivity factor")
    ax.set_ylabel("Residual (z)")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)


def plot_residual_histogram(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.hist(lines_df["residual_norm"], bins=25, alpha=0.75)
    ax.set_xlabel("Normalized residual")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)


def plot_residual_qq(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    stats.probplot(lines_df["residual_norm"], dist="norm", plot=ax)
    ax.set_title("QQ plot")
    _save(fig, out_path)


def plot_influence(influence_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.bar(np.arange(len(influence_df)), influence_df["delta_alpha"] * 1e6)
    ax.set_xlabel("Line index")
    ax.set_ylabel("Leave-one-out Δα/α (ppm)")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)


def plot_q_vs_shift(lines_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.errorbar(lines_df["q_cm1"], lines_df["delta_lambda"], yerr=lines_df["sigma_lambda_obs"], fmt="o")
    ax.set_xlabel("q (cm$^{-1}$)")
    ax.set_ylabel("Observed shift (Å)")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)


def plot_residuals(residuals, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.hist(residuals, bins=20, alpha=0.7)
    ax.set_xlabel("Residual (cm$^{-1}$)")
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.3)
    _save(fig, out_path)
