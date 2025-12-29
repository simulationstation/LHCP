from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from wdalpha.inference.model import infer_delta_alpha


@dataclass
class StressRunResult:
    run_name: str
    species_set: str
    distortion_model: str
    continuum_setting: str
    n_lines_gold: int
    delta_alpha: float
    sigma: float
    jitter: float
    fit_metric: float
    notes: str = ""


def _apply_species_split(lines: pd.DataFrame, species_set: str) -> pd.DataFrame:
    if species_set == "Fe":
        return lines[lines["species"].str.contains("Fe", case=False)]
    if species_set == "Ni":
        return lines[lines["species"].str.contains("Ni", case=False)]
    return lines


def _apply_wavelength_split(lines: pd.DataFrame, split: str) -> pd.DataFrame:
    if split == "full":
        return lines
    median = lines["lambda0_ang"].median()
    if split == "blue":
        return lines[lines["lambda0_ang"] <= median]
    if split == "red":
        return lines[lines["lambda0_ang"] > median]
    return lines


def run_stress_suite(
    line_sets: Dict[str, pd.DataFrame],
    atomic: pd.DataFrame,
    base_config: Dict[str, object],
    distortion_degrees: List[int],
    robust_losses: List[str],
    species_sets: List[str],
    wavelength_splits: List[str],
    jobs: int,
) -> pd.DataFrame:
    runs: List[Dict[str, object]] = []
    for cont_name, lines in line_sets.items():
        for loss in robust_losses:
            for degree in distortion_degrees:
                for species_set in species_sets:
                    for split in wavelength_splits:
                        run_name = f"{cont_name}_{loss}_deg{degree}_{species_set}_{split}"
                        runs.append(
                            {
                                "run_name": run_name,
                                "continuum": cont_name,
                                "loss": loss,
                                "degree": degree,
                                "species_set": species_set,
                                "split": split,
                                "lines": lines,
                            }
                        )

    def _run_one(payload: Dict[str, object]) -> Dict[str, object]:
        lines = _apply_species_split(payload["lines"], payload["species_set"])
        lines = _apply_wavelength_split(lines, payload["split"])
        if len(lines) < 2:
            return {
                "run_name": payload["run_name"],
                "species_set": payload["species_set"],
                "distortion_model": base_config["distortion_model"],
                "continuum_setting": payload["continuum"],
                "N_lines_gold": len(lines),
                "delta_alpha_over_alpha": np.nan,
                "sigma": np.nan,
                "jitter": np.nan,
                "fit_metric": np.nan,
                "notes": "insufficient lines",
            }
        config = dict(base_config)
        config["robust_loss"] = payload["loss"]
        config["distortion_degree"] = payload["degree"]
        res = infer_delta_alpha(lines, atomic, **config)
        return {
            "run_name": payload["run_name"],
            "species_set": payload["species_set"],
            "distortion_model": f"{config['distortion_model']}_deg{config['distortion_degree']}",
            "continuum_setting": payload["continuum"],
            "N_lines_gold": len(lines),
            "delta_alpha_over_alpha": res.delta_alpha,
            "sigma": res.delta_alpha_err,
            "jitter": res.jitter,
            "fit_metric": res.chi2 / res.dof if res.dof > 0 else np.nan,
            "notes": "",
        }

    rows: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for row in executor.map(_run_one, runs):
            rows.append(row)
    return pd.DataFrame(rows)
