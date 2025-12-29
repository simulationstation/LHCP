from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Field


class ContinuumConfig(BaseModel):
    spline_s: float = Field(0.001, description="Smoothing factor for spline continuum.")
    mask_width_aa: float = Field(0.4, description="Half-width to mask around each line.")


class LineFitConfig(BaseModel):
    window_aa: float = Field(0.6, description="Half-width around line for fitting.")
    max_components: int = Field(2, description="Max Gaussian components for blends.")


class InferenceConfig(BaseModel):
    include_lab_uncertainty: bool = Field(True, description="Include lab wavelength uncertainty.")


class RunConfig(BaseModel):
    target: str = "G191-B2B"
    seed: int = 42
    data_root: Path = Path("data")
    results_root: Path = Path("results")
    continuum: ContinuumConfig = Field(default_factory=ContinuumConfig)
    linefit: LineFitConfig = Field(default_factory=LineFitConfig)
    inference: InferenceConfig = Field(default_factory=InferenceConfig)
    n_workers: int = 4
    line_list: Optional[Path] = None


class AtomicTableConfig(BaseModel):
    nist_path: Path
    supplemental_paths: List[Path] = Field(default_factory=list)
    out_path: Path
