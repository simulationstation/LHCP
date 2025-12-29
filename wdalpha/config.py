from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class ContinuumConfig(BaseModel):
    spline_s: float = Field(0.001, description="Smoothing factor for spline continuum.")
    mask_width_aa: float = Field(0.4, description="Half-width to mask around each line.")


class LineFitConfig(BaseModel):
    window_aa: float = Field(0.6, description="Half-width around line for fitting.")
    max_components: int = Field(2, description="Max Gaussian components for blends.")
    min_snr: float = Field(5.0, description="Minimum depth SNR proxy for keeping lines.")
    max_reduced_chi2: float = Field(4.0, description="Maximum reduced chi2 for acceptable fits.")
    max_asymmetry: float = Field(0.35, description="Maximum asymmetry proxy for acceptable fits.")
    saturation_threshold: float = Field(0.1, description="Flux floor for saturation flag.")
    edge_buffer_aa: float = Field(0.2, description="Edge buffer to flag line fits.")


class InferenceConfig(BaseModel):
    include_lab_uncertainty: bool = Field(True, description="Include lab wavelength uncertainty.")
    distortion_model: str = Field("poly", description="Distortion model: poly or binned.")
    distortion_degree: int = Field(1, description="Polynomial degree for wavelength distortion.")
    distortion_bins: int = Field(4, description="Number of bins for binned distortion model.")
    robust_loss: str = Field("cauchy", description="Robust loss for least squares.")
    huber_scale: float = Field(1.5, description="Huber scale for robust loss.")
    jitter_init: float = Field(3e-7, description="Initial jitter in redshift units.")
    jitter_per_species: bool = Field(False, description="Fit separate jitter per species.")
    max_iter: int = Field(2000, description="Maximum optimizer iterations.")


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


class ContinuumSetting(BaseModel):
    name: str = Field("default", description="Continuum setting name.")
    spline_s: float = Field(0.001, description="Spline smoothing factor.")
    mask_width_aa: float = Field(0.4, description="Half-width to mask around each line.")


class GoldSelectionConfig(BaseModel):
    sigma_lambda_max: float = Field(0.002, description="Maximum centroid uncertainty in Angstrom.")
    reduced_chi2_max: float = Field(3.0, description="Maximum reduced chi2.")
    snr_min: float = Field(5.0, description="Minimum SNR proxy.")
    max_asymmetry: float = Field(0.35, description="Maximum asymmetry proxy.")
    max_window_contam: float = Field(0.4, description="Maximum contamination proxy.")
    excluded_flags: List[str] = Field(default_factory=lambda: ["BLEND_SUSPECT", "POOR_FIT", "EDGE"])


class StressConfig(BaseModel):
    distortion_degrees: List[int] = Field(default_factory=lambda: [0, 1, 2])
    robust_losses: List[str] = Field(default_factory=lambda: ["cauchy", "linear"])
    species_sets: List[str] = Field(default_factory=lambda: ["Fe", "Ni", "both"])
    wavelength_splits: List[str] = Field(default_factory=lambda: ["full", "blue", "red"])


class HolisticRunConfig(BaseModel):
    seed: int = 42
    target: str = "G191-B2B"
    data_root: Path = Path("data")
    results_root: Path = Path("results")
    spectrum_glob: str = "hlsp/wd-linelist/*g191-b2b*e140h*coadd-spec*.fits"
    atomic_path: Path = Path("data/atomic/analysis_lines_clean.csv")
    use_hlsp_linelists: Optional[Path] = None
    do_e230h: bool = False
    jobs: int = 4
    continuum_settings: List[ContinuumSetting] = Field(default_factory=lambda: [ContinuumSetting()])
    linefit: LineFitConfig = Field(default_factory=LineFitConfig)
    gold: GoldSelectionConfig = Field(default_factory=GoldSelectionConfig)
    inference: InferenceConfig = Field(default_factory=InferenceConfig)
    stress: StressConfig = Field(default_factory=StressConfig)
    max_lines: Optional[int] = Field(None, description="Maximum number of lines for quick smoke runs.")
    extra: Dict[str, str] = Field(default_factory=dict)
