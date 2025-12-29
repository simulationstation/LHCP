from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import yaml

from wdalpha.config import ContinuumSetting, HolisticRunConfig
from wdalpha.io.fits import check_spectrum_sanity, read_spectrum
from wdalpha.inference.model import build_design, build_influence, infer_delta_alpha
from wdalpha.inference.stress import run_stress_suite
from wdalpha.lines.fit import apply_quality_flags, estimate_z_guess, fit_lines
from wdalpha.preprocess.continuum import build_continuum
from wdalpha.report.plots import (
    plot_influence,
    plot_residual_histogram,
    plot_residual_qq,
    plot_residual_vs_q,
    plot_residual_vs_sensitivity,
    plot_residual_vs_wavelength,
)
from wdalpha.report.summary import build_exec_summary
from wdalpha.utils.units import omega0_from_lambda, sensitivity_factor


def load_holistic_config(default_path: Path, override_path: Path | None, cli_overrides: Dict[str, object]) -> HolisticRunConfig:
    data = yaml.safe_load(default_path.read_text())
    if override_path:
        override = yaml.safe_load(override_path.read_text())
        data.update(override)
    data.update({k: v for k, v in cli_overrides.items() if v is not None})
    return HolisticRunConfig(**data)


def _load_hlsp_linelists(path: Path) -> List[float]:
    wavelengths = []
    for line in path.read_text().splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        try:
            wavelengths.append(float(parts[0]))
        except ValueError:
            continue
    return wavelengths


def _estimate_z_from_hlsp(hlsp_lines: List[float], atomic: pd.DataFrame) -> float:
    if not hlsp_lines:
        return 0.0
    atomic_lambda = atomic["lambda0_ang"].to_numpy()
    z_vals = []
    for obs in hlsp_lines:
        idx = int(np.argmin(np.abs(atomic_lambda - obs)))
        z_vals.append(obs / atomic_lambda[idx] - 1.0)
    return float(np.median(z_vals)) if z_vals else 0.0


def _normalize_atomic_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "lambda0_ang" not in out.columns:
        if "wavelength_aa" in out.columns:
            out["lambda0_ang"] = out["wavelength_aa"]
        elif "wavelength" in out.columns:
            out["lambda0_ang"] = out["wavelength"]
        else:
            raise ValueError("Atomic table missing lambda0_ang/wavelength_aa/wavelength.")
    if "q_cm1" not in out.columns:
        if "q_cm-1" in out.columns:
            out["q_cm1"] = out["q_cm-1"]
        elif "q" in out.columns:
            out["q_cm1"] = out["q"]
    if "line_id" not in out.columns:
        out["line_id"] = [f"line_{i:04d}" for i in range(len(out))]
    if "lambda0_unc_ang" not in out.columns:
        if "sigma_lambda0_ang" in out.columns:
            out["lambda0_unc_ang"] = out["sigma_lambda0_ang"]
        elif "wavelength_unc_aa" in out.columns:
            out["lambda0_unc_ang"] = out["wavelength_unc_aa"]
    if "join_flag" not in out.columns:
        out["join_flag"] = "OK"
    out = out.sort_values("lambda0_ang").reset_index(drop=True)
    return out


def _build_line_table(
    atomic: pd.DataFrame,
    wave_min: float,
    wave_max: float,
    z_guess: float,
    edge_buffer_aa: float,
) -> List[Dict[str, object]]:
    rows = []
    for _, row in atomic.iterrows():
        lambda0 = float(row["lambda0_ang"])
        lambda_expected = lambda0 * (1.0 + z_guess)
        if lambda_expected < wave_min or lambda_expected > wave_max:
            continue
        edge_flag = lambda_expected < wave_min + edge_buffer_aa or lambda_expected > wave_max - edge_buffer_aa
        rows.append(
            {
                "line_id": row["line_id"],
                "species": row["species"],
                "lambda0": lambda0,
                "lambda_expected": lambda_expected,
                "edge_flag": edge_flag,
            }
        )
    return rows


def _flags_to_str(flags: List[str]) -> str:
    return ";".join(sorted(set(flags))) if flags else ""


def _score_gold(
    df: pd.DataFrame,
    config: HolisticRunConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if df.empty:
        return df, df
    reasons = []
    for _, row in df.iterrows():
        line_reasons = []
        if not np.isfinite(row["lambda_obs"]):
            line_reasons.append("fit_failed")
        if row.get("join_flag", "OK") != "OK":
            line_reasons.append("join_flag")
        if np.isfinite(row["sigma_lambda_obs"]) and row["sigma_lambda_obs"] > config.gold.sigma_lambda_max:
            line_reasons.append("sigma_lambda_obs")
        if np.isfinite(row["reduced_chi2_local"]) and row["reduced_chi2_local"] > config.gold.reduced_chi2_max:
            line_reasons.append("reduced_chi2")
        if np.isfinite(row["snr"]) and row["snr"] < config.gold.snr_min:
            line_reasons.append("low_snr")
        if np.isfinite(row["asymmetry"]) and abs(row["asymmetry"]) > config.gold.max_asymmetry:
            line_reasons.append("asymmetry")
        if np.isfinite(row["window_contam"]) and row["window_contam"] > config.gold.max_window_contam:
            line_reasons.append("window_contam")
        flags = set(str(row.get("flags", "")).split(";")) if row.get("flags") else set()
        if flags.intersection(config.gold.excluded_flags):
            line_reasons.append("flags")
        reasons.append(";".join(line_reasons))
    out = df.copy()
    out["exclusion_reason"] = reasons
    gold = out[out["exclusion_reason"] == ""].copy()
    return out, gold


def _write_sanity_report(run_dir: Path, sanity: Dict[str, str]) -> None:
    lines = [
        "# HLSP Sanity Report",
        "",
        f"- Coverage: {sanity['coverage']}",
        f"- Pixels: {sanity['n_pixels']}",
    ]
    if sanity["errors"]:
        lines.append(f"- Errors: {sanity['errors']}")
    if sanity["warnings"]:
        lines.append(f"- Warnings: {sanity['warnings']}")
    (run_dir / "hlsp_sanity_report.md").write_text("\n".join(lines))


def _write_logs(run_dir: Path, argv: List[str]) -> None:
    import numpy
    import pandas
    import scipy

    lines = [
        f"command: {' '.join(argv)}",
        f"python: {sys.version}",
        f"numpy: {numpy.__version__}",
        f"pandas: {pandas.__version__}",
        f"scipy: {scipy.__version__}",
    ]
    (run_dir / "logs.txt").write_text("\n".join(lines))


def _prepare_lines(
    spectrum,
    atomic: pd.DataFrame,
    config: HolisticRunConfig,
    setting: ContinuumSetting,
    z_guess: float,
) -> pd.DataFrame:
    line_table = _build_line_table(
        atomic,
        spectrum.wavelength.min(),
        spectrum.wavelength.max(),
        z_guess,
        config.linefit.edge_buffer_aa,
    )
    if not line_table:
        return pd.DataFrame()
    line_centers = [row["lambda_expected"] for row in line_table]
    norm_flux, norm_error, continuum, mask = build_continuum(
        spectrum.wavelength,
        spectrum.flux,
        spectrum.error,
        line_centers,
        setting.mask_width_aa,
        setting.spline_s,
    )
    results = fit_lines(
        spectrum.wavelength,
        norm_flux,
        norm_error,
        line_table,
        config.linefit.window_aa,
        config.linefit.max_components,
        config.linefit.edge_buffer_aa,
        jobs=config.jobs,
    )
    fitted_ids = {res.line_id for res in results}
    rows: List[Dict[str, object]] = []
    for res in results:
        apply_quality_flags(
            res,
            config.linefit.min_snr,
            config.linefit.max_reduced_chi2,
            config.linefit.max_asymmetry,
            config.linefit.saturation_threshold,
        )
        rows.append(
            {
                "line_id": res.line_id,
                "species": res.species,
                "lambda0_ang": res.lambda0,
                "lambda_obs": res.lambda_obs,
                "sigma_lambda_obs": res.sigma_lambda_obs,
                "depth": res.depth,
                "width": res.width,
                "continuum_c0": res.continuum_c0,
                "continuum_c1": res.continuum_c1,
                "chi2_local": res.chi2_local,
                "reduced_chi2_local": res.reduced_chi2_local,
                "asymmetry": res.asymmetry,
                "snr": res.snr,
                "window_contam": res.window_contam,
                "flags": _flags_to_str(res.flags),
                "n_components": res.n_components,
            }
        )
    missing = [row for row in line_table if row["line_id"] not in fitted_ids]
    for row in missing:
        rows.append(
            {
                "line_id": row["line_id"],
                "species": row["species"],
                "lambda0_ang": row["lambda0"],
                "lambda_obs": np.nan,
                "sigma_lambda_obs": np.nan,
                "depth": np.nan,
                "width": np.nan,
                "continuum_c0": np.nan,
                "continuum_c1": np.nan,
                "chi2_local": np.nan,
                "reduced_chi2_local": np.nan,
                "asymmetry": np.nan,
                "snr": np.nan,
                "window_contam": np.nan,
                "flags": "FIT_FAILED",
                "n_components": 0,
            }
        )
    return pd.DataFrame(rows)


def run_holistic(config: HolisticRunConfig, run_dir: Path, argv: List[str]) -> Dict[str, object]:
    np.random.seed(config.seed)
    run_dir.mkdir(parents=True, exist_ok=True)
    _write_logs(run_dir, argv)

    atomic_path = Path(config.atomic_path)
    if not atomic_path.is_absolute() and not atomic_path.exists():
        atomic_path = config.data_root / atomic_path
    atomic = pd.read_csv(atomic_path)
    atomic = _normalize_atomic_table(atomic)
    if config.max_lines:
        atomic = atomic.head(config.max_lines).copy()
    atomic.to_csv(run_dir / "atomic_used.csv", index=False)

    spectrum_path = None
    if config.extra.get("spectrum"):
        spectrum_path = Path(config.extra["spectrum"])
    else:
        matches = sorted(config.data_root.glob(config.spectrum_glob))
        if matches:
            spectrum_path = matches[0]
    if spectrum_path is None or not spectrum_path.exists():
        raise FileNotFoundError("Spectrum not found. Provide --spectrum or set spectrum_glob.")

    spectrum = read_spectrum(spectrum_path)
    sanity = check_spectrum_sanity(spectrum)
    _write_sanity_report(run_dir, sanity)
    if sanity["errors"]:
        raise RuntimeError(f"Spectrum sanity check failed: {sanity['errors']}")

    z_guess = estimate_z_guess(spectrum.wavelength, spectrum.flux, atomic["lambda0_ang"], config.linefit.window_aa)
    if config.use_hlsp_linelists and Path(config.use_hlsp_linelists).exists():
        hlsp_lines = _load_hlsp_linelists(Path(config.use_hlsp_linelists))
        if hlsp_lines:
            z_guess = _estimate_z_from_hlsp(hlsp_lines, atomic)
    line_sets: Dict[str, pd.DataFrame] = {}
    all_lines = pd.DataFrame()
    baseline = pd.DataFrame()
    for setting in config.continuum_settings:
        line_df = _prepare_lines(spectrum, atomic, config, setting, z_guess)
        # Don't include lambda0_ang in merge - it's already in line_df from _prepare_lines
        line_df = line_df.merge(atomic[["line_id", "join_flag", "q_cm1"]], on="line_id", how="left")
        line_df = line_df.sort_values("lambda0_ang")
        scored, gold = _score_gold(line_df, config)
        line_sets[setting.name] = gold
        if setting.name == config.continuum_settings[0].name:
            all_lines = scored
            baseline = gold

    all_lines.to_csv(run_dir / "line_measurements_all.csv", index=False)
    baseline.to_csv(run_dir / "line_measurements_gold.csv", index=False)
    exclusions = all_lines[all_lines["exclusion_reason"] != ""].copy()
    if not exclusions.empty:
        exclusions.sort_values(["exclusion_reason", "reduced_chi2_local"], ascending=False).to_csv(
            run_dir / "line_exclusions.csv", index=False
        )

    if baseline.empty or len(baseline) < 2:
        raise RuntimeError("Insufficient gold lines for inference.")

    infer_config = {
        "distortion_model": config.inference.distortion_model,
        "distortion_degree": config.inference.distortion_degree,
        "distortion_bins": config.inference.distortion_bins,
        "robust_loss": config.inference.robust_loss,
        "huber_scale": config.inference.huber_scale,
        "include_lab_uncertainty": config.inference.include_lab_uncertainty,
        "jitter_init": config.inference.jitter_init,
        "jitter_per_species": config.inference.jitter_per_species,
        "max_iter": config.inference.max_iter,
    }
    result = infer_delta_alpha(gold, atomic, **infer_config)
    result_payload = {
        "delta_alpha": result.delta_alpha,
        "delta_alpha_err": result.delta_alpha_err,
        "z0": result.z0,
        "z0_err": result.z0_err,
        "jitter": result.jitter,
        "chi2": result.chi2,
        "dof": result.dof,
        "loss": result.loss,
    }
    (run_dir / "inferred_alpha.json").write_text(json.dumps(result_payload, indent=2))

    omega0 = omega0_from_lambda(gold["lambda0_ang"].to_numpy())
    sensitivity = sensitivity_factor(gold["q_cm1"].to_numpy(), omega0)
    _, distortion = build_design(
        gold["lambda0_ang"].to_numpy(),
        gold["q_cm1"].to_numpy(),
        config.inference.distortion_model,
        config.inference.distortion_degree,
        config.inference.distortion_bins,
    )
    distortion_term = distortion @ result.distortion_coeffs if distortion.size else 0.0
    z_model = result.z0 + distortion_term + sensitivity * result.delta_alpha
    z_obs = gold["lambda_obs"].to_numpy() / gold["lambda0_ang"].to_numpy() - 1.0
    sigma_z = gold["sigma_lambda_obs"].to_numpy() / gold["lambda0_ang"].to_numpy()
    residual_z = z_obs - z_model
    residual_norm = residual_z / np.maximum(1e-9, sigma_z)

    diag = gold.copy()
    diag["sensitivity"] = sensitivity
    diag["residual_z"] = residual_z
    diag["residual_norm"] = residual_norm
    diag.to_csv(run_dir / "line_diagnostics.csv", index=False)

    plot_dir = run_dir / "diagnostic_plots"
    plot_residual_vs_wavelength(diag, plot_dir / "residual_vs_wavelength.png")
    plot_residual_vs_q(diag, plot_dir / "residual_vs_q.png")
    plot_residual_vs_sensitivity(diag, plot_dir / "residual_vs_sensitivity.png")
    plot_residual_histogram(diag, plot_dir / "residual_hist.png")
    plot_residual_qq(diag, plot_dir / "residual_qq.png")

    influence = build_influence(gold, atomic, infer_config, jobs=config.jobs)
    influence.to_csv(run_dir / "influence.csv", index=False)
    plot_influence(influence, plot_dir / "influence.png")

    summary = run_stress_suite(
        line_sets,
        atomic,
        infer_config,
        config.stress.distortion_degrees,
        config.stress.robust_losses,
        config.stress.species_sets,
        config.stress.wavelength_splits,
        config.jobs,
    )
    summary.to_csv(run_dir / "summary_delta_alpha.csv", index=False)

    gating = {
        "species_tension": "PASS",
        "chi2_target": "PASS",
        "continuum_stability": "PASS",
    }
    fe = summary[summary["species_set"] == "Fe"]["delta_alpha_over_alpha"].dropna()
    ni = summary[summary["species_set"] == "Ni"]["delta_alpha_over_alpha"].dropna()
    if fe.size and ni.size:
        if np.abs(np.median(fe) - np.median(ni)) > 3 * np.sqrt(
            np.nanmedian(summary["sigma"].to_numpy()) ** 2
        ):
            gating["species_tension"] = "FAIL: species tension"
    if result.dof > 0 and result.chi2 / result.dof > 2.0:
        gating["chi2_target"] = "FAIL: error model / distortion insufficient"
    cont_spread = summary.groupby("continuum_setting")["delta_alpha_over_alpha"].median()
    if cont_spread.size > 1 and cont_spread.max() - cont_spread.min() > 3 * result.delta_alpha_err:
        gating["continuum_stability"] = "FAIL: continuum sensitivity dominates"

    build_exec_summary(run_dir, result_payload, gating, summary)
    # Convert config to dict with string paths for yaml serialization
    config_dict = config.model_dump()
    for key, value in config_dict.items():
        if isinstance(value, Path):
            config_dict[key] = str(value)
    (run_dir / "config_used.yaml").write_text(yaml.safe_dump(config_dict))
    return result_payload
