from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import yaml
from rich.console import Console

from wdalpha.config import RunConfig
from wdalpha.holistic import load_holistic_config, run_holistic
from wdalpha.io.fits import read_spectrum
from wdalpha.io.mast import download_hlsp, download_mast
from wdalpha.io.tables import build_atomic_table
from wdalpha.lines.fit import fit_lines
from wdalpha.preprocess.continuum import build_mask, normalize_spectrum
from wdalpha.inference.mm_regression import infer_delta_alpha, save_result
from wdalpha.report.plots import plot_q_vs_shift, plot_residuals

console = Console()


def _timestamp() -> str:
    return datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")


def cmd_download(args: argparse.Namespace) -> None:
    out_dir = Path(args.out)
    if args.mode == "hlsp":
        download_hlsp(args.target, out_dir)
    else:
        download_mast(args.target, out_dir)


def cmd_build_atomic(args: argparse.Namespace) -> None:
    supplemental = [Path(p) for p in (args.supp or [])]
    build_atomic_table(Path(args.nistsource), supplemental, Path(args.out))


def cmd_fit_lines(args: argparse.Namespace) -> None:
    spectrum = read_spectrum(Path(args.spectrum))
    atomic = pd.read_csv(args.atomic)
    mask = build_mask(spectrum.wavelength, atomic["wavelength_aa"], args.mask_width)
    norm_flux, norm_error, continuum = normalize_spectrum(
        spectrum.wavelength, spectrum.flux, spectrum.error, mask, args.spline_s
    )
    line_table = [
        {
            "line_id": row.get("line_id", row.get("species")),
            "species": row["species"],
            "lambda0": row["wavelength_aa"],
            "lambda_expected": row["wavelength_aa"],
            "edge_flag": False,
        }
        for _, row in atomic.iterrows()
    ]
    results = fit_lines(
        spectrum.wavelength,
        norm_flux,
        norm_error,
        line_table,
        args.window,
        args.max_components,
        edge_buffer_aa=0.2,
        jobs=4,
    )
    rows = [
        {
            "line_id": res.line_id,
            "species": res.species,
            "lambda_obs": res.lambda_obs,
            "sigma_lambda_obs": res.sigma_lambda_obs,
            "chi2": res.chi2_local,
            "n_components": res.n_components,
        }
        for res in results
    ]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)


def cmd_infer(args: argparse.Namespace) -> None:
    lines = pd.read_csv(args.lines)
    atomic = pd.read_csv(args.atomic)
    result = infer_delta_alpha(lines, atomic, include_lab_uncertainty=not args.no_lab_unc)
    save_result(result, Path(args.out))


def cmd_report(args: argparse.Namespace) -> None:
    run_dir = Path(args.run)
    lines = pd.read_csv(run_dir / "line_measurements.csv")
    atomic = pd.read_csv(run_dir / "atomic_used.csv")
    merged = lines.merge(atomic, on="species", suffixes=("_obs", "_lab"))
    merged["delta_lambda"] = merged["lambda_obs"] - merged["wavelength_aa"]
    plot_q_vs_shift(merged, run_dir / "diagnostic_plots" / "q_vs_shift.png")
    residuals = merged["delta_lambda"].to_numpy()
    plot_residuals(residuals, run_dir / "diagnostic_plots" / "residuals.png")


def cmd_run(args: argparse.Namespace) -> None:
    config = RunConfig()
    np.random.seed(config.seed)
    run_dir = Path(args.out) / f"run_{_timestamp()}"
    run_dir.mkdir(parents=True, exist_ok=True)

    data_root = Path(args.data_root)
    spectra_dir = data_root / "spectra"
    spectra_dir.mkdir(parents=True, exist_ok=True)

    console.log("Downloading HLSP data")
    spectra = download_hlsp(args.target, spectra_dir)
    if not spectra:
        raise RuntimeError("No spectra downloaded")

    console.log("Loading atomic table")
    atomic = pd.read_csv(args.atomic)
    atomic.to_csv(run_dir / "atomic_used.csv", index=False)

    spectrum = read_spectrum(Path(spectra[0]))
    windows = build_line_windows(atomic, config.linefit.window_aa)
    mask = build_mask(spectrum.wavelength, atomic["wavelength_aa"], config.continuum.mask_width_aa)
    norm_flux, norm_error, continuum = normalize_spectrum(
        spectrum.wavelength, spectrum.flux, spectrum.error, mask, config.continuum.spline_s
    )
    results = fit_lines(
        spectrum.wavelength,
        norm_flux,
        norm_error,
        windows,
        max_components=config.linefit.max_components,
    )
    rows = [
        {
            "line_id": res.line_id,
            "species": res.species,
            "lambda_obs": res.lambda_obs,
            "sigma_lambda_obs": res.sigma_lambda_obs,
            "chi2": res.chi2,
            "n_components": res.n_components,
        }
        for res in results
    ]
    lines_path = run_dir / "line_measurements.csv"
    pd.DataFrame(rows).to_csv(lines_path, index=False)

    result = infer_delta_alpha(pd.DataFrame(rows), atomic, include_lab_uncertainty=config.inference.include_lab_uncertainty)
    save_result(result, run_dir / "inferred_alpha.json", metadata={"target": args.target})

    plot_q_vs_shift(
        pd.DataFrame(rows).merge(atomic, on="species", suffixes=("_obs", "_lab")),
        run_dir / "diagnostic_plots" / "q_vs_shift.png",
    )
    plot_residuals(
        pd.DataFrame(rows)["lambda_obs"].to_numpy(),
        run_dir / "diagnostic_plots" / "residuals.png",
    )

    config_path = run_dir / "config_used.yaml"
    config_path.write_text(yaml.safe_dump(config.model_dump()))
    (run_dir / "logs.txt").write_text("Run completed")


def cmd_holistic_run(args: argparse.Namespace) -> None:
    default_path = Path("default_config.yaml")
    config = load_holistic_config(
        default_path,
        Path(args.config) if args.config else None,
        {
            "target": args.target,
            "data_root": args.data_root,
            "results_root": args.out,
            "atomic_path": args.atomic,
            "do_e230h": args.do_e230h,
            "jobs": args.jobs,
        },
    )
    if args.spectrum:
        config.extra["spectrum"] = args.spectrum
    if args.use_hlsp_linelists:
        config.use_hlsp_linelists = Path(args.use_hlsp_linelists)
    run_dir = Path(args.out)
    run_holistic(config, run_dir, sys.argv)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="wdalpha")
    sub = parser.add_subparsers(dest="command", required=True)

    download = sub.add_parser("download", help="Download spectra")
    download.add_argument("--target", required=True)
    download.add_argument("--mode", choices=["hlsp", "mast"], default="hlsp")
    download.add_argument("--out", default="data/spectra")
    download.set_defaults(func=cmd_download)

    atomic = sub.add_parser("build-atomic", help="Build combined atomic table")
    atomic.add_argument("--nistsource", required=True)
    atomic.add_argument("--supp", nargs="*")
    atomic.add_argument("--out", required=True)
    atomic.set_defaults(func=cmd_build_atomic)

    fit = sub.add_parser("fit-lines", help="Fit spectral lines")
    fit.add_argument("--spectrum", required=True)
    fit.add_argument("--atomic", required=True)
    fit.add_argument("--out", required=True)
    fit.add_argument("--window", type=float, default=0.6)
    fit.add_argument("--mask-width", type=float, default=0.4)
    fit.add_argument("--spline-s", type=float, default=0.001)
    fit.add_argument("--max-components", type=int, default=2)
    fit.set_defaults(func=cmd_fit_lines)

    infer = sub.add_parser("infer", help="Infer delta alpha/alpha")
    infer.add_argument("--lines", required=True)
    infer.add_argument("--atomic", required=True)
    infer.add_argument("--out", required=True)
    infer.add_argument("--no-lab-unc", action="store_true")
    infer.set_defaults(func=cmd_infer)

    report = sub.add_parser("report", help="Generate diagnostic plots")
    report.add_argument("--run", required=True)
    report.set_defaults(func=cmd_report)

    run = sub.add_parser("run", help="Fast path pipeline")
    run.add_argument("--target", required=True)
    run.add_argument("--data-root", default="data")
    run.add_argument("--atomic", required=True)
    run.add_argument("--out", default="results")
    run.set_defaults(func=cmd_run)

    holistic = sub.add_parser("holistic-run", help="Holistic end-to-end pipeline")
    holistic.add_argument("--target", required=True)
    holistic.add_argument("--data-root", default="data")
    holistic.add_argument("--out", required=True)
    holistic.add_argument("--spectrum")
    holistic.add_argument("--atomic")
    holistic.add_argument("--use-hlsp-linelist")
    holistic.add_argument("--do-e230h", action="store_true")
    holistic.add_argument("--jobs", type=int, default=4)
    holistic.add_argument("--config")
    holistic.set_defaults(func=cmd_holistic_run)

    return parser


def main(argv: List[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
