#!/usr/bin/env python3
"""
Run Holistic V3 Pipeline - Blend-aware line-group modeling

This script runs the complete V3 pipeline:
1. Load spectrum and atomic data
2. Estimate global shift z_guess
3. Run blend-aware group fitting
4. Select analysis/gold sets
5. Calibrate uncertainties via injection/recovery
6. Run baseline inference with calibrated sigma
7. Run stress harness
8. Generate diagnostic plots and executive summary
"""
import sys
import json
import time
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from wdalpha.holistic_v3 import (
    V3Config,
    setup_logging, log_progress,
    estimate_z_guess_xcorr, estimate_z_guess_hlsp,
    fit_all_groups,
    select_analysis_set,
    calibrate_uncertainties,
    infer_delta_alpha,
    run_stress_harness,
    leave_one_out_analysis,
    generate_attrition_audit,
    evaluate_gates,
    generate_exec_summary,
    generate_identifiability_report,
)
from wdalpha.io.fits import read_spectrum


def load_hlsp_linelist(path: Path) -> pd.DataFrame:
    """Load HLSP linelist from text file."""
    # Skip header lines starting with #
    with open(path, 'r') as f:
        lines = f.readlines()

    # Find data lines
    data_lines = []
    for line in lines:
        if not line.startswith('#') and line.strip():
            data_lines.append(line.strip())

    if not data_lines:
        return pd.DataFrame()

    # Parse columns based on HLSP format
    # WAVEOBS WOBSUNC EQWIDTH EQWUNC METAL ION WAVEREST WRESTUNC VELOCITY VELUNC COMMS ORIGIN
    records = []
    for line in data_lines:
        parts = line.split()
        if len(parts) >= 7:
            try:
                records.append({
                    'WAVEOBS': float(parts[0]),
                    'WOBSUNC': float(parts[1]) / 1000,  # mA to A
                    'EQWIDTH': float(parts[2]),
                    'METAL': parts[4],
                    'ION': parts[5],
                    'WAVEREST': float(parts[6]) if parts[6] != '0.0000' else 0,
                })
            except:
                continue

    return pd.DataFrame(records)


def plot_diagnostics(
    measurements_df: pd.DataFrame,
    analysis_df: pd.DataFrame,
    stress_df: pd.DataFrame,
    injection_df: pd.DataFrame,
    loo_df: pd.DataFrame,
    baseline_result: dict,
    sigma_inflation: float,
    output_dir: Path
):
    """Generate diagnostic plots."""
    plot_dir = output_dir / 'diagnostic_plots'
    plot_dir.mkdir(exist_ok=True)

    # 1. Centroid uncertainty distribution
    if 'sigma_lambda_obs' in analysis_df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(analysis_df['sigma_lambda_obs'] * 1000, bins=30, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Centroid uncertainty (mA)')
        ax.set_ylabel('Count')
        ax.set_title('V3 Analysis Set: Centroid Uncertainties')
        ax.axvline(50, color='red', linestyle='--', label='Old threshold')
        plt.tight_layout()
        plt.savefig(plot_dir / 'centroid_uncertainties.png', dpi=150)
        plt.close()

    # 2. Chi2 distribution
    if 'reduced_chi2_local' in analysis_df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        chi2_vals = analysis_df['reduced_chi2_local']
        ax.hist(chi2_vals[chi2_vals < 20], bins=30, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Reduced chi2 (local)')
        ax.set_ylabel('Count')
        ax.set_title('V3 Analysis Set: Local Fit Quality')
        plt.tight_layout()
        plt.savefig(plot_dir / 'chi2_distribution.png', dpi=150)
        plt.close()

    # 3. Number of components distribution
    if 'n_components' in analysis_df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(analysis_df['n_components'], bins=range(1, 6), edgecolor='black', alpha=0.7, align='left')
        ax.set_xlabel('Number of components fitted')
        ax.set_ylabel('Count')
        ax.set_title('V3 Analysis Set: Blend Model Complexity')
        ax.set_xticks(range(1, 5))
        plt.tight_layout()
        plt.savefig(plot_dir / 'component_counts.png', dpi=150)
        plt.close()

    # 4. Species distribution
    if 'species' in analysis_df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        species_counts = analysis_df['species'].value_counts()
        ax.bar(species_counts.index, species_counts.values, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Species')
        ax.set_ylabel('Count')
        ax.set_title('V3 Analysis Set: Species Distribution')
        plt.tight_layout()
        plt.savefig(plot_dir / 'species_distribution.png', dpi=150)
        plt.close()

    # 5. Pull distribution from injection/recovery
    if len(injection_df) > 0 and 'pull' in injection_df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        pulls = injection_df['pull'].values
        ax.hist(pulls, bins=30, edgecolor='black', alpha=0.7, density=True)
        # Overlay expected Gaussian
        x = np.linspace(-4, 4, 100)
        ax.plot(x, np.exp(-0.5 * x**2) / np.sqrt(2 * np.pi), 'r-', lw=2, label='N(0,1)')
        ax.set_xlabel('Pull')
        ax.set_ylabel('Density')
        ax.set_title(f'V3 Injection/Recovery: Pull Distribution (inflation={sigma_inflation:.2f})')
        ax.legend()
        plt.tight_layout()
        plt.savefig(plot_dir / 'pulls_injection.png', dpi=150)
        plt.close()

    # 6. Stress test results
    if len(stress_df) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot delta_alpha for different configurations
        for species in ['Fe', 'Ni', 'both']:
            subset = stress_df[
                (stress_df['species'] == species) &
                (stress_df['wavelength_split'] == 'full') &
                (stress_df['likelihood'] == 'student_t')
            ]
            if len(subset) > 0:
                x = subset['distortion_degree'].values
                y = subset['delta_alpha'].values
                yerr = subset['delta_alpha_err'].values * sigma_inflation
                ax.errorbar(x + {'Fe': -0.1, 'Ni': 0, 'both': 0.1}[species],
                           y, yerr=yerr, fmt='o-', label=species, capsize=3)

        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Distortion degree')
        ax.set_ylabel('delta_alpha/alpha')
        ax.set_title('V3 Stress Test: delta_alpha vs Distortion Degree')
        ax.legend()
        ax.set_xticks([0, 1, 2])
        plt.tight_layout()
        plt.savefig(plot_dir / 'stress_distortion.png', dpi=150)
        plt.close()

    # 7. Wavelength coverage
    if 'lambda0_ang' in analysis_df.columns:
        fig, ax = plt.subplots(figsize=(12, 5))

        fe_df = analysis_df[analysis_df['species'].str.contains('Fe')]
        ni_df = analysis_df[analysis_df['species'].str.contains('Ni')]

        ax.hist(fe_df['lambda0_ang'], bins=20, alpha=0.6, label=f'Fe V ({len(fe_df)})', color='blue')
        ax.hist(ni_df['lambda0_ang'], bins=20, alpha=0.6, label=f'Ni V ({len(ni_df)})', color='red')
        ax.set_xlabel('Rest wavelength (A)')
        ax.set_ylabel('Count')
        ax.set_title('V3 Analysis Set: Wavelength Coverage')
        ax.legend()
        plt.tight_layout()
        plt.savefig(plot_dir / 'wavelength_coverage.png', dpi=150)
        plt.close()

    # 8. LOO influence
    if len(loo_df) > 0:
        fig, ax = plt.subplots(figsize=(10, 5))
        top_loo = loo_df.head(20)
        ax.barh(range(len(top_loo)), top_loo['influence'].values[::-1])
        ax.set_yticks(range(len(top_loo)))
        ax.set_yticklabels([f"{r['line_id']} ({r['species']})" for _, r in top_loo.iloc[::-1].iterrows()])
        ax.set_xlabel('|delta from baseline|')
        ax.set_title('V3: Most Influential Lines (LOO Analysis)')
        plt.tight_layout()
        plt.savefig(plot_dir / 'loo_influence.png', dpi=150)
        plt.close()

    log_progress(f"Diagnostic plots saved to {plot_dir}")


def main():
    start_time = time.time()

    # Configuration
    config = V3Config(
        seed=42,
        target="G191-B2B",
        data_root=Path("data"),
        results_root=Path("results/holistic_v3"),
        jobs=4,

        # Group fitting
        window_half_width_aa=0.4,
        max_components=3,
        min_component_depth=0.02,
        shared_width=True,

        # Relaxed analysis thresholds
        analysis_sigma_max=0.10,
        analysis_chi2_max=50.0,
        analysis_snr_min=2.0,
        analysis_asymmetry_max=2.0,
        analysis_contam_max=0.9,
        analysis_excluded_flags=['EDGE'],

        # Inference
        use_student_t=True,
        student_t_df=4.0,
        distortion_degree=1,
        condition_number_threshold=100.0,
        auto_orthogonalize=True,
        ridge_lambda=1e-6,

        # Calibration
        calibration_n_trials=50,
    )

    # Create output directory
    run_dir = config.results_root
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "diagnostic_plots").mkdir(exist_ok=True)
    (run_dir / "line_groups").mkdir(exist_ok=True)

    # Setup logging
    setup_logging(run_dir)
    log_progress("=" * 60)
    log_progress("Holistic V3 Pipeline Started")
    log_progress("=" * 60)

    np.random.seed(config.seed)

    # Paths
    spectrum_path = config.data_root / "hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_coadd-spec.fits"
    atomic_path = config.data_root / "atomic/analysis_lines_clean.csv"
    hlsp_linelist_path = config.data_root / "hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_linelist.txt"

    log_progress(f"Spectrum: {spectrum_path}")
    log_progress(f"Atomic table: {atomic_path}")
    log_progress(f"HLSP linelist: {hlsp_linelist_path}")

    # =========================================================================
    # Load data
    # =========================================================================
    log_progress("Loading spectrum...")
    spectrum = read_spectrum(spectrum_path)
    wavelength = spectrum.wavelength
    flux = spectrum.flux
    error = spectrum.error
    log_progress(f"Spectrum: {len(wavelength)} points, {wavelength.min():.2f}-{wavelength.max():.2f} A")

    log_progress("Loading atomic table...")
    atomic_df = pd.read_csv(atomic_path)
    # Standardize column names
    if 'q_cm-1' in atomic_df.columns and 'q_cm1' not in atomic_df.columns:
        atomic_df['q_cm1'] = atomic_df['q_cm-1']
    log_progress(f"Atomic table: {len(atomic_df)} lines")

    log_progress("Loading HLSP linelist...")
    hlsp_df = None
    if hlsp_linelist_path.exists():
        hlsp_df = load_hlsp_linelist(hlsp_linelist_path)
        log_progress(f"HLSP linelist: {len(hlsp_df)} lines")
    else:
        log_progress("HLSP linelist not found, proceeding without it")

    # =========================================================================
    # Estimate z_guess
    # =========================================================================
    log_progress("Estimating global shift z_guess...")

    z_guess_info = {}
    try:
        # Try HLSP method first if available
        if hlsp_df is not None and len(hlsp_df) > 0:
            z_guess, method, z_unc = estimate_z_guess_hlsp(
                hlsp_df, atomic_df['lambda0_ang'].values
            )
            z_guess_info = {'z_guess': z_guess, 'method': method, 'uncertainty': z_unc}
            log_progress(f"z_guess from HLSP: {z_guess:.6f} +/- {z_unc:.6f}")
    except Exception as e:
        log_progress(f"HLSP z_guess failed: {e}")

    if not z_guess_info:
        # Fall back to cross-correlation
        z_guess, method, z_unc = estimate_z_guess_xcorr(
            wavelength, flux, atomic_df['lambda0_ang'].values
        )
        z_guess_info = {'z_guess': z_guess, 'method': method, 'uncertainty': z_unc}
        log_progress(f"z_guess from xcorr: {z_guess:.6f} +/- {z_unc:.6f}")

    # Save z_guess
    with open(run_dir / 'z_guess.json', 'w') as f:
        json.dump(z_guess_info, f, indent=2)

    z_guess = z_guess_info['z_guess']

    # =========================================================================
    # Normalize spectrum
    # =========================================================================
    log_progress("Normalizing spectrum...")
    from wdalpha.preprocess.continuum import build_continuum

    line_centers = [row['lambda0_ang'] * (1 + z_guess) for _, row in atomic_df.iterrows()]
    norm_flux, norm_error, continuum, mask = build_continuum(
        wavelength, flux, error, line_centers,
        mask_width_aa=0.4, spline_s=0.001
    )
    log_progress("Spectrum normalized")

    # =========================================================================
    # Group fitting
    # =========================================================================
    log_progress("Starting group fitting...")
    measurements_df = fit_all_groups(
        wavelength, norm_flux, norm_error,
        atomic_df, hlsp_df,
        z_guess, config
    )

    # Save all measurements
    measurements_df.to_csv(run_dir / 'line_measurements_all.csv', index=False)
    log_progress(f"Saved {len(measurements_df)} measurements")

    # =========================================================================
    # Select analysis and gold sets
    # =========================================================================
    log_progress("Selecting analysis and gold sets...")
    analysis_df, gold_df = select_analysis_set(measurements_df, config)

    analysis_df.to_csv(run_dir / 'line_measurements_analysis.csv', index=False)
    gold_df.to_csv(run_dir / 'line_measurements_gold.csv', index=False)

    n_analysis = len(analysis_df)
    n_gold = len(gold_df)
    log_progress(f"Analysis set: {n_analysis} lines")
    log_progress(f"Gold set: {n_gold} lines")

    if n_analysis < 10:
        log_progress("WARNING: Very few analysis lines, results may be unreliable")

    # =========================================================================
    # Generate attrition audit
    # =========================================================================
    log_progress("Generating attrition audit...")
    generate_attrition_audit(
        atomic_df, measurements_df, analysis_df, gold_df,
        config, run_dir
    )

    # =========================================================================
    # Calibrate uncertainties
    # =========================================================================
    log_progress("Calibrating uncertainties...")
    if n_analysis >= 10:
        sigma_inflation, injection_df = calibrate_uncertainties(
            analysis_df, config, n_trials=config.calibration_n_trials
        )
        injection_df.to_csv(run_dir / 'injection_recovery.csv', index=False)

        # Generate injection/recovery summary
        with open(run_dir / 'injection_recovery.md', 'w') as f:
            f.write("# Injection/Recovery Calibration\n\n")
            f.write(f"- N trials: {len(injection_df)}\n")
            f.write(f"- Median |pull|: {np.median(np.abs(injection_df['pull'])):.3f}\n")
            f.write(f"- Target median |pull|: {config.target_median_pull:.3f}\n")
            f.write(f"- **Inflation factor: {sigma_inflation:.3f}**\n")
    else:
        sigma_inflation = 5.0  # Conservative default
        injection_df = pd.DataFrame()
        log_progress(f"Too few lines for calibration, using default inflation={sigma_inflation}")

    # =========================================================================
    # Baseline inference
    # =========================================================================
    log_progress("Running baseline inference...")
    baseline_result = infer_delta_alpha(
        analysis_df['lambda_obs'].values,
        analysis_df['sigma_lambda_obs'].values,
        analysis_df['q_cm1'].values,
        analysis_df['lambda0_ang'].values,
        config,
        sigma_inflation=sigma_inflation
    )

    if baseline_result['converged']:
        da = baseline_result['delta_alpha']
        da_err = baseline_result['delta_alpha_err'] * sigma_inflation
        log_progress(f"Baseline: delta_alpha/alpha = {da:.4e} +/- {da_err:.4e}")
        log_progress(f"  chi2/dof = {baseline_result['chi2']:.1f}/{baseline_result['dof']} = {baseline_result['reduced_chi2']:.2f}")
        log_progress(f"  Condition number = {baseline_result['condition_number']:.1f}")
        log_progress(f"  Orthogonalized = {baseline_result['orthogonalized']}")
    else:
        log_progress(f"Baseline inference failed: {baseline_result.get('message', 'unknown')}")

    # Save baseline result
    with open(run_dir / 'inferred_alpha.json', 'w') as f:
        # Convert numpy types for JSON serialization
        result_dict = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
                       for k, v in baseline_result.items()}
        result_dict['sigma_inflation'] = float(sigma_inflation)
        json.dump(result_dict, f, indent=2)

    # Generate identifiability report
    generate_identifiability_report(baseline_result, config, run_dir)

    # =========================================================================
    # Stress harness
    # =========================================================================
    log_progress("Running stress harness...")
    stress_df = run_stress_harness(analysis_df, config, sigma_inflation)
    stress_df.to_csv(run_dir / 'summary_delta_alpha.csv', index=False)
    log_progress(f"Stress harness: {len(stress_df)} configurations")

    # =========================================================================
    # Leave-one-out analysis
    # =========================================================================
    log_progress("Running leave-one-out analysis...")
    if n_analysis <= 100:  # Only if manageable size
        loo_df = leave_one_out_analysis(analysis_df, config, sigma_inflation, baseline_result)
        loo_df.to_csv(run_dir / 'loo_analysis.csv', index=False)
    else:
        loo_df = pd.DataFrame()
        log_progress("Skipping LOO (too many lines)")

    # =========================================================================
    # Evaluate gates
    # =========================================================================
    log_progress("Evaluating gates...")
    gates = evaluate_gates(analysis_df, baseline_result, stress_df, sigma_inflation, config)

    for gate_name, gate_info in gates.items():
        log_progress(f"  {gate_name}: {gate_info['status']} - {gate_info.get('note', '')}")

    # =========================================================================
    # Generate plots
    # =========================================================================
    log_progress("Generating diagnostic plots...")
    plot_diagnostics(
        measurements_df, analysis_df, stress_df,
        injection_df, loo_df, baseline_result,
        sigma_inflation, run_dir
    )

    # =========================================================================
    # Generate executive summary
    # =========================================================================
    log_progress("Generating executive summary...")
    generate_exec_summary(
        config, baseline_result, stress_df, loo_df,
        gates, sigma_inflation, z_guess_info,
        n_analysis, n_gold, run_dir
    )

    # =========================================================================
    # Final summary
    # =========================================================================
    elapsed = time.time() - start_time
    log_progress("=" * 60)
    log_progress(f"Pipeline complete in {elapsed:.1f}s")
    log_progress(f"Output directory: {run_dir}")
    log_progress("=" * 60)

    # Print summary to stdout
    print("\n" + "=" * 60)
    print("HOLISTIC V3 PIPELINE COMPLETE")
    print("=" * 60)

    if baseline_result['converged']:
        da = baseline_result['delta_alpha']
        da_err = baseline_result['delta_alpha_err'] * sigma_inflation
        print(f"\nBaseline Result ({n_analysis} lines, calibrated sigma):")
        print(f"  delta_alpha/alpha = {da:.4e} +/- {da_err:.4e}")
        print(f"  chi2/dof = {baseline_result['reduced_chi2']:.2f}")

    print("\nGates:")
    for gate_name, gate_info in gates.items():
        print(f"  {gate_name}: {gate_info['status']}")

    print(f"\nOutputs: {run_dir}")
    print("  - EXEC_SUMMARY.md")
    print("  - summary_delta_alpha.csv")
    print("  - diagnostic_plots/")


if __name__ == "__main__":
    main()
