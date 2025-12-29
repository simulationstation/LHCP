#!/usr/bin/env python3
"""
Run Holistic V2 Pipeline - Single integrated execution

This script runs the complete redesigned pipeline:
1. Load spectrum and atomic data
2. Fit lines with existing fitter
3. Select analysis set (relaxed) and gold set (strict)
4. Run baseline inference with Student-t likelihood
5. Run stress harness
6. Run injection/recovery test
7. Generate diagnostic plots
8. Produce EXEC_SUMMARY with gates
"""
import sys
import time
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from wdalpha.holistic_v2 import (
    V2Config, InferenceResultV2,
    setup_logging, log_progress,
    infer_delta_alpha_v2,
    select_analysis_set,
    run_stress_harness,
    run_injection_recovery,
    generate_exec_summary,
)
from wdalpha.io.fits import read_spectrum
from wdalpha.lines.fit import fit_lines, apply_quality_flags, estimate_z_guess
from wdalpha.preprocess.continuum import build_continuum


def main():
    start_time = time.time()

    # Configuration
    config = V2Config(
        seed=42,
        target="G191-B2B",
        data_root=Path("data"),
        results_root=Path("results/holistic_v2"),
        jobs=4,
        # Relaxed analysis thresholds
        analysis_sigma_max=0.05,
        analysis_chi2_max=20.0,
        analysis_snr_min=2.0,
        analysis_asymmetry_max=1.0,
        analysis_contam_max=0.8,
        analysis_excluded_flags=['EDGE'],
        # Inference settings
        use_student_t=True,
        student_t_df=4.0,
        distortion_degree=1,
        auto_orthogonalize=True,
    )

    # Create output directory
    run_dir = config.results_root
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "diagnostic_plots").mkdir(exist_ok=True)

    # Setup logging
    log_path = setup_logging(run_dir)
    log_progress(log_path, "=" * 60)
    log_progress(log_path, "Holistic V2 Pipeline Started")
    log_progress(log_path, "=" * 60)

    np.random.seed(config.seed)

    # Locate files
    spectrum_path = config.data_root / "hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_coadd-spec.fits"
    atomic_path = config.data_root / "atomic/analysis_lines_clean.csv"

    log_progress(log_path, f"Spectrum: {spectrum_path}")
    log_progress(log_path, f"Atomic table: {atomic_path}")

    # Load data
    log_progress(log_path, "Loading spectrum...")
    spectrum = read_spectrum(spectrum_path)
    log_progress(log_path, f"Spectrum loaded: {len(spectrum.wavelength)} points, {spectrum.wavelength.min():.2f}-{spectrum.wavelength.max():.2f} Å")

    log_progress(log_path, "Loading atomic table...")
    atomic = pd.read_csv(atomic_path)
    if 'q_cm-1' in atomic.columns and 'q_cm1' not in atomic.columns:
        atomic['q_cm1'] = atomic['q_cm-1']
    if 'sigma_lambda0_ang' in atomic.columns and 'lambda0_unc_ang' not in atomic.columns:
        atomic['lambda0_unc_ang'] = atomic['sigma_lambda0_ang']
    if 'line_id' not in atomic.columns:
        atomic['line_id'] = [f'line_{i:04d}' for i in range(len(atomic))]

    log_progress(log_path, f"Atomic table loaded: {len(atomic)} lines")

    # Save atomic used
    atomic.to_csv(run_dir / "atomic_used.csv", index=False)

    # Check wavelength coverage
    in_range = atomic[(atomic['lambda0_ang'] >= spectrum.wavelength.min()) &
                       (atomic['lambda0_ang'] <= spectrum.wavelength.max())]
    log_progress(log_path, f"Lines in E140H range: {len(in_range)}")

    # Estimate z_guess from spectrum
    z_guess = estimate_z_guess(
        spectrum.wavelength,
        spectrum.flux,
        atomic['lambda0_ang'],
        config.window_aa
    )
    log_progress(log_path, f"Estimated z_guess: {z_guess:.6f}")

    # Build line table
    log_progress(log_path, "Building line table...")
    line_table = []
    for _, row in atomic.iterrows():
        lambda0 = float(row['lambda0_ang'])
        lambda_expected = lambda0 * (1.0 + z_guess)
        if lambda_expected < spectrum.wavelength.min() or lambda_expected > spectrum.wavelength.max():
            continue
        edge_flag = (lambda_expected < spectrum.wavelength.min() + config.edge_buffer_aa or
                     lambda_expected > spectrum.wavelength.max() - config.edge_buffer_aa)
        line_table.append({
            'line_id': row['line_id'],
            'species': row['species'],
            'lambda0': lambda0,
            'lambda_expected': lambda_expected,
            'edge_flag': edge_flag,
        })

    log_progress(log_path, f"Lines to fit: {len(line_table)}")

    # Normalize spectrum
    log_progress(log_path, "Normalizing spectrum...")
    line_centers = [row['lambda_expected'] for row in line_table]
    norm_flux, norm_error, continuum, mask = build_continuum(
        spectrum.wavelength,
        spectrum.flux,
        spectrum.error,
        line_centers,
        mask_width_aa=0.4,
        spline_s=0.001,
    )
    log_progress(log_path, "Spectrum normalized")

    # Fit lines
    log_progress(log_path, "Fitting lines...")
    fit_start = time.time()
    results = fit_lines(
        spectrum.wavelength,
        norm_flux,
        norm_error,
        line_table,
        config.window_aa,
        config.max_components,
        config.edge_buffer_aa,
        jobs=config.jobs,
    )
    fit_time = time.time() - fit_start
    log_progress(log_path, f"Line fitting complete: {len(results)} fits in {fit_time:.1f}s")

    # Apply quality flags
    for res in results:
        apply_quality_flags(
            res,
            min_snr=config.analysis_snr_min,
            max_reduced_chi2=config.analysis_chi2_max,
            max_asymmetry=config.analysis_asymmetry_max,
            saturation_threshold=0.1,
        )

    # Build all measurements dataframe
    fitted_ids = {res.line_id for res in results}
    all_rows = []

    for res in results:
        flags_str = ';'.join(sorted(set(res.flags))) if res.flags else ''
        all_rows.append({
            'line_id': res.line_id,
            'species': res.species,
            'lambda0_ang': res.lambda0,
            'lambda_obs': res.lambda_obs,
            'sigma_lambda_obs': res.sigma_lambda_obs,
            'depth': res.depth,
            'width': res.width,
            'continuum_c0': res.continuum_c0,
            'continuum_c1': res.continuum_c1,
            'chi2_local': res.chi2_local,
            'reduced_chi2_local': res.reduced_chi2_local,
            'asymmetry': res.asymmetry,
            'snr': res.snr,
            'window_contam': res.window_contam,
            'flags': flags_str,
            'n_components': res.n_components,
        })

    # Add missing lines
    for row in line_table:
        if row['line_id'] not in fitted_ids:
            all_rows.append({
                'line_id': row['line_id'],
                'species': row['species'],
                'lambda0_ang': row['lambda0'],
                'lambda_obs': np.nan,
                'sigma_lambda_obs': np.nan,
                'depth': np.nan,
                'width': np.nan,
                'continuum_c0': np.nan,
                'continuum_c1': np.nan,
                'chi2_local': np.nan,
                'reduced_chi2_local': np.nan,
                'asymmetry': np.nan,
                'snr': np.nan,
                'window_contam': np.nan,
                'flags': 'FIT_FAILED',
                'n_components': 0,
            })

    all_lines = pd.DataFrame(all_rows)
    all_lines = all_lines.sort_values('lambda0_ang').reset_index(drop=True)

    # Merge with atomic data for q values
    all_lines = all_lines.merge(
        atomic[['line_id', 'q_cm1', 'lambda0_unc_ang']],
        on='line_id',
        how='left'
    )

    # Save all measurements
    all_lines.to_csv(run_dir / "line_measurements_all.csv", index=False)

    n_converged = all_lines['lambda_obs'].notna().sum()
    n_failed = all_lines['lambda_obs'].isna().sum()
    log_progress(log_path, f"Converged: {n_converged}, Failed: {n_failed}")

    # Select analysis and gold sets
    log_progress(log_path, "Selecting analysis and gold sets...")
    analysis_set, gold_set = select_analysis_set(all_lines, config)

    n_analysis = len(analysis_set)
    n_gold = len(gold_set)
    log_progress(log_path, f"Analysis set: {n_analysis} lines")
    log_progress(log_path, f"Gold set: {n_gold} lines")

    # Save sets
    analysis_set.to_csv(run_dir / "line_measurements_analysis.csv", index=False)
    gold_set.to_csv(run_dir / "line_measurements_gold.csv", index=False)

    # Check minimum lines gate
    if n_analysis < config.min_analysis_lines:
        log_progress(log_path, f"WARNING: N_analysis ({n_analysis}) < min ({config.min_analysis_lines})")

    # Run baseline inference on analysis set
    log_progress(log_path, "Running baseline inference (Student-t, both species)...")
    baseline_result = infer_delta_alpha_v2(analysis_set, atomic, config, use_student_t=True)

    log_progress(log_path, f"Baseline: Δα/α = {baseline_result.delta_alpha:.4e} ± {baseline_result.delta_alpha_err:.4e}")
    log_progress(log_path, f"  z0 = {baseline_result.z0:.4e}, jitter = {baseline_result.jitter:.4e}")
    log_progress(log_path, f"  χ² = {baseline_result.chi2:.2f}, dof = {baseline_result.dof}")
    log_progress(log_path, f"  Condition number = {baseline_result.condition_number:.1f}")
    log_progress(log_path, f"  Orthogonalized = {baseline_result.orthogonalized}")

    # Save baseline result
    result_dict = {
        'delta_alpha': baseline_result.delta_alpha,
        'delta_alpha_err': baseline_result.delta_alpha_err,
        'z0': baseline_result.z0,
        'z0_err': baseline_result.z0_err,
        'jitter': baseline_result.jitter,
        'chi2': baseline_result.chi2,
        'dof': baseline_result.dof,
        'n_lines': baseline_result.n_lines,
        'condition_number': baseline_result.condition_number,
        'orthogonalized': baseline_result.orthogonalized,
        'likelihood': baseline_result.likelihood_type,
    }

    import json
    with open(run_dir / "inferred_alpha.json", 'w') as f:
        json.dump(result_dict, f, indent=2)

    # Run stress harness
    log_progress(log_path, "Running stress harness...")
    stress_results = run_stress_harness(analysis_set, atomic, config, log_path)
    stress_results.to_csv(run_dir / "summary_delta_alpha.csv", index=False)
    log_progress(log_path, f"Stress harness complete: {len(stress_results)} configurations")

    # Run injection/recovery test
    log_progress(log_path, "Running injection/recovery test...")
    injection_results = run_injection_recovery(analysis_set, atomic, config)
    injection_results.to_csv(run_dir / "injection_recovery.csv", index=False)
    log_progress(log_path, "Injection/recovery complete")

    # Generate diagnostic plots
    log_progress(log_path, "Generating diagnostic plots...")

    # Plot 1: Residuals vs wavelength
    fig, ax = plt.subplots(figsize=(10, 6))
    valid = analysis_set[analysis_set['lambda_obs'].notna()].copy()
    if len(valid) > 0:
        # Compute residuals using baseline model
        from wdalpha.holistic_v2 import omega0_from_lambda, sensitivity_factor

        lambda0 = valid['lambda0_ang'].values
        lambda_obs = valid['lambda_obs'].values
        q_cm1 = valid['q_cm1'].values
        z_obs = lambda_obs / lambda0 - 1.0

        omega0 = omega0_from_lambda(lambda0)
        sensitivity = sensitivity_factor(q_cm1, omega0)
        z_model = baseline_result.z0 + sensitivity * baseline_result.delta_alpha
        residual = z_obs - z_model

        ax.scatter(lambda0, residual * 1e6, alpha=0.6, s=20)
        ax.axhline(0, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('Wavelength (Å)')
        ax.set_ylabel('Residual (ppm)')
        ax.set_title('Residuals vs Wavelength')

    plt.tight_layout()
    plt.savefig(run_dir / "diagnostic_plots/residual_vs_wavelength.png", dpi=150)
    plt.close()

    # Plot 2: Residuals vs q
    fig, ax = plt.subplots(figsize=(10, 6))
    if len(valid) > 0:
        ax.scatter(q_cm1, residual * 1e6, alpha=0.6, s=20)
        ax.axhline(0, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('q (cm⁻¹)')
        ax.set_ylabel('Residual (ppm)')
        ax.set_title('Residuals vs q-coefficient')

    plt.tight_layout()
    plt.savefig(run_dir / "diagnostic_plots/residual_vs_q.png", dpi=150)
    plt.close()

    # Plot 3: Histogram of residuals
    fig, ax = plt.subplots(figsize=(8, 6))
    if len(valid) > 0:
        ax.hist(residual * 1e6, bins=30, alpha=0.7, edgecolor='black')
        ax.axvline(0, color='red', linestyle='--')
        ax.set_xlabel('Residual (ppm)')
        ax.set_ylabel('Count')
        ax.set_title('Residual Distribution')

    plt.tight_layout()
    plt.savefig(run_dir / "diagnostic_plots/residual_hist.png", dpi=150)
    plt.close()

    # Plot 4: Species comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    fe_stress = stress_results[(stress_results['species'] == 'Fe') &
                                (stress_results['wavelength_split'] == 'full') &
                                stress_results['delta_alpha'].notna()]
    ni_stress = stress_results[(stress_results['species'] == 'Ni') &
                                (stress_results['wavelength_split'] == 'full') &
                                stress_results['delta_alpha'].notna()]
    both_stress = stress_results[(stress_results['species'] == 'both') &
                                  (stress_results['wavelength_split'] == 'full') &
                                  stress_results['delta_alpha'].notna()]

    y_pos = []
    labels = []
    values = []
    errors = []

    if len(fe_stress) > 0:
        y_pos.append(0)
        labels.append('Fe only')
        values.append(fe_stress['delta_alpha'].iloc[0])
        errors.append(fe_stress['delta_alpha_err'].iloc[0])

    if len(ni_stress) > 0:
        y_pos.append(1)
        labels.append('Ni only')
        values.append(ni_stress['delta_alpha'].iloc[0])
        errors.append(ni_stress['delta_alpha_err'].iloc[0])

    if len(both_stress) > 0:
        y_pos.append(2)
        labels.append('Both')
        values.append(both_stress['delta_alpha'].iloc[0])
        errors.append(both_stress['delta_alpha_err'].iloc[0])

    if len(y_pos) > 0:
        ax.errorbar(values, y_pos, xerr=errors, fmt='o', capsize=5, markersize=8)
        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel('Δα/α')
        ax.set_title('Species Comparison')

    plt.tight_layout()
    plt.savefig(run_dir / "diagnostic_plots/species_comparison.png", dpi=150)
    plt.close()

    log_progress(log_path, "Diagnostic plots saved")

    # Generate exec summary with gates
    log_progress(log_path, "Generating executive summary...")
    gates = generate_exec_summary(
        run_dir, config, baseline_result, stress_results,
        n_analysis, n_gold, injection_results
    )

    log_progress(log_path, "Gates:")
    for gate, status in gates.items():
        log_progress(log_path, f"  {gate}: {status}")

    # Final timing
    total_time = time.time() - start_time
    log_progress(log_path, "=" * 60)
    log_progress(log_path, f"Pipeline complete in {total_time:.1f}s")
    log_progress(log_path, f"Output directory: {run_dir}")
    log_progress(log_path, "=" * 60)

    # Print summary for user
    print("\n" + "=" * 60)
    print("HOLISTIC V2 PIPELINE COMPLETE")
    print("=" * 60)
    print(f"\nBaseline Result (Student-t, {n_analysis} lines):")
    print(f"  Δα/α = {baseline_result.delta_alpha:.4e} ± {baseline_result.delta_alpha_err:.4e}")
    print(f"\nGates:")
    for gate, status in gates.items():
        print(f"  {gate}: {status}")
    print(f"\nOutputs: {run_dir}")
    print(f"  - EXEC_SUMMARY.md")
    print(f"  - summary_delta_alpha.csv")
    print(f"  - inferred_alpha.json")


if __name__ == "__main__":
    main()
