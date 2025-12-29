# Holistic Pipeline Run Notes

## Overview

This run executed the Codex-generated wdalpha holistic pipeline on G191-B2B E140H HST/STIS
spectroscopy to measure potential variation in the fine-structure constant (Δα/α).

## Key Results

**Primary Result:**
- Δα/α = **(-6.52 ± 7.0) × 10⁻⁴**
- This is consistent with zero (|Δα/α| < 1σ)
- **No statistically significant variation detected**

**Quality Gates:** All PASS
- species_tension: PASS
- chi2_target: PASS (χ² = 18.0, dof = 14, reduced χ² = 1.29)
- continuum_stability: PASS

## Data Summary

- Target: G191-B2B white dwarf
- Spectrum: HST/STIS E140H coadd (HLSP WD-LINELIST)
- Atomic lines analyzed: 349 Fe V + Ni V transitions
- Gold sample (passing quality cuts): 10 lines
- Continuum settings tested: 3 (default, smooth, tight)

## Bug Fixes Applied

Several bugs in the Codex-generated code were fixed:

1. **Jitter objective function (model.py)**
   - Issue: least_squares with jitter in denominator has degenerate solution (jitter → ∞)
   - Fix: Switched to proper Gaussian NLL with log-variance penalty term
   - This is the mathematically correct approach per ChatGPT feedback

2. **Column naming after pandas merge (model.py, holistic.py)**
   - Issue: Columns like `lambda0_ang`, `q_cm1`, `species` get suffixed after merge
   - Fix: Added conditional column selection to handle `_obs`/`_lab` suffixes

3. **Atomic table column normalization (holistic.py)**
   - Issue: Column names `q_cm-1` and `sigma_lambda0_ang` not recognized
   - Fix: Added column name mapping in `_normalize_atomic_table()`

4. **Test design (test_inference_model.py)**
   - Issue: Original test had collinear design matrix (q correlated with λ) and SNR < 1
   - Fix: Used realistic non-monotonic q values and higher delta_alpha for adequate SNR

5. **build_influence index handling (model.py)**
   - Issue: Gold dataframe has non-contiguous indices after filtering
   - Fix: Reset index at start of function

6. **Config YAML serialization (holistic.py)**
   - Issue: Path objects can't be serialized by yaml.safe_dump
   - Fix: Convert Path objects to strings before serialization

## Artifacts Produced

```
results/holistic_run/
├── EXEC_SUMMARY.md           # Human-readable summary
├── atomic_used.csv           # Atomic data used
├── diagnostic_plots/         # Diagnostic visualizations
│   ├── influence.png
│   ├── residual_hist.png
│   ├── residual_qq.png
│   ├── residual_vs_q.png
│   ├── residual_vs_sensitivity.png
│   └── residual_vs_wavelength.png
├── hlsp_sanity_report.md     # Spectrum sanity check
├── inferred_alpha.json       # Primary result
├── influence.csv             # Leave-one-out influence analysis
├── line_diagnostics.csv      # Per-line diagnostic info
├── line_exclusions.csv       # Lines excluded from gold sample
├── line_measurements_all.csv # All line measurements
├── line_measurements_gold.csv # Gold sample lines
├── logs.txt                  # Environment info
└── summary_delta_alpha.csv   # Stress test results
```

## Interpretation

The measured Δα/α = (-6.52 ± 7.0) × 10⁻⁴ is:
- Consistent with zero at the 1σ level
- Provides an upper limit of |Δα/α| < 1.4 × 10⁻³ at 95% confidence
- Less constraining than cosmological measurements but consistent with α being constant

The relatively large error bars are due to:
- Small number of gold lines (10) after quality cuts
- High crowding in the UV spectrum leading to many blended lines
- Conservative quality criteria excluding potentially useful lines

## Environment

- Python 3.12.3
- numpy 1.26.4, scipy 1.12.0, pandas 2.2.2, astropy 6.0.1
- Robust loss function: Cauchy
- Distortion model: polynomial degree 1

## Notes for Future Work

1. Relaxing gold selection criteria could increase sample size
2. E230H analysis was not performed (skipped per user request)
3. The pipeline now uses proper NLL for jitter estimation, avoiding degeneracy issues
4. Consider adding regularization for distortion vs sensitivity collinearity (per ChatGPT)
