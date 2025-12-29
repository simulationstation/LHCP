# EXEC_SUMMARY - Holistic V3 Pipeline

Generated: 2025-12-29 07:45:19

## Configuration

- Target: G191-B2B
- Seed: 42
- Group fitting: max_components=3, shared_width=True
- Distortion degree: 1
- Likelihood: Student-t (df=4.0)
- Uncertainty inflation: 0.500

## Global Shift Estimation

- Method: hlsp
- z_guess: 0.000080
- Uncertainty: 0.000000

## Baseline Result (calibrated uncertainties)

- **N_analysis = 248**, N_gold = 106
- Baseline fit did not converge!

## Species Breakdown

| Species | N | delta_alpha | sigma (calibrated) |
|---------|---|-------------|--------------------|
| Fe | 185 | 7.0908e-05 | 1.4887e-01 |
| Ni | 63 | 4.7316e-06 | 1.6783e-01 |

## Stress Test Summary

- Max drift (both, full): 8.6505e-04
- Total configurations tested: 43

## Most Influential Lines (LOO)

## Gates

| Gate | Status | Note |
|------|--------|------|
| gate_n | **PASS** | N_analysis=248, target>=100 |
| gate_chi2 | **FAIL** | chi2/dof=inf (calibrated sigma) |
| gate_stability | **PASS** | drift=8.65e-04, 3*sigma=1.50e+00 |
| gate_species | **PASS** | Fe=7.09e-05+/-1.49e-01, Ni=4.73e-06+/-1.68e-01, tension=0.0sigma |

## Artifacts

- line_measurements_all.csv
- line_measurements_analysis.csv
- line_measurements_gold.csv
- summary_delta_alpha.csv
- injection_recovery.csv
- loo_analysis.csv
- audit_attrition.csv
- audit_exclusion_reasons.csv
- identifiability_report.md
- diagnostic_plots/
