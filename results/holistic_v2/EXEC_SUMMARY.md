# EXEC_SUMMARY - Holistic V2 Pipeline

## Baseline Result (Student-t, Both species, Full wavelength)

- **Δα/α** = -4.5640e-04 ± 3.1392e-05
- z0 = 5.0791e-05 ± 2.3023e-06
- jitter = 3.0000e-07
- χ² = 30509.81, dof = 65
- N_analysis = 69, N_gold = 0
- Condition number = 4.0
- Orthogonalized = False
- Likelihood = student_t

## Gates

| Gate | Status | Note |
|------|--------|------|
| N_analysis ≥ 50 | **PASS** | N_analysis=69 >= 50 |
| χ²/dof < 3 | **FAIL** | reduced_chi2=469.38 >= 3.0 |
| Stability | **FAIL** | drift=7.40e-03 > 3*sigma=9.42e-05 |
| Species | **FAIL** | Fe=-8.45e-04, Ni=1.85e-04, tension |

## Injection Recovery Summary

- Inject 0.0e+00: bias = 1.31e-04 ± 4.63e-04, mean_pull = 3.86
- Inject 1.0e-04: bias = -1.97e-04 ± 3.00e-04, mean_pull = -6.05
- Inject 3.0e-04: bias = -3.91e-05 ± 3.55e-04, mean_pull = -1.21

## Stress Test Summary (first 10 rows)

| Dist | Likelihood | Species | Split | N | Δα/α | σ | Cond |
|------|------------|---------|-------|---|------|---|------|
| 0 | student_t | Fe | full | 53 | -1.26e-03 | 4.04e-05 | 1.0 |
| 0 | student_t | Fe | blue | 26 | -1.44e-03 | 6.36e-05 | 1.0 |
| 0 | student_t | Fe | red | 27 | -1.43e-04 | 5.41e-05 | 1.0 |
| 0 | student_t | Ni | full | 16 | -7.94e-05 | 6.60e-05 | 1.0 |
| 0 | student_t | Ni | blue | 8 | -1.09e-03 | 3.19e-04 | 1.0 |
| 0 | student_t | Ni | red | 8 | -1.66e-04 | 6.91e-05 | 1.0 |
| 0 | student_t | both | full | 69 | -6.07e-04 | 3.12e-05 | 1.0 |
| 0 | student_t | both | blue | 34 | -1.96e-04 | 4.49e-05 | 1.0 |
| 0 | student_t | both | red | 35 | -1.72e-04 | 4.61e-05 | 1.0 |
| 0 | gaussian | Fe | full | 53 | -8.55e-04 | 4.27e-04 | 1.0 |

## Artifacts

- line_measurements_all.csv
- line_measurements_analysis.csv
- line_measurements_gold.csv
- summary_delta_alpha.csv
- inferred_alpha.json
- injection_recovery.csv
- diagnostic_plots/
