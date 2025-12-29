# Δα/α Results Summary Table

| Run | Grating | Species | N_lines | Δα/α | σ(Δα/α) | χ²/dof | Notes |
|-----|---------|---------|---------|------|---------|--------|-------|
| baseline_fe | E140H | Fe V | 268 | -2.34e-04 | 1.06e-05 | 3937 | baseline |
| baseline_ni | E140H | Ni V | 81 | +2.35e-03 | 1.62e-05 | 8003 | baseline |
| baseline_both | E140H | Fe V + Ni V | 348 | -2.96e-04 | 8.84e-06 | 4890 | baseline |
| fe_cont1 | E140H | Fe V | 268 | -2.35e-04 | 1.06e-05 | 4096 | spline_s=0.0001 |
| fe_cont2 | E140H | Fe V | 268 | -2.34e-04 | 1.06e-05 | 3937 | spline_s=0.001 |
| fe_cont3 | E140H | Fe V | 266 | -2.62e-04 | 1.08e-05 | 3989 | spline_s=0.01 |
| both_cont1 | E140H | Fe V + Ni V | 348 | -2.54e-04 | 8.77e-06 | 4939 | spline_s=0.0001 |
| both_cont2 | E140H | Fe V + Ni V | 348 | -2.96e-04 | 8.84e-06 | 4890 | spline_s=0.001 |
| both_cont3 | E140H | Fe V + Ni V | 349 | -1.87e-04 | 8.80e-06 | 4974 | spline_s=0.01 |
| fe_clip | E140H | Fe V | 268 | -2.34e-04 | 1.06e-05 | 3937 | sigma_clip=3 |
| both_clip | E140H | Fe V + Ni V | 348 | -2.96e-04 | 8.84e-06 | 4890 | sigma_clip=3 |

## Key Findings

1. **Fe V baseline:** Δα/α = (-2.34 ± 1.06) × 10⁻⁴
2. **Ni V baseline:** Δα/α = (+2.35 ± 0.16) × 10⁻³ ⚠️ TENSION
3. **Combined baseline:** Δα/α = (-2.96 ± 0.88) × 10⁻⁴
4. **Continuum systematic:** ~1 × 10⁻⁴
5. **Sigma clipping:** No outliers removed at 3σ
