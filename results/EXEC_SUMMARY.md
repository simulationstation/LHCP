# Executive Summary: Δα/α Analysis of G191-B2B

**Generated:** 2025-12-29
**Target:** G191-B2B (WD 0501+527)
**Analysis Pipeline:** wdalpha v0.1.0

---

## Data Files Used

| File | Path | Description |
|------|------|-------------|
| E140H Spectrum | `data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_coadd-spec.fits` | 69,079 points, 1160-1685 Å, VACUUM |
| Atomic Lines | `data/atomic/analysis_lines_clean.csv` | 351 clean lines with q-coefficients |

---

## Baseline Results (E140H)

| Species | Δα/α | σ(Δα/α) | N_lines | χ²/dof |
|---------|------|---------|---------|--------|
| **Fe V only** | -2.34 × 10⁻⁴ | 1.06 × 10⁻⁵ | 268 | 3937 |
| **Ni V only** | +2.35 × 10⁻³ | 1.62 × 10⁻⁵ | 81 | 8003 |
| **Combined** | -2.96 × 10⁻⁴ | 8.84 × 10⁻⁶ | 348 | 4890 |

### Critical Observation: Fe/Ni Tension

**Fe V and Ni V yield inconsistent results:**
- Fe V: Δα/α = -2.34 × 10⁻⁴ (α smaller on WD)
- Ni V: Δα/α = +2.35 × 10⁻³ (α larger on WD)

The disagreement is **~10× in magnitude** and **opposite in sign**. This is a major red flag indicating:
1. Systematic errors in q-coefficients
2. Wavelength calibration issues
3. Line blending not properly handled
4. Different systematic effects for different species

---

## Stress Test Results

### Continuum Sensitivity (Fe V)

| spline_s | Δα/α | Change from baseline |
|----------|------|---------------------|
| 0.0001 | -2.35 × 10⁻⁴ | +0.01 × 10⁻⁴ |
| 0.001 (baseline) | -2.34 × 10⁻⁴ | — |
| 0.01 | -2.62 × 10⁻⁴ | -0.28 × 10⁻⁴ |

### Continuum Sensitivity (Combined)

| spline_s | Δα/α | Change from baseline |
|----------|------|---------------------|
| 0.0001 | -2.54 × 10⁻⁴ | +0.42 × 10⁻⁴ |
| 0.001 (baseline) | -2.96 × 10⁻⁴ | — |
| 0.01 | -1.87 × 10⁻⁴ | +1.09 × 10⁻⁴ |

**Systematic uncertainty from continuum:** ~1 × 10⁻⁴ (comparable to statistical error)

### Sigma Clipping (3σ)

No outliers removed - all lines within 3σ of median residual.

---

## Quality Assessment

### χ²/dof Analysis

All fits show **extremely high χ²/dof** (3900-8000 instead of ~1):

| Possible Cause | Evidence |
|----------------|----------|
| Underestimated errors | Likely - formal errors don't include systematics |
| Model inadequacy | Likely - single Gaussian may not capture line shapes |
| Wavelength distortions | Possible - STIS may have calibration drifts |
| Line blending | Likely - many lines overlap in crowded regions |

### Recommendations for Refinement

1. **Inflate errors** by √(χ²/dof) to get realistic uncertainties
2. **Implement robust regression** (Huber loss) to downweight outliers
3. **Add wavelength distortion term** (low-order polynomial vs λ)
4. **Investigate Fe/Ni tension** - check q-coefficient sources
5. **Use line-by-line cross-validation** to identify problematic lines

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| Total runs completed | 11 |
| E140H lines fit (combined) | 348/349 |
| Fe V lines in E140H | 268 |
| Ni V lines in E140H | 81 |
| Max Δα/α deviation (continuum) | 1.09 × 10⁻⁴ |
| Baseline result (combined) | (-2.96 ± 0.88) × 10⁻⁴ |

---

## Conclusion

The pipeline successfully processes real HST/STIS data and produces Δα/α measurements. However, the **Fe/Ni tension** and **high χ²/dof** indicate that the current analysis is **not publication-ready**. The infrastructure is in place for:

1. More sophisticated continuum modeling
2. Robust statistical methods
3. Systematic error quantification
4. Multi-grating analysis

**Next steps:** Address Fe/Ni inconsistency, implement robust regression, add wavelength distortion terms.

---

*Analysis performed with wdalpha pipeline on WD-LINELIST HLSP data (Preval et al. 2013)*
