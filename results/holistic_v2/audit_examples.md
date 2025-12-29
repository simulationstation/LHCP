# Attrition Audit Report

## Waterfall Summary

| Stage | Description | Count | Dropped | Reason |
|-------|-------------|-------|---------|--------|
| S0 | Total atomic lines | 351 | 0 |  |
| S1 | In E140H range | 349 | 2 | OUT_OF_RANGE |
| S2 | Attempted fits | 349 | 0 | NOT_IN_FIT_TABLE |
| S3 | Fit converged | 275 | 74 | FIT_FAILED |
| S4 | Basic quality | 10 | 265 | QUALITY_CUTS |
| S5 | Gold (current) | 10 | 0 | ADDITIONAL_CUTS |

## Root Cause Analysis

### Why only 10 gold lines?

Starting from 351 atomic lines:

1. **Out of range**: 2 lines outside E140H coverage
2. **Fit failed**: 74 lines (21.2%) - no absorption detected in window
3. **Quality cuts (converged)**: 265 lines failed quality thresholds:
   - has_blend_suspect_flag: 275
   - has_poor_fit_flag: 205
   - reduced_chi2_exceeded: 187
   - asymmetry_exceeded: 187
   - sigma_lambda_obs_exceeded: 100

### Key Findings

- The **dominant loss** is from strict quality cuts on converged fits
- 275 lines have BLEND_SUSPECT flag (but this is expected in crowded WD spectra)
- 187 lines exceed asymmetry threshold (may be real blends or fitting artifacts)
- 187 lines have high reduced chi2 (poor local fit quality)

### Recommended Actions

1. **Relax guillotine cuts** - use robust weighting instead of hard exclusion
2. **Keep more lines** - use Student-t likelihood to downweight outliers
3. **Threshold adjustments**: sigma_max=0.02, chi2_max=10.0 would yield ~50+ lines
