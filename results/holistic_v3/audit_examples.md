# V3 Attrition Audit Report

## Waterfall Summary

| Stage | Description | Count | Dropped | Reason |
|-------|-------------|-------|---------|--------|
| S0 | Total atomic lines | 351 | 0 |  |
| S1 | Attempted fits | 351 | 0 | OUT_OF_RANGE_OR_EDGE |
| S2 | Fit converged | 336 | 15 | FIT_FAILED |
| S3 | Analysis set | 248 | 88 | QUALITY_CUTS |
| S4 | Gold set | 106 | 142 | STRICT_CUTS |

## Key Improvements in V3

- **Blend-aware fitting**: 85 lines with blends successfully modeled
- **Analysis set size**: 248 lines (target >= 50)
- **Relaxed thresholds**: sigma_max=0.1, chi2_max=50.0

## Top Exclusion Reasons (from converged fits)

- **blend_resolved_kept**: 125 (37.2%)
- **blend_modeled_kept**: 85 (25.3%)
- **low_snr**: 42 (12.5%)
- **chi2_exceeded**: 39 (11.6%)
- **sigma_exceeded**: 16 (4.8%)
