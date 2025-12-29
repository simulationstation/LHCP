# HLSP Sanity Report
Generated: 2025-12-29T09:34:57Z

## Summary

Both E140H and E230H coadd FITS files pass all sanity checks.

## File Details

### E140H Coadd

| Property | Value |
|----------|-------|
| Filename | `hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_coadd-spec.fits` |
| Size | 2,770,560 bytes |
| HDU with data | 1 (SPECDATA) |
| Columns | WAVE, FLUX, ERROR, TIM_EL, SIGNOI |
| N points | 69,079 |
| Wavelength range | 1160.00 - 1684.99 Å |
| Wavelength convention | **VAC** (vacuum) |
| Median flux | 8.955e-12 erg/s/cm²/Å |
| Median S/N | 71.9 |
| Min/Max flux | 7.331e-14 / 1.976e-11 |
| NaN in flux | 0 |
| NaN in error | 0 |
| Zero errors | 0 |

### E230H Coadd

| Property | Value |
|----------|-------|
| Filename | `hlsp_wd-linelist_hst_stis_g191-b2b_e230h_v1_coadd-spec.fits` |
| Size | 4,294,080 bytes |
| HDU with data | 1 (SPECDATA) |
| Columns | WAVE, FLUX, ERROR, TIM_EL, SIGNOI |
| N points | 107,143 |
| Wavelength range | 1650.00 - 3149.99 Å |
| Wavelength convention | **VAC** (vacuum) |
| Median flux | 1.436e-12 erg/s/cm²/Å |
| Median S/N | 42.7 |
| Min/Max flux | 9.262e-15 / 6.683e-12 |
| NaN in flux | 0 |
| NaN in error | 0 |
| Zero errors | 0 |

## Wavelength Convention Verification

Both files have `AIRORVAC = 'VAC'` in the primary header, confirming VACUUM wavelengths.
This is consistent with the WD-LINELIST documentation (Preval et al. 2013).

## Checks Passed

- ✓ WAVE, FLUX, ERROR columns present
- ✓ AIRORVAC = VAC (vacuum wavelengths)
- ✓ No NaN values in flux or error
- ✓ No zero errors
- ✓ Wavelength coverage is plausible for STIS E140H/E230H
- ✓ S/N > 40 (high quality data)

## Data Ready for Analysis

Both coadds are suitable for line fitting and Δα/α inference.
