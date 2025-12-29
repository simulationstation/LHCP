# LHCP Data Collection for Δα/α Measurements in White Dwarf Spectra

## Overview

This dataset supports measurement of fine-structure constant (α) variation
in strong gravitational fields using Fe V and Ni V absorption lines in the
hot DA white dwarf G191-B2B.

**Collection Date:** 2025-12-28
**Target:** G191-B2B (WD 0501+527)

---

## Data Status Summary

| Category | Status | Files | Description |
|----------|--------|-------|-------------|
| HLSP Spectra | Complete | 4 | WD-LINELIST coadds + CALSPEC |
| Raw STIS x1d | Sample | 3 | E140H sample for validation |
| NIST Atomic Data | Complete | 1126 lines | Fe V + Ni V wavelengths |
| Fe V Q-Coefficients | Complete | 346 lines | From Hu et al. 2021 |
| Ni V Q-Coefficients | Complete | 184 lines | From Webb et al. 2025 |
| Combined Line Table | Complete | 578 w/q | Merged NIST + q-coefficients |
| Stellar Parameters | Complete | 1 file | Mass, radius, Teff, log g |

---

## Directory Structure

```
LHCP/
├── data/
│   ├── hlsp/
│   │   ├── calspec/
│   │   │   ├── g191b2b_mod_011.fits      # Model spectrum
│   │   │   └── g191b2b_fos_003.fits      # FOS spectrum
│   │   └── wd-linelist/
│   │       ├── hlsp_wd-linelist_..._cspec.fits  # Coadded STIS E140H+E230H
│   │       └── hlsp_wd-linelist_..._sed.fits    # Multi-instrument SED
│   ├── raw/
│   │   └── mastDownload/HST/             # Sample raw x1d files (3)
│   ├── atomic/
│   │   ├── combined_lines_q.csv          # Combined line table with q-coefficients
│   │   ├── nist/
│   │   │   ├── nist_parsed_Fe_V_*.csv    # 462 Fe V lines
│   │   │   ├── nist_parsed_Ni_V_*.csv    # 664 Ni V lines
│   │   │   └── nist_combined_*.csv       # Combined NIST table
│   │   └── papers/
│   │       ├── FeV_q_coefficients_Hu2021.csv   # Fe V q-values (346 lines)
│   │       ├── NiV_q_coefficients_G191_*.csv   # Ni V q-values (4 continuum levels)
│   │       ├── source_2007.10905/        # Hu et al. 2021 LaTeX source
│   │       └── source_2410.01849/        # Webb et al. 2025 source
│   ├── mast_query_results.csv            # 65 STIS observation inventory
│   └── stellar_parameters_g191b2b.json   # G191-B2B stellar parameters
├── scripts/
│   ├── parse_fev_q_coefficients.py       # Parse Fe V q from LaTeX
│   ├── download_hlsp_direct.py           # HLSP download script
│   ├── build_combined_table.py           # Build combined line table
│   ├── update_manifests.py               # Update manifest files
│   ├── fetch_nist_data.py                # NIST retrieval
│   └── query_mast_g191b2b.py             # MAST observation query
├── manifests/
│   ├── manifest_files.csv                # File inventory (42 files)
│   ├── manifest_obs.csv                  # Observation inventory (65 obs)
│   └── manifest_lines.csv                # Line inventory (1126 lines)
└── README_DATA.md
```

---

## Key Data Files

### 1. Combined Atomic Line Table

**File:** `data/atomic/combined_lines_q.csv`

This is the primary analysis file containing:
- 1126 total lines (462 Fe V + 664 Ni V)
- 578 lines with q-coefficients for Δα/α analysis
- Columns: species, wavelength, uncertainty, energy levels, Aki, fik, q-coefficient, sources

| Species | Total Lines | With Q | Wavelength Range (Å) | Q Range (cm⁻¹) |
|---------|------------|--------|---------------------|----------------|
| Fe V | 462 | 347 (75%) | 1100.14 - 1845.04 | 545 to 4562 |
| Ni V | 664 | 231 (35%) | 1105.58 - 1397.55 | 612 to 5009 |

### 2. HLSP Spectra

**WD-LINELIST Products (Preval et al.):**
- `hlsp_wd-linelist_hst_stis_g191-b2b_e140h-e230h_v1_cspec.fits` - Coadded STIS
- `hlsp_wd-linelist_multi_multi_g191-b2b_multi_v1_sed.fits` - Multi-instrument SED

**CALSPEC Products:**
- `g191b2b_mod_011.fits` - Model spectrum (3.9 MB)
- `g191b2b_fos_003.fits` - FOS spectrum (98 KB)

### 3. Stellar Parameters

**File:** `data/stellar_parameters_g191b2b.json`

| Parameter | Value | Uncertainty | Source |
|-----------|-------|-------------|--------|
| Mass | 0.54 M☉ | ±0.03 | Rauch+ 2013 |
| Radius | 0.0204 R☉ | ±0.001 | Preval+ 2013 |
| Teff | 60,000 K | ±2,000 | Rauch+ 2013 |
| log g | 7.58 (cgs) | ±0.05 | Rauch+ 2013 |
| z_grav | 5.04 × 10⁻⁵ | ±4 × 10⁻⁶ | Computed |

---

## Q-Coefficient Sources

### Fe V: Hu et al. (2021)
- **Paper:** MNRAS 500, 1466-1475
- **arXiv:** 2007.10905
- **Method:** New calculations superseding Ong et al. (2013)
- **File:** `FeV_q_coefficients_Hu2021.csv`

### Ni V: Webb et al. (2025)
- **Paper:** Open Journal of Astrophysics, Vol 8
- **arXiv:** 2410.01849
- **Method:** Multiple continuum placement scenarios (0.5, 0.75, 1.0, 1.25)
- **Files:** `NiV_q_coefficients_G191_*.csv` (4 files)

---

## MAST Observations

**Total:** 65 HST/STIS observations
**Gratings:** E140H (25 obs) + E230H (40 obs)
**Date Range:** MJD 51164 - 52172 (Dec 1998 - Sep 2001)
**Total Exposure:** ~90,000 seconds

Full observation list in `manifests/manifest_obs.csv`.

---

## Quick Start for Δα/α Analysis

1. **Load combined line table:**
   ```python
   import pandas as pd
   lines = pd.read_csv('data/atomic/combined_lines_q.csv')

   # Filter to lines with q-coefficients
   analysis_lines = lines[lines['q_cm-1'].notna()]
   print(f"Lines available for analysis: {len(analysis_lines)}")
   ```

2. **Load coadded spectrum:**
   ```python
   from astropy.io import fits
   hdul = fits.open('data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h-e230h_v1_cspec.fits')
   wavelength = hdul[1].data['WAVELENGTH']
   flux = hdul[1].data['FLUX']
   ```

3. **Apply gravitational redshift:**
   ```python
   import json
   with open('data/stellar_parameters_g191b2b.json') as f:
       params = json.load(f)
   z_grav = params['gravitational_redshift']['z_grav']
   ```

---

## Downloading Additional Data

### Full Raw STIS Dataset

To download all 65 observations:

```python
from astroquery.mast import Observations
import pandas as pd

obs_df = pd.read_csv('manifests/manifest_obs.csv')
for obsid in obs_df['obsid']:
    products = Observations.get_product_list(str(obsid))
    x1d = Observations.filter_products(products, productSubGroupDescription=['X1D'])
    Observations.download_products(x1d, download_dir='data/raw/')
```

---

## References

### Primary Publications

1. **Hu et al. (2021)** - MNRAS 500, 1466
   "Measuring the fine structure constant on a white dwarf surface; a detailed analysis of Fe V absorption in G191-B2B"
   arXiv:2007.10905

2. **Webb et al. (2025)** - Open J. Astrophys. 8
   "Searching for new physics using high precision absorption spectroscopy: Ni V in G191-B2B"
   arXiv:2410.01849

3. **Preval et al. (2013)** - MNRAS 436, 659
   "A comprehensive near- and far-ultraviolet spectroscopic study of G191-B2B"

### Data Sources

- NIST ASD v5.12: https://physics.nist.gov/asd
- MAST HLSP: https://archive.stsci.edu/hlsp/
- WD-LINELIST: https://archive.stsci.edu/prepds/wd-linelist/

---

## Changelog

- **2025-12-28:** Dataset complete
  - Downloaded HLSP products (WD-LINELIST + CALSPEC)
  - Parsed Fe V q-coefficients from Hu et al. LaTeX tables
  - Built combined atomic table with 578 lines having q-coefficients
  - Updated all manifests
  - Sample raw STIS x1d files downloaded for validation

---

*Generated by LHCP data collection pipeline*
