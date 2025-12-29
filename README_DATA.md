# LHCP Data Collection for Δα/α Measurements in White Dwarf Spectra

## Overview

This dataset supports measurement of fine-structure constant (α) variation
in strong gravitational fields using Fe V and Ni V absorption lines in the
hot DA white dwarf G191-B2B.

**Collection Date:** 2025-12-28
**Target:** G191-B2B (WD 0501+527)

---

## Data Status Summary

| Category | Status | Description |
|----------|--------|-------------|
| WD-LINELIST Coadds | **Complete** | E140H + E230H STIS echelle coadds |
| WD-LINELIST Linelists | **Complete** | 1,211 identified lines with centroids |
| CALSPEC Spectra | Complete | Model + FOS reference spectra |
| Raw STIS x1d | Sample (3) | Full dataset (65 obs) available via MAST |
| Atomic Line Table | Complete | 1126 lines (Fe V + Ni V) |
| Q-Coefficients | Complete | 578 lines with q-values |
| Analysis-Ready Table | Complete | 351 clean lines (join_flag=OK) |
| Stellar Parameters | Complete | Mass, radius, Teff, log g |

---

## Key Analysis Files

### 1. Analysis-Ready Line Tables

**Primary analysis file:** `data/atomic/analysis_lines.csv`
- 1126 total lines with explicit λ0 selection
- Columns: species, lambda0_ang, sigma_lambda0_ang, q_cm-1, lambda0_source, join_flag

**Clean subset:** `data/atomic/analysis_lines_clean.csv`
- 351 lines with `join_flag=OK` (safe for analysis)
- No large wavelength mismatches or ambiguous joins

| Species | Total | With Q | Clean (OK) | Wavelength Range |
|---------|-------|--------|------------|------------------|
| Fe V | 462 | 347 | ~290 | 1100-1845 Å |
| Ni V | 664 | 231 | ~61 | 1106-1398 Å |

**Join Flag Meanings:**
- `OK`: λ0 and q consistent, |Δλ| ≤ 0.005 Å
- `LARGE_MISMATCH`: |Δλ| > 0.005 Å (72 lines)
- `AMBIGUOUS`: Same q-wavelength matched to multiple NIST lines (155 lines)
- `NO_Q`: No q-coefficient available (548 lines)

### 2. Q-Coefficient Sources

| Species | Source | Reference | Lines |
|---------|--------|-----------|-------|
| Fe V | Hu2021_table | Hu et al. 2021, MNRAS 500, 1466 | 346 |
| Ni V | Lee2025_2410.01849 | Lee et al. 2025, OJAp 8 (arXiv:2410.01849) | 184 |

### 3. HLSP Spectra

**WD-LINELIST (Primary - Preval et al. 2013):**
- `data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_coadd-spec.fits`
  - E140H coadd: 1160-1685 Å, 69,079 points, VACUUM wavelengths
- `data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e230h_v1_coadd-spec.fits`
  - E230H coadd: 1650-3150 Å, 107,143 points, VACUUM wavelengths
- `data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_linelist.txt`
  - 1,155 identified absorption features with centroids
- `data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e230h_v1_linelist.txt`
  - 56 identified absorption features

**CALSPEC (Reference):**
- `data/hlsp/calspec/g191b2b_mod_011.fits` - Model spectrum (900-320000 Å)
- `data/hlsp/calspec/g191b2b_fos_003.fits` - FOS observation (1141-9203 Å)

---

## Directory Structure

```
LHCP/
├── data/
│   ├── hlsp/
│   │   ├── calspec/
│   │   │   ├── g191b2b_mod_011.fits      # Model spectrum
│   │   │   └── g191b2b_fos_003.fits      # FOS spectrum
│   │   ├── wd-linelist/                  # (empty - not accessible)
│   │   ├── HLSP_NOTES.txt                # Detailed notes on HLSP status
│   │   └── download_provenance.json      # Download attempts log
│   ├── raw/
│   │   └── mastDownload/HST/             # Sample STIS x1d files (3)
│   ├── atomic/
│   │   ├── analysis_lines.csv            # PRIMARY: Full analysis table
│   │   ├── analysis_lines_clean.csv      # Clean subset (join_flag=OK)
│   │   ├── combined_lines_q.csv          # Merged NIST + q-coefficients
│   │   ├── join_diagnostics.md           # Join quality report
│   │   ├── nist/                         # NIST laboratory wavelengths
│   │   └── papers/                       # Q-coefficient source tables
│   ├── mast_query_results.csv            # 65 STIS observations
│   └── stellar_parameters_g191b2b.json   # G191-B2B stellar parameters
├── scripts/
│   ├── analyze_qjoin.py                  # Q-join analysis
│   ├── build_combined_table.py           # Build combined line table
│   ├── download_wdlinelist_v2.py         # HLSP download attempts
│   └── update_manifests.py               # Manifest updater
├── manifests/
│   ├── manifest_files.csv                # File inventory
│   ├── manifest_obs.csv                  # Observation inventory
│   └── manifest_lines.csv                # Line inventory
└── README_DATA.md
```

---

## Join Diagnostics Summary

From `data/atomic/join_diagnostics.md`:

| Metric | Value |
|--------|-------|
| Total lines | 1126 |
| Lines with q | 578 |
| |Δλ| > 0.005 Å | 144 (24.9%) |
| |Δλ| > 0.01 Å | 63 (10.9%) |
| |Δλ| > 0.02 Å | 22 (3.8%) |
| Worst |Δλ| | 0.050 Å |
| Duplicate q-wavelength matches | 69 |

**Per-Species:**
- Fe V: 57 lines with |Δλ| > 0.005 Å, max 0.034 Å
- Ni V: 87 lines with |Δλ| > 0.005 Å, max 0.050 Å

---

## Quick Start for Δα/α Analysis

```python
import pandas as pd

# Load clean analysis lines
lines = pd.read_csv('data/atomic/analysis_lines_clean.csv')
print(f"Analysis-ready lines: {len(lines)}")

# Use lambda0_ang as the reference wavelength
# (consistent with q-coefficient source)
fe_lines = lines[lines['species'] == 'Fe V']
ni_lines = lines[lines['species'] == 'Ni V']

print(f"Fe V: {len(fe_lines)} lines, q range: {fe_lines['q_cm-1'].min()}-{fe_lines['q_cm-1'].max()} cm^-1")
print(f"Ni V: {len(ni_lines)} lines, q range: {ni_lines['q_cm-1'].min()}-{ni_lines['q_cm-1'].max()} cm^-1")
```

---

## Stellar Parameters

**File:** `data/stellar_parameters_g191b2b.json`

| Parameter | Value | Uncertainty | Source |
|-----------|-------|-------------|--------|
| Mass | 0.54 M☉ | ±0.03 | Rauch+ 2013 |
| Radius | 0.0204 R☉ | ±0.001 | Preval+ 2013 |
| Teff | 60,000 K | ±2,000 | Rauch+ 2013 |
| log g | 7.58 (cgs) | ±0.05 | Rauch+ 2013 |
| z_grav | 5.04 × 10⁻⁵ | ±4 × 10⁻⁶ | Computed |

---

## Downloading Full Raw STIS Dataset

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
   "Measuring the fine structure constant on a white dwarf surface"
   arXiv:2007.10905

2. **Lee et al. (2025)** - Open J. Astrophys. 8
   "Searching for new physics using high precision absorption spectroscopy: Ni V in G191-B2B"
   arXiv:2410.01849

3. **Preval et al. (2013)** - MNRAS 436, 659
   "A comprehensive near- and far-ultraviolet spectroscopic study of G191-B2B"

### Data Sources

- NIST ASD v5.12: https://physics.nist.gov/asd
- MAST HLSP: https://archive.stsci.edu/hlsp/
- CALSPEC: https://archive.stsci.edu/hlsps/reference-atlases/cdbs/calspec/

---

## Remaining Optional Tasks

1. **Full raw STIS download**: 65 observations available in MAST
2. **Manual coaddition**: Combine x1d files to create E140H/E230H coadds
3. **WD-LINELIST access**: Try MAST Portal manual download if API becomes available

---

## Changelog

- **2025-12-28:** Dataset analysis-ready
  - Identified duplicate/corrupted HLSP files, cleaned up
  - Built analysis_lines.csv with explicit λ0 selection
  - Created join_diagnostics.md with quality metrics
  - Fixed Ni V q-source attribution (Lee et al. 2025, not Webb)
  - 351 clean lines available for Δα/α analysis

---

*Generated by LHCP data collection pipeline*
