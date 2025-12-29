# LHCP Data Collection for Δα/α Measurements in White Dwarf Spectra

## Overview

This dataset collection supports measurement of possible fine-structure constant (α) variation
in strong gravitational fields using Fe V and Ni V absorption lines in the hot DA white dwarf
G191-B2B.

**Collection Date:** 2025-12-28
**Target:** G191-B2B (WD 0501+527)

---

## Summary of Collected Data

### ✅ Successfully Retrieved

| Category | Description | Status |
|----------|-------------|--------|
| **NIST Atomic Data** | Fe V and Ni V laboratory wavelengths (1100-1900 Å) | Complete |
| **Q-Coefficients** | Ni V sensitivity coefficients from Webb et al. 2025 | Complete |
| **Q-Coefficients** | Fe V sensitivity coefficients from Hu et al. 2021 (in LaTeX table) | Complete |
| **Stellar Parameters** | Mass, radius, log g, Teff for G191-B2B | Complete |
| **MAST Observation List** | 65 HST/STIS observations (25 E140H, 40 E230H) | Complete |
| **Source Documentation** | arXiv sources with full atomic data tables | Complete |

### ⚠️ Requires Additional Steps

| Category | Description | Action Required |
|----------|-------------|-----------------|
| **HLSP Spectra** | Preval et al. coadded STIS/FUSE spectra | Use MAST portal for download |
| **Raw STIS Data** | Individual x1d/wavecal products | Use astroquery download script |
| **Fe V q-table** | Machine-readable extraction from LaTeX | Parse from source_2007.10905/table/summary_atomic_data.tex |

---

## Directory Structure

```
LHCP/
├── data/
│   ├── raw/                    # (empty) For MAST observation downloads
│   ├── hlsp/                   # (empty) For HLSP coadded spectra
│   ├── atomic/
│   │   ├── nist/               # NIST ASD laboratory wavelengths
│   │   │   ├── nist_raw_Fe_V_*.csv         # Raw NIST download
│   │   │   ├── nist_raw_Ni_V_*.csv         # Raw NIST download
│   │   │   ├── nist_parsed_Fe_V_*.csv      # Parsed Fe V lines
│   │   │   ├── nist_parsed_Ni_V_*.csv      # Parsed Ni V lines
│   │   │   ├── nist_combined_FeV_NiV_*.csv # Combined table
│   │   │   └── nist_metadata_*.txt         # Retrieval metadata
│   │   └── papers/             # Literature atomic data
│   │       ├── q_coefficient_sources.md    # Source documentation
│   │       ├── NiV_q_coefficients_G191_*.csv  # Ni V q values
│   │       ├── arXiv_2007.10905_source.tar.gz # Hu et al. 2021 (Fe V)
│   │       ├── arXiv_2410.01849_source.tar.gz # Webb et al. 2025 (Ni V)
│   │       ├── source_2007.10905/          # Extracted Fe V tables
│   │       └── source_2410.01849/          # Extracted Ni V tables
│   ├── mast_query_results.csv  # STIS observation inventory
│   └── stellar_parameters_g191b2b.json
├── scripts/
│   ├── query_mast_g191b2b.py   # MAST query script
│   ├── fetch_nist_data.py      # NIST retrieval script
│   ├── parse_fev_q_coefficients.py  # Fe V q extraction
│   └── utils.py                # Wavelength conversion utilities
├── manifests/
│   ├── manifest_files.csv      # File inventory
│   ├── manifest_obs.csv        # Observation inventory
│   └── manifest_lines.csv      # Combined line table
└── README_DATA.md              # This file
```

---

## Data Sources and Provenance

### 1. NIST Atomic Spectra Database (ASD)

**Version:** 5.12
**Access Date:** 2025-12-28
**URL:** https://physics.nist.gov/asd

| Species | Lines | Wavelength Range (Å) | File |
|---------|-------|---------------------|------|
| Fe V | 462 | 1100.14 - 1845.04 | nist_parsed_Fe_V_*.csv |
| Ni V | 664 | 1105.58 - 1397.55 | nist_parsed_Ni_V_*.csv |

**Notes:**
- Wavelengths are in VACUUM (NIST uses vacuum below 2000 Å)
- Uncertainties included where available
- Oscillator strengths (f values) and transition probabilities (Aki) included

### 2. Q-Coefficient Tables

#### Fe V (Primary: Hu et al. 2021)
- **Paper:** MNRAS 500, 1466-1475 (2021)
- **arXiv:** 2007.10905
- **Location:** `data/atomic/papers/source_2007.10905/table/summary_atomic_data.tex`
- **Format:** LaTeX longtable (machine-readable parsing required)
- **Notes:** These supersede Ong et al. (2013) calculations

#### Ni V (Primary: Webb et al. 2025)
- **Paper:** Open Journal of Astrophysics, Vol 8 (2025)
- **arXiv:** 2410.01849
- **Location:** `data/atomic/papers/NiV_q_coefficients_G191_*.csv`
- **Format:** Tab-separated text converted to CSV
- **Notes:** Multiple continuum placement scenarios (0.5, 0.75, 1.0, 1.25)

### 3. Stellar Parameters (G191-B2B)

| Parameter | Adopted Value | Uncertainty | Source |
|-----------|---------------|-------------|--------|
| Mass | 0.54 M☉ | ±0.03 | Rauch+ 2013 / Preval+ 2013 |
| Radius | 0.0204 R☉ | ±0.001 | Preval+ 2013 |
| Teff | 60,000 K | ±2,000 | Rauch+ 2013 |
| log g | 7.58 (cgs) | ±0.05 | Rauch+ 2013 / Gianninas+ 2011 |
| Parallax | 19.05 mas | ±0.06 | Gaia DR3 |
| Distance | 52.49 pc | ±0.15 | From Gaia |

**Computed Gravitational Redshift:**
- z_grav = 5.04 × 10⁻⁵ ± 4 × 10⁻⁶
- v_grav = 15.1 ± 1.1 km/s
- Observed: 19 ± 3 km/s (Reid & Wegner 1988)

### 4. HST/STIS Observations

**Total Observations:** 65 (E140H: 25, E230H: 40)
**Date Range:** 1998-12-17 to 2001-09-19
**HST Proposals:** 8067, 8421, 8915

| Grating | Wavelength (nm) | Total Exposure (s) | Products |
|---------|-----------------|-------------------|----------|
| E140H | 114-170 | 40,230 | x1d, wav, raw |
| E230H | 162-315 | 50,118 | x1d, wav, raw |

All data is PUBLIC (no proprietary restrictions).

---

## Missing Data and Alternatives

### 1. HLSP Coadded Spectra (Not Downloaded)

The Preval et al. HLSP products were not automatically downloaded due to URL changes at MAST.

**To Download:**
```bash
# Option 1: MAST Portal
# Visit: https://archive.stsci.edu/prepds/wd-linelist/
# Navigate to G191-B2B products and download manually

# Option 2: astroquery
python3 << 'EOF'
from astroquery.mast import Observations
obs = Observations.query_criteria(provenance_name="wd-linelist", target_name="G191B2B")
products = Observations.get_product_list(obs)
Observations.download_products(products)
EOF
```

### 2. Raw STIS Observations (Not Downloaded)

The x1d, wavecal, and raw products are available from MAST.

**To Download:**
```python
from astroquery.mast import Observations
import pandas as pd

# Load observation list
obs_df = pd.read_csv('data/mast_query_results.csv')

# Download all products for E140H observations
for obsid in obs_df[obs_df['filters'].str.contains('E140H')]['obsid']:
    products = Observations.get_product_list(obsid)
    Observations.download_products(products,
                                   productSubGroupDescription=['X1D', 'WAV'],
                                   download_dir='data/raw/')
```

### 3. Fe V Q-Coefficients (In LaTeX Table)

The Fe V q-coefficients are embedded in a LaTeX longtable format. To extract:

```bash
# The table is in: data/atomic/papers/source_2007.10905/table/summary_atomic_data.tex
# Column layout: ID, config_lower, term_lower, J_lower, config_upper, term_upper, J_upper,
#                E_lower, E_upper, wavelength_K14, flag, unc_K14, wavelength_W19, flag, unc_W19,
#                f, Gamma, q
# Last column (q) contains q-coefficients in cm^-1
```

---

## Minimal Viable Analysis Set

For a first Δα/α measurement run, you need:

### Required Files

1. **Spectrum:** Download E140H coadded spectrum from MAST HLSP
2. **Fe V Line Table:** `data/atomic/nist/nist_parsed_Fe_V_*.csv` (462 lines)
3. **Fe V Q-Coefficients:** Extract from `source_2007.10905/table/summary_atomic_data.tex`
4. **Stellar Parameters:** `data/stellar_parameters_g191b2b.json`

### Analysis Steps

1. **Line Identification:** Match observed absorption features to Fe V lab wavelengths
2. **Wavelength Measurement:** Fit Gaussian profiles to isolated Fe V lines
3. **Velocity Calculation:** v = c × (λ_obs - λ_lab) / λ_lab
4. **Δα/α Extraction:**
   - Δv = 2c × (Δα/α) × q × λ / E
   - Or use full VPFIT modeling with q-coefficients

### Expected Precision

- Individual line: ~1 km/s (STIS resolution)
- Combined (N~100 lines): ~0.1 km/s
- Δα/α sensitivity: ~10⁻⁵ (for q ~ 1000-4000 cm⁻¹)

---

## References

### Primary Publications

1. **Hu et al. (2021)** - MNRAS 500, 1466
   "Measuring the fine structure constant on a white dwarf surface"
   arXiv:2007.10905

2. **Webb et al. (2025)** - Open J. Astrophys. 8
   "Searching for new physics using high precision absorption spectroscopy"
   arXiv:2410.01849

3. **Preval et al. (2013)** - MNRAS 436, 659
   "A comprehensive near- and far-ultraviolet spectroscopic study of G191-B2B"

4. **Berengut et al. (2013)** - Phys. Rev. Lett. 111, 010801
   arXiv:1305.1337

### Data Sources

- NIST ASD: https://physics.nist.gov/asd
- MAST HLSP: https://archive.stsci.edu/prepds/wd-linelist/
- Gaia DR3: https://gea.esac.esa.int/archive/

### Contact for Theoretical Data

- Vladimir Dzuba: v.dzuba@unsw.edu.au (q-coefficients)
- John Webb: j.webb@unsw.edu.au (observational analysis)

---

## Changelog

- **2025-12-28:** Initial data collection
  - Retrieved 1126 Fe V + Ni V lines from NIST
  - Downloaded arXiv sources with q-coefficient tables
  - Compiled stellar parameters from literature
  - Generated MAST observation inventory

---

*Generated by LHCP data collection pipeline*
