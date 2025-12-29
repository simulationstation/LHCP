#!/usr/bin/env python3
"""
Fetch Fe V and Ni V spectral line data from NIST Atomic Spectra Database (ASD)
for white dwarf Delta alpha/alpha measurements in the UV range (1100-1900 Angstrom).

NIST ASD Version: 5.12
Data retrieved from: https://physics.nist.gov/cgi-bin/ASD/lines1.pl

Note: NIST uses VACUUM wavelengths below 2000 Angstrom.
"""

import os
import re
import csv
import urllib.request
import urllib.parse
from datetime import datetime
from pathlib import Path


# NIST ASD endpoint
NIST_ASD_URL = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"

# NIST ASD Version (as of late 2024/2025)
NIST_ASD_VERSION = "5.12"

# Wavelength range in nm (1100-1900 Angstrom = 110-190 nm)
WAVELENGTH_LOW_NM = 110.0
WAVELENGTH_HIGH_NM = 190.0

# Output directories
DATA_DIR = Path("/home/primary/LHCP/data/atomic/nist")
SCRIPT_DIR = Path("/home/primary/LHCP/scripts")


def build_nist_query_url(species: str, low_w: float, upp_w: float) -> str:
    """
    Build NIST ASD query URL for spectral lines.

    Parameters:
    -----------
    species : str
        Atomic species (e.g., 'Fe V', 'Ni V')
    low_w : float
        Lower wavelength bound in nm
    upp_w : float
        Upper wavelength bound in nm

    Returns:
    --------
    str : Complete URL for NIST ASD query
    """
    params = {
        'spectra': species,
        'output_type': '0',
        'low_w': str(low_w),
        'upp_w': str(upp_w),
        'unit': '1',          # nm
        'submit': 'Retrieve Data',
        'de': '0',
        'format': '2',        # CSV format
        'line_out': '0',      # All lines
        'en_unit': '0',       # cm^-1
        'output': '0',
        'bibrefs': '1',
        'page_size': '15',
        'show_obs_wl': '1',   # Observed wavelengths
        'show_calc_wl': '1',  # Ritz wavelengths
        'unc_out': '1',       # Uncertainties
        'order_out': '0',     # Order by wavelength
        'show_av': '2',       # Show vacuum wavelengths
        'A_out': '1',         # Transition probabilities
        'intens_out': 'on',   # Relative intensities
        'allowed_out': '1',   # Allowed transitions
        'forbid_out': '1',    # Forbidden transitions
        'conf_out': 'on',     # Configurations
        'term_out': 'on',     # Terms
        'enrg_out': 'on',     # Energy levels
        'J_out': 'on',        # J values
        'f_out': 'on',        # Oscillator strengths
    }

    query_string = urllib.parse.urlencode(params)
    return f"{NIST_ASD_URL}?{query_string}"


def fetch_nist_data(species: str, low_w: float, upp_w: float) -> str:
    """
    Fetch spectral line data from NIST ASD.

    Parameters:
    -----------
    species : str
        Atomic species (e.g., 'Fe V', 'Ni V')
    low_w : float
        Lower wavelength bound in nm
    upp_w : float
        Upper wavelength bound in nm

    Returns:
    --------
    str : Raw CSV data from NIST
    """
    url = build_nist_query_url(species, low_w, upp_w)
    print(f"Fetching data from NIST ASD for {species}...")
    print(f"URL: {url[:100]}...")

    headers = {
        'User-Agent': 'Mozilla/5.0 (Python scientific data retrieval)'
    }

    request = urllib.request.Request(url, headers=headers)

    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            data = response.read().decode('utf-8', errors='replace')
            return data
    except urllib.error.URLError as e:
        print(f"Error fetching data for {species}: {e}")
        raise


def clean_csv_value(val: str) -> str:
    """
    Clean NIST CSV value by removing Excel-style formatting.
    NIST uses ="value" format for CSV export.
    """
    if not val:
        return ""
    val = val.strip()
    # Remove ="..." wrapper
    if val.startswith('="') and val.endswith('"'):
        val = val[2:-1]
    elif val.startswith('"') and val.endswith('"'):
        val = val[1:-1]
    return val.strip()


def parse_nist_csv(raw_data: str, species: str) -> list:
    """
    Parse NIST ASD CSV data into structured records.

    Parameters:
    -----------
    raw_data : str
        Raw CSV data from NIST
    species : str
        Species name for record

    Returns:
    --------
    list : List of dictionaries with parsed line data
    """
    lines = raw_data.strip().split('\n')

    if len(lines) < 2:
        print(f"Warning: No data found for {species}")
        return []

    # Find header line (first line with column names)
    header_line = lines[0]
    headers = [h.strip() for h in header_line.split(',')]

    records = []

    for line_num, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue

        # Parse CSV line (handle quoted values with commas)
        try:
            # Use csv module to handle complex CSV parsing
            reader = csv.reader([line])
            values = next(reader)
        except Exception as e:
            print(f"Warning: Could not parse line {line_num}: {e}")
            continue

        # Create record dictionary
        record = {}
        for i, header in enumerate(headers):
            if i < len(values):
                record[header] = clean_csv_value(values[i])
            else:
                record[header] = ""

        # Add species
        record['species'] = species

        # Skip records without wavelength data
        obs_wl = record.get('obs_wl_vac(nm)', '')
        ritz_wl = record.get('ritz_wl_vac(nm)', '')

        if not obs_wl and not ritz_wl:
            continue

        records.append(record)

    return records


def convert_nm_to_angstrom(nm_str: str) -> str:
    """Convert wavelength from nm to Angstrom."""
    if not nm_str:
        return ""
    try:
        nm = float(nm_str)
        return f"{nm * 10:.4f}"
    except ValueError:
        return ""


def format_level_info(conf: str, term: str, j: str) -> str:
    """Format level information as 'config term J'."""
    parts = []
    if conf:
        parts.append(conf)
    if term:
        parts.append(term)
    if j:
        parts.append(f"J={j}")
    return " ".join(parts)


def create_parsed_table(records: list) -> list:
    """
    Create standardized parsed table from NIST records.

    Output columns:
    - species: Fe V or Ni V
    - wavelength_vac: Vacuum wavelength in Angstrom
    - wavelength_unc: Wavelength uncertainty in Angstrom
    - lower_level: Lower level configuration, term, J
    - upper_level: Upper level configuration, term, J
    - Aki: Transition probability (s^-1)
    - fik: Oscillator strength
    - source: NIST ASD reference
    """
    parsed = []

    for rec in records:
        # Use observed wavelength if available, otherwise Ritz
        obs_wl = rec.get('obs_wl_vac(nm)', '')
        ritz_wl = rec.get('ritz_wl_vac(nm)', '')
        obs_unc = rec.get('unc_obs_wl', '')
        ritz_unc = rec.get('unc_ritz_wl', '')

        if obs_wl:
            wl_nm = obs_wl
            unc_nm = obs_unc
        elif ritz_wl:
            wl_nm = ritz_wl
            unc_nm = ritz_unc
        else:
            continue

        # Convert to Angstrom
        wl_ang = convert_nm_to_angstrom(wl_nm)
        unc_ang = convert_nm_to_angstrom(unc_nm)

        if not wl_ang:
            continue

        # Format level info
        lower_level = format_level_info(
            rec.get('conf_i', ''),
            rec.get('term_i', ''),
            rec.get('J_i', '')
        )
        upper_level = format_level_info(
            rec.get('conf_k', ''),
            rec.get('term_k', ''),
            rec.get('J_k', '')
        )

        # Get transition probability and oscillator strength
        Aki = rec.get('gA(s^-1)', '')
        fik = rec.get('fik', '')

        # Get reference
        line_ref = rec.get('line_ref', '')
        tp_ref = rec.get('tp_ref', '')
        source = f"NIST_ASD_v{NIST_ASD_VERSION}"
        if line_ref:
            source += f";{line_ref}"
        if tp_ref and tp_ref != line_ref:
            source += f";{tp_ref}"

        parsed.append({
            'species': rec.get('species', ''),
            'wavelength_vac': wl_ang,
            'wavelength_unc': unc_ang,
            'lower_level': lower_level,
            'upper_level': upper_level,
            'Aki': Aki,
            'fik': fik,
            'source': source
        })

    return parsed


def save_raw_data(data: str, species: str, timestamp: str) -> Path:
    """Save raw NIST data to file."""
    species_clean = species.replace(' ', '_')
    filename = f"nist_raw_{species_clean}_{timestamp}.csv"
    filepath = DATA_DIR / filename

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(data)

    print(f"Saved raw data to: {filepath}")
    return filepath


def save_parsed_table(records: list, species: str, timestamp: str) -> Path:
    """Save parsed table to CSV."""
    species_clean = species.replace(' ', '_')
    filename = f"nist_parsed_{species_clean}_{timestamp}.csv"
    filepath = DATA_DIR / filename

    fieldnames = ['species', 'wavelength_vac', 'wavelength_unc',
                  'lower_level', 'upper_level', 'Aki', 'fik', 'source']

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Saved parsed data to: {filepath}")
    return filepath


def save_combined_table(all_records: list, timestamp: str) -> Path:
    """Save combined table with all species."""
    filename = f"nist_combined_FeV_NiV_{timestamp}.csv"
    filepath = DATA_DIR / filename

    fieldnames = ['species', 'wavelength_vac', 'wavelength_unc',
                  'lower_level', 'upper_level', 'Aki', 'fik', 'source']

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_records)

    print(f"Saved combined data to: {filepath}")
    return filepath


def save_metadata(timestamp: str, species_stats: dict, wavelength_stats: dict) -> Path:
    """Save retrieval metadata."""
    filename = f"nist_metadata_{timestamp}.txt"
    filepath = DATA_DIR / filename

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write("NIST Atomic Spectra Database - Data Retrieval Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"NIST ASD Version: {NIST_ASD_VERSION}\n")
        f.write(f"Retrieval Timestamp: {timestamp}\n")
        f.write(f"Retrieval Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n\n")

        f.write("Wavelength Information:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Requested Range: {WAVELENGTH_LOW_NM*10:.1f} - {WAVELENGTH_HIGH_NM*10:.1f} Angstrom\n")
        f.write(f"                 ({WAVELENGTH_LOW_NM:.1f} - {WAVELENGTH_HIGH_NM:.1f} nm)\n")
        f.write("Wavelength Type: VACUUM (NIST uses VAC below 2000 Angstrom)\n\n")

        f.write("Species Statistics:\n")
        f.write("-" * 40 + "\n")
        for species, count in species_stats.items():
            f.write(f"{species}: {count} lines\n")
        f.write(f"Total: {sum(species_stats.values())} lines\n\n")

        f.write("Wavelength Coverage:\n")
        f.write("-" * 40 + "\n")
        for species, stats in wavelength_stats.items():
            f.write(f"{species}:\n")
            f.write(f"  Min wavelength: {stats['min']:.4f} Angstrom\n")
            f.write(f"  Max wavelength: {stats['max']:.4f} Angstrom\n")

        f.write("\n" + "=" * 60 + "\n")
        f.write("Data Source: NIST Atomic Spectra Database\n")
        f.write("URL: https://physics.nist.gov/asd\n")
        f.write("Reference: Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team.\n")
        f.write("           NIST Atomic Spectra Database (ver. 5.12), [Online].\n")
        f.write("           Available: https://physics.nist.gov/asd [Accessed: " +
                datetime.now().strftime('%Y-%m-%d') + "]\n")

    print(f"Saved metadata to: {filepath}")
    return filepath


def main():
    """Main function to fetch and process NIST ASD data."""

    print("=" * 70)
    print("NIST Atomic Spectra Database - Data Retrieval for Fe V and Ni V")
    print("Wavelength Range: 1100-1900 Angstrom (110-190 nm) - VACUUM")
    print("=" * 70)
    print()

    # Create output directory
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Generate timestamp for filenames
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # Species to query
    species_list = ['Fe V', 'Ni V']

    all_parsed_records = []
    species_stats = {}
    wavelength_stats = {}
    access_issues = []

    for species in species_list:
        print(f"\n{'='*50}")
        print(f"Processing {species}")
        print('='*50)

        try:
            # Fetch raw data
            raw_data = fetch_nist_data(species, WAVELENGTH_LOW_NM, WAVELENGTH_HIGH_NM)

            # Save raw data
            save_raw_data(raw_data, species, timestamp)

            # Parse data
            records = parse_nist_csv(raw_data, species)
            print(f"Parsed {len(records)} raw records for {species}")

            # Create parsed table
            parsed_records = create_parsed_table(records)
            print(f"Created {len(parsed_records)} parsed records for {species}")

            # Save individual parsed table
            save_parsed_table(parsed_records, species, timestamp)

            # Add to combined records
            all_parsed_records.extend(parsed_records)

            # Calculate statistics
            species_stats[species] = len(parsed_records)

            if parsed_records:
                wavelengths = [float(r['wavelength_vac']) for r in parsed_records if r['wavelength_vac']]
                if wavelengths:
                    wavelength_stats[species] = {
                        'min': min(wavelengths),
                        'max': max(wavelengths)
                    }

        except Exception as e:
            print(f"ERROR processing {species}: {e}")
            access_issues.append(f"{species}: {str(e)}")
            species_stats[species] = 0

    # Save combined table
    if all_parsed_records:
        save_combined_table(all_parsed_records, timestamp)

    # Save metadata
    save_metadata(timestamp, species_stats, wavelength_stats)

    # Print summary
    print("\n" + "=" * 70)
    print("RETRIEVAL SUMMARY")
    print("=" * 70)
    print(f"\nNIST ASD Version: {NIST_ASD_VERSION}")
    print(f"Timestamp: {timestamp}")
    print(f"\nWavelength Range: {WAVELENGTH_LOW_NM*10:.1f} - {WAVELENGTH_HIGH_NM*10:.1f} Angstrom (VACUUM)")
    print(f"                  {WAVELENGTH_LOW_NM:.1f} - {WAVELENGTH_HIGH_NM:.1f} nm (VACUUM)")

    print("\nLines Retrieved:")
    print("-" * 40)
    for species, count in species_stats.items():
        print(f"  {species}: {count} lines")
    print(f"  TOTAL: {sum(species_stats.values())} lines")

    print("\nWavelength Coverage:")
    print("-" * 40)
    for species, stats in wavelength_stats.items():
        print(f"  {species}: {stats['min']:.4f} - {stats['max']:.4f} Angstrom")

    if access_issues:
        print("\nAccess Issues Encountered:")
        print("-" * 40)
        for issue in access_issues:
            print(f"  - {issue}")
    else:
        print("\nNo access issues encountered.")

    print("\nOutput Files:")
    print("-" * 40)
    print(f"  Raw data: {DATA_DIR}/nist_raw_*_{timestamp}.csv")
    print(f"  Parsed data: {DATA_DIR}/nist_parsed_*_{timestamp}.csv")
    print(f"  Combined data: {DATA_DIR}/nist_combined_FeV_NiV_{timestamp}.csv")
    print(f"  Metadata: {DATA_DIR}/nist_metadata_{timestamp}.txt")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

    return species_stats, wavelength_stats, access_issues


if __name__ == "__main__":
    main()
