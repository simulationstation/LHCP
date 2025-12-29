#!/usr/bin/env python3
"""
Build combined atomic line table for Δα/α analysis.

Combines:
- NIST Fe V and Ni V laboratory wavelengths
- Fe V q-coefficients from Hu et al. (2021) arXiv:2007.10905
- Ni V q-coefficients from Webb et al. (2025) arXiv:2410.01849
"""

import csv
from pathlib import Path
from datetime import datetime
import glob

DATA_DIR = Path("/home/primary/LHCP/data")
OUTPUT_DIR = Path("/home/primary/LHCP/data/atomic")

# Wavelength matching tolerance in Angstroms
MATCH_TOLERANCE = 0.05


def load_csv(filepath):
    """Load CSV file into list of dicts."""
    records = []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)
    return records


def find_nist_files():
    """Find the most recent NIST parsed files."""
    fe_pattern = str(DATA_DIR / "atomic/nist/nist_parsed_Fe_V_*.csv")
    ni_pattern = str(DATA_DIR / "atomic/nist/nist_parsed_Ni_V_*.csv")

    fe_files = sorted(glob.glob(fe_pattern))
    ni_files = sorted(glob.glob(ni_pattern))

    return fe_files[-1] if fe_files else None, ni_files[-1] if ni_files else None


def load_fev_q():
    """Load Fe V q-coefficients from Hu et al. (2021)."""
    filepath = DATA_DIR / "atomic/papers/FeV_q_coefficients_Hu2021.csv"
    if not filepath.exists():
        print(f"Warning: Fe V q-coefficient file not found: {filepath}")
        return {}

    records = load_csv(filepath)

    # Build lookup by wavelength
    q_lookup = {}
    for r in records:
        wl = r.get('wavelength_ang', '')
        q = r.get('q_cm-1', '')
        if wl and q:
            try:
                wl_float = float(wl)
                q_lookup[wl_float] = {
                    'q_cm-1': q,
                    'q_source': 'Hu_et_al_2021',
                    'wl_source': r.get('wavelength_source', ''),
                    'f_hu': r.get('f', ''),
                    'gamma_s-1': r.get('gamma_s-1', ''),
                }
            except ValueError:
                pass

    print(f"  Loaded {len(q_lookup)} Fe V q-coefficients from Hu et al. (2021)")
    return q_lookup


def load_niv_q():
    """Load Ni V q-coefficients from Webb et al. (2025)."""
    # Use the 1.0 continuum placement file as default
    filepath = DATA_DIR / "atomic/papers/NiV_q_coefficients_G191_1.0.csv"
    if not filepath.exists():
        print(f"Warning: Ni V q-coefficient file not found: {filepath}")
        return {}

    records = []
    with open(filepath, 'r', encoding='utf-8') as f:
        # Custom parsing for this format
        for line in f:
            line = line.strip()
            if not line or line.startswith('lam_K'):
                continue

            # Parse comma-separated but handle embedded commas in numbers
            parts = line.split(',')
            if len(parts) >= 5:
                try:
                    lam_k = float(parts[0].strip())
                    q = int(parts[4].strip())
                    records.append({'wavelength': lam_k, 'q': q})
                except (ValueError, IndexError):
                    pass

    # Build lookup
    q_lookup = {}
    for r in records:
        q_lookup[r['wavelength']] = {
            'q_cm-1': str(r['q']),
            'q_source': 'Webb_et_al_2025',
            'wl_source': 'K_exp',
        }

    print(f"  Loaded {len(q_lookup)} Ni V q-coefficients from Webb et al. (2025)")
    return q_lookup


def match_wavelength(target_wl, lookup, tolerance=MATCH_TOLERANCE):
    """Find closest matching wavelength in lookup within tolerance."""
    best_match = None
    best_diff = tolerance + 1

    for wl in lookup.keys():
        diff = abs(wl - target_wl)
        if diff < best_diff and diff <= tolerance:
            best_diff = diff
            best_match = wl

    return best_match


def build_combined_table():
    """Build the combined atomic line table."""
    print("Building combined atomic line table...")
    print("=" * 60)

    # Find NIST files
    fe_file, ni_file = find_nist_files()
    if not fe_file or not ni_file:
        print("Error: NIST files not found")
        return None

    print(f"Loading NIST data...")
    print(f"  Fe V: {fe_file}")
    print(f"  Ni V: {ni_file}")

    fe_nist = load_csv(fe_file)
    ni_nist = load_csv(ni_file)

    print(f"  Loaded {len(fe_nist)} Fe V lines from NIST")
    print(f"  Loaded {len(ni_nist)} Ni V lines from NIST")

    # Load q-coefficients
    print(f"\nLoading q-coefficients...")
    fev_q = load_fev_q()
    niv_q = load_niv_q()

    # Build combined table
    combined = []

    # Process Fe V
    print(f"\nMatching Fe V lines...")
    fe_matched = 0
    for line in fe_nist:
        try:
            wl = float(line['wavelength_vac'])
        except (ValueError, KeyError):
            continue

        record = {
            'species': 'Fe V',
            'wavelength_vac_ang': line.get('wavelength_vac', ''),
            'wavelength_unc_ang': line.get('wavelength_unc', ''),
            'lower_level': line.get('lower_level', ''),
            'upper_level': line.get('upper_level', ''),
            'Aki_s-1': line.get('Aki', ''),
            'fik': line.get('fik', ''),
            'nist_source': line.get('source', ''),
            'q_cm-1': '',
            'q_source': '',
            'wl_match_ang': '',
            'f_q': '',
            'gamma_s-1': '',
        }

        # Try to match with q-coefficient
        match_wl = match_wavelength(wl, fev_q, tolerance=MATCH_TOLERANCE)
        if match_wl is not None:
            q_data = fev_q[match_wl]
            record['q_cm-1'] = q_data['q_cm-1']
            record['q_source'] = q_data['q_source']
            record['wl_match_ang'] = f"{match_wl:.4f}"
            record['f_q'] = q_data.get('f_hu', '')
            record['gamma_s-1'] = q_data.get('gamma_s-1', '')
            fe_matched += 1

        combined.append(record)

    print(f"  Fe V: {fe_matched}/{len(fe_nist)} lines matched with q-coefficients")

    # Process Ni V
    print(f"\nMatching Ni V lines...")
    ni_matched = 0
    for line in ni_nist:
        try:
            wl = float(line['wavelength_vac'])
        except (ValueError, KeyError):
            continue

        record = {
            'species': 'Ni V',
            'wavelength_vac_ang': line.get('wavelength_vac', ''),
            'wavelength_unc_ang': line.get('wavelength_unc', ''),
            'lower_level': line.get('lower_level', ''),
            'upper_level': line.get('upper_level', ''),
            'Aki_s-1': line.get('Aki', ''),
            'fik': line.get('fik', ''),
            'nist_source': line.get('source', ''),
            'q_cm-1': '',
            'q_source': '',
            'wl_match_ang': '',
            'f_q': '',
            'gamma_s-1': '',
        }

        # Try to match with q-coefficient
        match_wl = match_wavelength(wl, niv_q, tolerance=MATCH_TOLERANCE)
        if match_wl is not None:
            q_data = niv_q[match_wl]
            record['q_cm-1'] = q_data['q_cm-1']
            record['q_source'] = q_data['q_source']
            record['wl_match_ang'] = f"{match_wl:.3f}"
            ni_matched += 1

        combined.append(record)

    print(f"  Ni V: {ni_matched}/{len(ni_nist)} lines matched with q-coefficients")

    # Sort by wavelength
    combined.sort(key=lambda x: float(x['wavelength_vac_ang']) if x['wavelength_vac_ang'] else 0)

    return combined


def save_combined_table(records, output_path):
    """Save combined table to CSV."""
    if not records:
        print("No records to save!")
        return

    fieldnames = list(records[0].keys())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"\nSaved {len(records)} records to {output_path}")


def print_summary(records):
    """Print summary statistics."""
    if not records:
        return

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    # Count by species
    fe_count = sum(1 for r in records if r['species'] == 'Fe V')
    ni_count = sum(1 for r in records if r['species'] == 'Ni V')

    # Count with q-coefficients
    fe_q = sum(1 for r in records if r['species'] == 'Fe V' and r['q_cm-1'])
    ni_q = sum(1 for r in records if r['species'] == 'Ni V' and r['q_cm-1'])

    # Wavelength ranges
    fe_wl = [float(r['wavelength_vac_ang']) for r in records if r['species'] == 'Fe V' and r['wavelength_vac_ang']]
    ni_wl = [float(r['wavelength_vac_ang']) for r in records if r['species'] == 'Ni V' and r['wavelength_vac_ang']]

    # Q-coefficient ranges for lines that have them
    fe_q_vals = [int(r['q_cm-1']) for r in records if r['species'] == 'Fe V' and r['q_cm-1'] and r['q_cm-1'].lstrip('-').isdigit()]
    ni_q_vals = [int(r['q_cm-1']) for r in records if r['species'] == 'Ni V' and r['q_cm-1'] and r['q_cm-1'].lstrip('-').isdigit()]

    print(f"\nFe V:")
    print(f"  Total lines: {fe_count}")
    print(f"  With q-coefficients: {fe_q} ({100*fe_q/fe_count:.1f}%)")
    if fe_wl:
        print(f"  Wavelength range: {min(fe_wl):.2f} - {max(fe_wl):.2f} Å")
    if fe_q_vals:
        print(f"  Q-coefficient range: {min(fe_q_vals)} to {max(fe_q_vals)} cm^-1")

    print(f"\nNi V:")
    print(f"  Total lines: {ni_count}")
    print(f"  With q-coefficients: {ni_q} ({100*ni_q/ni_count:.1f}%)")
    if ni_wl:
        print(f"  Wavelength range: {min(ni_wl):.2f} - {max(ni_wl):.2f} Å")
    if ni_q_vals:
        print(f"  Q-coefficient range: {min(ni_q_vals)} to {max(ni_q_vals)} cm^-1")

    print(f"\nTotal lines: {len(records)}")
    print(f"Total with q-coefficients: {fe_q + ni_q}")


def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = OUTPUT_DIR / f"combined_lines_q_{timestamp}.csv"

    records = build_combined_table()

    if records:
        save_combined_table(records, output_file)
        print_summary(records)

        # Also create a symlink/copy with fixed name
        fixed_name = OUTPUT_DIR / "combined_lines_q.csv"
        if fixed_name.exists():
            fixed_name.unlink()
        # Copy instead of symlink for portability
        with open(output_file, 'r') as src, open(fixed_name, 'w') as dst:
            dst.write(src.read())
        print(f"\nAlso saved as: {fixed_name}")


if __name__ == "__main__":
    main()
