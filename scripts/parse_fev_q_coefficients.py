#!/usr/bin/env python3
"""
Parse Fe V q-coefficients from Hu et al. (2021) arXiv:2007.10905
LaTeX source table (summary_atomic_data.tex).

This script extracts q-coefficients that supersede Ong et al. (2013).
"""

import re
import csv
from pathlib import Path

LATEX_FILE = Path("/home/primary/LHCP/data/atomic/papers/source_2007.10905/table/summary_atomic_data.tex")
OUTPUT_FILE = Path("/home/primary/LHCP/data/atomic/papers/FeV_q_coefficients_Hu2021.csv")

def parse_latex_table(filepath):
    """Parse Fe V atomic data from LaTeX longtable."""
    records = []

    with open(filepath, 'r') as f:
        content = f.read()

    # Find all data rows (after header, before endlastfoot)
    # Pattern: ID& config1 & term1 & J1 & config2 & term2 & J2 & E1 & E2 & wavelength...
    lines = content.split('\n')

    for line in lines:
        # Skip non-data lines
        if not re.match(r'^\d+&', line.strip()):
            continue

        # Clean up the line
        line = line.strip()
        if line.endswith('\\\\'):
            line = line[:-2]

        # Split by &
        parts = [p.strip() for p in line.split('&')]

        if len(parts) < 17:
            continue

        try:
            # Parse fields
            line_id = int(parts[0])
            config_lower = parts[1].strip()
            term_lower = parts[2].strip()
            j_lower = parts[3].strip()
            config_upper = parts[4].strip()
            term_upper = parts[5].strip()
            j_upper = parts[6].strip()
            E_lower = parts[7].strip()
            E_upper = parts[8].strip()

            # Wavelengths (K14 or W19)
            wl_k14 = parts[9].strip()
            wl_k14_unc = parts[11].strip() if len(parts) > 11 else ''
            wl_w19 = parts[12].strip() if len(parts) > 12 else ''
            wl_w19_unc = parts[14].strip() if len(parts) > 14 else ''

            # Oscillator strength
            f_val = parts[15].strip() if len(parts) > 15 else ''

            # Damping constant
            gamma = parts[16].strip() if len(parts) > 16 else ''

            # Q-coefficient (last column)
            q_val = parts[17].strip() if len(parts) > 17 else ''

            # Clean up wavelength values (remove $ markers, stars, diamonds)
            wl_k14 = re.sub(r'[\$\*]', '', wl_k14).strip()
            wl_w19 = re.sub(r'[\$\*\diamond]', '', wl_w19).strip()
            wl_k14_unc = re.sub(r'[\$\*]', '', wl_k14_unc).strip()
            wl_w19_unc = re.sub(r'[\$\*\diamond]', '', wl_w19_unc).strip()

            # Clean gamma (remove footnote markers)
            gamma = re.sub(r'\$[^\$]*\$', '', gamma).strip()
            gamma = re.sub(r'\^\mathsection', '', gamma).strip()

            # Clean f value (remove footnotes)
            f_val = re.sub(r'\$[^\$]*\$', '', f_val).strip()
            f_val = re.sub(r'\^\dagger', '', f_val).strip()

            # Clean q value
            q_val = re.sub(r'\s+', '', q_val).strip()

            # Skip if no q-value
            if not q_val:
                continue

            # Use K14 wavelength preferentially, else W19
            wavelength = wl_k14 if wl_k14 else wl_w19
            wavelength_unc = wl_k14_unc if wl_k14 else wl_w19_unc
            wavelength_source = 'K14' if wl_k14 else 'W19'

            record = {
                'id': line_id,
                'species': 'Fe V',
                'config_lower': config_lower,
                'term_lower': term_lower,
                'J_lower': j_lower,
                'config_upper': config_upper,
                'term_upper': term_upper,
                'J_upper': j_upper,
                'E_lower_cm-1': E_lower,
                'E_upper_cm-1': E_upper,
                'wavelength_ang': wavelength,
                'wavelength_unc_ang': wavelength_unc,
                'wavelength_source': wavelength_source,
                'f': f_val,
                'gamma_s-1': gamma,
                'q_cm-1': q_val,
                'source': 'Hu_et_al_2021_MNRAS_500_1466'
            }
            records.append(record)

        except Exception as e:
            print(f"Warning: Could not parse line: {line[:80]}... Error: {e}")
            continue

    return records


def save_csv(records, output_path):
    """Save records to CSV."""
    if not records:
        print("No records to save!")
        return

    fieldnames = list(records[0].keys())

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Saved {len(records)} records to {output_path}")


def main():
    print("Parsing Fe V q-coefficients from Hu et al. (2021)")
    print(f"Input: {LATEX_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print()

    records = parse_latex_table(LATEX_FILE)

    print(f"Parsed {len(records)} Fe V lines with q-coefficients")

    # Filter to only those with valid q values
    valid_records = [r for r in records if r['q_cm-1']]
    print(f"Lines with valid q-coefficients: {len(valid_records)}")

    if valid_records:
        save_csv(valid_records, OUTPUT_FILE)

        # Print sample
        print("\nSample records:")
        for r in valid_records[:5]:
            print(f"  λ={r['wavelength_ang']} Å, q={r['q_cm-1']} cm^-1")

        # Statistics
        wavelengths = [float(r['wavelength_ang']) for r in valid_records if r['wavelength_ang']]
        q_values = [int(r['q_cm-1']) for r in valid_records if r['q_cm-1'].isdigit() or (r['q_cm-1'].startswith('-') and r['q_cm-1'][1:].isdigit())]

        print(f"\nWavelength range: {min(wavelengths):.4f} - {max(wavelengths):.4f} Å")
        print(f"Q-coefficient range: {min(q_values)} to {max(q_values)} cm^-1")


if __name__ == "__main__":
    main()
