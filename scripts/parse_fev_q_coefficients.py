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

    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    lines = content.split('\n')

    for line in lines:
        # Skip non-data lines - data lines start with a number followed by &
        stripped = line.strip()
        if not re.match(r'^\d+\s*&', stripped):
            continue

        # Remove trailing \\ if present
        if stripped.endswith('\\\\'):
            stripped = stripped[:-2].strip()

        # Split by &
        parts = [p.strip() for p in stripped.split('&')]

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

            # Energy levels
            E_lower = parts[7].strip()
            E_upper = parts[8].strip()

            # K14 wavelength and uncertainty (columns 9, 10, 11)
            wl_k14_raw = parts[9].strip()
            k14_flag = parts[10].strip() if len(parts) > 10 else ''
            wl_k14_unc_raw = parts[11].strip() if len(parts) > 11 else ''

            # W19 wavelength and uncertainty (columns 12, 13, 14)
            wl_w19_raw = parts[12].strip() if len(parts) > 12 else ''
            w19_flag = parts[13].strip() if len(parts) > 13 else ''
            wl_w19_unc_raw = parts[14].strip() if len(parts) > 14 else ''

            # Oscillator strength (column 15)
            f_val_raw = parts[15].strip() if len(parts) > 15 else ''

            # Damping constant (column 16)
            gamma_raw = parts[16].strip() if len(parts) > 16 else ''

            # Q-coefficient (column 17)
            q_val_raw = parts[17].strip() if len(parts) > 17 else ''

            # Clean wavelength values - remove LaTeX markers and flags
            def clean_wavelength(wl):
                if not wl:
                    return ''
                # Remove $\star$, $\diamond$, and other markers
                wl = re.sub(r'\$[^$]*\$', '', wl)
                wl = re.sub(r'\\star', '', wl)
                wl = re.sub(r'\\diamond', '', wl)
                return wl.strip()

            def clean_numeric(val):
                if not val:
                    return ''
                # Remove LaTeX markers, footnotes
                val = re.sub(r'\$[^$]*\$', '', val)
                val = re.sub(r'\^\s*\\?mathsection', '', val)
                val = re.sub(r'\^\s*\\?dagger', '', val)
                val = re.sub(r'\\mathsection', '', val)
                val = re.sub(r'\\dagger', '', val)
                return val.strip()

            wl_k14 = clean_wavelength(wl_k14_raw)
            wl_k14_unc = clean_wavelength(wl_k14_unc_raw)
            wl_w19 = clean_wavelength(wl_w19_raw)
            wl_w19_unc = clean_wavelength(wl_w19_unc_raw)
            f_val = clean_numeric(f_val_raw)
            gamma = clean_numeric(gamma_raw)
            q_val = clean_numeric(q_val_raw)

            # Check flags for marked lines
            is_k14_used = '$\\star$' in k14_flag or '\\star' in k14_flag or '*' in wl_k14_raw
            is_w19_used = '$\\diamond$' in w19_flag or '\\diamond' in w19_flag

            # Use K14 wavelength preferentially, else W19
            if wl_k14:
                wavelength = wl_k14
                wavelength_unc = wl_k14_unc
                wavelength_source = 'K14'
            elif wl_w19:
                wavelength = wl_w19
                wavelength_unc = wl_w19_unc
                wavelength_source = 'W19'
            else:
                wavelength = ''
                wavelength_unc = ''
                wavelength_source = ''

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
                'wavelength_K14_ang': wl_k14,
                'wavelength_K14_unc_ang': wl_k14_unc,
                'wavelength_W19_ang': wl_w19,
                'wavelength_W19_unc_ang': wl_w19_unc,
                'wavelength_ang': wavelength,
                'wavelength_unc_ang': wavelength_unc,
                'wavelength_source': wavelength_source,
                'used_K14': is_k14_used,
                'used_W19': is_w19_used,
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Saved {len(records)} records to {output_path}")


def main():
    print("Parsing Fe V q-coefficients from Hu et al. (2021)")
    print(f"Input: {LATEX_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print()

    if not LATEX_FILE.exists():
        print(f"ERROR: Input file not found: {LATEX_FILE}")
        return

    records = parse_latex_table(LATEX_FILE)
    print(f"Parsed {len(records)} Fe V lines total")

    # Filter to those with valid q values
    valid_q_records = [r for r in records if r['q_cm-1']]
    print(f"Lines with q-coefficients: {len(valid_q_records)}")

    # Filter to those with valid wavelengths
    valid_records = [r for r in records if r['wavelength_ang']]
    print(f"Lines with wavelengths: {len(valid_records)}")

    if records:
        save_csv(records, OUTPUT_FILE)

        # Print sample
        print("\nSample records (first 5):")
        for r in records[:5]:
            print(f"  ID={r['id']}, λ={r['wavelength_ang']} Å ({r['wavelength_source']}), q={r['q_cm-1']} cm^-1")

        # Statistics
        wavelengths = []
        q_values = []
        for r in records:
            if r['wavelength_ang']:
                try:
                    wavelengths.append(float(r['wavelength_ang']))
                except ValueError:
                    pass
            if r['q_cm-1']:
                try:
                    q_values.append(int(r['q_cm-1']))
                except ValueError:
                    pass

        if wavelengths:
            print(f"\nWavelength range: {min(wavelengths):.4f} - {max(wavelengths):.4f} Å")
        if q_values:
            print(f"Q-coefficient range: {min(q_values)} to {max(q_values)} cm^-1")

        # Count marked lines
        k14_used = sum(1 for r in records if r['used_K14'])
        w19_used = sum(1 for r in records if r['used_W19'])
        print(f"\nLines marked with star (K14 used): {k14_used}")
        print(f"Lines marked with diamond (W19 used): {w19_used}")


if __name__ == "__main__":
    main()
