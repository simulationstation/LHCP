#!/usr/bin/env python3
"""
Analyze q-coefficient join quality and build analysis-ready line tables.
"""

import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

DATA_DIR = Path("/home/primary/LHCP/data/atomic")


def load_combined_table():
    """Load the combined lines table."""
    filepath = DATA_DIR / "combined_lines_q.csv"
    records = []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)
    return records


def analyze_joins(records):
    """Analyze wavelength match quality."""
    stats = {
        'total': len(records),
        'with_q': 0,
        'fe_total': 0,
        'fe_with_q': 0,
        'ni_total': 0,
        'ni_with_q': 0,
        'delta_gt_005': 0,
        'delta_gt_01': 0,
        'delta_gt_02': 0,
        'deltas': [],
        'duplicate_wl_matches': defaultdict(list),
        'issues': [],
        'worst_delta': 0,
        'worst_delta_line': None,
    }

    for r in records:
        species = r['species']
        q = r.get('q_cm-1', '').strip()
        wl_nist = r.get('wavelength_vac_ang', '').strip()
        wl_match = r.get('wl_match_ang', '').strip()

        if species == 'Fe V':
            stats['fe_total'] += 1
        elif species == 'Ni V':
            stats['ni_total'] += 1

        if q:
            stats['with_q'] += 1
            if species == 'Fe V':
                stats['fe_with_q'] += 1
            elif species == 'Ni V':
                stats['ni_with_q'] += 1

            # Track duplicate wl_match usage
            if wl_match:
                stats['duplicate_wl_matches'][wl_match].append(r)

            # Compute delta
            if wl_nist and wl_match:
                try:
                    delta = abs(float(wl_nist) - float(wl_match))
                    stats['deltas'].append({
                        'delta': delta,
                        'species': species,
                        'wl_nist': wl_nist,
                        'wl_match': wl_match,
                        'q': q
                    })

                    if delta > stats['worst_delta']:
                        stats['worst_delta'] = delta
                        stats['worst_delta_line'] = r

                    if delta > 0.005:
                        stats['delta_gt_005'] += 1
                    if delta > 0.01:
                        stats['delta_gt_01'] += 1
                    if delta > 0.02:
                        stats['delta_gt_02'] += 1

                except ValueError:
                    pass

    # Find actual duplicates (same wl_match mapped to multiple NIST lines)
    stats['actual_duplicates'] = {
        k: v for k, v in stats['duplicate_wl_matches'].items()
        if len(v) > 1
    }

    return stats


def write_diagnostics(stats, output_path):
    """Write join diagnostics report."""
    with open(output_path, 'w') as f:
        f.write("# Q-Coefficient Join Diagnostics\n\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n\n")

        f.write("## Summary\n\n")
        f.write(f"| Metric | Value |\n")
        f.write(f"|--------|-------|\n")
        f.write(f"| Total lines | {stats['total']} |\n")
        f.write(f"| Lines with q | {stats['with_q']} |\n")
        f.write(f"| Fe V total | {stats['fe_total']} |\n")
        f.write(f"| Fe V with q | {stats['fe_with_q']} ({100*stats['fe_with_q']/stats['fe_total']:.1f}%) |\n")
        f.write(f"| Ni V total | {stats['ni_total']} |\n")
        f.write(f"| Ni V with q | {stats['ni_with_q']} ({100*stats['ni_with_q']/stats['ni_total']:.1f}%) |\n")
        f.write("\n")

        f.write("## Wavelength Match Quality\n\n")
        f.write(f"| Threshold | Count | Percentage |\n")
        f.write(f"|-----------|-------|------------|\n")
        f.write(f"| |Δλ| > 0.005 Å | {stats['delta_gt_005']} | {100*stats['delta_gt_005']/stats['with_q']:.1f}% |\n")
        f.write(f"| |Δλ| > 0.01 Å | {stats['delta_gt_01']} | {100*stats['delta_gt_01']/stats['with_q']:.1f}% |\n")
        f.write(f"| |Δλ| > 0.02 Å | {stats['delta_gt_02']} | {100*stats['delta_gt_02']/stats['with_q']:.1f}% |\n")
        f.write("\n")

        f.write(f"**Worst match:** Δλ = {stats['worst_delta']:.4f} Å\n")
        if stats['worst_delta_line']:
            wl = stats['worst_delta_line']
            f.write(f"  - Species: {wl['species']}\n")
            f.write(f"  - NIST λ: {wl['wavelength_vac_ang']} Å\n")
            f.write(f"  - Match λ: {wl['wl_match_ang']} Å\n")
            f.write(f"  - q: {wl['q_cm-1']} cm⁻¹\n")
        f.write("\n")

        f.write("## Duplicate Wavelength Matches\n\n")
        f.write(f"Number of q-wavelengths matched to multiple NIST lines: {len(stats['actual_duplicates'])}\n\n")

        if stats['actual_duplicates']:
            f.write("| q-λ (Å) | # NIST lines | NIST wavelengths |\n")
            f.write("|---------|--------------|------------------|\n")
            for wl, lines in sorted(stats['actual_duplicates'].items())[:20]:
                nist_wls = [l['wavelength_vac_ang'] for l in lines]
                f.write(f"| {wl} | {len(lines)} | {', '.join(nist_wls[:5])} |\n")
            if len(stats['actual_duplicates']) > 20:
                f.write(f"\n... and {len(stats['actual_duplicates']) - 20} more\n")
        f.write("\n")

        f.write("## Per-Species Summary\n\n")
        f.write("### Fe V\n")
        fe_deltas = [d for d in stats['deltas'] if d['species'] == 'Fe V']
        if fe_deltas:
            f.write(f"- Lines with q: {stats['fe_with_q']}\n")
            f.write(f"- |Δλ| > 0.005 Å: {sum(1 for d in fe_deltas if d['delta'] > 0.005)}\n")
            f.write(f"- Max |Δλ|: {max(d['delta'] for d in fe_deltas):.4f} Å\n")
        f.write("\n")

        f.write("### Ni V\n")
        ni_deltas = [d for d in stats['deltas'] if d['species'] == 'Ni V']
        if ni_deltas:
            f.write(f"- Lines with q: {stats['ni_with_q']}\n")
            f.write(f"- |Δλ| > 0.005 Å: {sum(1 for d in ni_deltas if d['delta'] > 0.005)}\n")
            f.write(f"- Max |Δλ|: {max(d['delta'] for d in ni_deltas):.4f} Å\n")
        f.write("\n")

    print(f"Diagnostics written to: {output_path}")
    return stats


def build_analysis_table(records, stats, output_path, clean_output_path):
    """Build analysis-ready line table with proper λ0 selection."""

    # Identify which wl_match values are ambiguous (used multiple times)
    ambiguous_matches = set(stats['actual_duplicates'].keys())

    analysis_records = []
    clean_records = []

    for r in records:
        q = r.get('q_cm-1', '').strip()
        wl_nist = r.get('wavelength_vac_ang', '').strip()
        wl_unc = r.get('wavelength_unc_ang', '').strip()
        wl_match = r.get('wl_match_ang', '').strip()
        q_source = r.get('q_source', '').strip()

        # Determine join_flag
        if not q:
            join_flag = 'NO_Q'
            lambda0 = wl_nist
            lambda0_source = 'NIST'
        else:
            # Use wl_match as λ0 when q is present (for consistency)
            lambda0 = wl_match if wl_match else wl_nist
            lambda0_source = q_source.replace('Hu_et_al_2021', 'Hu2021_table')

            # Correct Ni V source naming: Webb_et_al_2025 -> Lee2025_2410.01849
            if 'Webb' in lambda0_source:
                lambda0_source = 'Lee2025_2410.01849'

            # Check for issues
            if wl_match in ambiguous_matches:
                join_flag = 'AMBIGUOUS'
            elif wl_nist and wl_match:
                try:
                    delta = abs(float(wl_nist) - float(wl_match))
                    if delta > 0.005:
                        join_flag = 'LARGE_MISMATCH'
                    else:
                        join_flag = 'OK'
                except ValueError:
                    join_flag = 'OK'
            else:
                join_flag = 'OK'

        analysis_record = {
            'species': r['species'],
            'lambda0_ang': lambda0,
            'sigma_lambda0_ang': wl_unc,
            'q_cm-1': q,
            'sigma_q_cm-1': '',  # Not available in source
            'lambda0_source': lambda0_source if q else 'NIST',
            'join_flag': join_flag,
            'nist_wavelength_ang': wl_nist,
            'lower_level': r.get('lower_level', ''),
            'upper_level': r.get('upper_level', ''),
            'Aki_s-1': r.get('Aki_s-1', ''),
            'fik': r.get('fik', ''),
        }
        analysis_records.append(analysis_record)

        if join_flag == 'OK':
            clean_records.append(analysis_record)

    # Write full analysis table
    fieldnames = list(analysis_records[0].keys())
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(analysis_records)

    print(f"Analysis table written: {output_path} ({len(analysis_records)} lines)")

    # Write clean subset
    with open(clean_output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(clean_records)

    print(f"Clean table written: {clean_output_path} ({len(clean_records)} lines)")

    # Return stats
    flag_counts = defaultdict(int)
    for r in analysis_records:
        flag_counts[r['join_flag']] += 1

    return {
        'total': len(analysis_records),
        'clean': len(clean_records),
        'flag_counts': dict(flag_counts),
        'with_q': sum(1 for r in analysis_records if r['q_cm-1']),
    }


def main():
    print("=" * 60)
    print("Q-Coefficient Join Analysis")
    print("=" * 60)

    # Load data
    print("\nLoading combined_lines_q.csv...")
    records = load_combined_table()
    print(f"Loaded {len(records)} records")

    # Analyze
    print("\nAnalyzing join quality...")
    stats = analyze_joins(records)

    # Write diagnostics
    diag_path = DATA_DIR / "join_diagnostics.md"
    write_diagnostics(stats, diag_path)

    # Build analysis tables
    print("\nBuilding analysis tables...")
    analysis_path = DATA_DIR / "analysis_lines.csv"
    clean_path = DATA_DIR / "analysis_lines_clean.csv"
    table_stats = build_analysis_table(records, stats, analysis_path, clean_path)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total lines: {table_stats['total']}")
    print(f"Lines with q: {table_stats['with_q']}")
    print(f"Clean lines (join_flag=OK): {table_stats['clean']}")
    print(f"\nJoin flag distribution:")
    for flag, count in sorted(table_stats['flag_counts'].items()):
        print(f"  {flag}: {count}")
    print(f"\nWorst |Δλ|: {stats['worst_delta']:.4f} Å")
    print(f"Lines with |Δλ| > 0.005 Å: {stats['delta_gt_005']}")


if __name__ == "__main__":
    main()
