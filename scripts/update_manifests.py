#!/usr/bin/env python3
"""
Update manifest files with current data inventory.
"""

import csv
import hashlib
import os
from pathlib import Path
from datetime import datetime
import glob

DATA_DIR = Path("/home/primary/LHCP/data")
MANIFEST_DIR = Path("/home/primary/LHCP/manifests")


def get_md5(filepath):
    """Compute MD5 hash of file."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_file_info(filepath):
    """Get file metadata."""
    stat = os.stat(filepath)
    return {
        'path': str(filepath),
        'relative_path': str(filepath.relative_to(DATA_DIR.parent)),
        'filename': filepath.name,
        'size_bytes': stat.st_size,
        'modified': datetime.fromtimestamp(stat.st_mtime).isoformat(),
        'md5': get_md5(filepath),
    }


def categorize_file(filepath):
    """Categorize file by type."""
    path_str = str(filepath)
    name = filepath.name.lower()

    if 'hlsp' in path_str:
        if 'calspec' in path_str:
            return 'CALSPEC spectrum'
        elif 'wd-linelist' in path_str:
            return 'WD-LINELIST HLSP'
        else:
            return 'HLSP product'
    elif 'raw' in path_str:
        if '_x1d' in name:
            return 'Raw STIS x1d'
        elif '_wav' in name:
            return 'Raw STIS wavecal'
        else:
            return 'Raw STIS product'
    elif 'nist' in path_str:
        if 'combined' in name:
            return 'Combined NIST table'
        elif 'parsed' in name:
            return 'Parsed NIST data'
        elif 'raw' in name:
            return 'Raw NIST download'
        else:
            return 'NIST data'
    elif 'papers' in path_str:
        if 'q_coefficient' in name.lower():
            return 'Q-coefficient table'
        elif '.tar.gz' in name:
            return 'arXiv source archive'
        elif 'source_' in path_str:
            return 'arXiv extracted file'
        else:
            return 'Literature data'
    elif 'stellar' in name:
        return 'Stellar parameters'
    elif 'combined_lines_q' in name:
        return 'Combined atomic table'
    elif 'mast' in name:
        return 'MAST query results'
    else:
        return 'Data file'


def update_file_manifest():
    """Update manifest_files.csv with all data files."""
    print("Updating file manifest...")

    # Find all files in data/
    all_files = []

    for pattern in ['**/*.fits', '**/*.csv', '**/*.json', '**/*.txt', '**/*.md', '**/*.tar.gz', '**/*.tex']:
        matches = DATA_DIR.glob(pattern)
        for f in matches:
            if f.is_file():
                all_files.append(f)

    # Get info for each file
    records = []
    for filepath in sorted(all_files):
        try:
            info = get_file_info(filepath)
            info['category'] = categorize_file(filepath)
            records.append(info)
        except Exception as e:
            print(f"  Warning: Could not process {filepath}: {e}")

    # Save manifest
    output_path = MANIFEST_DIR / "manifest_files.csv"
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

    fieldnames = ['relative_path', 'filename', 'category', 'size_bytes', 'modified', 'md5']

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)

    print(f"  Saved {len(records)} files to {output_path}")

    # Print summary by category
    categories = {}
    for r in records:
        cat = r['category']
        categories[cat] = categories.get(cat, 0) + 1

    print("\n  File counts by category:")
    for cat, count in sorted(categories.items()):
        print(f"    {cat}: {count}")

    return records


def update_lines_manifest():
    """Update manifest_lines.csv with combined line table summary."""
    print("\nUpdating lines manifest...")

    # Find the combined lines file
    combined_file = DATA_DIR / "atomic" / "combined_lines_q.csv"
    if not combined_file.exists():
        print(f"  Warning: Combined lines file not found: {combined_file}")
        return

    # Read and summarize
    with open(combined_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        lines = list(reader)

    # Create summary manifest
    output_path = MANIFEST_DIR / "manifest_lines.csv"

    # Just copy the combined file as the lines manifest
    with open(combined_file, 'r') as src, open(output_path, 'w') as dst:
        dst.write(src.read())

    print(f"  Saved {len(lines)} lines to {output_path}")

    # Summary stats
    fe_total = sum(1 for l in lines if l['species'] == 'Fe V')
    ni_total = sum(1 for l in lines if l['species'] == 'Ni V')
    fe_q = sum(1 for l in lines if l['species'] == 'Fe V' and l['q_cm-1'])
    ni_q = sum(1 for l in lines if l['species'] == 'Ni V' and l['q_cm-1'])

    print(f"\n  Summary:")
    print(f"    Fe V: {fe_total} lines ({fe_q} with q)")
    print(f"    Ni V: {ni_total} lines ({ni_q} with q)")
    print(f"    Total: {len(lines)} lines ({fe_q + ni_q} with q)")


def update_obs_manifest():
    """Check observation manifest status."""
    print("\nChecking observation manifest...")

    obs_file = MANIFEST_DIR / "manifest_obs.csv"
    if obs_file.exists():
        with open(obs_file, 'r') as f:
            reader = csv.DictReader(f)
            obs = list(reader)
        print(f"  manifest_obs.csv exists with {len(obs)} observations")

        # Count by instrument
        instruments = {}
        for o in obs:
            inst = o.get('instrument', 'Unknown')
            instruments[inst] = instruments.get(inst, 0) + 1

        print("  Observations by instrument:")
        for inst, count in sorted(instruments.items()):
            print(f"    {inst}: {count}")
    else:
        print("  Warning: manifest_obs.csv not found")


def main():
    print("=" * 60)
    print("Updating LHCP Manifests")
    print("=" * 60)

    update_file_manifest()
    update_lines_manifest()
    update_obs_manifest()

    print("\n" + "=" * 60)
    print("Manifest update complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
