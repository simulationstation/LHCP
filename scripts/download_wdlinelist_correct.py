#!/usr/bin/env python3
"""
Download WD-LINELIST HLSP products correctly from MAST.
Previous downloads were HTML error pages, not actual FITS files.
"""

import os
import hashlib
from datetime import datetime
from pathlib import Path

try:
    from astroquery.mast import Observations
except ImportError:
    print("ERROR: astroquery not installed")
    exit(1)

HLSP_DIR = Path("/home/primary/LHCP/data/hlsp/wd-linelist")
PROVENANCE_FILE = Path("/home/primary/LHCP/data/hlsp/download_provenance.json")


def get_md5(filepath):
    """Compute MD5 hash."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def is_valid_fits(filepath):
    """Check if file is a valid FITS file."""
    with open(filepath, 'rb') as f:
        header = f.read(80)
    return header.startswith(b'SIMPLE')


def main():
    print("=" * 60)
    print("WD-LINELIST HLSP Download (Corrected)")
    print("=" * 60)

    # Remove corrupted files
    corrupted = []
    for f in HLSP_DIR.glob("*.fits"):
        if not is_valid_fits(f):
            print(f"Removing corrupted file: {f.name}")
            corrupted.append(str(f))
            f.unlink()

    # Query MAST for WD-LINELIST products
    print("\nQuerying MAST for WD-LINELIST products...")

    # Method 1: Search by HLSP collection
    try:
        obs = Observations.query_criteria(
            obs_collection="HLSP",
            provenance_name="WD-LINELIST"
        )
        print(f"Found {len(obs)} WD-LINELIST observations in HLSP collection")

        # Filter for G191-B2B
        g191_obs = []
        for row in obs:
            target = str(row.get('target_name', '')).upper().replace(' ', '').replace('-', '')
            if 'G191' in target or 'WD0501' in target:
                g191_obs.append(row)
                print(f"  Match: {row.get('target_name', '')} - {row.get('obs_id', '')}")

        if g191_obs:
            print(f"\nDownloading {len(g191_obs)} G191-B2B products...")
            for row in g191_obs:
                obsid = row['obsid']
                products = Observations.get_product_list(obsid)

                # Download all products
                manifest = Observations.download_products(
                    products,
                    download_dir=str(HLSP_DIR.parent),
                    cache=False  # Don't use cache, force fresh download
                )

                if manifest is not None:
                    for r in manifest:
                        if r['Status'] == 'COMPLETE':
                            local = Path(r['Local Path'])
                            print(f"  Downloaded: {local.name}")
                            if is_valid_fits(local):
                                print(f"    Valid FITS: YES")
                            else:
                                print(f"    WARNING: Not a valid FITS file!")

    except Exception as e:
        print(f"Error with HLSP query: {e}")

    # Method 2: Try direct product search
    print("\nTrying direct product search for G191-B2B HLSP...")
    try:
        for target in ["G191-B2B", "G191B2B", "g191-b2b", "g191b2b"]:
            obs = Observations.query_criteria(
                target_name=target,
                obs_collection="HLSP"
            )
            if len(obs) > 0:
                print(f"Found {len(obs)} HLSP observations for {target}")
                for row in obs:
                    print(f"  {row.get('obs_id', 'N/A')}: {row.get('dataproduct_type', 'N/A')}")

                    # Get and download products
                    products = Observations.get_product_list(row['obsid'])
                    print(f"    {len(products)} products available")

                    if len(products) > 0:
                        manifest = Observations.download_products(
                            products,
                            download_dir=str(HLSP_DIR),
                            cache=False
                        )
                        if manifest is not None:
                            for r in manifest:
                                if r['Status'] == 'COMPLETE':
                                    print(f"    Downloaded: {r['Local Path']}")
                break  # Found results, stop trying other target names

    except Exception as e:
        print(f"Error with direct search: {e}")

    # Check what we have now
    print("\n" + "=" * 60)
    print("CURRENT HLSP DIRECTORY STATUS")
    print("=" * 60)

    fits_files = list(HLSP_DIR.glob("**/*.fits"))
    valid_fits = []
    invalid_fits = []

    for f in fits_files:
        if is_valid_fits(f):
            valid_fits.append(f)
            md5 = get_md5(f)
            print(f"  [VALID] {f.name}: {f.stat().st_size} bytes, MD5: {md5[:8]}...")
        else:
            invalid_fits.append(f)
            print(f"  [INVALID] {f.name}: NOT A VALID FITS FILE")

    # Write provenance
    import json
    provenance = {
        "retrieval_time": datetime.now().isoformat(),
        "method": "astroquery.mast.Observations",
        "query": "obs_collection=HLSP, target_name=G191-B2B",
        "corrupted_removed": corrupted,
        "valid_fits_downloaded": [str(f) for f in valid_fits],
        "invalid_fits": [str(f) for f in invalid_fits]
    }

    with open(PROVENANCE_FILE, 'w') as f:
        json.dump(provenance, f, indent=2)
    print(f"\nProvenance written to: {PROVENANCE_FILE}")


if __name__ == "__main__":
    main()
