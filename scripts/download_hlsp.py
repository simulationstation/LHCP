#!/usr/bin/env python3
"""
Download HLSP (High-Level Science Products) for G191-B2B from MAST.

Uses astroquery to search for and download the Preval et al. products
from the WD-LINELIST HLSP collection.
"""

import os
from pathlib import Path

try:
    from astroquery.mast import Observations
except ImportError:
    print("ERROR: astroquery not installed. Run: pip3 install astroquery --break-system-packages")
    exit(1)

# Output directory
HLSP_DIR = Path("/home/primary/LHCP/data/hlsp")
TARGET_NAMES = ["G191-B2B", "G191B2B", "G191 B2B", "WD0501+527", "WD 0501+527"]


def search_hlsp():
    """Search for HLSP products related to G191-B2B."""
    print("Searching MAST for G191-B2B HLSP products...")

    all_products = []

    # Try different search strategies

    # Strategy 1: Search by provenance_name for known HLSPs
    hlsp_provenances = [
        "WD-LINELIST",
        "MUSCLES",
        "CALSPEC",
        "STISCLEAN",
    ]

    for prov in hlsp_provenances:
        print(f"  Searching provenance: {prov}")
        try:
            obs = Observations.query_criteria(
                provenance_name=prov,
                dataproduct_type="spectrum"
            )
            if len(obs) > 0:
                print(f"    Found {len(obs)} observations")
                # Filter for G191-B2B
                for row in obs:
                    target = str(row.get('target_name', '')).upper().replace(' ', '').replace('-', '')
                    if 'G191' in target or 'WD0501' in target:
                        all_products.append({
                            'obsid': row['obsid'],
                            'target': row.get('target_name', ''),
                            'provenance': prov,
                            'obs_collection': row.get('obs_collection', ''),
                            'instrument': row.get('instrument_name', ''),
                            'filters': row.get('filters', ''),
                        })
                        print(f"      Match: {row.get('target_name', '')} ({row['obsid']})")
        except Exception as e:
            print(f"    Error: {e}")

    # Strategy 2: Direct target search for high-level products
    for target_name in TARGET_NAMES:
        print(f"  Searching target: {target_name}")
        try:
            obs = Observations.query_criteria(
                target_name=target_name,
                dataproduct_type="spectrum",
                intentType="science"
            )
            if len(obs) > 0:
                # Filter for HLSP-like products
                for row in obs:
                    prov = str(row.get('provenance_name', '')).upper()
                    coll = str(row.get('obs_collection', '')).upper()
                    if 'HLSP' in coll or prov in ['WD-LINELIST', 'CALSPEC', 'MUSCLES']:
                        all_products.append({
                            'obsid': row['obsid'],
                            'target': row.get('target_name', ''),
                            'provenance': row.get('provenance_name', ''),
                            'obs_collection': row.get('obs_collection', ''),
                            'instrument': row.get('instrument_name', ''),
                            'filters': row.get('filters', ''),
                        })
        except Exception as e:
            print(f"    Error: {e}")

    # Deduplicate
    seen = set()
    unique_products = []
    for p in all_products:
        if p['obsid'] not in seen:
            seen.add(p['obsid'])
            unique_products.append(p)

    return unique_products


def download_products(observations, output_dir):
    """Download products for given observations."""
    output_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []

    for obs in observations:
        obsid = obs['obsid']
        print(f"\nProcessing obsid: {obsid} ({obs.get('provenance', '')})")

        try:
            # Get product list
            products = Observations.get_product_list(obsid)

            if len(products) == 0:
                print(f"  No products found for {obsid}")
                continue

            print(f"  Found {len(products)} products")

            # Filter for science products (spectra)
            science_products = Observations.filter_products(
                products,
                productType=["SCIENCE"],
                extension=["fits", "FITS"]
            )

            if len(science_products) == 0:
                # Try without extension filter
                science_products = Observations.filter_products(
                    products,
                    productType=["SCIENCE"]
                )

            if len(science_products) == 0:
                print(f"  No science products found, downloading all")
                science_products = products

            print(f"  Downloading {len(science_products)} products...")

            # Download
            manifest = Observations.download_products(
                science_products,
                download_dir=str(output_dir),
                cache=True
            )

            if manifest is not None and len(manifest) > 0:
                for row in manifest:
                    local_path = row['Local Path']
                    status = row['Status']
                    if status == 'COMPLETE':
                        downloaded.append(local_path)
                        print(f"    Downloaded: {os.path.basename(local_path)}")
                    else:
                        print(f"    {status}: {os.path.basename(local_path)}")

        except Exception as e:
            print(f"  Error downloading {obsid}: {e}")

    return downloaded


def search_calspec():
    """Search specifically for CALSPEC standard star data."""
    print("\nSearching CALSPEC for G191-B2B...")

    try:
        # CALSPEC is in the HST obs_collection
        obs = Observations.query_criteria(
            obs_collection="HST",
            provenance_name="CALSPEC"
        )

        if len(obs) > 0:
            print(f"Found {len(obs)} CALSPEC observations")
            g191_obs = []
            for row in obs:
                target = str(row.get('target_name', '')).upper().replace(' ', '').replace('-', '')
                if 'G191' in target:
                    g191_obs.append({
                        'obsid': row['obsid'],
                        'target': row.get('target_name', ''),
                        'provenance': 'CALSPEC',
                        'obs_collection': row.get('obs_collection', ''),
                        'instrument': row.get('instrument_name', ''),
                        'filters': row.get('filters', ''),
                    })
                    print(f"  Found: {row.get('target_name', '')} ({row['obsid']})")
            return g191_obs
    except Exception as e:
        print(f"Error searching CALSPEC: {e}")

    return []


def main():
    print("=" * 60)
    print("HLSP Download for G191-B2B")
    print("=" * 60)

    # Search for products
    products = search_hlsp()

    # Also search CALSPEC
    calspec = search_calspec()
    products.extend(calspec)

    if not products:
        print("\nNo HLSP products found via astroquery.")
        print("Alternative: Try manual download from MAST portal:")
        print("  https://archive.stsci.edu/prepds/wd-linelist/")
        print("  https://archive.stsci.edu/hlsp/calspec")
        return

    print(f"\nFound {len(products)} unique HLSP observations")

    # Download
    downloaded = download_products(products, HLSP_DIR)

    print("\n" + "=" * 60)
    print(f"Download complete. {len(downloaded)} files downloaded to {HLSP_DIR}")
    print("=" * 60)

    # List downloaded files
    if downloaded:
        print("\nDownloaded files:")
        for f in downloaded:
            print(f"  {f}")


if __name__ == "__main__":
    main()
