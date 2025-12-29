#!/usr/bin/env python3
"""
Download HLSP products for G191-B2B using direct MAST URLs and alternative methods.
"""

import os
import urllib.request
import urllib.error
from pathlib import Path
import hashlib
import json

HLSP_DIR = Path("/home/primary/LHCP/data/hlsp")
HLSP_DIR.mkdir(parents=True, exist_ok=True)

# Known CALSPEC files for G191-B2B (multiple versions/formats exist)
CALSPEC_BASE_URL = "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/"
CALSPEC_FILES = [
    "g191b2b_stisnic_008.fits",     # STIS+NICMOS composite
    "g191b2b_mod_011.fits",         # Model spectrum
    "g191b2b_fos_003.fits",         # FOS spectrum
    "g191b2b_stis_006.fits",        # STIS only
]

# Alternative CALSPEC URLs (older CDN)
CALSPEC_ALT_URLS = [
    "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/calspec/",
    "ftp://ftp.stsci.edu/cdbs/current_calspec/",
]

# WD-LINELIST HLSP URLs (if available)
WD_LINELIST_BASE = "https://archive.stsci.edu/hlsp/wd-linelist/"
WD_LINELIST_FILES = [
    "g191b2b/hlsp_wd-linelist_multi_multi_g191-b2b_multi_v1_sed.fits",
    "g191b2b/hlsp_wd-linelist_hst_stis_g191-b2b_e140h-e230h_v1_cspec.fits",
]

# Alternative: MAST CAOM API query
MAST_API_URL = "https://mast.stsci.edu/api/v0.1/Download/file"


def download_file(url, dest_path, description=""):
    """Download a file with progress indication."""
    print(f"  Downloading: {os.path.basename(dest_path)}")
    print(f"    URL: {url}")

    try:
        # Add headers to avoid 403 errors
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0 (LHCP Data Pipeline)')

        with urllib.request.urlopen(req, timeout=60) as response:
            content = response.read()

        dest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(dest_path, 'wb') as f:
            f.write(content)

        # Compute MD5
        md5 = hashlib.md5(content).hexdigest()
        print(f"    Success: {len(content)} bytes, MD5: {md5[:8]}...")
        return True, md5

    except urllib.error.HTTPError as e:
        print(f"    HTTP Error {e.code}: {e.reason}")
        return False, None
    except urllib.error.URLError as e:
        print(f"    URL Error: {e.reason}")
        return False, None
    except Exception as e:
        print(f"    Error: {e}")
        return False, None


def try_download_calspec():
    """Try to download CALSPEC files from various sources."""
    print("\n" + "=" * 60)
    print("Attempting CALSPEC downloads")
    print("=" * 60)

    downloaded = []

    for filename in CALSPEC_FILES:
        dest = HLSP_DIR / "calspec" / filename

        if dest.exists():
            print(f"  Already exists: {filename}")
            downloaded.append(str(dest))
            continue

        # Try primary URL
        url = CALSPEC_BASE_URL + filename
        success, md5 = download_file(url, dest)

        if not success:
            # Try alternative URLs
            for alt_base in CALSPEC_ALT_URLS:
                url = alt_base + filename
                success, md5 = download_file(url, dest)
                if success:
                    break

        if success:
            downloaded.append(str(dest))

    return downloaded


def try_download_wdlinelist():
    """Try to download WD-LINELIST HLSP files."""
    print("\n" + "=" * 60)
    print("Attempting WD-LINELIST downloads")
    print("=" * 60)

    downloaded = []

    for filepath in WD_LINELIST_FILES:
        filename = os.path.basename(filepath)
        dest = HLSP_DIR / "wd-linelist" / filename

        if dest.exists():
            print(f"  Already exists: {filename}")
            downloaded.append(str(dest))
            continue

        url = WD_LINELIST_BASE + filepath
        success, md5 = download_file(url, dest)

        if success:
            downloaded.append(str(dest))

    return downloaded


def try_mast_portal_download():
    """Try to find and download via MAST portal API."""
    print("\n" + "=" * 60)
    print("Attempting MAST Portal API search")
    print("=" * 60)

    try:
        from astroquery.mast import Mast
        import requests

        # Search for HLSP via the Portal API
        service = "Mast.Caom.Filtered"
        params = {
            "columns": "*",
            "filters": [
                {"paramName": "obs_collection", "values": ["HLSP"]},
                {"paramName": "target_name", "values": ["G191-B2B", "G191B2B", "g191b2b"]},
            ]
        }

        # This is a more direct API call
        api_url = "https://mast.stsci.edu/api/v0.1/invoke"
        headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain"}

        data = {
            "request": json.dumps({
                "service": service,
                "params": params,
                "format": "json"
            })
        }

        response = requests.post(api_url, data=data, headers=headers, timeout=30)

        if response.status_code == 200:
            result = response.json()
            if result.get('data'):
                print(f"  Found {len(result['data'])} HLSP observations")
                for row in result['data'][:5]:
                    print(f"    {row.get('obs_id', 'N/A')}: {row.get('dataproduct_type', 'N/A')}")
                return True
        else:
            print(f"  API returned status {response.status_code}")

    except ImportError:
        print("  requests not available for MAST API call")
    except Exception as e:
        print(f"  MAST API error: {e}")

    return False


def download_stis_raw_sample():
    """Download a sample of raw STIS x1d products to demonstrate the pipeline."""
    print("\n" + "=" * 60)
    print("Attempting raw STIS x1d sample download")
    print("=" * 60)

    try:
        from astroquery.mast import Observations
        import pandas as pd

        # Load our observation manifest
        manifest_path = Path("/home/primary/LHCP/manifests/manifest_obs.csv")
        if not manifest_path.exists():
            print("  Observation manifest not found")
            return []

        obs_df = pd.read_csv(manifest_path)

        # Get first 3 E140H observations as a sample
        e140h_obs = obs_df[obs_df['instrument'].str.contains('FUV', na=False)].head(3)

        if len(e140h_obs) == 0:
            print("  No E140H observations found in manifest")
            return []

        downloaded = []
        raw_dir = Path("/home/primary/LHCP/data/raw")
        raw_dir.mkdir(parents=True, exist_ok=True)

        for _, row in e140h_obs.iterrows():
            obsid = str(row['obsid'])
            print(f"\n  Processing obsid: {obsid}")

            try:
                # Get product list
                products = Observations.get_product_list(obsid)

                if len(products) == 0:
                    print(f"    No products found")
                    continue

                # Filter for x1d files
                x1d_products = Observations.filter_products(
                    products,
                    productSubGroupDescription=["X1D"]
                )

                if len(x1d_products) == 0:
                    # Try science products
                    x1d_products = Observations.filter_products(
                        products,
                        productType=["SCIENCE"],
                        extension=["fits"]
                    )

                if len(x1d_products) > 0:
                    print(f"    Found {len(x1d_products)} products to download")
                    manifest = Observations.download_products(
                        x1d_products[:5],  # Limit to 5 per observation
                        download_dir=str(raw_dir),
                        cache=True
                    )

                    if manifest is not None:
                        for row in manifest:
                            if row['Status'] == 'COMPLETE':
                                downloaded.append(row['Local Path'])
                                print(f"    Downloaded: {os.path.basename(row['Local Path'])}")

            except Exception as e:
                print(f"    Error: {e}")

        return downloaded

    except ImportError:
        print("  astroquery or pandas not available")
        return []
    except Exception as e:
        print(f"  Error: {e}")
        return []


def main():
    print("=" * 60)
    print("HLSP and Raw Data Download for G191-B2B")
    print("=" * 60)

    all_downloaded = []

    # Try CALSPEC
    calspec = try_download_calspec()
    all_downloaded.extend(calspec)

    # Try WD-LINELIST
    wdlinelist = try_download_wdlinelist()
    all_downloaded.extend(wdlinelist)

    # Try MAST Portal API
    try_mast_portal_download()

    # Try raw STIS sample
    stis_raw = download_stis_raw_sample()
    all_downloaded.extend(stis_raw)

    # Summary
    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)
    print(f"Total files downloaded: {len(all_downloaded)}")

    if all_downloaded:
        print("\nDownloaded files:")
        for f in all_downloaded:
            print(f"  {f}")

    # Check what we have
    print("\n" + "=" * 60)
    print("CURRENT DATA DIRECTORY STATUS")
    print("=" * 60)

    for subdir in ['hlsp/calspec', 'hlsp/wd-linelist', 'raw']:
        path = Path("/home/primary/LHCP/data") / subdir
        if path.exists():
            files = list(path.rglob("*"))
            data_files = [f for f in files if f.is_file()]
            print(f"  {subdir}/: {len(data_files)} files")
            for f in data_files[:3]:
                print(f"    {f.name}")
            if len(data_files) > 3:
                print(f"    ... and {len(data_files) - 3} more")
        else:
            print(f"  {subdir}/: (not created)")


if __name__ == "__main__":
    main()
