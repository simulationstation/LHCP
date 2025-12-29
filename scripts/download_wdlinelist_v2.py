#!/usr/bin/env python3
"""
Download WD-LINELIST HLSP products using direct URL patterns and MAST Portal API.
"""

import os
import json
import hashlib
import urllib.request
import urllib.error
from datetime import datetime
from pathlib import Path

HLSP_DIR = Path("/home/primary/LHCP/data/hlsp/wd-linelist")
HLSP_DIR.mkdir(parents=True, exist_ok=True)

# Known WD-LINELIST URL patterns from MAST HLSP
# Format: https://archive.stsci.edu/hlsp/wd-linelist/<target>/<filename>
WD_LINELIST_BASE = "https://archive.stsci.edu/hlsp/wd-linelist/"

# Known file patterns for G191-B2B
KNOWN_FILES = [
    # Coadded spectra
    "g191-b2b/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_cspec.fits",
    "g191-b2b/hlsp_wd-linelist_hst_stis_g191-b2b_e230h_v1_cspec.fits",
    "g191-b2b/hlsp_wd-linelist_multi_multi_g191-b2b_multi_v1_sed.fits",
    "g191-b2b/hlsp_wd-linelist_fuse_lwrs_g191-b2b_fuse_v1_cspec.fits",
    # Line lists
    "g191-b2b/hlsp_wd-linelist_multi_multi_g191-b2b_multi_v1_linelist.txt",
    "g191-b2b/hlsp_wd-linelist_multi_multi_g191-b2b_multi_v1_linelist.ecsv",
    # Alternative target naming
    "g191b2b/hlsp_wd-linelist_hst_stis_g191b2b_e140h_v1_cspec.fits",
    "g191b2b/hlsp_wd-linelist_hst_stis_g191b2b_e230h_v1_cspec.fits",
]


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


def download_file(url, dest_path):
    """Download a file."""
    try:
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0 (LHCP Pipeline)')

        with urllib.request.urlopen(req, timeout=60) as response:
            content = response.read()

        # Check if we got HTML instead of FITS
        if content[:100].startswith(b'<!DOCTYPE') or content[:100].startswith(b'<html'):
            print(f"    Got HTML instead of data file")
            return False, None

        dest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(dest_path, 'wb') as f:
            f.write(content)

        md5 = hashlib.md5(content).hexdigest()
        return True, md5

    except urllib.error.HTTPError as e:
        if e.code != 404:
            print(f"    HTTP Error {e.code}: {e.reason}")
        return False, None
    except urllib.error.URLError as e:
        print(f"    URL Error: {e.reason}")
        return False, None
    except Exception as e:
        print(f"    Error: {e}")
        return False, None


def try_mast_portal_api():
    """Try MAST Portal API for HLSP products."""
    print("\nTrying MAST Portal API...")

    try:
        import requests

        # MAST Portal API
        api_url = "https://mast.stsci.edu/api/v0.1/invoke"

        # Search for HLSP observations
        payload = {
            "service": "Mast.Hlsp.Search",
            "params": {"search": "wd-linelist"},
            "format": "json"
        }

        headers = {
            "Content-type": "application/x-www-form-urlencoded",
            "Accept": "application/json"
        }

        response = requests.post(
            api_url,
            data={"request": json.dumps(payload)},
            headers=headers,
            timeout=30
        )

        if response.status_code == 200:
            result = response.json()
            if 'data' in result and len(result['data']) > 0:
                print(f"  Found {len(result['data'])} WD-LINELIST entries")
                for item in result['data'][:10]:
                    print(f"    {item.get('target', 'N/A')}: {item.get('mission', 'N/A')}")
                return result['data']
            else:
                print("  No data in response")
        else:
            print(f"  API error: {response.status_code}")

    except ImportError:
        print("  requests not installed")
    except Exception as e:
        print(f"  Error: {e}")

    return []


def try_caom_search():
    """Try CAOM API search for HLSP."""
    print("\nTrying CAOM API search...")

    try:
        import requests

        api_url = "https://mast.stsci.edu/api/v0.1/invoke"

        # Search CAOM for HLSP
        payload = {
            "service": "Mast.Caom.Filtered",
            "params": {
                "columns": "obsid,obs_id,target_name,obs_collection,dataproduct_type,calib_level,dataURL",
                "filters": [
                    {"paramName": "obs_collection", "values": ["HLSP"]},
                    {"paramName": "project", "values": ["WD-LINELIST", "wd-linelist"]}
                ]
            },
            "format": "json"
        }

        headers = {"Content-type": "application/x-www-form-urlencoded"}

        response = requests.post(
            api_url,
            data={"request": json.dumps(payload)},
            headers=headers,
            timeout=30
        )

        if response.status_code == 200:
            result = response.json()
            if 'data' in result:
                print(f"  Found {len(result['data'])} HLSP observations")
                # Filter for G191-B2B
                g191 = [d for d in result['data']
                        if 'G191' in str(d.get('target_name', '')).upper()]
                print(f"  G191-B2B matches: {len(g191)}")
                for item in g191[:5]:
                    print(f"    {item.get('obs_id', 'N/A')}: {item.get('dataURL', 'N/A')}")
                return g191

    except ImportError:
        print("  requests not installed")
    except Exception as e:
        print(f"  Error: {e}")

    return []


def main():
    print("=" * 60)
    print("WD-LINELIST HLSP Download (v2)")
    print("=" * 60)

    downloaded = []
    failed = []

    # Method 1: Try direct URLs
    print("\nMethod 1: Direct URL downloads...")
    for filepath in KNOWN_FILES:
        filename = os.path.basename(filepath)
        dest = HLSP_DIR / filename
        url = WD_LINELIST_BASE + filepath

        print(f"  Trying: {filename}")
        success, md5 = download_file(url, dest)

        if success:
            # Verify it's valid
            if filename.endswith('.fits') and not is_valid_fits(dest):
                print(f"    Invalid FITS file, removing")
                dest.unlink()
                failed.append(filename)
            else:
                downloaded.append({
                    'filename': filename,
                    'path': str(dest),
                    'url': url,
                    'md5': md5,
                    'size': dest.stat().st_size
                })
                print(f"    SUCCESS: {dest.stat().st_size} bytes")
        else:
            failed.append(filename)

    # Method 2: Try MAST Portal API
    portal_results = try_mast_portal_api()

    # Method 3: Try CAOM search
    caom_results = try_caom_search()

    # Summary
    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)

    if downloaded:
        print(f"\nSuccessfully downloaded {len(downloaded)} files:")
        for d in downloaded:
            print(f"  {d['filename']}: {d['size']} bytes, MD5: {d['md5'][:8]}...")
    else:
        print("\nNo files downloaded via direct URLs.")
        print("The WD-LINELIST HLSP may have moved or require different access.")
        print("\nAlternative access methods:")
        print("  1. MAST Portal: https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html")
        print("  2. Search for: Collection=HLSP, Project=WD-LINELIST")
        print("  3. Or use CALSPEC spectra as alternative (already downloaded)")

    # Write provenance
    provenance = {
        "retrieval_time": datetime.now().isoformat(),
        "method": "direct_url_download",
        "base_url": WD_LINELIST_BASE,
        "files_attempted": KNOWN_FILES,
        "files_downloaded": downloaded,
        "files_failed": failed,
        "portal_results_count": len(portal_results),
        "caom_results_count": len(caom_results)
    }

    prov_file = HLSP_DIR.parent / "download_provenance.json"
    with open(prov_file, 'w') as f:
        json.dump(provenance, f, indent=2)
    print(f"\nProvenance written to: {prov_file}")

    # Check current HLSP status
    print("\n" + "=" * 60)
    print("CURRENT HLSP DIRECTORY STATUS")
    print("=" * 60)

    for subdir in ['calspec', 'wd-linelist']:
        path = HLSP_DIR.parent / subdir
        if path.exists():
            fits_files = list(path.glob("*.fits"))
            other_files = [f for f in path.glob("*") if f.is_file() and not f.suffix == '.fits']
            valid = [f for f in fits_files if is_valid_fits(f)]
            print(f"\n{subdir}/:")
            print(f"  FITS files: {len(fits_files)} ({len(valid)} valid)")
            print(f"  Other files: {len(other_files)}")
            for f in valid:
                print(f"    {f.name}: {f.stat().st_size} bytes")


if __name__ == "__main__":
    main()
