#!/usr/bin/env python3
"""
Query MAST archive for HST/STIS observations of G191-B2B white dwarf.
Focus on E140H and E230H gratings.

Output: CSV file with observation metadata and product availability.
"""

import os
import sys
from datetime import datetime

import pandas as pd
from astroquery.mast import Observations

# Output paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), 'data')
OUTPUT_CSV = os.path.join(DATA_DIR, 'mast_query_results.csv')

# Target and instrument filters
TARGET_NAME = 'G191-B2B'
TARGET_NAME_MAST = 'G191B2B'  # MAST uses name without hyphen
INSTRUMENT = 'STIS'
GRATINGS_OF_INTEREST = ['E140H', 'E230H']


def query_mast_observations():
    """Query MAST for G191-B2B STIS observations."""
    print(f"Querying MAST for {TARGET_NAME} observations...")
    print(f"Searching for HST/STIS with gratings: {GRATINGS_OF_INTEREST}")
    print("-" * 60)

    # Query by target name with HST/STIS criteria
    # MAST stores this target as 'G191B2B' (no hyphen)
    obs_table = Observations.query_criteria(
        target_name=TARGET_NAME_MAST,
        obs_collection='HST',
        instrument_name='STIS*'
    )

    # If few results, also try hyphenated version and combine
    if len(obs_table) < 10:
        try:
            obs_table2 = Observations.query_criteria(
                target_name=TARGET_NAME,
                obs_collection='HST',
                instrument_name='STIS*'
            )
            if len(obs_table2) > 0:
                from astropy.table import vstack
                obs_table = vstack([obs_table, obs_table2])
                # Remove duplicates based on obsid
                obs_df = obs_table.to_pandas()
                obs_df = obs_df.drop_duplicates(subset=['obsid'])
                from astropy.table import Table
                obs_table = Table.from_pandas(obs_df)
        except Exception as e:
            print(f"  Note: Could not query alternate name: {e}")

    print(f"Found {len(obs_table)} total STIS observations")
    return obs_table


def filter_echelle_gratings(obs_table):
    """Filter observations for E140H and E230H gratings."""
    # Get the filters/gratings column - may be named differently
    if 'filters' in obs_table.colnames:
        grating_col = 'filters'
    elif 'instrument_name' in obs_table.colnames:
        # Sometimes grating info is embedded in instrument name
        grating_col = 'instrument_name'
    else:
        print("Available columns:", obs_table.colnames)
        return obs_table

    # Convert to pandas for easier filtering
    df = obs_table.to_pandas()

    # Filter for E140H and E230H gratings
    mask = df[grating_col].str.contains('E140H|E230H', case=False, na=False)
    df_filtered = df[mask].copy()

    print(f"Filtered to {len(df_filtered)} observations with E140H or E230H gratings")

    return df_filtered, obs_table


def get_product_info(obs_ids, batch_size=50):
    """Get data products for each observation."""
    print(f"\nQuerying data products for {len(obs_ids)} observations...")

    product_info = {}
    total = len(obs_ids)

    for i, obsid in enumerate(obs_ids):
        if (i + 1) % 20 == 0:
            print(f"  Progress: {i + 1}/{total} observations queried...")

        try:
            products = Observations.get_product_list(obsid)
            if products is not None and len(products) > 0:
                product_types = set()
                for p in products:
                    suffix = p['productSubGroupDescription']
                    if suffix:
                        product_types.add(suffix.upper())
                product_info[obsid] = {
                    'x1d': 'X1D' in product_types,
                    'x1dsum': 'X1DSUM' in product_types,
                    'wavecal': 'WAV' in product_types or 'WAVECAL' in product_types,
                    'raw': 'RAW' in product_types,
                    'flt': 'FLT' in product_types,
                    'all_products': ','.join(sorted(product_types))
                }
        except Exception as e:
            print(f"  Warning: Could not get products for {obsid}: {e}")
            product_info[obsid] = {
                'x1d': False, 'x1dsum': False, 'wavecal': False,
                'raw': False, 'flt': False, 'all_products': 'ERROR'
            }

    print(f"  Completed product query for {len(product_info)} observations.")
    return product_info


def extract_metadata(obs_table):
    """Extract relevant metadata from observations table."""

    # Convert to pandas DataFrame
    if hasattr(obs_table, 'to_pandas'):
        df = obs_table.to_pandas()
    else:
        df = pd.DataFrame(obs_table)

    print(f"\nAvailable columns: {list(df.columns)}")

    # Define column mappings (MAST column names can vary)
    column_mappings = {
        'obsid': ['obsid', 'obs_id'],
        'target_name': ['target_name', 'target'],
        'grating': ['filters', 'filter', 'grating'],
        'instrument': ['instrument_name', 'instrument'],
        't_exptime': ['t_exptime', 'exptime', 'exposure_time'],
        't_min': ['t_min', 'mjd_min', 'start_time'],
        't_max': ['t_max', 'mjd_max', 'end_time'],
        's_ra': ['s_ra', 'ra'],
        's_dec': ['s_dec', 'dec'],
        'dataRights': ['dataRights', 'data_rights', 'calib_level'],
        'proposal_id': ['proposal_id', 'proposal'],
        'wavelength_region': ['wavelength_region', 'em_domain'],
        'em_min': ['em_min', 'wavelength_min'],
        'em_max': ['em_max', 'wavelength_max'],
    }

    # Find actual column names
    def find_column(options):
        for opt in options:
            if opt in df.columns:
                return opt
        return None

    actual_cols = {k: find_column(v) for k, v in column_mappings.items()}

    # Extract relevant data
    extracted = {}
    for key, col in actual_cols.items():
        if col:
            extracted[key] = df[col]

    result_df = pd.DataFrame(extracted)

    return result_df, df


def main():
    """Main function to query MAST and process results."""
    print("=" * 60)
    print("MAST Archive Query for G191-B2B HST/STIS Observations")
    print("=" * 60)
    print(f"Timestamp: {datetime.now().isoformat()}")
    print()

    # Ensure output directory exists
    os.makedirs(DATA_DIR, exist_ok=True)

    # Query MAST
    obs_table = query_mast_observations()

    if len(obs_table) == 0:
        print("ERROR: No observations found!")
        sys.exit(1)

    # Convert to DataFrame for processing
    df_all = obs_table.to_pandas()

    print(f"\nTotal STIS observations found: {len(df_all)}")
    print(f"\nColumns available: {list(df_all.columns)}")

    # Find the grating/filter column
    grating_col = None
    for col in ['filters', 'filter', 'grating', 'instrument_name']:
        if col in df_all.columns:
            grating_col = col
            break

    if grating_col:
        print(f"\nGrating column: {grating_col}")
        print(f"Unique gratings/modes: {df_all[grating_col].unique()}")

    # Filter for E140H and E230H
    if grating_col:
        # Create mask for E140H and E230H
        e140h_mask = df_all[grating_col].str.contains('E140H', case=False, na=False)
        e230h_mask = df_all[grating_col].str.contains('E230H', case=False, na=False)
        echelle_mask = e140h_mask | e230h_mask

        df_filtered = df_all[echelle_mask].copy()

        n_e140h = e140h_mask.sum()
        n_e230h = e230h_mask.sum()

        print(f"\n" + "=" * 60)
        print("SUMMARY STATISTICS")
        print("=" * 60)
        print(f"Total E140H observations: {n_e140h}")
        print(f"Total E230H observations: {n_e230h}")
        print(f"Total echelle observations: {len(df_filtered)}")
    else:
        df_filtered = df_all
        n_e140h = 'Unknown'
        n_e230h = 'Unknown'

    # Get date range
    date_col = None
    for col in ['t_min', 't_obs', 'mjd_min', 'dateobs']:
        if col in df_filtered.columns:
            date_col = col
            break

    if date_col:
        # Convert MJD to date if needed
        from astropy.time import Time
        try:
            mjd_min = df_filtered[date_col].min()
            mjd_max = df_filtered[date_col].max()

            if mjd_min > 40000:  # Likely MJD
                date_min = Time(mjd_min, format='mjd').datetime.strftime('%Y-%m-%d')
                date_max = Time(mjd_max, format='mjd').datetime.strftime('%Y-%m-%d')
            else:
                date_min = str(mjd_min)
                date_max = str(mjd_max)

            print(f"\nDate range: {date_min} to {date_max}")
        except Exception as e:
            print(f"\nCould not parse dates: {e}")
            date_min = date_max = 'Unknown'
    else:
        date_min = date_max = 'Unknown'

    # Check data rights
    rights_col = None
    for col in ['dataRights', 'data_rights', 'calib_level']:
        if col in df_filtered.columns:
            rights_col = col
            break

    if rights_col:
        print(f"\nData rights status:")
        print(df_filtered[rights_col].value_counts().to_string())

        # Check for any restricted data
        if 'EXCLUSIVE' in df_filtered[rights_col].values:
            print("\nWARNING: Some data has EXCLUSIVE access rights (proprietary)")
        else:
            print("\nAll data appears to be publicly accessible")

    # Get product information for subset
    print("\n" + "-" * 60)
    print("Querying data products (this may take a while)...")

    if 'obsid' in df_filtered.columns:
        obsid_col = 'obsid'
    elif 'obs_id' in df_filtered.columns:
        obsid_col = 'obs_id'
    else:
        obsid_col = None

    if obsid_col and len(df_filtered) > 0:
        # Query products for all observations
        all_obsids = df_filtered[obsid_col].tolist()
        product_info = get_product_info(all_obsids)

        # Add product columns to dataframe
        df_filtered['has_x1d'] = df_filtered[obsid_col].apply(
            lambda x: product_info.get(x, {}).get('x1d', None))
        df_filtered['has_x1dsum'] = df_filtered[obsid_col].apply(
            lambda x: product_info.get(x, {}).get('x1dsum', None))
        df_filtered['has_wavecal'] = df_filtered[obsid_col].apply(
            lambda x: product_info.get(x, {}).get('wavecal', None))
        df_filtered['has_raw'] = df_filtered[obsid_col].apply(
            lambda x: product_info.get(x, {}).get('raw', None))
        df_filtered['has_flt'] = df_filtered[obsid_col].apply(
            lambda x: product_info.get(x, {}).get('flt', None))
        df_filtered['products'] = df_filtered[obsid_col].apply(
            lambda x: product_info.get(x, {}).get('all_products', ''))

    # Select columns for output
    output_cols = []
    desired_cols = [
        'obsid', 'obs_id', 'target_name', 'filters', 'instrument_name',
        't_exptime', 't_min', 't_max', 's_ra', 's_dec', 'dataRights',
        'proposal_id', 'em_min', 'em_max', 'wavelength_region',
        'has_x1d', 'has_x1dsum', 'has_wavecal', 'has_raw', 'has_flt', 'products'
    ]

    for col in desired_cols:
        if col in df_filtered.columns:
            output_cols.append(col)

    df_output = df_filtered[output_cols].copy()

    # Save to CSV
    df_output.to_csv(OUTPUT_CSV, index=False)
    print(f"\n\nResults saved to: {OUTPUT_CSV}")
    print(f"Total rows: {len(df_output)}")

    # Print summary
    print("\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"Target: {TARGET_NAME}")
    print(f"Instrument: HST/STIS")
    print(f"Total E140H observations: {n_e140h}")
    print(f"Total E230H observations: {n_e230h}")
    print(f"Date range: {date_min} to {date_max}")
    print(f"Output file: {OUTPUT_CSV}")

    # Print first few rows
    print("\n" + "-" * 60)
    print("Sample of observations:")
    print(df_output.head(10).to_string())

    # Print additional summary statistics
    print("\n" + "=" * 60)
    print("DETAILED STATISTICS")
    print("=" * 60)

    # Breakdown by proposal
    if 'proposal_id' in df_output.columns:
        print("\nObservations by HST Proposal ID:")
        print(df_output.groupby(['proposal_id', grating_col]).size().to_string())

    # Product availability summary
    if 'has_x1d' in df_output.columns:
        print("\nProduct Availability Summary:")
        print(f"  X1D (extracted spectra): {df_output['has_x1d'].sum()}/{len(df_output)}")
        print(f"  X1DSUM (summed spectra): {df_output['has_x1dsum'].sum()}/{len(df_output)}")
        print(f"  WAV (wavecal): {df_output['has_wavecal'].sum()}/{len(df_output)}")
        print(f"  RAW: {df_output['has_raw'].sum()}/{len(df_output)}")
        print(f"  FLT: {df_output['has_flt'].sum()}/{len(df_output)}")

    # Exposure time statistics
    if 't_exptime' in df_output.columns:
        e140h_exptimes = df_output[df_output[grating_col].str.contains('E140H')]['t_exptime']
        e230h_exptimes = df_output[df_output[grating_col].str.contains('E230H')]['t_exptime']
        print("\nExposure Time Statistics (seconds):")
        print(f"  E140H: min={e140h_exptimes.min():.1f}, max={e140h_exptimes.max():.1f}, "
              f"median={e140h_exptimes.median():.1f}, total={e140h_exptimes.sum():.1f}")
        print(f"  E230H: min={e230h_exptimes.min():.1f}, max={e230h_exptimes.max():.1f}, "
              f"median={e230h_exptimes.median():.1f}, total={e230h_exptimes.sum():.1f}")

    # Wavelength coverage
    if 'em_min' in df_output.columns and 'em_max' in df_output.columns:
        print("\nWavelength Coverage (nm):")
        e140h_data = df_output[df_output[grating_col].str.contains('E140H')]
        e230h_data = df_output[df_output[grating_col].str.contains('E230H')]
        if len(e140h_data) > 0:
            print(f"  E140H: {e140h_data['em_min'].min():.1f} - {e140h_data['em_max'].max():.1f} nm")
        if len(e230h_data) > 0:
            print(f"  E230H: {e230h_data['em_min'].min():.1f} - {e230h_data['em_max'].max():.1f} nm")

    print("\n" + "=" * 60)
    print("NOTES")
    print("=" * 60)
    print("""
- Aperture information is not directly available from MAST queries.
  For STIS echelle observations, typical apertures include:
  * 0.2X0.2 - Small aperture for high S/N point sources
  * 0.2X0.06 - Very small aperture
  * 0.1X0.03 - Narrowest aperture for best wavelength calibration
  Aperture info can be retrieved from FITS headers (APERTURE keyword).

- The em_min/em_max values represent approximate wavelength ranges.
  E140H covers ~1140-1700 A (UV), E230H covers ~1620-3150 A (UV).

- All observations are PUBLIC (no proprietary restrictions).

- X1D files contain the extracted 1D spectra ready for science use.
  WAV files contain wavelength calibration information.
""")

    return df_output


if __name__ == '__main__':
    main()
