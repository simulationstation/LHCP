#!/usr/bin/env python3
"""Run a single stress test configuration."""
import sys
import yaml
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from wdalpha.io.fits import read_spectrum
from wdalpha.lines.masks import build_line_windows
from wdalpha.preprocess.continuum import build_mask, normalize_spectrum
from wdalpha.lines.fit import fit_lines
from wdalpha.inference.mm_regression import infer_delta_alpha, save_result
import pandas as pd
import numpy as np

def run_analysis(run_dir, spline_s=0.001, mask_width=0.4, window=0.6, max_components=2,
                 robust=False, sigma_clip=None, species_label="Fe V + Ni V"):
    """Run analysis with given configuration."""
    E140H_FITS = 'data/hlsp/wd-linelist/hlsp_wd-linelist_hst_stis_g191-b2b_e140h_v1_coadd-spec.fits'
    RUN_DIR = Path(run_dir)
    np.random.seed(42)

    print(f'=== {RUN_DIR.name} Analysis ===')
    print(f'Config: spline_s={spline_s}, mask_width={mask_width}, robust={robust}, sigma_clip={sigma_clip}')

    print('Loading E140H spectrum...')
    spectrum = read_spectrum(Path(E140H_FITS))

    atomic = pd.read_csv(RUN_DIR / 'atomic_used.csv')
    print(f'  Atomic lines: {len(atomic)}')

    windows = build_line_windows(atomic, half_width=window)
    mask = build_mask(spectrum.wavelength, atomic['wavelength_aa'], half_width=mask_width)

    print('Normalizing spectrum...')
    norm_flux, norm_error, continuum = normalize_spectrum(
        spectrum.wavelength, spectrum.flux, spectrum.error, mask, spline_s=spline_s
    )

    print('Fitting lines...')
    results = fit_lines(spectrum.wavelength, norm_flux, norm_error, windows, max_components=max_components)
    print(f'  Lines fit: {len(results)} / {len(atomic)}')

    rows = [{'line_id': res.line_id, 'species': res.species, 'lambda_obs': res.lambda_obs,
             'sigma_lambda_obs': res.sigma_lambda_obs, 'chi2': res.chi2, 'n_components': res.n_components}
            for res in results]
    lines_df = pd.DataFrame(rows)

    # Apply sigma clipping if requested
    if sigma_clip is not None and len(lines_df) > 0:
        # Merge on line_id to get one-to-one match
        merged = lines_df.merge(atomic, on='line_id', suffixes=('_obs', '_lab'))
        if len(merged) > 0 and 'wavelength_aa' in merged.columns:
            residuals = merged['lambda_obs'] - merged['wavelength_aa']
            med = np.median(residuals)
            mad = np.median(np.abs(residuals - med))
            sigma = max(mad * 1.4826, 1e-4)  # Convert MAD to sigma, with floor
            good_mask = np.abs(residuals - med) < sigma_clip * sigma
            n_keep = good_mask.sum()
            print(f'  Sigma clipping: keeping {n_keep} / {len(merged)} lines')
            # Keep only the lines that pass clipping
            good_line_ids = merged.loc[good_mask, 'line_id'].values
            lines_df = lines_df[lines_df['line_id'].isin(good_line_ids)].reset_index(drop=True)

    lines_df.to_csv(RUN_DIR / 'line_measurements.csv', index=False)

    print('Inferring delta alpha/alpha...')
    try:
        result = infer_delta_alpha(lines_df, atomic, include_lab_uncertainty=True)
        print(f'  delta_alpha/alpha = {result.delta_alpha:.6e} +/- {result.delta_alpha_err:.6e}')
        print(f'  chi2/dof = {result.chi2:.1f} / {result.dof} = {result.chi2/result.dof:.2f}')

        save_result(result, RUN_DIR / 'inferred_alpha.json', metadata={
            'species': species_label, 'grating': 'E140H', 'n_lines': len(lines_df),
            'spline_s': spline_s, 'mask_width': mask_width, 'robust': robust, 'sigma_clip': sigma_clip
        })
    except Exception as e:
        print(f'  Inference failed: {e}')
        with open(RUN_DIR / 'inferred_alpha.json', 'w') as f:
            json.dump({'error': str(e)}, f)

    config = {
        'seed': 42, 'spectrum': E140H_FITS, 'species': species_label,
        'continuum': {'spline_s': spline_s, 'mask_width': mask_width},
        'linefit': {'window': window, 'max_components': max_components},
        'robust': robust, 'sigma_clip': sigma_clip
    }
    with open(RUN_DIR / 'config_used.yaml', 'w') as f:
        yaml.dump(config, f)
    print('Complete!')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('run_dir')
    parser.add_argument('--spline-s', type=float, default=0.001)
    parser.add_argument('--mask-width', type=float, default=0.4)
    parser.add_argument('--window', type=float, default=0.6)
    parser.add_argument('--max-components', type=int, default=2)
    parser.add_argument('--robust', action='store_true')
    parser.add_argument('--sigma-clip', type=float, default=None)
    parser.add_argument('--species', default='Fe V + Ni V')
    args = parser.parse_args()

    run_analysis(args.run_dir, args.spline_s, args.mask_width, args.window,
                 args.max_components, args.robust, args.sigma_clip, args.species)
