from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from astropy.io import fits
from astropy import units as u


@dataclass
class Spectrum1D:
    wavelength: np.ndarray
    flux: np.ndarray
    error: np.ndarray
    meta: Dict[str, str]


def _find_spectrum_columns(hdu: fits.BinTableHDU) -> Tuple[str, str, str]:
    columns = {name.upper(): name for name in hdu.columns.names}
    wave_key = None
    flux_key = None
    err_key = None
    for candidate in ("WAVELENGTH", "LAMBDA", "WAVE"):
        if candidate in columns:
            wave_key = columns[candidate]
            break
    for candidate in ("FLUX", "NET", "COUNT_RATE"):
        if candidate in columns:
            flux_key = columns[candidate]
            break
    for candidate in ("ERROR", "ERR", "FLUXERR", "STATERROR"):
        if candidate in columns:
            err_key = columns[candidate]
            break
    if not (wave_key and flux_key and err_key):
        raise ValueError("Spectrum columns not found in FITS table.")
    return wave_key, flux_key, err_key


def read_spectrum(path: Path) -> Spectrum1D:
    """Read a 1D spectrum from STIS/FUSE FITS files.

    Returns wavelength in Angstrom.
    """
    with fits.open(path) as hdul:
        primary_header = hdul[0].header if hdul else {}
        table_hdu = None
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                try:
                    _find_spectrum_columns(hdu)
                    table_hdu = hdu
                    break
                except ValueError:
                    continue
        if table_hdu is None:
            raise ValueError(f"No spectral table found in {path}")
        wave_key, flux_key, err_key = _find_spectrum_columns(table_hdu)
        wave = np.array(table_hdu.data[wave_key], dtype=float)
        flux = np.array(table_hdu.data[flux_key], dtype=float)
        error = np.array(table_hdu.data[err_key], dtype=float)
        wave_unit = table_hdu.columns[wave_key].unit or "Angstrom"

    wavelength = (wave * u.Unit(wave_unit)).to(u.AA).value
    meta = {
        "path": str(path),
        "wave_unit": wave_unit,
        "obsid": str(primary_header.get("OBSID", "")),
        "instrument": str(primary_header.get("INSTRUME", "")),
        "telescope": str(primary_header.get("TELESCOP", "")),
        "wcs_vacuum": str(primary_header.get("VACUUM", "")),
    }
    return Spectrum1D(wavelength=wavelength, flux=flux, error=error, meta=meta)


def check_spectrum_sanity(spectrum: Spectrum1D) -> Dict[str, str]:
    errors = []
    warnings = []
    wave = spectrum.wavelength
    flux = spectrum.flux
    error = spectrum.error
    if wave.ndim != 1 or flux.ndim != 1 or error.ndim != 1:
        errors.append("Spectrum arrays must be 1D.")
    if wave.size == 0:
        errors.append("Spectrum has zero length.")
    if not np.all(np.isfinite(wave)):
        errors.append("Non-finite wavelength values detected.")
    if not np.all(np.isfinite(flux)):
        errors.append("Non-finite flux values detected.")
    if not np.all(np.isfinite(error)):
        errors.append("Non-finite error values detected.")
    if wave.size > 1 and not np.all(np.diff(wave) > 0):
        errors.append("Wavelength array is not strictly increasing.")
    if np.any(error <= 0):
        errors.append("Non-positive error values detected.")
    if wave.size > 0:
        coverage = f"{wave.min():.2f}-{wave.max():.2f} AA"
    else:
        coverage = "unknown"
    if spectrum.meta.get("wcs_vacuum", "").strip() == "":
        warnings.append("VACUUM keyword not set in header; verify wavelength convention.")
    return {
        "errors": "; ".join(errors) if errors else "",
        "warnings": "; ".join(warnings) if warnings else "",
        "coverage": coverage,
        "n_pixels": str(wave.size),
    }
