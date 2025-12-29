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

    Returns wavelength in vacuum Angstrom.
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
    }
    return Spectrum1D(wavelength=wavelength, flux=flux, error=error, meta=meta)
