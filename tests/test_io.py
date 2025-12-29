import numpy as np
from astropy.io import fits
from wdalpha.io.fits import check_spectrum_sanity, read_spectrum


def test_read_spectrum(tmp_path):
    wave = np.linspace(1000, 1100, 10)
    flux = np.ones_like(wave)
    err = np.full_like(wave, 0.1)
    col_wave = fits.Column(name="WAVELENGTH", array=wave, format="D", unit="Angstrom")
    col_flux = fits.Column(name="FLUX", array=flux, format="D")
    col_err = fits.Column(name="ERROR", array=err, format="D")
    hdu = fits.BinTableHDU.from_columns([col_wave, col_flux, col_err])
    path = tmp_path / "spec.fits"
    hdu.writeto(path)

    spec = read_spectrum(path)
    assert spec.wavelength.size == 10
    assert np.allclose(spec.flux, flux)
    sanity = check_spectrum_sanity(spec)
    assert sanity["errors"] == ""
