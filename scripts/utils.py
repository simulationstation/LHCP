#!/usr/bin/env python3
"""
Utility functions for LHCP data collection and processing.
"""

import hashlib
import os
from datetime import datetime

# Physical constants for wavelength conversion
# IAU standard: https://www.iau.org/static/resolutions/IAU2012_English.pdf
C_VACUUM = 299792458.0  # m/s (exact)

# Morton (1991) formula coefficients for air-to-vacuum conversion
# Valid for λ > 2000 Å in air
def air_to_vacuum(wavelength_air_ang):
    """
    Convert air wavelength to vacuum wavelength.
    Uses the IAU standard formula (Morton 1991, ApJS 77, 119).

    Parameters:
    -----------
    wavelength_air_ang : float
        Wavelength in air in Angstroms

    Returns:
    --------
    float : Wavelength in vacuum in Angstroms

    Note: Below 2000 Å, NIST provides vacuum wavelengths directly.
    This function is for wavelengths > 2000 Å.
    """
    # Edlen (1966) formula, as parameterized by Morton
    sigma2 = (1e4 / wavelength_air_ang) ** 2
    n_minus_1 = 6.4328e-5 + 2.94981e-2 / (146 - sigma2) + 2.5540e-4 / (41 - sigma2)
    return wavelength_air_ang * (1 + n_minus_1)

def vacuum_to_air(wavelength_vac_ang):
    """
    Convert vacuum wavelength to air wavelength.
    Uses inverted Edlen formula.

    Parameters:
    -----------
    wavelength_vac_ang : float
        Wavelength in vacuum in Angstroms

    Returns:
    --------
    float : Wavelength in air in Angstroms
    """
    sigma2 = (1e4 / wavelength_vac_ang) ** 2
    n_minus_1 = 6.4328e-5 + 2.94981e-2 / (146 - sigma2) + 2.5540e-4 / (41 - sigma2)
    return wavelength_vac_ang / (1 + n_minus_1)

def compute_md5(filepath):
    """Compute MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def compute_sha256(filepath):
    """Compute SHA256 checksum of a file."""
    hash_sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()

def get_file_info(filepath):
    """Get file size and modification time."""
    stat = os.stat(filepath)
    return {
        'size_bytes': stat.st_size,
        'mtime': datetime.fromtimestamp(stat.st_mtime).isoformat()
    }

def add_to_manifest(manifest_path, entry_dict):
    """Append an entry to a CSV manifest file."""
    import csv

    # Read header to get field order
    with open(manifest_path, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

    # Append new row
    with open(manifest_path, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writerow(entry_dict)

def timestamp_now():
    """Return current UTC timestamp in ISO format."""
    return datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

# Gravitational redshift calculation
G = 6.67430e-11  # m^3 kg^-1 s^-2 (CODATA 2018)
C = 299792458.0  # m/s (exact)
M_SUN = 1.98841e30  # kg (IAU nominal solar mass)
R_SUN = 6.957e8  # m (IAU nominal solar radius)

def gravitational_redshift(mass_msun, radius_rsun):
    """
    Compute gravitational redshift z_grav = GM/(Rc^2) for a star.

    Parameters:
    -----------
    mass_msun : float
        Stellar mass in solar masses
    radius_rsun : float
        Stellar radius in solar radii

    Returns:
    --------
    float : Gravitational redshift z_grav (dimensionless)
    """
    M = mass_msun * M_SUN
    R = radius_rsun * R_SUN
    return G * M / (R * C**2)

def gravitational_potential(mass_msun, radius_rsun):
    """
    Compute gravitational potential φ = GM/R (in m^2/s^2).
    """
    M = mass_msun * M_SUN
    R = radius_rsun * R_SUN
    return G * M / R


if __name__ == "__main__":
    # Test wavelength conversion
    print("Testing air-to-vacuum conversion:")
    print(f"  5000 Å (air) -> {air_to_vacuum(5000):.4f} Å (vacuum)")
    print(f"  Expected ~5001.39 Å")

    # Test gravitational redshift for typical WD
    print("\nTesting gravitational redshift:")
    print(f"  M=0.6 M_sun, R=0.012 R_sun:")
    z = gravitational_redshift(0.6, 0.012)
    print(f"  z_grav = {z:.6e}")
    print(f"  v_grav = {z * C / 1000:.2f} km/s")
