from __future__ import annotations

from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import units as u


def barycentric_correction(ra: float, dec: float, obs_time: str, location: EarthLocation) -> float:
    """Return barycentric velocity correction (m/s).

    Placeholder for future expansion; for now raises if used.
    """
    _ = (ra, dec, obs_time, location)
    raise NotImplementedError("Barycentric correction not implemented; absorbed into velocity offset.")
