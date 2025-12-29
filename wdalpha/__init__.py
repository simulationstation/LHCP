"""wdalpha: Measure fractional fine-structure constant variation from WD spectra."""
from importlib.metadata import version

__all__ = ["__version__"]

try:
    __version__ = version("wdalpha")
except Exception:  # pragma: no cover - package metadata not available in dev
    __version__ = "0.1.0"
