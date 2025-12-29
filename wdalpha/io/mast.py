from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

from astroquery.mast import Observations
from rich.console import Console
from rich.progress import track

console = Console()


def _download_product(row, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    filename = Path(row["productFilename"])
    dest = out_dir / filename.name
    if dest.exists():
        return dest
    manifest = Observations.download_products([row], download_dir=str(out_dir), cache=True)
    local_path = Path(manifest[0]["Local Path"]).expanduser()
    if local_path.exists() and local_path != dest:
        local_path.rename(dest)
    return dest


def download_hlsp(target: str, out_dir: Path, products: Iterable[str] | None = None) -> List[Path]:
    """Download HLSP coadds for a target via MAST.

    Parameters
    ----------
    target : str
        Target name (e.g., G191-B2B).
    out_dir : Path
        Output directory for downloaded files.
    products : Iterable[str] | None
        Optional substring filters for productFilename.
    """
    console.log(f"Querying MAST HLSP for target {target}")
    obs = Observations.query_object(target, dataproduct_type="spectrum")
    hlsp = obs[obs["obs_collection"] == "HLSP"]
    products_table = Observations.get_product_list(hlsp)
    if products:
        mask = [any(p in name for p in products) for name in products_table["productFilename"]]
        products_table = products_table[mask]
    products_table = Observations.filter_products(products_table, productType="SCIENCE")

    downloaded: List[Path] = []
    for row in track(products_table, description="Downloading HLSP products"):
        downloaded.append(_download_product(row, out_dir))
    return downloaded


def download_mast(target: str, out_dir: Path, instrument: str | None = None) -> List[Path]:
    console.log(f"Querying MAST for target {target}")
    obs = Observations.query_object(target, dataproduct_type="spectrum")
    if instrument:
        obs = obs[obs["instrument_name"] == instrument]
    products_table = Observations.get_product_list(obs)
    products_table = Observations.filter_products(products_table, productType="SCIENCE")

    downloaded: List[Path] = []
    for row in track(products_table, description="Downloading MAST products"):
        downloaded.append(_download_product(row, out_dir))
    return downloaded
