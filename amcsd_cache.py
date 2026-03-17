from __future__ import annotations

import json
import urllib.request
import zipfile
from pathlib import Path

from cif_matcher import parse_cif_file


AMCSD_CIF_ZIP_URL = "https://www.rruff.net/AMS/zipped_files/cif.zip"


def default_cache_dir() -> Path:
    return Path(__file__).resolve().parent / ".cache" / "amcsd"


def download_amcsd_zip(zip_path: Path, url: str = AMCSD_CIF_ZIP_URL) -> None:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, zip_path.open("wb") as destination:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            destination.write(chunk)


def extract_amcsd_zip(zip_path: Path, extract_dir: Path) -> None:
    extract_dir.mkdir(parents=True, exist_ok=True)
    marker = extract_dir / ".complete"
    if marker.exists():
        return
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(extract_dir)
    marker.write_text("ok\n")


def build_amcsd_index(cif_dir: Path, index_path: Path) -> Path:
    index_path.parent.mkdir(parents=True, exist_ok=True)
    entries = []
    for cif_path in sorted(cif_dir.rglob("*.cif")):
        try:
            structure = parse_cif_file(cif_path)
        except Exception:
            continue
        entries.append(
            {
                "relative_path": str(cif_path.relative_to(cif_dir)),
                "name": structure.name,
                "formula": structure.formula,
                "space_group": structure.space_group,
                "cell_a": structure.cell_a,
                "is_cubic": structure.is_cubic,
                "element_count": len({site.element for site in structure.atom_sites}),
            }
        )

    index_path.write_text(json.dumps(entries))
    return index_path


def ensure_amcsd_library(cache_dir: Path | None = None, force_refresh: bool = False) -> tuple[Path, Path]:
    cache_root = cache_dir or default_cache_dir()
    zip_path = cache_root / "amcsd_cif.zip"
    cif_dir = cache_root / "cif"
    index_path = cache_root / "amcsd_index.json"

    if force_refresh or not zip_path.exists():
        download_amcsd_zip(zip_path)
    extract_amcsd_zip(zip_path, cif_dir)
    if force_refresh or not index_path.exists():
        build_amcsd_index(cif_dir, index_path)

    return cif_dir, index_path


def load_amcsd_index(index_path: Path) -> list[dict[str, object]]:
    return json.loads(index_path.read_text())


def select_amcsd_candidate_files(
    cif_dir: Path,
    index_path: Path,
    observed_a_angstrom: float,
    max_elements: int = 3,
    max_relative_a_diff: float = 0.25,
    limit: int = 400,
) -> list[Path]:
    entries = load_amcsd_index(index_path)
    filtered = []
    for entry in entries:
        if not entry.get("is_cubic"):
            continue
        if int(entry.get("element_count", 99)) > max_elements:
            continue
        cell_a = float(entry["cell_a"])
        if abs(cell_a - observed_a_angstrom) / observed_a_angstrom > max_relative_a_diff:
            continue
        filtered.append((abs(cell_a - observed_a_angstrom), entry))

    filtered.sort(key=lambda item: item[0])
    if limit <= 0:
        return [cif_dir / str(entry["relative_path"]) for _, entry in filtered]
    return [cif_dir / str(entry["relative_path"]) for _, entry in filtered[:limit]]
