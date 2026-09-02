#!/usr/bin/env python3
"""Build the compact GSE129308 vignette fixture from the GEO raw archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as error:  # pragma: no cover - regeneration dependency
    raise SystemExit(
        "Install the regeneration dependencies with `python -m pip install h5py numpy`."
    ) from error


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129308/"
    "suppl/GSE129308_RAW.tar"
)
SOURCE_SHA256 = "dcd8ffee06d71c50f30e04b1aa7a6658dfc9a359e3a53b4e59f0f8afeabe9e0a"
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT / "inst" / "extdata" / "GSE129308_AD_pseudobulk_counts.csv"
)
SOURCE_FILES = (
    "GSM3704357_1-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704358_1-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704359_2-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704360_2-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704361_3-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704362_3-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704363_4-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704364_4-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704365_5-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704366_5-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704367_6-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704368_6-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704369_7-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704370_7-AT8_filtered_feature_bc_matrix.h5",
    "GSM3704371_8-MAP2_filtered_feature_bc_matrix.h5",
    "GSM3704372_8-AT8_filtered_feature_bc_matrix.h5",
)
GENES = (
    "NEFM",
    "APP",
    "SQSTM1",
    "HSP90AA1",
    "YWHAE",
    "WASF1",
    "CNTNAP1",
    "GOT2",
    "HSP90AB1",
    "HSPA8",
    "CCT2",
    "CCT3",
    "CCT4",
    "CCT5",
    "CCT6A",
    "CCT7",
    "CCT8",
    "TCP1",
    "DLG4",
    "SHANK1",
    "SHANK2",
    "SHANK3",
    "ACTB",
    "TUBB",
    "TUBB3",
    "SYN1",
    "SYP",
    "SNAP25",
    "STX1A",
    "VAMP2",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire_archive(path: Path | None, temporary_directory: Path) -> Path:
    archive = path.resolve() if path else temporary_directory / "GSE129308_RAW.tar"
    if path is None:
        print(f"Downloading {SOURCE_URL}")
        urllib.request.urlretrieve(SOURCE_URL, archive)
    observed = sha256(archive)
    if observed != SOURCE_SHA256:
        raise SystemExit(
            f"GSE129308 archive checksum mismatch: expected {SOURCE_SHA256}, got {observed}."
        )
    return archive


def extract_member(archive: tarfile.TarFile, name: str, destination: Path) -> None:
    source = archive.extractfile(name)
    if source is None:
        raise SystemExit(f"Required archive member is missing: {name}")
    with source, destination.open("wb") as target:
        shutil.copyfileobj(source, target)


def aggregate_sample(path: Path) -> tuple[int, int, dict[str, int]]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        feature_names = [value.decode("utf-8") for value in matrix["features/name"][:]]
        indices = matrix["indices"][:]
        values = matrix["data"][:]
        totals = np.bincount(indices, weights=values, minlength=len(feature_names))
        library_size = int(np.asarray(values, dtype=np.int64).sum())
        cell_count = int(len(matrix["barcodes"]))

    positions: dict[str, list[int]] = {}
    for index, name in enumerate(feature_names):
        positions.setdefault(name, []).append(index)
    counts = {
        gene: int(round(float(totals[positions.get(gene, [])].sum())))
        for gene in GENES
    }
    return cell_count, library_size, counts


def source_metadata(name: str) -> tuple[str, int, str, str]:
    geo_sample, label, *_ = name.split("_")
    donor_text, marker = label.split("-")
    donor = int(donor_text)
    population = "NFT-free" if marker == "MAP2" else "NFT-bearing"
    return geo_sample, donor, marker, population


def build_fixture(archive_path: Path, output: Path, temporary_directory: Path) -> None:
    rows: list[dict[str, str | int]] = []
    with tarfile.open(archive_path, "r") as archive:
        for source_name in SOURCE_FILES:
            local_h5 = temporary_directory / source_name
            extract_member(archive, source_name, local_h5)
            cell_count, library_size, counts = aggregate_sample(local_h5)
            local_h5.unlink()
            geo_sample, donor, marker, population = source_metadata(source_name)
            rows.append(
                {
                    "sample": f"donor_{donor}_{marker}",
                    "donor": donor,
                    "population": population,
                    "sort_marker": marker,
                    "geo_sample": geo_sample,
                    "cell_count": cell_count,
                    "library_size": library_size,
                    **counts,
                }
            )

    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "donor",
        "population",
        "sort_marker",
        "geo_sample",
        "cell_count",
        "library_size",
        *GENES,
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if len(rows) != 16:
        raise SystemExit(f"Expected 16 paired samples, wrote {len(rows)}.")
    print(f"Wrote {output} ({output.stat().st_size} bytes).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive",
        type=Path,
        help="Existing checksum-matched GSE129308_RAW.tar; download when omitted.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="genefunnel-gse129308-") as temporary:
        temporary_directory = Path(temporary)
        archive = acquire_archive(arguments.archive, temporary_directory)
        build_fixture(archive, arguments.output.resolve(), temporary_directory)


if __name__ == "__main__":
    main()
