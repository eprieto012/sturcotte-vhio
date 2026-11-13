from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Report whether sample-sheet CELL barcodes are present in Data R1."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(
            "/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data"
        ),
        help="Root folder containing patient folders.",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        default=["pt600280"],
        help="Patient folders to scan.",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="Max reads per FASTQ file to inspect (0 = all reads).",
    )
    parser.add_argument(
        "--barcode-len",
        type=int,
        default=16,
        help="Cell barcode length to extract from R1 (default: 16).",
    )
    parser.add_argument(
        "--allow-r1-fallback",
        action="store_true",
        help=(
            "If Data/*_R1.fastq.gz is missing, allow lane/general R1 fallback. "
            "Default checks only Data R1."
        ),
    )
    return parser.parse_args()


def find_condition_dirs(root: Path, patients: list[str]) -> list[Path]:
    dirs: list[Path] = []
    for patient in patients:
        patient_dir = root / patient
        if not patient_dir.exists():
            continue
        for subdir in sorted(patient_dir.iterdir()):
            if subdir.is_dir() and (subdir / "10x-sample-sheet.tsv").exists():
                dirs.append(subdir)
    return dirs


def load_cells(sample_sheet: Path) -> set[str]:
    cells: set[str] = set()
    with sample_sheet.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"No header in {sample_sheet}")
        cell_col = next(
            (name for name in reader.fieldnames if name and name.strip().startswith("CELL")),
            None,
        )
        if cell_col is None:
            raise RuntimeError(
                f"No CELL* column in {sample_sheet}; headers={reader.fieldnames}"
            )
        for row in reader:
            barcode = (row.get(cell_col) or "").strip()
            if barcode:
                cells.add(barcode)
    return cells


def select_r1_files(cond_dir: Path, allow_fallback: bool) -> list[Path]:
    data_files = sorted((cond_dir / "Data").glob("*_R1.fastq.gz"))
    if data_files or not allow_fallback:
        return data_files

    lane_files = sorted(cond_dir.rglob("*_R1_001.fastq.gz"))
    if lane_files:
        return lane_files

    return sorted(cond_dir.rglob("*R1*.fastq.gz"))


def iter_fastq_sequences(path: Path, max_reads: int = 0) -> Iterable[str]:
    reads = 0
    with gzip.open(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:
                yield line.strip()
                reads += 1
                if max_reads and reads >= max_reads:
                    break


def collect_r1_prefixes(
    r1_fastq_files: list[Path], max_reads: int, barcode_len: int
) -> set[str]:
    prefixes: set[str] = set()
    for fastq in r1_fastq_files:
        for seq in iter_fastq_sequences(fastq, max_reads=max_reads):
            prefix = seq[:barcode_len]
            if prefix:
                prefixes.add(prefix)
    return prefixes


def folder_label(cond_dir: Path, root: Path) -> str:
    return str(cond_dir.relative_to(root))


def main() -> None:
    args = parse_args()
    cond_dirs = find_condition_dirs(args.root, args.patients)
    if not cond_dirs:
        raise SystemExit("No condition folders with 10x-sample-sheet.tsv were found.")

    print("=== CELL barcodes found in Data R1 (per folder) ===")
    print("folder\tCELL_total\tfound_in_R1\tmissing_from_R1\tall_found\tmissing_examples")

    for cond_dir in sorted(cond_dirs):
        label = folder_label(cond_dir, args.root)
        cells = load_cells(cond_dir / "10x-sample-sheet.tsv")
        r1_files = select_r1_files(cond_dir, allow_fallback=args.allow_r1_fallback)
        r1_prefixes = collect_r1_prefixes(
            r1_files, max_reads=args.max_reads, barcode_len=args.barcode_len
        )

        found = len(cells & r1_prefixes)
        missing_values = sorted(cells - r1_prefixes)
        missing_count = len(missing_values)
        all_found = "YES" if missing_count == 0 else "NO"
        preview = ",".join(missing_values[:5]) if missing_values else "-"

        print(
            f"{label}\t{len(cells)}\t{found}\t{missing_count}\t{all_found}\t{preview}"
        )


if __name__ == "__main__":
    main()
