#!/usr/bin/env python3
"""Create a deterministic inventory of FASTA bins from labeled directories."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO, Dict, List, Optional, Sequence, Tuple


FASTA_SUFFIXES = (".fa", ".fasta", ".fna")


class SummaryError(ValueError):
    """A bin-discovery or FASTA-validation error."""


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError as exc:
        raise SummaryError(f"cannot read {path}: {exc}") from exc


def _open_maybe_gzip(path: Path) -> BinaryIO:
    try:
        if _is_gzip(path):
            return gzip.open(path, "rb")
        return path.open("rb")
    except (OSError, gzip.BadGzipFile) as exc:
        raise SummaryError(f"cannot open {path}: {exc}") from exc


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise SummaryError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def _is_fasta_path(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def _bin_id(path: Path) -> str:
    name = path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
        lower = lower[:-3]
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    raise SummaryError(f"not a recognized FASTA path: {path}")


def _n50(lengths: Sequence[int]) -> int:
    total = sum(lengths)
    if total == 0:
        return 0
    threshold = (total + 1) // 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    raise AssertionError("unreachable N50 state")


def _fasta_stats(path: Path) -> Tuple[int, int, int, int]:
    lengths: List[int] = []
    current_length: Optional[int] = None
    try:
        with _open_maybe_gzip(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if raw_line.startswith(b">"):
                    if not raw_line[1:].strip():
                        raise SummaryError(
                            f"empty FASTA header in {path} at line {line_number}"
                        )
                    if current_length is not None:
                        if current_length == 0:
                            raise SummaryError(
                                f"empty FASTA sequence in {path} before "
                                f"line {line_number}"
                            )
                        lengths.append(current_length)
                    current_length = 0
                else:
                    if current_length is None:
                        raise SummaryError(
                            f"content before first FASTA header in {path} at "
                            f"line {line_number}"
                        )
                    sequence_line = raw_line.rstrip(b"\r\n")
                    current_length += sum(
                        len(piece) for piece in sequence_line.split()
                    )
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise SummaryError(f"cannot read bin FASTA {path}: {exc}") from exc
    if current_length is not None:
        if current_length == 0:
            raise SummaryError(f"empty FASTA sequence at end of {path}")
        lengths.append(current_length)
    if not lengths:
        raise SummaryError(f"bin FASTA is empty: {path}")
    bp = sum(lengths)
    return len(lengths), bp, _n50(lengths), max(lengths)


def _parse_labeled_dirs(values: Sequence[str]) -> List[Tuple[str, Path]]:
    parsed: List[Tuple[str, Path]] = []
    seen_labels: set[str] = set()
    for value in values:
        if "=" not in value:
            raise SummaryError(
                f"--bins-label must have the form label=dir, got: {value!r}"
            )
        label, directory_text = value.split("=", 1)
        if not label or not directory_text:
            raise SummaryError(
                f"--bins-label must have non-empty label and dir: {value!r}"
            )
        if any(character in label for character in "\t\r\n"):
            raise SummaryError(f"label contains a tab or newline: {label!r}")
        if label in seen_labels:
            raise SummaryError(f"duplicate --bins-label label: {label}")
        seen_labels.add(label)
        directory = Path(directory_text).resolve()
        if not directory.exists():
            raise SummaryError(f"bin directory does not exist: {directory}")
        if not directory.is_dir():
            raise SummaryError(f"bin path is not a directory: {directory}")
        parsed.append((label, directory))
    return parsed


def _atomic_write(path: Path, data: bytes) -> None:
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(file_descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass


def summarize_bins(
    bins_labels: Sequence[str],
    out_dir: Path,
) -> Tuple[Path, Path]:
    labeled_dirs = _parse_labeled_dirs(bins_labels)
    records: List[dict] = []
    seen_within_label: Dict[Tuple[str, str], Path] = {}

    for label, directory in labeled_dirs:
        paths = sorted(
            (
                path
                for path in directory.rglob("*")
                if path.is_file() and _is_fasta_path(path)
            ),
            key=lambda path: str(path.relative_to(directory)),
        )
        for path in paths:
            bin_id = _bin_id(path)
            key = (label, bin_id)
            if key in seen_within_label:
                raise SummaryError(
                    f"duplicate bin ID {bin_id!r} within label {label!r}: "
                    f"{seen_within_label[key]} and {path}"
                )
            seen_within_label[key] = path
            contigs, bp, n50, longest = _fasta_stats(path)
            records.append(
                {
                    "label": label,
                    "bin_id": bin_id,
                    "bin_key": f"{label}:{bin_id}",
                    "path": str(path.resolve()),
                    "contigs": contigs,
                    "bp": bp,
                    "n50": n50,
                    "longest_contig": longest,
                    "sha256": _sha256_file(path),
                }
            )

    labels_by_bin_id: Dict[str, List[str]] = defaultdict(list)
    for record in records:
        labels_by_bin_id[record["bin_id"]].append(record["label"])
    duplicate_bin_ids = [
        {
            "bin_id": bin_id,
            "labels": sorted(labels),
            "bin_keys": sorted(
                record["bin_key"]
                for record in records
                if record["bin_id"] == bin_id
            ),
        }
        for bin_id, labels in sorted(labels_by_bin_id.items())
        if len(labels) > 1
    ]

    label_summaries: Dict[str, dict] = {}
    for label, _directory in labeled_dirs:
        label_records = [
            record for record in records if record["label"] == label
        ]
        label_summaries[label] = {
            "bins": len(label_records),
            "contigs": sum(record["contigs"] for record in label_records),
            "bp": sum(record["bp"] for record in label_records),
        }
    summary = {
        "schema_version": 1,
        "totals": {
            "bins": len(records),
            "contigs": sum(record["contigs"] for record in records),
            "bp": sum(record["bp"] for record in records),
        },
        "labels": label_summaries,
        "duplicate_bin_ids_across_labels": duplicate_bin_ids,
    }

    header = (
        "label\tbin_id\tbin_key\tpath\tcontigs\tbp\tn50\t"
        "longest_contig\tsha256\n"
    )
    inventory_lines = [header]
    for record in records:
        inventory_lines.append(
            "\t".join(
                str(record[column])
                for column in (
                    "label",
                    "bin_id",
                    "bin_key",
                    "path",
                    "contigs",
                    "bp",
                    "n50",
                    "longest_contig",
                    "sha256",
                )
            )
            + "\n"
        )

    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise SummaryError(
            f"cannot create summary output directory {out_dir}: {exc}"
        ) from exc
    if not out_dir.is_dir():
        raise SummaryError(f"summary output path is not a directory: {out_dir}")
    inventory_path = out_dir / "bin_inventory.tsv"
    summary_path = out_dir / "summary.json"
    _atomic_write(
        inventory_path, "".join(inventory_lines).encode("utf-8")
    )
    _atomic_write(
        summary_path,
        (
            json.dumps(
                summary,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8"),
    )
    completion_path = out_dir / "summary.complete.json"
    completion = {
        "schema_version": 1,
        "bin_inventory_sha256": _sha256_file(inventory_path),
        "summary_sha256": _sha256_file(summary_path),
    }
    _atomic_write(
        completion_path,
        (
            json.dumps(
                completion,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8"),
    )
    return inventory_path, summary_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Recursively inventory FASTA bins from labeled directories."
    )
    parser.add_argument(
        "--bins-label",
        action="append",
        required=True,
        metavar="LABEL=DIR",
        help="repeat for each binning method/input directory",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("."),
        help="output directory (default: current directory)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        inventory_path, summary_path = summarize_bins(
            args.bins_label, args.out_dir
        )
    except SummaryError as exc:
        print(f"summarize_bins.py: error: {exc}", file=sys.stderr)
        return 2
    print(f"summarize_bins.py: wrote {inventory_path} and {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
