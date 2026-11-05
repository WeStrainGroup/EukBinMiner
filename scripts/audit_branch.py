#!/usr/bin/env python3
"""Audit one binner branch against the exact Tiara input FASTA."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, TextIO


FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
)


@dataclass(frozen=True)
class Record:
    identifier: str
    header: str
    sequence: str


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def read_fasta(path: Path) -> Iterator[Record]:
    header: Optional[str] = None
    sequence_parts: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield make_record(path, header, sequence_parts)
                header = line[1:]
                sequence_parts = []
            elif not line:
                continue
            elif header is None:
                raise ValueError(
                    f"{path}: sequence before first FASTA header at line {line_number}"
                )
            else:
                sequence_parts.append("".join(line.split()))
    if header is not None:
        yield make_record(path, header, sequence_parts)


def make_record(path: Path, header: str, sequence_parts: list[str]) -> Record:
    identifier_fields = header.split()
    if not identifier_fields:
        raise ValueError(f"{path}: empty FASTA identifier")
    sequence = "".join(sequence_parts)
    if not sequence:
        raise ValueError(f"{path}: empty sequence for {identifier_fields[0]}")
    return Record(identifier_fields[0], header, sequence)


def load_source(path: Path) -> OrderedDict[str, Record]:
    records: OrderedDict[str, Record] = OrderedDict()
    for record in read_fasta(path):
        if record.identifier in records:
            raise ValueError(f"{path}: duplicate source ID {record.identifier}")
        records[record.identifier] = record
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return records


def is_fasta(path: Path) -> bool:
    return path.is_file() and path.name.lower().endswith(FASTA_SUFFIXES)


def strip_fasta_suffix(name: str) -> str:
    lowered = name.lower()
    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if lowered.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, content: str) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite audit output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def format_fasta(records: Iterator[Record]) -> str:
    output: list[str] = []
    for record in records:
        output.append(f">{record.header}\n")
        sequence = record.sequence
        for start in range(0, len(sequence), 80):
            output.append(sequence[start : start + 80] + "\n")
    return "".join(output)


def audit(source_fasta: Path, bins_dir: Path, label: str, out_dir: Path) -> None:
    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(f"Refusing to reuse non-empty audit directory: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    source = load_source(source_fasta)
    assignments: dict[str, tuple[str, str]] = {}
    bin_files = sorted(path for path in bins_dir.rglob("*") if is_fasta(path))
    per_bin: list[dict[str, object]] = []

    for bin_path in bin_files:
        relative_path = bin_path.relative_to(bins_dir)
        bin_id = strip_fasta_suffix(relative_path.as_posix())
        contig_count = 0
        bin_bp = 0
        seen_in_file: set[str] = set()
        for record in read_fasta(bin_path):
            if record.identifier in seen_in_file:
                raise ValueError(f"{bin_path}: duplicate contig {record.identifier}")
            seen_in_file.add(record.identifier)
            if record.identifier not in source:
                raise ValueError(
                    f"{bin_path}: contig {record.identifier} is absent from source FASTA"
                )
            if record.identifier in assignments:
                previous_bin = assignments[record.identifier][0]
                raise ValueError(
                    f"Contig {record.identifier} occurs in multiple bins: "
                    f"{previous_bin}, {bin_id}"
                )
            source_sequence = source[record.identifier].sequence.upper()
            if record.sequence.upper() != source_sequence:
                raise ValueError(
                    f"{bin_path}: sequence mismatch for {record.identifier}"
                )
            assignments[record.identifier] = (bin_id, relative_path.as_posix())
            contig_count += 1
            bin_bp += len(record.sequence)
        if contig_count == 0:
            raise ValueError(f"{bin_path}: bin FASTA has no records")
        per_bin.append(
            {
                "bin_id": bin_id,
                "bin_file": relative_path.as_posix(),
                "contigs": contig_count,
                "bp": bin_bp,
            }
        )

    unbinned_ids = [
        identifier for identifier in source.keys() if identifier not in assignments
    ]
    source_bp = sum(len(record.sequence) for record in source.values())
    assigned_bp = sum(len(source[identifier].sequence) for identifier in assignments)
    unbinned_bp = source_bp - assigned_bp

    assignment_lines = [
        "contig_id\tlength_bp\tstatus\tbin_id\tbin_file\n"
    ]
    for identifier, record in source.items():
        if identifier in assignments:
            bin_id, bin_file = assignments[identifier]
            assignment_lines.append(
                f"{identifier}\t{len(record.sequence)}\tbinned\t{bin_id}\t{bin_file}\n"
            )
        else:
            assignment_lines.append(
                f"{identifier}\t{len(record.sequence)}\tunbinned\t\t\n"
            )

    unbinned_lengths = [len(source[identifier].sequence) for identifier in unbinned_ids]
    length_thresholds = [500, 1000, 1500, 3000, 10000]
    audit_result = {
        "label": label,
        "source_fasta": str(source_fasta),
        "source_sha256": sha256_file(source_fasta),
        "bins_dir": str(bins_dir),
        "bin_files": len(bin_files),
        "nonempty_bins": sum(1 for item in per_bin if int(item["contigs"]) > 0),
        "source_contigs": len(source),
        "source_bp": source_bp,
        "assigned_contigs": len(assignments),
        "assigned_bp": assigned_bp,
        "assigned_contig_fraction": len(assignments) / len(source),
        "assigned_bp_fraction": assigned_bp / source_bp,
        "unbinned_contigs": len(unbinned_ids),
        "unbinned_bp": unbinned_bp,
        "unbinned_at_or_above_threshold": {
            str(threshold): {
                "contigs": sum(length >= threshold for length in unbinned_lengths),
                "bp": sum(
                    length for length in unbinned_lengths if length >= threshold
                ),
            }
            for threshold in length_thresholds
        },
        "bins": per_bin,
    }

    for numeric_key in (
        "assigned_contig_fraction",
        "assigned_bp_fraction",
    ):
        if not math.isfinite(float(audit_result[numeric_key])):
            raise ValueError(f"Non-finite audit metric: {numeric_key}")

    atomic_text(out_dir / "assignments.tsv", "".join(assignment_lines))
    atomic_text(
        out_dir / "unbinned.list",
        "".join(f"{identifier}\n" for identifier in unbinned_ids),
    )
    atomic_text(
        out_dir / "unbinned.fa",
        format_fasta(source[identifier] for identifier in unbinned_ids),
    )
    atomic_text(
        out_dir / "audit.json",
        json.dumps(audit_result, indent=2, sort_keys=True) + "\n",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-fasta", required=True, type=Path)
    parser.add_argument("--bins-dir", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.source_fasta.is_file():
        raise FileNotFoundError(args.source_fasta)
    if not args.bins_dir.is_dir():
        raise NotADirectoryError(args.bins_dir)
    audit(args.source_fasta, args.bins_dir, args.label, args.out_dir)


if __name__ == "__main__":
    main()
