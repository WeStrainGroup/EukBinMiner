#!/usr/bin/env python3
"""Filter one prepared FASTA/depth set by an exclusive contig-length cutoff."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Optional, Sequence

from depth_schema import DepthSchemaError, identify_depth_schema


OUTPUT_NAMES = (
    "input.fa",
    "metabat.depth.tsv",
    "remag.coverage.tsv",
    "concoct.coverage.tsv",
    "filter.json",
)
REQUIRED_DEPTH_COLUMNS = ("contigName", "contigLen", "totalAvgDepth")


class FilterError(ValueError):
    """An invalid input, filtering, or output-safety condition."""


@dataclass(frozen=True)
class Record:
    identifier: str
    header: str
    sequence: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise FilterError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def read_fasta(path: Path) -> list[Record]:
    records: list[Record] = []
    seen: set[str] = set()
    header: Optional[str] = None
    sequence_parts: list[str] = []

    def finish_record() -> None:
        if header is None:
            return
        fields = header.split()
        if not fields:
            raise FilterError(f"{path}: empty FASTA header")
        identifier = fields[0]
        sequence = "".join(sequence_parts)
        if not sequence:
            raise FilterError(f"{path}: empty sequence for {identifier}")
        if identifier in seen:
            raise FilterError(f"{path}: duplicate FASTA ID {identifier}")
        seen.add(identifier)
        records.append(Record(identifier, header, sequence))

    try:
        with path.open("rt", encoding="utf-8", newline=None) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if line.startswith(">"):
                    finish_record()
                    header = line[1:]
                    sequence_parts = []
                elif not line:
                    continue
                elif header is None:
                    raise FilterError(
                        f"{path}: sequence before first header at line {line_number}"
                    )
                else:
                    sequence_parts.append("".join(line.split()))
    except (OSError, UnicodeError) as exc:
        raise FilterError(f"cannot read FASTA {path}: {exc}") from exc
    finish_record()
    if not records:
        raise FilterError(f"{path}: no FASTA records")
    return records


def format_fasta(records: Sequence[Record]) -> str:
    output: list[str] = []
    for record in records:
        output.append(f">{record.header}\n")
        for start in range(0, len(record.sequence), 80):
            output.append(record.sequence[start : start + 80] + "\n")
    return "".join(output)


def parse_exact_nonnegative_length(
    value: str,
    identifier: str,
    line_number: int,
) -> int:
    try:
        parsed = Decimal(value)
    except (InvalidOperation, ValueError) as exc:
        raise FilterError(
            f"invalid length for {identifier} at row {line_number}: {value!r}"
        ) from exc
    if (
        not parsed.is_finite()
        or parsed < 0
        or parsed != parsed.to_integral_value()
    ):
        raise FilterError(
            f"length must be a finite non-negative integer for {identifier} "
            f"at row {line_number}: {value!r}"
        )
    return int(parsed)


def parse_depth(
    path: Path,
    selected_lengths: dict[str, int],
    sample: str,
) -> tuple[
    str,
    list[str],
    OrderedDict[str, str],
    OrderedDict[str, tuple[str, ...]],
    tuple[str, ...],
    int,
]:
    selected_rows_by_id: dict[str, str] = {}
    selected_depth_by_id: dict[str, str] = {}
    selected_coverages_by_id: dict[str, tuple[str, ...]] = {}
    seen: set[str] = set()
    row_count = 0
    try:
        with path.open("rt", encoding="utf-8", newline=None) as handle:
            raw_header = handle.readline()
            if not raw_header:
                raise FilterError(f"depth table is empty: {path}")
            header = raw_header.rstrip("\r\n")
            columns = header.split("\t")
            try:
                schema = identify_depth_schema(
                    columns,
                    sample,
                    require_jgi=True,
                )
            except DepthSchemaError as exc:
                raise FilterError(str(exc)) from exc
            index = {column: columns.index(column) for column in columns}
            for line_number, raw_line in enumerate(handle, start=2):
                row = raw_line.rstrip("\r\n")
                if not row:
                    continue
                fields = row.split("\t")
                if len(fields) != len(columns):
                    raise FilterError(
                        f"depth row {line_number} has {len(fields)} fields; "
                        f"expected {len(columns)}"
                    )
                identifier = fields[index["contigName"]]
                if identifier in seen:
                    raise FilterError(f"duplicate depth ID: {identifier}")
                seen.add(identifier)
                row_count += 1
                if identifier not in selected_lengths:
                    continue
                length = parse_exact_nonnegative_length(
                    fields[index["contigLen"]],
                    identifier,
                    line_number,
                )
                try:
                    depth_value = float(fields[index["totalAvgDepth"]])
                except ValueError as exc:
                    raise FilterError(
                        f"invalid depth for {identifier} at row {line_number}"
                    ) from exc
                if length != selected_lengths[identifier]:
                    raise FilterError(
                        f"length mismatch for {identifier}: FASTA has "
                        f"{selected_lengths[identifier]}, depth has {length}"
                    )
                if not math.isfinite(depth_value) or depth_value < 0:
                    raise FilterError(
                        f"depth must be finite and non-negative for {identifier}"
                    )
                coverage_values: list[str] = []
                for name, coverage_index, variance_index in zip(
                    schema.coverage_names,
                    schema.coverage_indices,
                    schema.variance_indices,
                ):
                    coverage_text = fields[coverage_index]
                    try:
                        coverage_value = float(coverage_text)
                    except ValueError as exc:
                        raise FilterError(
                            f"invalid coverage {name!r} for {identifier} at "
                            f"row {line_number}"
                        ) from exc
                    if not math.isfinite(coverage_value) or coverage_value < 0:
                        raise FilterError(
                            f"coverage {name!r} must be finite and "
                            f"non-negative for {identifier}"
                        )
                    if variance_index is not None:
                        variance_text = fields[variance_index]
                        try:
                            variance_value = float(variance_text)
                        except ValueError as exc:
                            raise FilterError(
                                f"invalid variance {name!r} for {identifier} "
                                f"at row {line_number}"
                            ) from exc
                        if (
                            not math.isfinite(variance_value)
                            or variance_value < 0
                        ):
                            raise FilterError(
                                f"variance {name!r} must be finite and "
                                f"non-negative for {identifier}"
                            )
                    coverage_values.append(coverage_text)
                selected_rows_by_id[identifier] = row + "\n"
                selected_depth_by_id[identifier] = fields[
                    index["totalAvgDepth"]
                ]
                selected_coverages_by_id[identifier] = tuple(coverage_values)
    except (OSError, UnicodeError) as exc:
        raise FilterError(f"cannot read depth table {path}: {exc}") from exc

    missing_ids = sorted(set(selected_lengths) - set(selected_depth_by_id))
    if missing_ids:
        preview = ", ".join(missing_ids[:10])
        raise FilterError(
            f"{len(missing_ids)} selected contig(s) missing from depth: {preview}"
        )
    # MetaBAT2 requires abundance rows to occur in exactly the same order as
    # the contigs in the assembly FASTA.
    selected_rows = [
        selected_rows_by_id[identifier] for identifier in selected_lengths
    ]
    selected_depth = OrderedDict(
        (
            identifier,
            selected_depth_by_id[identifier],
        )
        for identifier in selected_lengths
    )
    selected_coverages = OrderedDict(
        (
            identifier,
            selected_coverages_by_id[identifier],
        )
        for identifier in selected_lengths
    )
    return (
        header + "\n",
        selected_rows,
        selected_depth,
        selected_coverages,
        schema.coverage_names,
        row_count,
    )


def write_stage(
    stage: Path,
    input_fasta: Path,
    metabat_depth: Path,
    sample: str,
    minimum_inclusive: int,
) -> None:
    records = read_fasta(input_fasta)
    selected = [
        record for record in records if len(record.sequence) >= minimum_inclusive
    ]
    if not selected:
        raise FilterError(
            f"no contigs remain after requiring length >= {minimum_inclusive} bp"
        )
    selected_lengths = {
        record.identifier: len(record.sequence) for record in selected
    }
    (
        header,
        selected_rows,
        selected_depth,
        selected_coverages,
        coverage_columns,
        depth_row_count,
    ) = parse_depth(
        metabat_depth,
        selected_lengths,
        sample,
    )

    (stage / "input.fa").write_text(
        format_fasta(selected), encoding="utf-8", newline="\n"
    )
    (stage / "metabat.depth.tsv").write_text(
        header + "".join(selected_rows), encoding="utf-8", newline="\n"
    )
    remag_lines = [
        f"{identifier}\t{depth}\n"
        for identifier, depth in selected_depth.items()
    ]
    (stage / "remag.coverage.tsv").write_text(
        "".join(remag_lines), encoding="utf-8", newline="\n"
    )
    (stage / "concoct.coverage.tsv").write_text(
        "contig_id\t"
        + "\t".join(coverage_columns)
        + "\n"
        + "".join(
            f"{identifier}\t" + "\t".join(values) + "\n"
            for identifier, values in selected_coverages.items()
        ),
        encoding="utf-8",
        newline="\n",
    )

    input_bp = sum(len(record.sequence) for record in records)
    selected_bp = sum(len(record.sequence) for record in selected)
    manifest = {
        "schema_version": 3,
        "sample": sample,
        "minimum_contig_bp": minimum_inclusive,
        "inputs": {
            "fasta": {
                "path": str(input_fasta.resolve()),
                "sha256": sha256_file(input_fasta),
            },
            "metabat_depth": {
                "path": str(metabat_depth.resolve()),
                "sha256": sha256_file(metabat_depth),
                "rows": depth_row_count,
            },
        },
        "depth": {
            "coverage_columns": list(coverage_columns),
            "coverage_column_count": len(coverage_columns),
            "remag_multi_coverage_policy": (
                "split_concoct_matrix_into_one_headerless_tsv_per_coverage_column"
            ),
        },
        "input": {"contigs": len(records), "bp": input_bp},
        "selected": {"contigs": len(selected), "bp": selected_bp},
        "dropped": {
            "contigs": len(records) - len(selected),
            "bp": input_bp - selected_bp,
        },
        "output_sha256": {
            name: sha256_file(stage / name)
            for name in OUTPUT_NAMES
            if name != "filter.json"
        },
    }
    (stage / "filter.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def filter_inputs(
    input_fasta: Path,
    metabat_depth: Path,
    sample: str,
    minimum_inclusive: int,
    out_dir: Path,
) -> str:
    if minimum_inclusive < 1:
        raise FilterError("--min-contig-bp must be positive")
    if not sample or any(character in sample for character in "\t\r\n"):
        raise FilterError("sample must be non-empty and contain no tab/newline")
    for path in (input_fasta, metabat_depth):
        if not path.is_file():
            raise FilterError(f"input file does not exist: {path}")

    out_dir.mkdir(parents=True, exist_ok=True)
    final_paths = {name: out_dir / name for name in OUTPUT_NAMES}
    present = [name for name, path in final_paths.items() if path.exists()]
    if present and len(present) != len(OUTPUT_NAMES):
        missing = [name for name in OUTPUT_NAMES if name not in present]
        raise FilterError(
            "refusing partial output set; present: "
            + ", ".join(present)
            + "; missing: "
            + ", ".join(missing)
        )
    all_present = len(present) == len(OUTPUT_NAMES)
    stage = Path(tempfile.mkdtemp(prefix=".filter_inputs.", dir=out_dir))
    try:
        write_stage(
            stage,
            input_fasta,
            metabat_depth,
            sample,
            minimum_inclusive,
        )
        if all_present:
            mismatches = [
                name
                for name in OUTPUT_NAMES
                if sha256_file(stage / name) != sha256_file(final_paths[name])
            ]
            if mismatches:
                raise FilterError(
                    "refusing to overwrite non-matching output(s): "
                    + ", ".join(mismatches)
                )
            return "already-current"
        for name in OUTPUT_NAMES:
            os.replace(stage / name, final_paths[name])
        return "created"
    finally:
        shutil.rmtree(stage, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--metabat-depth", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--min-contig-bp", required=True, type=int)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        status = filter_inputs(
            args.input_fasta,
            args.metabat_depth,
            args.sample,
            args.min_contig_bp,
            args.out_dir,
        )
    except (FilterError, OSError) as exc:
        print(f"filter_prepared_inputs.py: error: {exc}", file=sys.stderr)
        return 2
    print(f"filter_prepared_inputs.py: {status}: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
