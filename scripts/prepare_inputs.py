#!/usr/bin/env python3
"""Validate Tiara FASTA/depth inputs and build deterministic binning inputs.

Only Tiara-selected FASTA records longer than the requested input cutoff are
required to occur in the depth table.  The depth table may contain any number
of other assembly contigs and additional columns.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from contextlib import ExitStack
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, List, Mapping, Optional, Sequence

from depth_schema import DepthSchemaError, identify_depth_schema


REQUIRED_DEPTH_COLUMNS = ("contigName", "contigLen", "totalAvgDepth")
OUTPUT_NAMES = (
    "euk.fa",
    "union.fa",
    "contig_origin.tsv",
    "euk.metabat.depth.tsv",
    "union.metabat.depth.tsv",
    "euk.remag.coverage.tsv",
    "union.remag.coverage.tsv",
    "euk.concoct.coverage.tsv",
    "union.concoct.coverage.tsv",
    "preflight.json",
)


class PreflightError(ValueError):
    """An input or output-safety validation failure."""


@dataclass(frozen=True)
class Contig:
    contig_id: str
    length: int


class RawWriter:
    """A small binary writer that can safely terminate concatenated records."""

    def __init__(self, handle: BinaryIO):
        self.handle = handle
        self.last_byte: Optional[bytes] = None

    def write(self, data: bytes) -> None:
        self.handle.write(data)
        if data:
            self.last_byte = data[-1:]

    def write_line(self, data: bytes) -> None:
        self.write(data)
        if not data.endswith((b"\n", b"\r")):
            self.write(b"\n")

    def ensure_newline(self) -> None:
        if self.last_byte is not None and self.last_byte not in (b"\n", b"\r"):
            self.write(b"\n")


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError as exc:
        raise PreflightError(f"cannot read {path}: {exc}") from exc


def _open_maybe_gzip(path: Path) -> BinaryIO:
    try:
        if _is_gzip(path):
            return gzip.open(path, "rb")
        return path.open("rb")
    except (OSError, gzip.BadGzipFile) as exc:
        raise PreflightError(f"cannot open {path}: {exc}") from exc


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise PreflightError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def _decode(value: bytes, description: str) -> str:
    try:
        return value.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PreflightError(f"{description} is not valid UTF-8") from exc


def _scan_fasta(
    path: Path,
    label: str,
    outputs: Iterable[RawWriter],
    minimum_length_inclusive: int,
) -> tuple[List[Contig], set[str], dict]:
    """Validate a FASTA and stream only length-eligible records to outputs."""

    writers = tuple(outputs)
    records: List[Contig] = []
    seen: set[str] = set()
    current_id: Optional[str] = None
    current_length = 0
    current_lines: list[bytes] = []
    input_contigs = 0
    input_bp = 0

    def finish_record() -> None:
        nonlocal current_id, current_length, input_contigs, input_bp
        if current_id is not None:
            if current_length == 0:
                raise PreflightError(
                    f"{label} FASTA has an empty sequence for {current_id}"
                )
            input_contigs += 1
            input_bp += current_length
            if current_length >= minimum_length_inclusive:
                records.append(Contig(current_id, current_length))
                for writer in writers:
                    for raw_record_line in current_lines:
                        writer.write(raw_record_line)

    try:
        with _open_maybe_gzip(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if raw_line.startswith(b">"):
                    finish_record()
                    header = raw_line[1:].rstrip(b"\r\n")
                    tokens = header.split(None, 1)
                    if not tokens:
                        raise PreflightError(
                            f"{label} FASTA has an empty header at line {line_number}"
                        )
                    current_id = _decode(
                        tokens[0], f"{label} FASTA ID at line {line_number}"
                    )
                    if current_id in seen:
                        raise PreflightError(
                            f"duplicate contig ID in {label} FASTA: {current_id}"
                        )
                    seen.add(current_id)
                    current_length = 0
                    current_lines = [raw_line]
                else:
                    if current_id is None:
                        raise PreflightError(
                            f"{label} FASTA has content before its first header "
                            f"at line {line_number}"
                        )
                    sequence_line = raw_line.rstrip(b"\r\n")
                    current_length += sum(
                        len(piece) for piece in sequence_line.split()
                    )
                    current_lines.append(raw_line)
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise PreflightError(f"cannot read {label} FASTA {path}: {exc}") from exc

    finish_record()
    if input_contigs == 0:
        raise PreflightError(f"{label} FASTA is empty: {path}")
    for writer in writers:
        writer.ensure_newline()
    selected_bp = sum(contig.length for contig in records)
    return records, seen, {
        "input_contigs": input_contigs,
        "input_bp": input_bp,
        "selected_contigs": len(records),
        "selected_bp": selected_bp,
        "dropped_contigs": input_contigs - len(records),
        "dropped_bp": input_bp - selected_bp,
    }


def _write_origin(
    path: Path, euk_contigs: Sequence[Contig], unknown_contigs: Sequence[Contig]
) -> None:
    with path.open("wb") as handle:
        handle.write(b"contig_id\torigin\n")
        for contig in euk_contigs:
            handle.write(f"{contig.contig_id}\teuk\n".encode("utf-8"))
        for contig in unknown_contigs:
            handle.write(f"{contig.contig_id}\tunknown\n".encode("utf-8"))


def _parse_nonnegative_value(
    value: str,
    contig_id: str,
    line_number: int,
    column: str,
) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise PreflightError(
            f"invalid {column} for {contig_id} at depth line "
            f"{line_number}: {value!r}"
        ) from exc
    if not math.isfinite(parsed) or parsed < 0:
        raise PreflightError(
            f"{column} must be finite and non-negative for {contig_id} "
            f"at depth line {line_number}: {value!r}"
        )
    return parsed


def _parse_nonnegative_depth(
    value: str, contig_id: str, line_number: int
) -> float:
    return _parse_nonnegative_value(
        value,
        contig_id,
        line_number,
        "totalAvgDepth",
    )


def _parse_length(value: str, contig_id: str, line_number: int) -> int:
    try:
        parsed = Decimal(value)
    except (InvalidOperation, ValueError) as exc:
        raise PreflightError(
            f"invalid contigLen for {contig_id} at depth line "
            f"{line_number}: {value!r}"
        ) from exc
    if (
        not parsed.is_finite()
        or parsed < 0
        or parsed != parsed.to_integral_value()
    ):
        raise PreflightError(
            f"contigLen must be a finite non-negative integer for {contig_id} "
            f"at depth line "
            f"{line_number}: {value!r}"
        )
    return int(parsed)


def _write_depth_products(
    depth_path: Path,
    sample: str,
    stage: Path,
    euk_contigs: Sequence[Contig],
    unknown_contigs: Sequence[Contig],
) -> dict:
    euk_lengths = {contig.contig_id: contig.length for contig in euk_contigs}
    unknown_lengths = {
        contig.contig_id: contig.length for contig in unknown_contigs
    }
    selected_lengths = dict(euk_lengths)
    selected_lengths.update(unknown_lengths)
    selected_ids = set(selected_lengths)
    euk_ids = set(euk_lengths)
    selected_depth_rows: dict[str, str] = {}
    selected_coverage: Dict[str, float] = {}
    selected_coverages: Dict[str, Dict[str, float]] = {}
    depth_row_count = 0
    depth_extra_count = 0
    malformed_extra_lines: list[int] = []
    identical_selected_duplicate_lines: list[int] = []

    output_paths = {
        "euk_meta": stage / "euk.metabat.depth.tsv",
        "union_meta": stage / "union.metabat.depth.tsv",
        "euk_remag": stage / "euk.remag.coverage.tsv",
        "union_remag": stage / "union.remag.coverage.tsv",
        "euk_concoct": stage / "euk.concoct.coverage.tsv",
        "union_concoct": stage / "union.concoct.coverage.tsv",
    }

    try:
        with ExitStack() as stack:
            handles = {
                key: RawWriter(stack.enter_context(path.open("wb")))
                for key, path in output_paths.items()
            }

            depth_handle = stack.enter_context(_open_maybe_gzip(depth_path))
            raw_header = depth_handle.readline()
            if not raw_header:
                raise PreflightError(f"depth TSV is empty: {depth_path}")
            header_text = _decode(
                raw_header.rstrip(b"\r\n"), "depth TSV header"
            )
            columns = header_text.split("\t")
            try:
                schema = identify_depth_schema(
                    columns,
                    sample,
                    require_jgi=True,
                )
            except DepthSchemaError as exc:
                raise PreflightError(str(exc)) from exc
            column_index = {name: columns.index(name) for name in columns}
            extra_columns = [
                column
                for column in columns
                if column not in REQUIRED_DEPTH_COLUMNS
            ]
            coverage_header = (
                "contig_id\t" + "\t".join(schema.coverage_names) + "\n"
            ).encode("utf-8")
            handles["euk_concoct"].write(coverage_header)
            handles["union_concoct"].write(coverage_header)
            handles["euk_meta"].write_line(raw_header)
            handles["union_meta"].write_line(raw_header)

            for line_number, raw_row in enumerate(depth_handle, start=2):
                row_without_newline = raw_row.rstrip(b"\r\n")
                if not row_without_newline:
                    continue
                row_text = _decode(
                    row_without_newline, f"depth TSV row {line_number}"
                )
                fields = row_text.split("\t")
                if len(fields) != len(columns):
                    contig_index = column_index["contigName"]
                    possible_contig_id = (
                        fields[contig_index]
                        if contig_index < len(fields)
                        else ""
                    )
                    if possible_contig_id in selected_ids:
                        raise PreflightError(
                            f"selected contig {possible_contig_id!r} has a "
                            f"malformed depth row at line {line_number}: "
                            f"{len(fields)} fields; expected {len(columns)}"
                        )
                    malformed_extra_lines.append(line_number)
                    continue
                contig_id = fields[column_index["contigName"]]
                if not contig_id:
                    raise PreflightError(
                        f"depth TSV has an empty contigName at line {line_number}"
                    )
                depth_row_count += 1
                if contig_id not in selected_ids:
                    depth_extra_count += 1
                    continue
                if contig_id in selected_depth_rows:
                    if selected_depth_rows[contig_id] == row_text:
                        identical_selected_duplicate_lines.append(line_number)
                        continue
                    raise PreflightError(
                        f"conflicting duplicate selected contig ID in depth "
                        f"TSV: {contig_id}"
                    )
                selected_depth_rows[contig_id] = row_text
                contig_length = _parse_length(
                    fields[column_index["contigLen"]], contig_id, line_number
                )
                depth_text = fields[column_index["totalAvgDepth"]]
                depth_value = _parse_nonnegative_depth(
                    depth_text, contig_id, line_number
                )
                coverage_texts = [
                    fields[index] for index in schema.coverage_indices
                ]
                coverage_values = {
                    name: _parse_nonnegative_value(
                        value,
                        contig_id,
                        line_number,
                        name,
                    )
                    for name, value in zip(
                        schema.coverage_names,
                        coverage_texts,
                    )
                }
                for name, variance_index in zip(
                    schema.coverage_names,
                    schema.variance_indices,
                ):
                    if variance_index is not None:
                        _parse_nonnegative_value(
                            fields[variance_index],
                            contig_id,
                            line_number,
                            f"{name}-var",
                        )

                expected_length = selected_lengths[contig_id]
                if contig_length != expected_length:
                    raise PreflightError(
                        f"length mismatch for {contig_id}: FASTA has "
                        f"{expected_length} bp, depth TSV has {contig_length} bp"
                    )

                selected_coverage[contig_id] = depth_value
                selected_coverages[contig_id] = coverage_values
                remag_row = (
                    f"{contig_id}\t{depth_text}\n".encode("utf-8")
                )
                concoct_row = (
                    f"{contig_id}\t" + "\t".join(coverage_texts) + "\n"
                ).encode("utf-8")
                handles["union_meta"].write_line(raw_row)
                handles["union_remag"].write(remag_row)
                handles["union_concoct"].write(concoct_row)
                if contig_id in euk_ids:
                    handles["euk_meta"].write_line(raw_row)
                    handles["euk_remag"].write(remag_row)
                    handles["euk_concoct"].write(concoct_row)
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise PreflightError(f"cannot read depth TSV {depth_path}: {exc}") from exc

    missing_depth = sorted(selected_ids - set(selected_coverage))
    if missing_depth:
        preview = ", ".join(missing_depth[:10])
        suffix = "" if len(missing_depth) <= 10 else ", ..."
        raise PreflightError(
            f"{len(missing_depth)} selected contig(s) are missing from depth TSV: "
            f"{preview}{suffix}"
        )

    return {
        "coverage_by_id": selected_coverage,
        "coverages_by_id": selected_coverages,
        "coverage_columns": list(schema.coverage_names),
        "coverage_column_count": schema.coverage_column_count,
        "coverage_source_columns": list(schema.source_coverage_names),
        "depth_source_format": schema.source_format,
        "total_avg_depth_policy": schema.total_avg_depth_policy,
        "depth_row_count": depth_row_count,
        "depth_extra_count": depth_extra_count,
        "depth_extra_columns": extra_columns,
        "malformed_extra_row_count": len(malformed_extra_lines),
        "malformed_extra_row_lines": malformed_extra_lines,
        "identical_selected_duplicate_row_count": len(
            identical_selected_duplicate_lines
        ),
        "identical_selected_duplicate_row_lines": (
            identical_selected_duplicate_lines
        ),
    }


def _coverage_summary(
    contigs: Sequence[Contig], coverage_by_id: Mapping[str, float]
) -> dict:
    values = [coverage_by_id[contig.contig_id] for contig in contigs]
    if not values:
        return {
            "sum": 0.0,
            "mean": None,
            "min": None,
            "max": None,
        }
    total = math.fsum(values)
    return {
        "sum": total,
        "mean": total / len(values),
        "min": min(values),
        "max": max(values),
    }


def _collection_summary(
    contigs: Sequence[Contig],
    coverage_by_id: Mapping[str, float],
    coverages_by_id: Mapping[str, Mapping[str, float]],
    coverage_columns: Sequence[str],
) -> dict:
    contig_count = len(contigs)
    bp = sum(contig.length for contig in contigs)
    return {
        "contigs": contig_count,
        "contig_count": contig_count,
        "bp": bp,
        "total_bp": bp,
        "coverage": _coverage_summary(contigs, coverage_by_id),
        "coverage_by_column": {
            column: _coverage_summary(
                contigs,
                {
                    contig.contig_id: coverages_by_id[contig.contig_id][column]
                    for contig in contigs
                },
            )
            for column in coverage_columns
        },
    }


def _validate_input_file(path: Path, label: str) -> None:
    if not path.exists():
        raise PreflightError(f"{label} does not exist: {path}")
    if not path.is_file():
        raise PreflightError(f"{label} is not a regular file: {path}")


def prepare_inputs(
    euk_fasta: Path,
    unknown_fasta: Optional[Path],
    depth_tsv: Path,
    sample: str,
    out_dir: Path,
    minimum_length_inclusive: int = 1,
) -> str:
    """Build all outputs, returning ``created`` or ``already-current``."""

    if not sample or "\t" in sample or "\r" in sample or "\n" in sample:
        raise PreflightError(
            "sample must be non-empty and contain no tab or newline"
        )
    if minimum_length_inclusive < 1:
        raise PreflightError("--min-contig-bp must be positive")
    input_paths = {
        "euk_fasta": euk_fasta.resolve(),
        "depth_tsv": depth_tsv.resolve(),
    }
    if unknown_fasta is not None:
        input_paths["unknown_fasta"] = unknown_fasta.resolve()
    for label, path in input_paths.items():
        _validate_input_file(path, label)

    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise PreflightError(f"cannot create output directory {out_dir}: {exc}") from exc
    if not out_dir.is_dir():
        raise PreflightError(f"output path is not a directory: {out_dir}")

    final_paths = {name: out_dir / name for name in OUTPUT_NAMES}
    present = [name for name, path in final_paths.items() if path.exists()]
    if present and len(present) != len(OUTPUT_NAMES):
        missing = [name for name in OUTPUT_NAMES if name not in present]
        raise PreflightError(
            "refusing to overwrite a partial output set; present: "
            + ", ".join(present)
            + "; missing: "
            + ", ".join(missing)
        )
    all_outputs_present = len(present) == len(OUTPUT_NAMES)

    try:
        stage = Path(tempfile.mkdtemp(prefix=".prepare_inputs.", dir=out_dir))
    except OSError as exc:
        raise PreflightError(
            f"cannot create temporary output directory in {out_dir}: {exc}"
        ) from exc

    try:
        with ExitStack() as stack:
            euk_writer = RawWriter(
                stack.enter_context((stage / "euk.fa").open("wb"))
            )
            union_writer = RawWriter(
                stack.enter_context((stage / "union.fa").open("wb"))
            )
            euk_contigs, euk_all_ids, euk_filter = _scan_fasta(
                input_paths["euk_fasta"],
                "euk",
                (euk_writer, union_writer),
                minimum_length_inclusive,
            )
            if unknown_fasta is None:
                unknown_contigs = []
                unknown_all_ids = set()
                unknown_filter = {
                    key: 0
                    for key in (
                        "input_contigs",
                        "input_bp",
                        "selected_contigs",
                        "selected_bp",
                        "dropped_contigs",
                        "dropped_bp",
                    )
                }
            else:
                unknown_contigs, unknown_all_ids, unknown_filter = _scan_fasta(
                    input_paths["unknown_fasta"],
                    "unknown",
                    (union_writer,),
                    minimum_length_inclusive,
                )

        overlap = sorted(euk_all_ids & unknown_all_ids)
        if overlap:
            preview = ", ".join(overlap[:10])
            suffix = "" if len(overlap) <= 10 else ", ..."
            raise PreflightError(
                f"euk and unknown FASTA IDs overlap ({len(overlap)}): "
                f"{preview}{suffix}"
            )

        union_contigs = [*euk_contigs, *unknown_contigs]
        _write_origin(stage / "contig_origin.tsv", euk_contigs, unknown_contigs)
        depth_result = _write_depth_products(
            input_paths["depth_tsv"],
            sample,
            stage,
            euk_contigs,
            unknown_contigs,
        )
        coverage_by_id = depth_result["coverage_by_id"]
        coverages_by_id = depth_result["coverages_by_id"]
        coverage_columns = depth_result["coverage_columns"]

        input_manifest = {
            label: {
                "path": str(path),
                "sha256": _sha256_file(path),
            }
            for label, path in input_paths.items()
        }
        output_sha256 = {
            name: _sha256_file(stage / name)
            for name in OUTPUT_NAMES
            if name != "preflight.json"
        }
        preflight = {
            "schema_version": 3,
            "sample": sample,
            "inputs": input_manifest,
            "filter": {
                "minimum_contig_bp": minimum_length_inclusive,
                "euk": euk_filter,
                "unknown": unknown_filter,
                "union": {
                    key: euk_filter[key] + unknown_filter[key]
                    for key in euk_filter
                },
            },
            "collections": {
                "euk": _collection_summary(
                    euk_contigs,
                    coverage_by_id,
                    coverages_by_id,
                    coverage_columns,
                ),
                "unknown": _collection_summary(
                    unknown_contigs,
                    coverage_by_id,
                    coverages_by_id,
                    coverage_columns,
                ),
                "union": _collection_summary(
                    union_contigs,
                    coverage_by_id,
                    coverages_by_id,
                    coverage_columns,
                ),
            },
            "depth": {
                "rows": depth_result["depth_row_count"],
                "selected_count": len(union_contigs),
                "depth_extra_count": depth_result["depth_extra_count"],
                "malformed_extra_row_count": depth_result[
                    "malformed_extra_row_count"
                ],
                "malformed_extra_row_lines": depth_result[
                    "malformed_extra_row_lines"
                ],
                "identical_selected_duplicate_row_count": depth_result[
                    "identical_selected_duplicate_row_count"
                ],
                "identical_selected_duplicate_row_lines": depth_result[
                    "identical_selected_duplicate_row_lines"
                ],
                "extra_columns": depth_result["depth_extra_columns"],
                "extra_column_count": len(depth_result["depth_extra_columns"]),
                "format": depth_result["depth_source_format"],
                "coverage_source_columns": depth_result[
                    "coverage_source_columns"
                ],
                "coverage_columns": coverage_columns,
                "coverage_column_count": depth_result[
                    "coverage_column_count"
                ],
                "total_avg_depth_policy": depth_result[
                    "total_avg_depth_policy"
                ],
                "remag_multi_coverage_policy": (
                    "split_concoct_matrix_into_one_headerless_tsv_per_coverage_column"
                ),
            },
            "depth_extra_count": depth_result["depth_extra_count"],
            "output_sha256": output_sha256,
        }
        preflight_bytes = (
            json.dumps(
                preflight,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
        with (stage / "preflight.json").open("wb") as handle:
            handle.write(preflight_bytes)

        if all_outputs_present:
            mismatches = [
                name
                for name in OUTPUT_NAMES
                if _sha256_file(stage / name) != _sha256_file(final_paths[name])
            ]
            if mismatches:
                raise PreflightError(
                    "refusing to overwrite existing output(s) whose checksum "
                    "does not match this run: "
                    + ", ".join(mismatches)
                )
            return "already-current"

        for name in OUTPUT_NAMES:
            os.replace(stage / name, final_paths[name])
        return "created"
    finally:
        shutil.rmtree(stage, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Tiara-selected FASTA/depth inputs and prepare euk and "
            "euk+unknown coverage files for binning."
        )
    )
    parser.add_argument(
        "--contigs",
        type=Path,
        help=(
            "one already-selected contig FASTA; mutually exclusive with "
            "--euk-fasta/--unknown-fasta"
        ),
    )
    parser.add_argument("--euk-fasta", type=Path)
    parser.add_argument("--unknown-fasta", type=Path)
    parser.add_argument("--depth-tsv", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--min-contig-bp",
        type=int,
        default=1,
        help=(
            "minimum selected contig length, inclusive, before requiring "
            "depth coverage (default: 1)"
        ),
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.contigs is not None:
        if args.euk_fasta is not None or args.unknown_fasta is not None:
            print(
                "prepare_inputs.py: error: --contigs cannot be combined with "
                "--euk-fasta/--unknown-fasta",
                file=sys.stderr,
            )
            return 2
        euk_fasta = args.contigs
        unknown_fasta = None
    else:
        if args.euk_fasta is None or args.unknown_fasta is None:
            print(
                "prepare_inputs.py: error: provide either --contigs or both "
                "--euk-fasta and --unknown-fasta",
                file=sys.stderr,
            )
            return 2
        euk_fasta = args.euk_fasta
        unknown_fasta = args.unknown_fasta
    try:
        status = prepare_inputs(
            euk_fasta=euk_fasta,
            unknown_fasta=unknown_fasta,
            depth_tsv=args.depth_tsv,
            sample=args.sample,
            out_dir=args.out_dir,
            minimum_length_inclusive=args.min_contig_bp,
        )
    except PreflightError as exc:
        print(f"prepare_inputs.py: error: {exc}", file=sys.stderr)
        return 2
    print(f"prepare_inputs.py: {status}: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
