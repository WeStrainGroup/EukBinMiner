#!/usr/bin/env python3
"""Materialize canonical FASTA and depth inputs for one CLI run."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import IO, Optional, Sequence


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path) -> IO[str]:
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if compressed:
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(
            descriptor,
            "wt",
            encoding="utf-8",
            newline="\n",
        ) as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def metadata_input_path(metadata: dict[str, object], name: str) -> Path:
    inputs = metadata.get("inputs")
    if not isinstance(inputs, dict) or name not in inputs:
        raise ValueError(f"input metadata lacks inputs.{name}")
    value = inputs[name]
    if isinstance(value, dict):
        value = value.get("path")
    if not isinstance(value, str) or not value:
        raise ValueError(f"input metadata has invalid inputs.{name}")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"source FASTA does not exist: {path}")
    return path


def read_depth_table(
    path: Path,
) -> tuple[
    list[str],
    dict[str, tuple[str, ...]],
    int,
    list[tuple[int, str]],
]:
    with open_text(path) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"depth table is empty: {path}")
        columns = header_line.rstrip("\r\n").split("\t")
        id_column = next(
            (
                name
                for name in (
                    "contigName",
                    "contig_id",
                    "contigname",
                    "contig",
                )
                if name in columns
            ),
            columns[0] if columns else "",
        )
        if not id_column:
            raise ValueError(f"depth table has no contig ID column: {path}")
        id_index = columns.index(id_column)
        rows: dict[str, tuple[str, ...]] = {}
        duplicate_rows_removed = 0
        malformed_rows: list[tuple[int, str]] = []
        for line_number, raw_line in enumerate(handle, start=2):
            fields = tuple(raw_line.rstrip("\r\n").split("\t"))
            if len(fields) != len(columns):
                identifier = (
                    fields[id_index]
                    if id_index < len(fields)
                    else ""
                )
                malformed_rows.append((line_number, identifier))
                continue
            if id_index >= len(fields) or not fields[id_index]:
                raise ValueError(
                    f"{path}: malformed contig ID at line {line_number}"
                )
            identifier = fields[id_index]
            previous = rows.get(identifier)
            if previous is None:
                rows[identifier] = fields
            elif fields == previous:
                duplicate_rows_removed += 1
            else:
                raise ValueError(
                    f"{path}: conflicting duplicate depth for "
                    f"{identifier!r} at line {line_number}"
                )
    return columns, rows, duplicate_rows_removed, malformed_rows


def materialize(
    input_metadata: Path,
    output: Path,
    summary: Path,
    minimum_contig_bp: int,
    depth_output: Optional[Path] = None,
) -> dict[str, object]:
    if minimum_contig_bp < 1:
        raise ValueError("--min-contig-bp must be positive")
    if output.exists():
        raise ValueError(f"refusing to overwrite output FASTA: {output}")
    if summary.exists():
        raise ValueError(f"refusing to overwrite summary: {summary}")
    if depth_output is not None:
        if depth_output.exists():
            raise ValueError(
                f"refusing to overwrite output depth table: {depth_output}"
            )
        resolved_outputs = {
            output.resolve(),
            summary.resolve(),
            depth_output.resolve(),
        }
        if len(resolved_outputs) != 3:
            raise ValueError(
                "output FASTA, depth table, and summary must be distinct"
            )
    try:
        metadata = json.loads(input_metadata.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read input metadata: {input_metadata}") from error
    if not isinstance(metadata, dict):
        raise ValueError(f"input metadata is not an object: {input_metadata}")

    sources = [
        metadata_input_path(metadata, "euk_fasta"),
        metadata_input_path(metadata, "unknown_fasta"),
    ]
    depth_source = metadata_input_path(metadata, "depth_tsv")
    (
        depth_columns,
        depth_rows,
        duplicate_depth_rows,
        malformed_depth_rows,
    ) = read_depth_table(depth_source)
    depth_ids = set(depth_rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    identifiers: set[str] = set()
    source_contigs = 0
    source_bp = 0
    dropped_no_depth = 0
    total_contigs = 0
    total_bp = 0
    selected_contigs = 0
    selected_bp = 0
    materialized_identifiers: list[str] = []
    try:
        with os.fdopen(
            descriptor,
            "wt",
            encoding="utf-8",
            newline="\n",
        ) as target:
            for source in sources:
                header: Optional[str] = None
                identifier: Optional[str] = None
                sequence: list[str] = []

                def finish() -> None:
                    nonlocal source_contigs, source_bp, dropped_no_depth
                    nonlocal total_contigs, total_bp
                    nonlocal selected_contigs, selected_bp
                    if identifier is None or header is None:
                        return
                    joined = "".join(sequence)
                    if not joined:
                        raise ValueError(
                            f"{source}: empty sequence for {identifier!r}"
                        )
                    length = len(joined)
                    source_contigs += 1
                    source_bp += length
                    if identifier not in depth_ids:
                        dropped_no_depth += 1
                        return
                    target.write(f">{header}\n")
                    for start in range(0, length, 80):
                        target.write(joined[start : start + 80] + "\n")
                    materialized_identifiers.append(identifier)
                    total_contigs += 1
                    total_bp += length
                    if length >= minimum_contig_bp:
                        selected_contigs += 1
                        selected_bp += length

                with open_text(source) as handle:
                    for line_number, raw_line in enumerate(handle, start=1):
                        line = raw_line.rstrip("\r\n")
                        if line.startswith(">"):
                            finish()
                            header = line[1:]
                            fields = header.split()
                            if not fields:
                                raise ValueError(
                                    f"{source}: empty FASTA header at "
                                    f"line {line_number}"
                                )
                            identifier = fields[0]
                            if identifier in identifiers:
                                raise ValueError(
                                    f"duplicate contig ID across source FASTAs: "
                                    f"{identifier!r}"
                                )
                            identifiers.add(identifier)
                            sequence = []
                        elif not line:
                            continue
                        elif identifier is None:
                            raise ValueError(
                                f"{source}: sequence before first header at "
                                f"line {line_number}"
                            )
                        else:
                            sequence.append("".join(line.split()))
                    finish()
            if source_contigs == 0:
                raise ValueError("original eukarya+unknown FASTA union is empty")
            malformed_input_rows = [
                (line_number, identifier)
                for line_number, identifier in malformed_depth_rows
                if identifier in identifiers
            ]
            if malformed_input_rows:
                line_number, identifier = malformed_input_rows[0]
                raise ValueError(
                    f"{depth_source}: malformed depth row for input contig "
                    f"{identifier!r} at line {line_number}"
                )
            target.flush()
            os.fsync(target.fileno())
        if depth_output is not None:
            depth_content = ["\t".join(depth_columns) + "\n"]
            depth_content.extend(
                "\t".join(depth_rows[identifier]) + "\n"
                for identifier in depth_rows
            )
            atomic_text(depth_output, "".join(depth_content))
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        if depth_output is not None and not output.exists():
            depth_output.unlink(missing_ok=True)
        raise

    result = {
        "schema_version": 1,
        "input_metadata": str(input_metadata.resolve()),
        "input_metadata_sha256": sha256_file(input_metadata),
        "sources": [
            {
                "path": str(source),
                "sha256": sha256_file(source),
                "bytes": source.stat().st_size,
            }
            for source in sources
        ],
        "depth_source": {
            "path": str(depth_source),
            "sha256": sha256_file(depth_source),
            "bytes": depth_source.stat().st_size,
            "contig_ids": len(depth_ids),
            "duplicate_rows_removed": duplicate_depth_rows,
            "malformed_rows_removed": len(malformed_depth_rows),
            "malformed_row_line_numbers": [
                line_number
                for line_number, _identifier in malformed_depth_rows[:100]
            ],
            "malformed_row_line_numbers_truncated": (
                len(malformed_depth_rows) > 100
            ),
        },
        "minimum_contig_bp_inclusive": minimum_contig_bp,
        "source_contigs": source_contigs,
        "source_bp": source_bp,
        "dropped_no_depth": dropped_no_depth,
        "total_contigs": total_contigs,
        "total_bp": total_bp,
        "selected_contigs": selected_contigs,
        "selected_bp": selected_bp,
        "output": str(output.resolve()),
        "output_sha256": sha256_file(output),
        "output_bytes": output.stat().st_size,
    }
    if depth_output is not None:
        result["depth_output"] = {
            "path": str(depth_output.resolve()),
            "sha256": sha256_file(depth_output),
            "bytes": depth_output.stat().st_size,
            "contig_ids": len(depth_rows),
            "duplicate_rows_removed": duplicate_depth_rows,
            "malformed_rows_removed": len(malformed_depth_rows),
        }
    atomic_text(
        summary,
        json.dumps(result, indent=2, sort_keys=True) + "\n",
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--min-contig-bp", required=True, type=int)
    parser.add_argument("--depth-output", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = materialize(
            args.input_metadata,
            args.output,
            args.summary,
            args.min_contig_bp,
            args.depth_output,
        )
    except (OSError, ValueError) as error:
        print(
            f"materialize_cli_contigs.py: error: {error}",
            file=os.sys.stderr,
        )
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
