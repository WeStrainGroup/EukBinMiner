#!/usr/bin/env python3
"""Split a contig-by-sample coverage matrix into ReMAG 0.4.4 TSV inputs.

ReMAG treats each ``-c`` TSV as one sample and reads its last column.  It does
not interpret several sample columns from one file, so this adapter writes one
headerless ``contig_id<TAB>coverage`` file per matrix column in stable order.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path
from typing import Optional, Sequence


class SplitError(ValueError):
    """The coverage matrix or requested output is unsafe."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_coverage(matrix: Path, out_dir: Path) -> list[dict[str, object]]:
    if not matrix.is_file() or matrix.stat().st_size == 0:
        raise SplitError(f"coverage matrix does not exist or is empty: {matrix}")
    if out_dir.exists():
        if any(out_dir.iterdir()):
            raise SplitError(f"refusing non-empty output directory: {out_dir}")
    else:
        out_dir.mkdir(parents=True)

    handles = []
    output_paths: list[Path] = []
    metadata: list[dict[str, object]] = []
    seen_ids: set[str] = set()
    try:
        with matrix.open("rt", encoding="utf-8", newline=None) as source:
            raw_header = source.readline()
            if not raw_header:
                raise SplitError(f"coverage matrix is empty: {matrix}")
            columns = raw_header.rstrip("\r\n").split("\t")
            if len(columns) < 2:
                raise SplitError("coverage matrix needs an ID and coverage column")
            if not columns[0]:
                raise SplitError("coverage matrix has an empty ID column name")
            coverage_names = columns[1:]
            if any(
                not name
                or name != name.strip()
                or any(character in name for character in "\r\n\t\x00")
                for name in coverage_names
            ):
                raise SplitError("coverage matrix has an empty or unsafe column name")
            duplicates = sorted(
                {name for name in coverage_names if coverage_names.count(name) > 1}
            )
            if duplicates:
                raise SplitError(
                    "coverage matrix has duplicate coverage name(s): "
                    + ", ".join(duplicates)
                )
            for index, name in enumerate(coverage_names, start=1):
                filename = f"coverage.{index:04d}.tsv"
                path = out_dir / filename
                output_paths.append(path)
                handles.append(path.open("wt", encoding="utf-8", newline="\n"))
                metadata.append(
                    {
                        "index": index,
                        "coverage_name": name,
                        "file": filename,
                    }
                )

            for line_number, raw_line in enumerate(source, start=2):
                row = raw_line.rstrip("\r\n")
                if not row:
                    continue
                fields = row.split("\t")
                if len(fields) != len(columns):
                    raise SplitError(
                        f"coverage row {line_number} has {len(fields)} fields; "
                        f"expected {len(columns)}"
                    )
                identifier = fields[0]
                if not identifier:
                    raise SplitError(f"empty contig ID at coverage row {line_number}")
                if identifier in seen_ids:
                    raise SplitError(f"duplicate contig ID in coverage: {identifier}")
                seen_ids.add(identifier)
                for name, value, handle in zip(
                    coverage_names,
                    fields[1:],
                    handles,
                ):
                    try:
                        parsed = float(value)
                    except ValueError as exc:
                        raise SplitError(
                            f"invalid coverage {name!r} for {identifier}: {value!r}"
                        ) from exc
                    if not math.isfinite(parsed) or parsed < 0:
                        raise SplitError(
                            f"coverage {name!r} must be finite and non-negative "
                            f"for {identifier}"
                        )
                    handle.write(f"{identifier}\t{value}\n")
    except (OSError, UnicodeError) as exc:
        raise SplitError(f"cannot split coverage matrix {matrix}: {exc}") from exc
    finally:
        for handle in handles:
            handle.close()

    for item, path in zip(metadata, output_paths):
        item["rows"] = len(seen_ids)
        item["sha256"] = sha256_file(path)
    manifest = out_dir / "manifest.tsv"
    with manifest.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write("index\tcoverage_name\tfile\trows\tsha256\n")
        for item in metadata:
            handle.write(
                f"{item['index']}\t{item['coverage_name']}\t{item['file']}\t"
                f"{item['rows']}\t{item['sha256']}\n"
            )
    return metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coverage-matrix", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        metadata = split_coverage(args.coverage_matrix, args.out_dir)
    except (SplitError, OSError) as exc:
        print(f"split_remag_coverage.py: error: {exc}", file=sys.stderr)
        return 2
    print(
        f"split_remag_coverage.py: wrote {len(metadata)} coverage file(s) "
        f"to {args.out_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
