#!/usr/bin/env python3
"""Merge one-row EukCC CSV results while retaining every raw-bin record."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Optional, Sequence


META_COLUMNS = (
    "bin_rank",
    "bin_id",
    "study",
    "sample",
    "sample_rank",
    "binner",
    "contigs",
    "bp",
    "n50",
    "longest_contig",
    "sha256",
    "bin_path",
)
CANONICAL_EUKCC_COLUMNS = (
    "completeness",
    "contamination",
    "ncbi_lng",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def clean(value: object) -> str:
    return str(value if value is not None else "").replace("\t", " ").replace(
        "\r", " "
    ).replace("\n", " ")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rb") as probe:
        is_gzip = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if is_gzip else Path.open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def explicit_status(result_dir: Path) -> str:
    status_path = result_dir / "status.tsv"
    if not status_path.is_file():
        return ""
    try:
        rows = read_tsv(status_path)
    except (OSError, UnicodeError, csv.Error):
        return "invalid_status_tsv"
    if len(rows) != 1:
        return f"invalid_status_rows_{len(rows)}"
    status = (rows[0].get("eukcc_status") or "").strip()
    if status in {
        "no_placement_markers",
        "no_suitable_marker_set",
        "placement_alignment_failed",
        "metaeuk_empty_proteins",
    }:
        return status
    return f"invalid_explicit_status:{clean(status)}"


def quality_row(csv_path: Path) -> tuple[Optional[dict[str, str]], str]:
    try:
        with csv_path.open("rt", encoding="utf-8-sig", newline="") as handle:
            header = handle.readline()
            if not header:
                return None, "empty_csv"
            delimiter = "\t" if header.count("\t") > header.count(",") else ","
            handle.seek(0)
            reader = csv.DictReader(handle, delimiter=delimiter)
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error):
        return None, "unreadable_csv"
    if len(rows) != 1:
        return None, f"csv_rows_{len(rows)}"
    lowered = {column.lower() for column in (reader.fieldnames or [])}
    if not any("complet" in column for column in lowered):
        return None, "not_quality_csv"
    if not any("contam" in column for column in lowered):
        return None, "not_quality_csv"
    return rows[0], "ok"


def extract_taxid(row: dict[str, str]) -> str:
    lineage_value = ""
    for key, value in row.items():
        if key.lower() in {"ncbi_lng", "ncbi_lineage"}:
            lineage_value = value or ""
            break
    identifiers = re.findall(r"(?<![A-Za-z])(\d+)(?![A-Za-z])", lineage_value)
    return identifiers[-1] if identifiers else ""


def merge(bin_manifest: Path, eukcc_root: Path, out_dir: Path) -> None:
    bins = read_tsv(bin_manifest)
    result_records: list[tuple[dict[str, str], str, str, dict[str, str]]] = []
    dynamic_columns: set[str] = set()
    status_counts: Counter[str] = Counter()
    for bin_row in bins:
        result_dir = eukcc_root / bin_row["bin_id"]
        candidates: list[tuple[Path, dict[str, str]]] = []
        diagnostic_states: list[str] = []
        if result_dir.is_dir():
            for csv_path in sorted(result_dir.rglob("*.csv")):
                row, state = quality_row(csv_path)
                diagnostic_states.append(state)
                if row is not None:
                    candidates.append((csv_path, row))
        if len(candidates) == 1:
            candidate_path, result_row = candidates[0]
            csv_path = str(candidate_path)
            state = "ok"
            dynamic_columns.update(result_row)
        elif len(candidates) > 1:
            csv_path = ""
            result_row = {}
            state = f"ambiguous_quality_csv_{len(candidates)}"
        elif result_dir.is_dir():
            csv_path = ""
            result_row = {}
            state = explicit_status(result_dir)
            if not state:
                state = (
                    "no_quality_csv:"
                    + ",".join(sorted(set(diagnostic_states)))
                    if diagnostic_states
                    else "no_csv"
                )
        else:
            csv_path = ""
            result_row = {}
            state = "missing_result"
        status_counts[state] += 1
        result_records.append((bin_row, state, csv_path, result_row))

    ordered_dynamic = [
        *CANONICAL_EUKCC_COLUMNS,
        *sorted(dynamic_columns.difference(CANONICAL_EUKCC_COLUMNS)),
    ]
    columns = (
        *META_COLUMNS,
        "eukcc_status",
        "eukcc_csv",
        "taxid",
        *ordered_dynamic,
    )
    lines = ["\t".join(columns) + "\n"]
    for bin_row, state, csv_path, result_row in result_records:
        values = {
            **{column: bin_row.get(column, "") for column in META_COLUMNS},
            "eukcc_status": state,
            "eukcc_csv": csv_path,
            "taxid": extract_taxid(result_row),
            **result_row,
        }
        lines.append(
            "\t".join(clean(values.get(column, "")) for column in columns) + "\n"
        )

    merged_path = out_dir / "eukcc_merged.tsv"
    atomic_write(merged_path, "".join(lines))
    summary_path = out_dir / "eukcc_merged.summary.json"
    atomic_write(
        summary_path,
        json.dumps(
            {
                "schema_version": 1,
                "bins": len(bins),
                "status": dict(sorted(status_counts.items())),
                "bin_manifest": str(bin_manifest.resolve()),
                "bin_manifest_sha256": sha256_file(bin_manifest),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    atomic_write(
        out_dir / "eukcc_merged.complete.json",
        json.dumps(
            {
                "merged_sha256": sha256_file(merged_path),
                "summary_sha256": sha256_file(summary_path),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-manifest", required=True, type=Path)
    parser.add_argument("--eukcc-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    merge(args.bin_manifest, args.eukcc_root, args.out_dir)
    print(f"merge_eukcc.py: wrote {args.out_dir / 'eukcc_merged.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
