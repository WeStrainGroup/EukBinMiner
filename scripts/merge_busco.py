#!/usr/bin/env python3
"""Merge BUSCO short summaries for all selected candidate bins."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Optional, Sequence


SUMMARY_PATTERN = re.compile(
    r"C:(?P<complete>[0-9.]+)%"
    r"\[S:(?P<single>[0-9.]+)%,D:(?P<duplicated>[0-9.]+)%\],"
    r"F:(?P<fragmented>[0-9.]+)%,"
    r"M:(?P<missing>[0-9.]+)%,"
    r"n:(?P<n>[0-9]+)"
)
COUNT_PATTERNS = {
    "complete": re.compile(r"(?m)^\s*(\d+)\s+Complete BUSCOs \(C\)\s*$"),
    "single": re.compile(
        r"(?m)^\s*(\d+)\s+Complete and single-copy BUSCOs \(S\)\s*$"
    ),
    "duplicated": re.compile(
        r"(?m)^\s*(\d+)\s+Complete and duplicated BUSCOs \(D\)\s*$"
    ),
    "fragmented": re.compile(r"(?m)^\s*(\d+)\s+Fragmented BUSCOs \(F\)\s*$"),
    "missing": re.compile(r"(?m)^\s*(\d+)\s+Missing BUSCOs \(M\)\s*$"),
    "n": re.compile(r"(?m)^\s*(\d+)\s+Total BUSCO groups searched\s*$"),
}
BUSCO_COLUMNS = (
    "busco_status",
    "busco_summary",
    "busco_complete",
    "busco_single",
    "busco_duplicated",
    "busco_fragmented",
    "busco_missing",
    "busco_n",
)
METRIC_COLUMNS = (
    "complete",
    "single",
    "duplicated",
    "fragmented",
    "missing",
    "n",
)


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


def parse_summary(path: Path) -> Optional[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = SUMMARY_PATTERN.search(text.replace(" ", ""))
    if match:
        values = match.groupdict()
        metrics = [float(values[key]) for key in METRIC_COLUMNS[:-1]]
        if int(values["n"]) <= 0 or any(
            not math.isfinite(value) or not 0 <= value <= 100 for value in metrics
        ):
            return None
        complete, single, duplicated, fragmented, missing = metrics
        if abs(complete - single - duplicated) > 0.21:
            return None
        if abs(complete + fragmented + missing - 100) > 0.31:
            return None
        return values

    # BUSCO 5/6 also writes a line-oriented count section.  This is the
    # format used by the user's existing Notion workflow and remains useful
    # when the compact C:[S,D],F,M,n summary is absent.
    counts: dict[str, int] = {}
    for name, pattern in COUNT_PATTERNS.items():
        count_match = pattern.search(text)
        if count_match is None:
            return None
        counts[name] = int(count_match.group(1))
    total = counts["n"]
    if (total <= 0
        or counts["complete"] != counts["single"] + counts["duplicated"]
        or counts["complete"] + counts["fragmented"] + counts["missing"] != total):
        return None

    def percentage(count: int) -> str:
        return f"{100.0 * count / total:.6f}".rstrip("0").rstrip(".")

    return {
        "complete": percentage(counts["complete"]),
        "single": percentage(counts["single"]),
        "duplicated": percentage(counts["duplicated"]),
        "fragmented": percentage(counts["fragmented"]),
        "missing": percentage(counts["missing"]),
        "n": str(total),
    }


def merge(candidates: Path, busco_root: Path, output: Path) -> None:
    with candidates.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = list(reader.fieldnames or [])
    output_columns = [*columns, *BUSCO_COLUMNS]
    lines = ["\t".join(output_columns) + "\n"]
    for row in rows:
        result_dir = busco_root / row["bin_id"]
        parsed_candidates: list[tuple[Path, dict[str, str]]] = []
        if result_dir.is_dir():
            for summary_path in sorted(result_dir.rglob("short_summary*.txt")):
                parsed = parse_summary(summary_path)
                if parsed is not None:
                    parsed_candidates.append((summary_path, parsed))
        unique_metrics = {
            tuple(parsed[column] for column in METRIC_COLUMNS)
            for _, parsed in parsed_candidates
        }
        if parsed_candidates and len(unique_metrics) == 1:
            # BUSCO 6 commonly writes both a top-level
            # short_summary.specific.*.txt and an identical copy inside the
            # lineage run directory.  Prefer the canonical top-level summary
            # but do not report identical duplicates as an ambiguity.
            summary_path, parsed = min(
                parsed_candidates,
                key=lambda candidate: (
                    0
                    if candidate[0].name.startswith("short_summary.specific.")
                    else 1,
                    len(candidate[0].relative_to(result_dir).parts),
                    str(candidate[0]),
                ),
            )
            busco = {
                "busco_status": "ok",
                "busco_summary": str(summary_path),
                "busco_complete": parsed["complete"],
                "busco_single": parsed["single"],
                "busco_duplicated": parsed["duplicated"],
                "busco_fragmented": parsed["fragmented"],
                "busco_missing": parsed["missing"],
                "busco_n": parsed["n"],
            }
        elif parsed_candidates:
            busco = {
                column: "" for column in BUSCO_COLUMNS
            }
            busco["busco_status"] = (
                f"ambiguous_summary_{len(parsed_candidates)}"
                f"_metrics_{len(unique_metrics)}"
            )
        else:
            busco = {column: "" for column in BUSCO_COLUMNS}
            busco["busco_status"] = (
                "unparseable_or_missing_summary"
                if result_dir.is_dir()
                else "missing_result"
            )
        values = {**row, **busco}
        lines.append(
            "\t".join(
                values.get(column, "")
                .replace("\t", " ")
                .replace("\r", " ")
                .replace("\n", " ")
                for column in output_columns
            )
            + "\n"
        )
    atomic_write(output, "".join(lines))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--busco-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    merge(args.candidates, args.busco_root, args.output)
    print(f"merge_busco.py: wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
