#!/usr/bin/env python3
"""Compare two Kraken-compatible classifiers by contig, taxid, and lineage."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import IO, Optional, Sequence


DOMAINS = {"Archaea", "Bacteria", "Eukaryota", "Viruses"}


def open_text(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def read_calls(path: Path) -> dict[str, tuple[str, str]]:
    calls: dict[str, tuple[str, str]] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}: malformed line {line_number}")
            status, identifier, taxid = fields[:3]
            if status not in {"C", "U"}:
                raise ValueError(
                    f"{path}: invalid status {status!r} at line {line_number}"
                )
            if identifier in calls:
                raise ValueError(
                    f"{path}: duplicate classification for {identifier!r}"
                )
            calls[identifier] = (status, taxid)
    return calls


def read_report_ancestry(path: Path) -> dict[str, tuple[str, ...]]:
    ancestry: dict[str, tuple[str, ...]] = {}
    stack: list[tuple[int, tuple[str, ...]]] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t", 5)
            if len(fields) != 6:
                raise ValueError(
                    f"{path}: malformed report line {line_number}"
                )
            taxid = fields[4]
            raw_name = fields[5]
            indent = len(raw_name) - len(raw_name.lstrip(" "))
            name = raw_name.strip()
            while stack and stack[-1][0] >= indent:
                stack.pop()
            parents = stack[-1][1] if stack else ()
            lineage = (*parents, name)
            ancestry[taxid] = lineage
            stack.append((indent, lineage))
    return ancestry


def domain(lineage: tuple[str, ...]) -> str:
    return next((name for name in lineage if name in DOMAINS), "")


def compare(
    baseline_output: Path,
    baseline_report: Path,
    candidate_output: Path,
    candidate_report: Path,
    output_directory: Path,
) -> dict[str, object]:
    if output_directory.exists():
        raise ValueError(
            f"comparison output already exists: {output_directory}"
        )
    baseline = read_calls(baseline_output)
    candidate = read_calls(candidate_output)
    if set(baseline) != set(candidate):
        raise ValueError(
            "classifier contig sets differ: "
            f"baseline_only={sorted(set(baseline) - set(candidate))[:20]}, "
            f"candidate_only={sorted(set(candidate) - set(baseline))[:20]}"
        )
    baseline_ancestry = read_report_ancestry(baseline_report)
    candidate_ancestry = read_report_ancestry(candidate_report)
    output_directory.mkdir(parents=True)
    rows: list[dict[str, str]] = []
    counters = {
        "status_agreement": 0,
        "taxid_agreement": 0,
        "domain_agreement": 0,
        "fungi_agreement": 0,
        "baseline_classified": 0,
        "candidate_classified": 0,
    }
    for identifier in sorted(baseline):
        baseline_status, baseline_taxid = baseline[identifier]
        candidate_status, candidate_taxid = candidate[identifier]
        baseline_lineage = baseline_ancestry.get(baseline_taxid, ())
        candidate_lineage = candidate_ancestry.get(candidate_taxid, ())
        baseline_domain = domain(baseline_lineage)
        candidate_domain = domain(candidate_lineage)
        baseline_fungi = "yes" if "Fungi" in baseline_lineage else "no"
        candidate_fungi = "yes" if "Fungi" in candidate_lineage else "no"
        status_equal = baseline_status == candidate_status
        taxid_equal = baseline_taxid == candidate_taxid
        domain_equal = baseline_domain == candidate_domain
        fungi_equal = baseline_fungi == candidate_fungi
        counters["status_agreement"] += int(status_equal)
        counters["taxid_agreement"] += int(taxid_equal)
        counters["domain_agreement"] += int(domain_equal)
        counters["fungi_agreement"] += int(fungi_equal)
        counters["baseline_classified"] += int(baseline_status == "C")
        counters["candidate_classified"] += int(candidate_status == "C")
        if not (status_equal and taxid_equal and domain_equal and fungi_equal):
            rows.append(
                {
                    "contig_id": identifier,
                    "baseline_status": baseline_status,
                    "candidate_status": candidate_status,
                    "baseline_taxid": baseline_taxid,
                    "candidate_taxid": candidate_taxid,
                    "baseline_domain": baseline_domain,
                    "candidate_domain": candidate_domain,
                    "baseline_fungi": baseline_fungi,
                    "candidate_fungi": candidate_fungi,
                }
            )

    columns = [
        "contig_id",
        "baseline_status",
        "candidate_status",
        "baseline_taxid",
        "candidate_taxid",
        "baseline_domain",
        "candidate_domain",
        "baseline_fungi",
        "candidate_fungi",
    ]
    with (output_directory / "differences.tsv").open(
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    total = len(baseline)
    summary: dict[str, object] = {
        "schema_version": 1,
        "baseline_output": str(baseline_output.resolve()),
        "baseline_report": str(baseline_report.resolve()),
        "candidate_output": str(candidate_output.resolve()),
        "candidate_report": str(candidate_report.resolve()),
        "contigs": total,
        "differences": len(rows),
        **counters,
        **{
            f"{name}_fraction": (
                counters[name] / total if total else 1.0
            )
            for name in (
                "status_agreement",
                "taxid_agreement",
                "domain_agreement",
                "fungi_agreement",
            )
        },
    }
    (output_directory / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-output", required=True, type=Path)
    parser.add_argument("--baseline-report", required=True, type=Path)
    parser.add_argument("--candidate-output", required=True, type=Path)
    parser.add_argument("--candidate-report", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--require-complete-agreement", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    summary = compare(
        args.baseline_output,
        args.baseline_report,
        args.candidate_output,
        args.candidate_report,
        args.output_directory,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.require_complete_agreement and summary["differences"] != 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
