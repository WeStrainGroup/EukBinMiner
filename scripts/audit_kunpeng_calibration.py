#!/usr/bin/env python3
"""Audit one Kun-peng calibration sample against the historical Kraken run."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import IO, Optional, Sequence


PREFIX_RE = re.compile(r"^r([0-9]{6})C(.+)$")


def open_text(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def read_classifier(
    path: Path,
    *,
    selected_rank: Optional[int] = None,
) -> dict[str, tuple[str, ...]]:
    rows: dict[str, tuple[str, ...]] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}: malformed line {line_number}")
            if fields[0] not in {"C", "U"}:
                raise ValueError(
                    f"{path}: invalid status at line {line_number}"
                )
            match = PREFIX_RE.fullmatch(fields[1])
            if match is None:
                raise ValueError(
                    f"{path}: unexpected contig ID at line {line_number}: "
                    f"{fields[1]!r}"
                )
            rank = int(match.group(1))
            if selected_rank is not None and rank != selected_rank:
                continue
            contig = match.group(2)
            if contig in rows:
                raise ValueError(f"{path}: duplicate contig {contig!r}")
            fields[1] = contig
            rows[contig] = tuple(fields)
    return rows


def read_taxonomy(path: Path) -> dict[str, str]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["contigs", "predictions"]:
            raise ValueError(f"{path}: unexpected taxonomy header")
        rows: dict[str, str] = {}
        for row in reader:
            contig = row["contigs"]
            if contig in rows:
                raise ValueError(f"{path}: duplicate contig {contig!r}")
            rows[contig] = row["predictions"]
    return rows


def canonical_groups(path: Path) -> tuple[tuple[str, ...], ...]:
    groups: dict[str, list[str]] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig_id", "status", "bin_id"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: unexpected assignments header")
        seen: set[str] = set()
        for row in reader:
            contig = row["contig_id"]
            if contig in seen:
                raise ValueError(f"{path}: duplicate contig {contig!r}")
            seen.add(contig)
            if row["status"] == "binned":
                groups.setdefault(row["bin_id"], []).append(contig)
    canonical = {tuple(sorted(contigs)) for contigs in groups.values()}
    if len(canonical) != len(groups):
        raise ValueError(f"{path}: duplicate canonical bin groups")
    return tuple(sorted(canonical))


def canonical_sha256(groups: tuple[tuple[str, ...], ...]) -> str:
    payload = "".join(
        "\t".join(group) + "\n"
        for group in groups
    ).encode("utf-8")
    return sha256_bytes(payload)


def audit(args: argparse.Namespace) -> dict[str, object]:
    baseline = read_classifier(
        args.baseline_output,
        selected_rank=args.rank,
    )
    candidate = read_classifier(args.candidate_output)
    baseline_taxonomy = read_taxonomy(args.baseline_taxonomy)
    candidate_taxonomy = read_taxonomy(args.candidate_taxonomy)
    baseline_groups = canonical_groups(args.baseline_assignments)
    candidate_groups = canonical_groups(args.candidate_assignments)

    baseline_ids = set(baseline)
    candidate_ids = set(candidate)
    common = sorted(baseline_ids & candidate_ids)
    classifier_differences = [
        {
            "contig_id": contig,
            "baseline": baseline[contig],
            "candidate": candidate[contig],
        }
        for contig in common
        if baseline[contig] != candidate[contig]
    ]
    status_agreement = sum(
        baseline[contig][0] == candidate[contig][0]
        for contig in common
    )
    taxid_agreement = sum(
        baseline[contig][2] == candidate[contig][2]
        for contig in common
    )
    full_line_agreement = sum(
        baseline[contig] == candidate[contig]
        for contig in common
    )

    taxonomy_ids = set(baseline_taxonomy) | set(candidate_taxonomy)
    taxonomy_differences = [
        {
            "contig_id": contig,
            "baseline": baseline_taxonomy.get(contig),
            "candidate": candidate_taxonomy.get(contig),
        }
        for contig in sorted(taxonomy_ids)
        if baseline_taxonomy.get(contig) != candidate_taxonomy.get(contig)
    ]
    baseline_taxonomy_bytes = args.baseline_taxonomy.read_bytes()
    candidate_taxonomy_bytes = args.candidate_taxonomy.read_bytes()
    groups_equal = baseline_groups == candidate_groups
    baseline_taxonomy_covers_classifier = (
        set(baseline_taxonomy) == baseline_ids
    )
    candidate_taxonomy_covers_classifier = (
        set(candidate_taxonomy) == candidate_ids
    )
    classifier_taxonomy_agreement = (
        baseline_ids == candidate_ids
        and not classifier_differences
        and not taxonomy_differences
        and baseline_taxonomy_covers_classifier
        and candidate_taxonomy_covers_classifier
    )
    complete = classifier_taxonomy_agreement and groups_equal
    return {
        "schema_version": 1,
        "rank": args.rank,
        "classifier_taxonomy_agreement": classifier_taxonomy_agreement,
        "complete_agreement": complete,
        "classifier": {
            "baseline_contigs": len(baseline),
            "candidate_contigs": len(candidate),
            "baseline_only": sorted(baseline_ids - candidate_ids),
            "candidate_only": sorted(candidate_ids - baseline_ids),
            "status_agreement": status_agreement,
            "taxid_agreement": taxid_agreement,
            "full_line_agreement": full_line_agreement,
            "differences": classifier_differences,
        },
        "taxonomy": {
            "baseline_contigs": len(baseline_taxonomy),
            "candidate_contigs": len(candidate_taxonomy),
            "byte_identical": (
                baseline_taxonomy_bytes == candidate_taxonomy_bytes
            ),
            "baseline_sha256": sha256_bytes(baseline_taxonomy_bytes),
            "candidate_sha256": sha256_bytes(candidate_taxonomy_bytes),
            "baseline_covers_classifier": (
                baseline_taxonomy_covers_classifier
            ),
            "candidate_covers_classifier": (
                candidate_taxonomy_covers_classifier
            ),
            "differences": taxonomy_differences,
        },
        "taxvamb": {
            "baseline_raw_bins": len(baseline_groups),
            "candidate_raw_bins": len(candidate_groups),
            "canonical_groups_equal": groups_equal,
            "baseline_canonical_sha256": canonical_sha256(baseline_groups),
            "candidate_canonical_sha256": canonical_sha256(candidate_groups),
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rank", required=True, type=int)
    parser.add_argument("--baseline-output", required=True, type=Path)
    parser.add_argument("--candidate-output", required=True, type=Path)
    parser.add_argument("--baseline-taxonomy", required=True, type=Path)
    parser.add_argument("--candidate-taxonomy", required=True, type=Path)
    parser.add_argument("--baseline-assignments", required=True, type=Path)
    parser.add_argument("--candidate-assignments", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--require-classifier-taxonomy-agreement",
        action="store_true",
    )
    parser.add_argument("--require-complete-agreement", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output.exists():
        raise ValueError(f"refusing to overwrite audit: {args.output}")
    result = audit(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if (
        args.require_classifier_taxonomy_agreement
        and not result["classifier_taxonomy_agreement"]
    ):
        return 2
    if args.require_complete_agreement and not result["complete_agreement"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
