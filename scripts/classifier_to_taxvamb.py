#!/usr/bin/env python3
"""Convert Kraken-compatible sequence/report output to TaxVAMB taxonomies."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import IO, Optional, Sequence


RANK_CODES = ("D", "P", "C", "O", "F", "G", "S")
RANK_PREFIXES = ("d", "p", "c", "o", "f", "g", "s")
PREFIX_RE = re.compile(r"^r([0-9]{6})C(.+)$")


def open_text(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def clean_name(value: str) -> str:
    cleaned = " ".join(value.replace(";", " ").replace("\t", " ").split())
    return cleaned or "unnamed"


def format_lineage(standard: dict[str, tuple[str, str]]) -> str:
    known = [standard.get(code) for code in RANK_CODES]
    deepest = max(
        (index for index, value in enumerate(known) if value is not None),
        default=-1,
    )
    if deepest < 0:
        return ""
    output: list[str] = []
    for index in range(deepest + 1):
        value = known[index]
        prefix = RANK_PREFIXES[index]
        if value is None:
            next_taxid = next(
                known[future][1]
                for future in range(index + 1, deepest + 1)
                if known[future] is not None
            )
            output.append(f"{prefix}__unclassified_for_{next_taxid}")
        else:
            name, taxid = value
            output.append(f"{prefix}__{clean_name(name)}__taxid_{taxid}")
    return ";".join(output)


def parse_report(path: Path) -> dict[str, str]:
    lineages: dict[str, str] = {}
    stack: list[tuple[int, dict[str, tuple[str, str]]]] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t", 5)
            if len(fields) != 6:
                raise ValueError(
                    f"{path}: report line {line_number} has {len(fields)} fields"
                )
            rank_code = fields[3]
            taxid = fields[4]
            raw_name = fields[5]
            indent = len(raw_name) - len(raw_name.lstrip(" "))
            name = raw_name.strip()
            while stack and stack[-1][0] >= indent:
                stack.pop()
            standard = dict(stack[-1][1]) if stack else {}
            if rank_code in RANK_CODES:
                standard[rank_code] = (name, taxid)
            lineages[taxid] = format_lineage(standard)
            stack.append((indent, standard))
    return lineages


def read_samples(path: Path) -> dict[int, dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = {int(row["rank"]): row for row in rows}
    if len(samples) != len(rows):
        raise ValueError("samples.tsv contains duplicate ranks")
    return samples


def convert(
    classification_output: Path,
    classification_report: Path,
    samples_tsv: Path,
    out_dir: Path,
    *,
    classification_taxonomy: Optional[Path] = None,
) -> None:
    if out_dir.exists():
        raise ValueError(f"TaxVAMB taxonomy output already exists: {out_dir}")
    samples = read_samples(samples_tsv)
    lineages = (
        parse_report(classification_report)
        if classification_taxonomy is None
        else None
    )
    taxonomy = None
    if classification_taxonomy is not None:
        from resolve_kunpeng_compat import K2Taxonomy

        taxonomy = K2Taxonomy(classification_taxonomy)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{out_dir.name}.", dir=out_dir.parent)
    )
    per_sample = stage / "per_sample"
    per_sample.mkdir()
    counts: Counter[str] = Counter()
    sample_counts: Counter[int] = Counter()
    observed_ids: dict[int, set[str]] = {
        rank: set() for rank in samples
    }
    current_rank: Optional[int] = None
    current_handle: Optional[IO[str]] = None
    try:
        with open_text(classification_output) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(
                        f"{classification_output}: malformed line {line_number}"
                    )
                status, combined_id, taxid = fields[:3]
                if status not in {"C", "U"}:
                    raise ValueError(
                        f"{classification_output}: invalid status {status!r} "
                        f"at line {line_number}"
                    )
                match = PREFIX_RE.match(combined_id)
                if match is None:
                    raise ValueError(
                        f"unexpected combined contig ID: {combined_id}"
                    )
                rank = int(match.group(1))
                contig_id = match.group(2)
                if rank not in samples:
                    raise ValueError(
                        f"classification output contains unknown rank {rank}"
                    )
                if contig_id in observed_ids[rank]:
                    raise ValueError(
                        f"duplicate classification for rank {rank} "
                        f"contig {contig_id!r}"
                    )
                observed_ids[rank].add(contig_id)
                if current_rank is not None and rank < current_rank:
                    raise ValueError(
                        "classification output is not grouped by sample rank"
                    )
                if rank != current_rank:
                    if current_handle is not None:
                        current_handle.close()
                    current_rank = rank
                    current_handle = (
                        per_sample / f"r{rank:06d}.taxonomy.tsv"
                    ).open("wt", encoding="utf-8", newline="\n")
                    current_handle.write("contigs\tpredictions\n")
                if status == "C" and taxonomy is not None:
                    if not taxid.isdigit():
                        raise ValueError(
                            f"{classification_output}: invalid taxid "
                            f"{taxid!r} at line {line_number}"
                        )
                    lineage = taxonomy.taxvamb_lineage(int(taxid))
                elif status == "C":
                    assert lineages is not None
                    lineage = lineages.get(taxid, "")
                else:
                    lineage = ""
                current_handle.write(f"{contig_id}\t{lineage}\n")
                sample_counts[rank] += 1
                counts["classified" if lineage else "unclassified"] += 1
                domain = lineage.split(";", 1)[0] if lineage else "unclassified"
                counts[f"domain:{domain}"] += 1
        if current_handle is not None:
            current_handle.close()
            current_handle = None
        for rank, row in samples.items():
            expected = int(row["contigs"])
            observed = sample_counts[rank]
            taxonomy_path = per_sample / f"r{rank:06d}.taxonomy.tsv"
            if expected == 0:
                taxonomy_path.write_text(
                    "contigs\tpredictions\n",
                    encoding="utf-8",
                )
            if observed != expected:
                raise ValueError(
                    f"rank {rank}: taxonomy rows {observed} != contigs {expected}"
                )
        summary = {
            "schema_version": 1,
            "samples": len(samples),
            "contigs": sum(sample_counts.values()),
            "classified": counts["classified"],
            "unclassified": counts["unclassified"],
            "domain_counts": {
                key.removeprefix("domain:"): value
                for key, value in sorted(counts.items())
                if key.startswith("domain:")
            },
            "classification_output": str(classification_output.resolve()),
            "classification_report": str(classification_report.resolve()),
            "classification_taxonomy": (
                str(classification_taxonomy.resolve())
                if classification_taxonomy is not None
                else None
            ),
            "samples_tsv": str(samples_tsv.resolve()),
        }
        (stage / "taxonomy.summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(stage, out_dir)
    except Exception:
        if current_handle is not None:
            current_handle.close()
        shutil.rmtree(stage, ignore_errors=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classification-output", required=True, type=Path)
    parser.add_argument("--classification-report", required=True, type=Path)
    parser.add_argument("--classification-taxonomy", type=Path)
    parser.add_argument("--samples-tsv", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    convert(
        args.classification_output,
        args.classification_report,
        args.samples_tsv,
        args.out_dir,
        classification_taxonomy=args.classification_taxonomy,
    )
    print(f"classifier_to_taxvamb.py: wrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
