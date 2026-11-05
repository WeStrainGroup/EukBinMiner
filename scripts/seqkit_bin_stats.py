#!/usr/bin/env python3
"""Inventory every per-sample bin with seqkit and merge QC candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import struct
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional, Sequence


CANDIDATE_COLUMNS = (
    "bin_rank",
    "bin_id",
    "study",
    "sample",
    "sample_rank",
    "binner",
    "original_bin_id",
    "contigs",
    "bp",
    "n50",
    "longest_contig",
    "sha256",
    "source_path",
    "bin_path",
    "materialization",
)
STATS_COLUMNS = CANDIDATE_COLUMNS + (
    "meets_min_bp",
    "seqkit_status",
    "seqkit_error",
)
FASTA_SUFFIXES = (".fa", ".fna", ".fasta")


def safe_component(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    if not cleaned:
        raise ValueError(f"cannot form a safe identifier from {value!r}")
    return cleaned


def original_bin_id(path: Path) -> str:
    name = path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
        lower = lower[:-3]
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"not a recognized FASTA path: {path}")


def is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return path.is_file() and any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(name)
    try:
        os.chmod(temporary, 0o660)
        with os.fdopen(
            descriptor, "wt", encoding="utf-8", newline="\n"
        ) as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_gzip_tsv(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[dict[str, object]],
) -> None:
    """Write deterministic gzip TSV atomically without relying on gzip.open fd."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    temporary = Path(name)
    try:
        os.chmod(temporary, 0o660)
        with temporary.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                mtime=0,
            ) as gzip_handle:
                header = "\t".join(columns) + "\n"
                gzip_handle.write(header.encode("utf-8"))
                for row in rows:
                    line = "\t".join(
                        str(row.get(column, ""))
                        .replace("\t", " ")
                        .replace("\r", " ")
                        .replace("\n", " ")
                        for column in columns
                    )
                    gzip_handle.write((line + "\n").encode("utf-8"))
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_manifest_index(manifest: Path, index_output: Path) -> None:
    """Write one fixed-width byte offset per sequential manifest rank."""
    offsets: list[int] = []
    with manifest.open("rb") as handle:
        header = handle.readline().decode("utf-8").rstrip("\r\n").split("\t")
        if "rank" not in header:
            raise ValueError(f"manifest lacks rank column: {manifest}")
        rank_column = header.index("rank")
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            fields = line.decode("utf-8").rstrip("\r\n").split("\t")
            expected = len(offsets) + 1
            if len(fields) != len(header) or int(fields[rank_column]) != expected:
                raise ValueError(
                    f"manifest rank is not sequential at row {expected}: {manifest}"
                )
            offsets.append(offset)
    payload = b"".join(struct.pack("<Q", offset) for offset in offsets)
    index_output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        dir=index_output.parent,
        prefix=f".{index_output.name}.",
        suffix=".tmp",
    )
    temporary = Path(name)
    try:
        os.chmod(temporary, 0o660)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, index_output)
    finally:
        temporary.unlink(missing_ok=True)
    atomic_write(
        index_output.with_suffix(index_output.suffix + ".json"),
        json.dumps(
            {
                "schema_version": 1,
                "manifest": str(manifest.resolve()),
                "manifest_sha256": sha256_file(manifest),
                "samples": len(offsets),
                "record_bytes": 8,
                "index_sha256": sha256_file(index_output),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )


def read_manifest_rank(
    manifest: Path,
    manifest_index: Path,
    rank: int,
) -> dict[str, str]:
    if rank < 1:
        raise ValueError("rank must be positive")
    with manifest_index.open("rb") as index_handle:
        index_handle.seek((rank - 1) * 8)
        packed = index_handle.read(8)
    if len(packed) != 8:
        raise ValueError(f"manifest rank {rank} is outside the index")
    offset = struct.unpack("<Q", packed)[0]
    with manifest.open("rb") as handle:
        header = handle.readline().decode("utf-8").rstrip("\r\n").split("\t")
        handle.seek(offset)
        values = handle.readline().decode("utf-8").rstrip("\r\n").split("\t")
    if len(values) != len(header):
        raise ValueError(f"invalid manifest row at rank {rank}")
    row = dict(zip(header, values))
    if int(row["rank"]) != rank:
        raise ValueError(f"manifest index rank mismatch: {row['rank']} != {rank}")
    return row


def parse_int(value: str) -> int:
    return int(value.replace(",", "").strip())


def _seqkit_once(
    seqkit: Path,
    paths: list[Path],
    threads: int,
    temporary_dir: Path,
    token: str,
) -> tuple[dict[Path, dict[str, str]], str]:
    list_path = temporary_dir / f"{token}.inputs.txt"
    list_path.write_text(
        "".join(f"{path.resolve()}\n" for path in paths),
        encoding="utf-8",
    )
    completed = subprocess.run(
        [
            str(seqkit),
            "stats",
            "-a",
            "-T",
            "-j",
            str(threads),
            "-X",
            str(list_path),
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        return {}, completed.stderr.strip() or f"seqkit exit {completed.returncode}"
    reader = csv.DictReader(completed.stdout.splitlines(), delimiter="\t")
    observed: dict[Path, dict[str, str]] = {}
    for row in reader:
        path = Path(row["file"]).resolve()
        if path in observed:
            return {}, f"duplicate seqkit result: {path}"
        observed[path] = row
    expected = {path.resolve() for path in paths}
    if set(observed) != expected:
        missing = sorted(str(path) for path in expected.difference(observed))
        extra = sorted(str(path) for path in set(observed).difference(expected))
        return (
            {},
            f"seqkit result mismatch; missing={missing[:3]} extra={extra[:3]}",
        )
    return observed, ""


def run_seqkit_resilient(
    seqkit: Path,
    paths: list[Path],
    threads: int,
    temporary_dir: Path,
) -> tuple[dict[Path, dict[str, str]], dict[Path, str]]:
    """Run all files together and bisect only failed batches."""
    observed: dict[Path, dict[str, str]] = {}
    errors: dict[Path, str] = {}
    serial = 0

    def visit(batch: list[Path]) -> None:
        nonlocal serial
        if not batch:
            return
        serial += 1
        results, error = _seqkit_once(
            seqkit,
            batch,
            threads,
            temporary_dir,
            f"batch{serial:06d}",
        )
        if not error:
            observed.update(results)
            return
        if len(batch) == 1:
            errors[batch[0].resolve()] = error[:1000]
            return
        midpoint = len(batch) // 2
        visit(batch[:midpoint])
        visit(batch[midpoint:])

    visit(paths)
    return observed, errors


def candidate_base(
    *,
    rank: int,
    study: str,
    sample: str,
    binner: str,
    path: Path,
) -> dict[str, object]:
    source = str(path.resolve())
    raw_id = original_bin_id(path)
    return {
        "bin_rank": 0,
        "bin_id": (
            f"r{rank:06d}__{safe_component(study)}__"
            f"{safe_component(sample)}__{safe_component(binner)}__"
            f"{safe_component(raw_id)}"
        ),
        "study": study,
        "sample": sample,
        "sample_rank": rank,
        "binner": binner,
        "original_bin_id": raw_id,
        "contigs": "",
        "bp": "",
        "n50": "",
        "longest_contig": "",
        "sha256": "",
        "source_path": source,
        "bin_path": source,
        "materialization": "source_reference",
    }


def inventory_sample(
    *,
    manifest: Path,
    manifest_index: Path,
    results_root: Path,
    rank: int,
    binners: Sequence[str],
    seqkit: Path,
    threads: int,
    minimum_bp: int,
    stats_output: Path,
    done_output: Path,
    candidate_output: Optional[Path] = None,
    temporary_root: Optional[Path] = None,
    require_complete_branches: bool = False,
    fail_on_seqkit_error: bool = False,
) -> None:
    sample_row = read_manifest_rank(manifest, manifest_index, rank)
    study = sample_row["study"]
    sample = sample_row["sample"]

    paths_by_binner: dict[str, list[Path]] = {}
    branch_state: dict[str, str] = {}
    for binner in binners:
        branch = (
            results_root
            / study
            / sample
            / "results"
            / "union"
            / binner
        )
        bins = branch / "bins"
        paths_by_binner[binner] = (
            sorted(path for path in bins.iterdir() if is_fasta(path))
            if bins.is_dir()
            else []
        )
        if not branch.exists():
            branch_state[binner] = "absent"
        elif not bins.is_dir():
            branch_state[binner] = "missing_bins_dir"
        elif (branch / "fingerprint.json").is_file():
            branch_state[binner] = "complete"
        else:
            branch_state[binner] = "unverified_output"

    if require_complete_branches:
        incomplete = {
            binner: state
            for binner, state in branch_state.items()
            if state != "complete"
        }
        if incomplete:
            raise ValueError(
                "sample does not have complete output from every requested "
                f"binner: {incomplete}"
            )

    all_paths = [
        path
        for binner in binners
        for path in paths_by_binner[binner]
    ]
    if temporary_root is None:
        temporary_root = stats_output.parent / ".seqkit-runtime"
    temporary_root.mkdir(parents=True, exist_ok=True)
    temp_parent = str(temporary_root.resolve())
    with tempfile.TemporaryDirectory(
        prefix=f"seqkit-bin-stats-r{rank:06d}.",
        dir=temp_parent,
    ) as temporary_name:
        observed, errors = run_seqkit_resilient(
            seqkit,
            all_paths,
            threads,
            Path(temporary_name),
        )

    stats_rows: list[dict[str, object]] = []
    by_binner: dict[str, Counter[str]] = defaultdict(Counter)
    for binner in binners:
        by_binner[binner]["discovered_bins"] += len(paths_by_binner[binner])
        for path in paths_by_binner[binner]:
            row = candidate_base(
                rank=rank,
                study=study,
                sample=sample,
                binner=binner,
                path=path,
            )
            resolved = path.resolve()
            if resolved in errors:
                row.update(
                    {
                        "meets_min_bp": "",
                        "seqkit_status": "error",
                        "seqkit_error": errors[resolved],
                    }
                )
                by_binner[binner]["seqkit_errors"] += 1
                stats_rows.append(row)
                continue
            values = observed[resolved]
            bp = parse_int(values["sum_len"])
            meets = int(bp >= minimum_bp)
            row.update(
                {
                    "contigs": parse_int(values["num_seqs"]),
                    "bp": bp,
                    "n50": parse_int(values["N50"]),
                    "longest_contig": parse_int(values["max_len"]),
                    "sha256": sha256_file(path) if meets else "",
                    "meets_min_bp": meets,
                    "seqkit_status": "ok",
                    "seqkit_error": "",
                }
            )
            by_binner[binner]["stats_ok"] += 1
            by_binner[binner]["raw_bp"] += bp
            if meets:
                by_binner[binner]["candidates"] += 1
                by_binner[binner]["candidate_bp"] += bp
            stats_rows.append(row)

    write_gzip_tsv(stats_output, STATS_COLUMNS, stats_rows)
    if errors and fail_on_seqkit_error:
        raise ValueError(
            "seqkit failed for one or more bin FASTA files: "
            + "; ".join(
                f"{path}: {message}"
                for path, message in sorted(
                    errors.items(),
                    key=lambda item: str(item[0]),
                )
            )[:2000]
        )

    candidate_rows = [
        {
            column: row.get(column, "")
            for column in CANDIDATE_COLUMNS
        }
        for row in stats_rows
        if row.get("meets_min_bp") == 1
    ]
    for bin_rank, row in enumerate(candidate_rows, start=1):
        row["bin_rank"] = bin_rank
    if candidate_output is not None:
        candidate_lines = ["\t".join(CANDIDATE_COLUMNS) + "\n"]
        candidate_lines.extend(
            "\t".join(str(row.get(column, "")) for column in CANDIDATE_COLUMNS)
            + "\n"
            for row in candidate_rows
        )
        atomic_write(candidate_output, "".join(candidate_lines))

    payload = {
        "schema_version": 1,
        "rank": rank,
        "study": study,
        "sample": sample,
        "minimum_bp_inclusive": minimum_bp,
        "binners": list(binners),
        "branch_state": branch_state,
        "discovered_bins": len(all_paths),
        "stats_ok": len(observed),
        "seqkit_errors": len(errors),
        "candidates": sum(
            values["candidates"] for values in by_binner.values()
        ),
        "by_binner": {
            binner: dict(by_binner[binner]) for binner in binners
        },
        "stats_path": str(stats_output.resolve()),
        "stats_sha256": sha256_file(stats_output),
    }
    if candidate_output is not None:
        payload.update(
            {
                "candidate_manifest": str(candidate_output.resolve()),
                "candidate_manifest_sha256": sha256_file(candidate_output),
            }
        )
    atomic_write(
        done_output,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )


def merge_stats(
    *,
    fragment_dir: Path,
    first_rank: int,
    last_rank: int,
    minimum_bp: int,
    inventory_output: Path,
    candidate_output: Path,
    summary_output: Path,
    done_output: Path,
    allow_missing: bool,
) -> None:
    if first_rank < 1 or last_rank < first_rank:
        raise ValueError("invalid rank range")
    candidate_rows: list[dict[str, object]] = []
    valid_fragments: list[Path] = []
    totals: Counter[str] = Counter()
    by_binner: dict[str, Counter[str]] = defaultdict(Counter)
    branch_state: dict[str, Counter[str]] = defaultdict(Counter)
    missing: list[int] = []
    invalid: list[int] = []

    for rank in range(first_rank, last_rank + 1):
        stem = f"r{rank:06d}"
        stats_path = fragment_dir / f"{stem}.stats.tsv.gz"
        done_path = fragment_dir / f"{stem}.done"
        if not stats_path.is_file() or not done_path.is_file():
            missing.append(rank)
            continue
        try:
            marker = json.loads(done_path.read_text(encoding="utf-8"))
            if marker["rank"] != rank:
                raise ValueError("rank mismatch")
            if marker["minimum_bp_inclusive"] != minimum_bp:
                raise ValueError("minimum bp mismatch")
            if marker["stats_sha256"] != sha256_file(stats_path):
                raise ValueError("stats checksum mismatch")
            rows = read_tsv(stats_path)
        except (OSError, KeyError, ValueError, json.JSONDecodeError):
            invalid.append(rank)
            continue

        valid_fragments.append(stats_path)
        totals["samples"] += 1
        for state_binner, state in marker.get("branch_state", {}).items():
            branch_state[state_binner][state] += 1
        for row in rows:
            totals["discovered_bins"] += 1
            binner = row["binner"]
            by_binner[binner]["discovered_bins"] += 1
            if row["seqkit_status"] == "ok":
                totals["stats_ok"] += 1
                by_binner[binner]["stats_ok"] += 1
                bp = int(row["bp"])
                totals["raw_bp"] += bp
                by_binner[binner]["raw_bp"] += bp
            else:
                totals["seqkit_errors"] += 1
                by_binner[binner]["seqkit_errors"] += 1
            if row["meets_min_bp"] == "1":
                candidate = {
                    column: row.get(column, "")
                    for column in CANDIDATE_COLUMNS
                }
                candidate_rows.append(candidate)
                by_binner[binner]["candidates"] += 1
                by_binner[binner]["candidate_bp"] += int(row["bp"])

    if (missing or invalid) and not allow_missing:
        raise ValueError(
            f"incomplete fragments: missing={missing[:20]} invalid={invalid[:20]}"
        )

    for candidate_rank, row in enumerate(candidate_rows, start=1):
        row["bin_rank"] = candidate_rank

    def inventory_records() -> Iterable[dict[str, object]]:
        raw_rank = 0
        for stats_path in valid_fragments:
            for row in read_tsv(stats_path):
                raw_rank += 1
                record: dict[str, object] = dict(row)
                record["bin_rank"] = raw_rank
                yield record

    write_gzip_tsv(inventory_output, STATS_COLUMNS, inventory_records())
    candidate_lines = ["\t".join(CANDIDATE_COLUMNS) + "\n"]
    candidate_lines.extend(
        "\t".join(str(row.get(column, "")) for column in CANDIDATE_COLUMNS)
        + "\n"
        for row in candidate_rows
    )
    atomic_write(candidate_output, "".join(candidate_lines))

    summary = {
        "schema_version": 1,
        "mode": "all_output_seqkit_stats",
        "first_rank": first_rank,
        "last_rank": last_rank,
        "expected_samples": last_rank - first_rank + 1,
        "completed_samples": totals["samples"],
        "missing_ranks": missing,
        "invalid_ranks": invalid,
        "minimum_bp_inclusive": minimum_bp,
        "discovered_bins": totals["discovered_bins"],
        "stats_ok": totals["stats_ok"],
        "seqkit_errors": totals["seqkit_errors"],
        "raw_bp": totals["raw_bp"],
        "candidates": len(candidate_rows),
        "candidate_bp": sum(int(row["bp"]) for row in candidate_rows),
        "by_binner": {
            binner: dict(values)
            for binner, values in sorted(by_binner.items())
        },
        "branch_state": {
            binner: dict(values)
            for binner, values in sorted(branch_state.items())
        },
        "inventory": str(inventory_output.resolve()),
        "inventory_sha256": sha256_file(inventory_output),
        "candidate_manifest": str(candidate_output.resolve()),
        "candidate_manifest_sha256": sha256_file(candidate_output),
    }
    atomic_write(
        summary_output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    atomic_write(
        done_output,
        json.dumps(
            {
                "inventory_sha256": sha256_file(inventory_output),
                "candidate_manifest_sha256": sha256_file(candidate_output),
                "summary_sha256": sha256_file(summary_output),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    index = subparsers.add_parser("index")
    index.add_argument("--manifest", required=True, type=Path)
    index.add_argument("--index-output", required=True, type=Path)

    sample = subparsers.add_parser("sample")
    sample.add_argument("--manifest", required=True, type=Path)
    sample.add_argument("--manifest-index", required=True, type=Path)
    sample.add_argument("--results-root", required=True, type=Path)
    sample.add_argument("--rank", required=True, type=int)
    sample.add_argument("--binner", action="append", required=True)
    sample.add_argument("--seqkit", required=True, type=Path)
    sample.add_argument("--threads", required=True, type=int)
    sample.add_argument("--min-bp", required=True, type=int)
    sample.add_argument("--stats-output", required=True, type=Path)
    sample.add_argument("--done-output", required=True, type=Path)
    sample.add_argument("--candidate-output", type=Path)
    sample.add_argument("--temporary-root", type=Path)
    sample.add_argument("--require-complete-branches", action="store_true")
    sample.add_argument("--fail-on-seqkit-error", action="store_true")

    merge = subparsers.add_parser("merge")
    merge.add_argument("--fragment-dir", required=True, type=Path)
    merge.add_argument("--first-rank", required=True, type=int)
    merge.add_argument("--last-rank", required=True, type=int)
    merge.add_argument("--min-bp", required=True, type=int)
    merge.add_argument("--inventory-output", required=True, type=Path)
    merge.add_argument("--candidate-output", required=True, type=Path)
    merge.add_argument("--summary-output", required=True, type=Path)
    merge.add_argument("--done-output", required=True, type=Path)
    merge.add_argument("--allow-missing", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "index":
        build_manifest_index(args.manifest, args.index_output)
        return 0
    if args.min_bp < 1:
        raise ValueError("--min-bp must be positive")
    if args.command == "sample":
        if args.threads < 1:
            raise ValueError("--threads must be positive")
        inventory_sample(
            manifest=args.manifest,
            manifest_index=args.manifest_index,
            results_root=args.results_root,
            rank=args.rank,
            binners=args.binner,
            seqkit=args.seqkit,
            threads=args.threads,
            minimum_bp=args.min_bp,
            stats_output=args.stats_output,
            done_output=args.done_output,
            candidate_output=args.candidate_output,
            temporary_root=args.temporary_root,
            require_complete_branches=args.require_complete_branches,
            fail_on_seqkit_error=args.fail_on_seqkit_error,
        )
    else:
        merge_stats(
            fragment_dir=args.fragment_dir,
            first_rank=args.first_rank,
            last_rank=args.last_rank,
            minimum_bp=args.min_bp,
            inventory_output=args.inventory_output,
            candidate_output=args.candidate_output,
            summary_output=args.summary_output,
            done_output=args.done_output,
            allow_missing=args.allow_missing,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
