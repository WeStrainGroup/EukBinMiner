#!/usr/bin/env python3
"""Build one rank-prefixed Kraken2 catalogue from prepared VAMB inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Iterator, Optional, Sequence


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    identifier: Optional[str] = None
    sequence: list[str] = []
    seen: set[str] = set()
    with path.open("rt", encoding="utf-8", newline=None) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence)
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"{path}: empty header at line {line_number}")
                identifier = fields[0]
                if identifier in seen:
                    raise ValueError(f"{path}: duplicate FASTA ID {identifier}")
                seen.add(identifier)
                sequence = []
            elif identifier is None:
                raise ValueError(
                    f"{path}: sequence before first header at line {line_number}"
                )
            else:
                sequence.append("".join(line.split()))
    if identifier is not None:
        yield identifier, "".join(sequence)


def read_manifest(path: Path, max_rank: Optional[int]) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if max_rank is not None:
        rows = [row for row in rows if int(row["rank"]) <= max_rank]
    if not rows:
        raise ValueError("no manifest rows selected")
    return rows


def build(
    manifest: Path,
    inputs_root: Path,
    out_dir: Path,
    max_rank: Optional[int],
) -> str:
    manifest_sha = sha256_file(manifest)
    if out_dir.exists():
        complete = out_dir / "catalogue.json"
        if not complete.is_file():
            raise ValueError(f"partial Kraken catalogue exists: {out_dir}")
        metadata = json.loads(complete.read_text(encoding="utf-8"))
        if (
            metadata.get("manifest_sha256") != manifest_sha
            or metadata.get("max_rank") != max_rank
        ):
            raise ValueError("existing Kraken catalogue selection differs")
        for name in ("catalogue.fa", "samples.tsv"):
            path = out_dir / name
            if (
                not path.is_file()
                or sha256_file(path) != metadata["outputs"][name]["sha256"]
            ):
                raise ValueError(f"existing Kraken catalogue differs: {path}")
        return "already-current"

    rows = read_manifest(manifest, max_rank)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{out_dir.name}.", dir=out_dir.parent)
    )
    total_contigs = 0
    total_bp = 0
    try:
        with (stage / "catalogue.fa").open(
            "wt", encoding="utf-8", newline="\n"
        ) as fasta_out, (stage / "samples.tsv").open(
            "wt", encoding="utf-8", newline="\n"
        ) as samples_out:
            samples_out.write(
                "rank\tstudy\tsample\tprefix\tcontigs\tbp\t"
                "input_path\tabundance_path\n"
            )
            for row in rows:
                rank = int(row["rank"])
                input_dir = inputs_root / f"r{rank:06d}"
                metadata_path = input_dir / "input.json"
                if not metadata_path.is_file():
                    raise FileNotFoundError(metadata_path)
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                prefix = f"r{rank:06d}C"
                observed_contigs = 0
                observed_bp = 0
                for identifier, sequence in fasta_records(input_dir / "input.fa"):
                    if not sequence:
                        raise ValueError(
                            f"{input_dir / 'input.fa'}: empty sequence {identifier}"
                        )
                    fasta_out.write(f">{prefix}{identifier}\n")
                    for start in range(0, len(sequence), 80):
                        fasta_out.write(sequence[start : start + 80] + "\n")
                    observed_contigs += 1
                    observed_bp += len(sequence)
                if (
                    observed_contigs != int(metadata["contigs"])
                    or observed_bp != int(metadata["bp"])
                ):
                    raise ValueError(
                        f"prepared input metadata mismatch for rank {rank}"
                    )
                samples_out.write(
                    "\t".join(
                        (
                            str(rank),
                            row["study"],
                            row["sample"],
                            prefix,
                            str(observed_contigs),
                            str(observed_bp),
                            str((input_dir / "input.fa").resolve()),
                            str((input_dir / "abundance.tsv").resolve()),
                        )
                    )
                    + "\n"
                )
                total_contigs += observed_contigs
                total_bp += observed_bp
        if total_contigs == 0:
            raise ValueError("Kraken catalogue has no contigs")
        metadata = {
            "schema_version": 1,
            "manifest": str(manifest.resolve()),
            "manifest_sha256": manifest_sha,
            "max_rank": max_rank,
            "samples": len(rows),
            "contigs": total_contigs,
            "bp": total_bp,
            "outputs": {
                name: {
                    "sha256": sha256_file(stage / name),
                    "bytes": (stage / name).stat().st_size,
                }
                for name in ("catalogue.fa", "samples.tsv")
            },
        }
        (stage / "catalogue.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(stage, out_dir)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return "created"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--inputs-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--max-rank", type=int)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    status = build(
        args.manifest, args.inputs_root, args.out_dir, args.max_rank
    )
    print(f"build_kraken_catalogue.py: {status}: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
