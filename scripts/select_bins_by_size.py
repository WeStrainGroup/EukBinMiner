#!/usr/bin/env python3
"""Select candidate bin FASTAs by total size while preserving raw bins."""

from __future__ import annotations

import argparse
import errno
import gzip
import os
import shutil
import sys
from pathlib import Path
from typing import BinaryIO, Optional, Sequence


FASTA_SUFFIXES = (".fa", ".fasta", ".fna")


class SelectionError(ValueError):
    """A bin discovery, FASTA validation, or output-safety failure."""


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return path.is_file() and any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def _open_maybe_gzip(path: Path) -> BinaryIO:
    try:
        with path.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"\x1f\x8b":
            return gzip.open(path, "rb")
        return path.open("rb")
    except OSError as exc:
        raise SelectionError(f"cannot open bin FASTA {path}: {exc}") from exc


def _fasta_stats(path: Path) -> tuple[int, int]:
    contigs = 0
    total_bp = 0
    current_id: Optional[bytes] = None
    current_length = 0

    def finish_record() -> None:
        nonlocal total_bp
        if current_id is None:
            return
        if current_length == 0:
            raise SelectionError(
                f"empty sequence for {current_id.decode('utf-8', 'replace')} in {path}"
            )
        total_bp += current_length

    try:
        with _open_maybe_gzip(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if raw_line.startswith(b">"):
                    finish_record()
                    fields = raw_line[1:].strip().split(None, 1)
                    if not fields:
                        raise SelectionError(
                            f"empty FASTA header in {path} at line {line_number}"
                        )
                    current_id = fields[0]
                    current_length = 0
                    contigs += 1
                else:
                    if current_id is None:
                        raise SelectionError(
                            f"content before first FASTA header in {path} at "
                            f"line {line_number}"
                        )
                    current_length += sum(
                        1
                        for byte in raw_line.rstrip(b"\r\n")
                        if not chr(byte).isspace()
                    )
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise SelectionError(f"cannot read bin FASTA {path}: {exc}") from exc

    finish_record()
    if contigs == 0:
        raise SelectionError(f"bin FASTA is empty: {path}")
    return contigs, total_bp


def _link_or_copy(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, destination)
        return "hardlink"
    except OSError as exc:
        if exc.errno not in {
            errno.EXDEV,
            errno.EPERM,
            errno.EACCES,
            errno.ENOTSUP,
        }:
            raise SelectionError(
                f"cannot link selected bin {source} to {destination}: {exc}"
            ) from exc
    try:
        shutil.copy2(source, destination)
    except OSError as exc:
        raise SelectionError(
            f"cannot copy selected bin {source} to {destination}: {exc}"
        ) from exc
    return "copy"


def select_bins(
    input_dir: Path,
    output_dir: Path,
    manifest: Path,
    min_bp: int,
) -> tuple[int, int]:
    if min_bp < 1:
        raise SelectionError(f"--min-bp must be positive, got {min_bp}")
    if not input_dir.is_dir():
        raise SelectionError(f"input bin directory does not exist: {input_dir}")

    input_resolved = input_dir.resolve()
    output_resolved = output_dir.resolve()
    if output_resolved == input_resolved or input_resolved in output_resolved.parents:
        raise SelectionError("--output-dir must not be inside --input-dir")
    if output_dir.exists() and not output_dir.is_dir():
        raise SelectionError(f"output path is not a directory: {output_dir}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SelectionError(f"refusing to reuse non-empty output directory: {output_dir}")
    if manifest.exists():
        raise SelectionError(f"refusing to overwrite manifest: {manifest}")

    paths = sorted(
        (path for path in input_dir.rglob("*") if _is_fasta(path)),
        key=lambda path: path.relative_to(input_dir).as_posix(),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "bin_file\tcontigs\tbp\tselected\tmaterialization\tselected_file\n"
    ]
    selected_count = 0

    for source in paths:
        relative = source.relative_to(input_dir)
        contigs, bp = _fasta_stats(source)
        if bp >= min_bp:
            destination = output_dir / relative
            method = _link_or_copy(source, destination)
            selected = "yes"
            selected_file = (Path("bins") / relative).as_posix()
            selected_count += 1
        else:
            method = ""
            selected = "no"
            selected_file = ""
        lines.append(
            f"{relative.as_posix()}\t{contigs}\t{bp}\t{selected}\t"
            f"{method}\t{selected_file}\n"
        )

    try:
        manifest.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise SelectionError(f"cannot write selection manifest {manifest}: {exc}") from exc
    return selected_count, len(paths)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--min-bp", required=True, type=int)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        selected, total = select_bins(
            args.input_dir,
            args.output_dir,
            args.manifest,
            args.min_bp,
        )
    except SelectionError as exc:
        print(f"select_bins_by_size.py: error: {exc}", file=sys.stderr)
        return 2
    print(
        f"select_bins_by_size.py: selected {selected}/{total} bins "
        f"at >= {args.min_bp} bp"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
