#!/usr/bin/env python3
"""Write a deterministic branch fingerprint from named values and input files."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence


class FingerprintError(ValueError):
    """An invalid field, file, or output-safety condition."""


def _parse_named_values(values: Sequence[str], option: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise FingerprintError(f"{option} requires NAME=VALUE, got: {value!r}")
        name, item = value.split("=", 1)
        if not name or any(character in name for character in "=\t\r\n"):
            raise FingerprintError(f"invalid name in {option}: {name!r}")
        if name in parsed:
            raise FingerprintError(f"duplicate name in {option}: {name}")
        parsed[name] = item
    return parsed


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise FingerprintError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def build_fingerprint(
    fields: Sequence[str],
    files: Sequence[str],
) -> dict:
    parsed_fields = _parse_named_values(fields, "--field")
    parsed_files = _parse_named_values(files, "--file")
    file_records: dict[str, dict[str, str]] = {}
    for name, path_text in parsed_files.items():
        path = Path(path_text).resolve()
        if not path.is_file():
            raise FingerprintError(f"fingerprint file does not exist: {path}")
        file_records[name] = {
            "path": str(path),
            "sha256": _sha256(path),
        }

    payload = {
        "schema_version": 1,
        "fields": parsed_fields,
        "files": file_records,
    }
    canonical = json.dumps(
        payload,
        sort_keys=True,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")
    payload["fingerprint_sha256"] = hashlib.sha256(canonical).hexdigest()
    return payload


def write_fingerprint(
    fields: Sequence[str],
    files: Sequence[str],
    output: Path,
) -> None:
    payload = build_fingerprint(fields, files)
    rendered = (
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    _atomic_write(output, rendered)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--field",
        action="append",
        default=[],
        metavar="NAME=VALUE",
    )
    parser.add_argument(
        "--file",
        action="append",
        default=[],
        metavar="NAME=PATH",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        write_fingerprint(args.field, args.file, args.output)
    except FingerprintError as exc:
        print(f"branch_fingerprint.py: error: {exc}", file=sys.stderr)
        return 2
    print(f"branch_fingerprint.py: wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
