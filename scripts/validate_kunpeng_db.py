#!/usr/bin/env python3
"""Validate Kun-peng hash_config and every shard header/size."""

from __future__ import annotations

import argparse
import json
import re
import struct
from pathlib import Path
from typing import Optional, Sequence

from resolve_kunpeng_compat import K2Options, K2Taxonomy


def validate(
    database: Path,
    *,
    portable: bool = False,
) -> dict[str, object]:
    database = database.expanduser().resolve()
    for name in ("hash_config.k2d", "opts.k2d", "taxo.k2d"):
        path = database / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing or empty Kun-peng database file: {path}")
    index_options = K2Options(database / "opts.k2d")
    taxonomy = K2Taxonomy(database / "taxo.k2d")

    config_path = database / "hash_config.k2d"
    config = config_path.read_bytes()
    if len(config) != 48:
        raise ValueError(
            f"{config_path}: expected 48 bytes, observed {len(config)}"
        )
    (
        version,
        partitions,
        hash_capacity,
        capacity,
        size,
        value_bits,
    ) = struct.unpack("<6Q", config)
    if (
        partitions < 1
        or hash_capacity < 1
        or capacity < 1
        or size > capacity
        or not 1 <= value_bits <= 63
    ):
        raise ValueError(f"{config_path}: invalid configuration values")

    expected_names = {
        f"hash_{index}.k2d" for index in range(1, partitions + 1)
    }
    observed_names = {
        path.name
        for path in database.glob("hash_*.k2d")
        if re.fullmatch(r"hash_[1-9][0-9]*\.k2d", path.name)
    }
    if observed_names != expected_names:
        raise ValueError(
            "shard set mismatch: "
            f"missing={sorted(expected_names - observed_names)}, "
            f"extra={sorted(observed_names - expected_names)}"
        )

    observed_capacity = 0
    shard_bytes = 0
    shards: list[dict[str, int | str]] = []
    for index in range(1, partitions + 1):
        shard = database / f"hash_{index}.k2d"
        with shard.open("rb") as handle:
            header = handle.read(16)
        if len(header) != 16:
            raise ValueError(f"{shard}: short header")
        header_index, shard_capacity = struct.unpack("<2Q", header)
        expected_size = 16 + shard_capacity * 4
        if header_index != index:
            raise ValueError(
                f"{shard}: header index {header_index} != {index}"
            )
        if not 1 <= shard_capacity <= hash_capacity:
            raise ValueError(
                f"{shard}: invalid shard capacity {shard_capacity}"
            )
        if shard.stat().st_size != expected_size:
            raise ValueError(
                f"{shard}: size {shard.stat().st_size} != {expected_size}"
            )
        observed_capacity += shard_capacity
        shard_bytes += expected_size
        shards.append(
            {
                "name": shard.name,
                "index": index,
                "capacity": shard_capacity,
                "bytes": expected_size,
            }
        )
    if observed_capacity != capacity:
        raise ValueError(
            f"shard capacity sum {observed_capacity} != {capacity}"
        )
    return {
        "schema_version": 1,
        "database": "." if portable else str(database),
        "version": version,
        "partitions": partitions,
        "hash_capacity": hash_capacity,
        "capacity": capacity,
        "size": size,
        "value_bits": value_bits,
        "shard_bytes": shard_bytes,
        "index_options": index_options.as_dict(),
        "taxonomy_nodes": taxonomy.node_count,
        "shards": shards,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--portable",
        action="store_true",
        help="record the validated database path as '.' in the JSON report",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    report = validate(args.db, portable=args.portable)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
