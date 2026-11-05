from __future__ import annotations

import struct
from pathlib import Path


def write_options(path: Path) -> None:
    path.write_bytes(
        struct.pack(
            "<QQQQ?7xQiii4x",
            35,
            31,
            0x3FFFFFFFF3333333,
            0xE37E28C4271B5A2D,
            True,
            0,
            1,
            0,
            0,
        )
    )


def write_taxonomy(path: Path) -> None:
    names = b"\0root\0Eukaryota\0"
    ranks = b"no rank\0superkingdom\0"
    nodes = (
        (0, 0, 0, 0, 0, 0, 0),
        (0, 2, 1, 1, 0, 1, 0),
        (1, 3, 0, 6, 8, 2759, 0),
    )
    payload = (
        b"K2TAXDAT"
        + struct.pack("<QQQ", len(nodes), len(names), len(ranks))
        + b"".join(struct.pack("<7Q", *node) for node in nodes)
        + names
        + ranks
    )
    path.write_bytes(payload)
