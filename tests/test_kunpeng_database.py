from __future__ import annotations

import struct
import tempfile
import unittest
from pathlib import Path

from mock_k2 import write_options, write_taxonomy


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"


class KunpengDatabaseValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.database = Path(self.temporary.name) / "database"
        self.database.mkdir()
        write_options(self.database / "opts.k2d")
        write_taxonomy(self.database / "taxo.k2d")
        (self.database / "hash_config.k2d").write_bytes(
            struct.pack("<6Q", 0, 2, 3, 5, 4, 16)
        )
        (self.database / "hash_1.k2d").write_bytes(
            struct.pack("<2Q", 1, 3) + b"\0" * 12
        )
        (self.database / "hash_2.k2d").write_bytes(
            struct.pack("<2Q", 2, 2) + b"\0" * 8
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_valid_database_reports_exact_layout(self) -> None:
        import sys

        sys.path.insert(0, str(SCRIPTS))
        try:
            from validate_kunpeng_db import validate

            report = validate(self.database)
        finally:
            sys.path.pop(0)
        self.assertEqual(report["partitions"], 2)
        self.assertEqual(report["capacity"], 5)
        self.assertEqual(report["shard_bytes"], 52)
        self.assertEqual(report["index_options"]["k"], 35)
        self.assertEqual(report["taxonomy_nodes"], 3)
        self.assertEqual(report["database"], str(self.database.resolve()))

    def test_portable_layout_does_not_retain_staging_path(self) -> None:
        import sys

        sys.path.insert(0, str(SCRIPTS))
        try:
            from validate_kunpeng_db import validate

            report = validate(self.database, portable=True)
        finally:
            sys.path.pop(0)
        self.assertEqual(report["database"], ".")

    def test_wrong_shard_size_is_rejected(self) -> None:
        import sys

        (self.database / "hash_2.k2d").write_bytes(
            struct.pack("<2Q", 2, 2) + b"\0" * 4
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from validate_kunpeng_db import validate

            with self.assertRaisesRegex(ValueError, "size"):
                validate(self.database)
        finally:
            sys.path.pop(0)


if __name__ == "__main__":
    unittest.main()
