from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

from mock_k2 import write_options, write_taxonomy


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"


class KunpengCompatibilityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.options = self.root / "opts.k2d"
        self.taxonomy = self.root / "taxo.k2d"
        write_options(self.options)
        write_taxonomy(self.taxonomy)
        self.sequence = (
            "ACGTTGCATGTCAGTACGATCGTAGCTAGCTACGATCGATGCTAGCATCGATCGTAGCTAGCAT"
        )
        self.fasta = self.root / "input.fa"
        self.fasta.write_text(
            f">r000001Cc1\n{self.sequence}\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_run_lengths_are_restored_exactly(self) -> None:
        raw = self.root / "raw.txt"
        raw.write_text(
            "C\tr000001Cc1\t2759\t64\t"
            "2759:2 0:2 2759:1 0:1 2759:2 0:1 2759:1\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from resolve_kunpeng_compat import resolve

            output = self.root / "compat.txt"
            summary = resolve(
                raw,
                self.fasta,
                self.options,
                self.taxonomy,
                output,
            )
        finally:
            sys.path.pop(0)
        self.assertEqual(summary["minimizer_events"], 10)
        self.assertEqual(summary["kmer_windows"], 30)
        self.assertEqual(summary["classified"], 1)
        self.assertEqual(summary["corrected_calls"], 0)
        fields = output.read_text(encoding="utf-8").rstrip().split("\t")
        self.assertEqual(fields[:4], ["C", "r000001Cc1", "2759", "64"])
        self.assertEqual(
            fields[4],
            "2759:4 0:10 2759:2 0:2 2759:7 0:2 2759:3",
        )

    def test_minimum_hit_groups_uses_nonzero_minimizer_events(self) -> None:
        raw = self.root / "raw-two-groups.txt"
        raw.write_text(
            "C\tr000001Cc1\t2759\t64\t2759:1 0:8 2759:1\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from resolve_kunpeng_compat import resolve

            output = self.root / "compat-two-groups.txt"
            summary = resolve(
                raw,
                self.fasta,
                self.options,
                self.taxonomy,
                output,
            )
        finally:
            sys.path.pop(0)
        self.assertEqual(summary["unclassified"], 1)
        self.assertEqual(summary["corrected_calls"], 1)
        self.assertTrue(
            output.read_text(encoding="utf-8").startswith(
                "U\tr000001Cc1\t0\t64\t"
            )
        )

    def test_event_count_mismatch_is_rejected(self) -> None:
        raw = self.root / "raw-short.txt"
        raw.write_text(
            "C\tr000001Cc1\t2759\t64\t2759:9\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from resolve_kunpeng_compat import resolve

            with self.assertRaisesRegex(ValueError, "hit events"):
                resolve(
                    raw,
                    self.fasta,
                    self.options,
                    self.taxonomy,
                    self.root / "should-not-exist.txt",
                )
        finally:
            sys.path.pop(0)

    def test_report_validated_zero_hit_omission_is_restored(self) -> None:
        second = "TGCA" * 20
        self.fasta.write_text(
            f">r000001Cc1\n{self.sequence}\n"
            f">r000001Cc2\n{second}\n",
            encoding="utf-8",
        )
        raw = self.root / "raw-omitted.txt"
        raw.write_text(
            "C\tr000001Cc1\t2759\t64\t"
            "2759:2 0:2 2759:1 0:1 2759:2 0:1 2759:1\n",
            encoding="utf-8",
        )
        report = self.root / "report.kreport2"
        report.write_text(
            " 50.00\t1\t1\tU\t0\tunclassified\n"
            " 50.00\t1\t0\tR\t1\troot\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from resolve_kunpeng_compat import resolve

            output = self.root / "compat-omitted.txt"
            summary = resolve(
                raw,
                self.fasta,
                self.options,
                self.taxonomy,
                output,
                classification_report=report,
            )
        finally:
            sys.path.pop(0)
        self.assertEqual(summary["raw_output_contigs"], 1)
        self.assertEqual(summary["raw_omitted_zero_hit_contigs"], 1)
        self.assertEqual(summary["unclassified"], 1)
        self.assertIn(
            "U\tr000001Cc2\t0\t80\t0:46\n",
            output.read_text(encoding="utf-8"),
        )

    def test_report_rejects_a_missing_classified_record(self) -> None:
        second = "TGCA" * 20
        self.fasta.write_text(
            f">r000001Cc1\n{self.sequence}\n"
            f">r000001Cc2\n{second}\n",
            encoding="utf-8",
        )
        raw = self.root / "raw-truncated.txt"
        raw.write_text(
            "C\tr000001Cc1\t2759\t64\t"
            "2759:2 0:2 2759:1 0:1 2759:2 0:1 2759:1\n",
            encoding="utf-8",
        )
        report = self.root / "truncated.kreport2"
        report.write_text(
            "100.00\t2\t0\tR\t1\troot\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from resolve_kunpeng_compat import resolve

            with self.assertRaisesRegex(ValueError, "classified counts differ"):
                resolve(
                    raw,
                    self.fasta,
                    self.options,
                    self.taxonomy,
                    self.root / "should-not-exist-truncated.txt",
                    classification_report=report,
                )
        finally:
            sys.path.pop(0)


if __name__ == "__main__":
    unittest.main()
