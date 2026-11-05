from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"


class ClassifierComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.report = self.root / "report.kreport2"
        self.report.write_text(
            "100.00\t2\t0\tR\t1\troot\n"
            "100.00\t2\t0\tD\t2759\t  Eukaryota\n"
            "50.00\t1\t1\tK\t4751\t    Fungi\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_exact_outputs_have_no_differences(self) -> None:
        output = self.root / "calls.txt"
        output.write_text(
            "C\tr000001Cc1\t4751\t100\tA:100\n"
            "C\tr000001Cc2\t2759\t100\tA:100\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from compare_classifier_outputs import compare

            summary = compare(
                output,
                self.report,
                output,
                self.report,
                self.root / "comparison",
            )
        finally:
            sys.path.pop(0)
        self.assertEqual(summary["differences"], 0)
        self.assertEqual(summary["fungi_agreement_fraction"], 1.0)

    def test_taxid_difference_is_reported(self) -> None:
        baseline = self.root / "baseline.txt"
        candidate = self.root / "candidate.txt"
        baseline.write_text(
            "C\tr000001Cc1\t4751\t100\tA:100\n",
            encoding="utf-8",
        )
        candidate.write_text(
            "C\tr000001Cc1\t2759\t100\tA:100\n",
            encoding="utf-8",
        )
        sys.path.insert(0, str(SCRIPTS))
        try:
            from compare_classifier_outputs import compare

            summary = compare(
                baseline,
                self.report,
                candidate,
                self.report,
                self.root / "comparison-different",
            )
        finally:
            sys.path.pop(0)
        self.assertEqual(summary["differences"], 1)
        self.assertEqual(summary["taxid_agreement"], 0)
        self.assertEqual(summary["domain_agreement"], 1)
        self.assertEqual(summary["fungi_agreement"], 0)


if __name__ == "__main__":
    unittest.main()
