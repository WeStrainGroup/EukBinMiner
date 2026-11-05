from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

from audit_kunpeng_calibration import audit  # noqa: E402


class KunpengCalibrationAuditTests(unittest.TestCase):
    def test_rank_prefix_and_bin_ids_are_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline_output = root / "kraken.output.gz"
            candidate_output = root / "kunpeng.compat.output.txt.gz"
            with gzip.open(baseline_output, "wt", encoding="utf-8") as handle:
                handle.write("C\tr000981Ca\t2\t100\t2:66\n")
                handle.write("U\tr000040Cother\t0\t100\t0:66\n")
            with gzip.open(candidate_output, "wt", encoding="utf-8") as handle:
                handle.write("C\tr000001Ca\t2\t100\t2:66\n")
            taxonomy = "contigs\tpredictions\na\td__Eukaryota\n"
            baseline_taxonomy = root / "baseline.taxonomy.tsv"
            candidate_taxonomy = root / "candidate.taxonomy.tsv"
            baseline_taxonomy.write_text(taxonomy, encoding="utf-8")
            candidate_taxonomy.write_text(taxonomy, encoding="utf-8")
            header = "contig_id\tlength_bp\tstatus\tbin_id\tbin_file\n"
            baseline_assignments = root / "baseline.assignments.tsv"
            candidate_assignments = root / "candidate.assignments.tsv"
            baseline_assignments.write_text(
                header + "a\t100\tbinned\t99\t99.fna\n",
                encoding="utf-8",
            )
            candidate_assignments.write_text(
                header + "a\t100\tbinned\t1\t1.fna\n",
                encoding="utf-8",
            )
            result = audit(
                argparse.Namespace(
                    rank=981,
                    baseline_output=baseline_output,
                    candidate_output=candidate_output,
                    baseline_taxonomy=baseline_taxonomy,
                    candidate_taxonomy=candidate_taxonomy,
                    baseline_assignments=baseline_assignments,
                    candidate_assignments=candidate_assignments,
                )
            )
            self.assertTrue(result["classifier_taxonomy_agreement"])
            self.assertTrue(result["complete_agreement"])
            self.assertEqual(result["classifier"]["full_line_agreement"], 1)
            self.assertTrue(result["taxvamb"]["canonical_groups_equal"])

    def test_classifier_gate_is_independent_of_stochastic_bins(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            classifier = root / "classifier.output.gz"
            with gzip.open(classifier, "wt", encoding="utf-8") as handle:
                handle.write("C\tr000001Ca\t2\t100\t2:66\n")
            taxonomy = root / "taxonomy.tsv"
            taxonomy.write_text(
                "contigs\tpredictions\na\td__Eukaryota\n",
                encoding="utf-8",
            )
            header = "contig_id\tlength_bp\tstatus\tbin_id\tbin_file\n"
            baseline_assignments = root / "baseline.assignments.tsv"
            candidate_assignments = root / "candidate.assignments.tsv"
            baseline_assignments.write_text(
                header + "a\t100\tbinned\t1\t1.fna\n",
                encoding="utf-8",
            )
            candidate_assignments.write_text(
                header + "a\t100\tunbinned\t\t\n",
                encoding="utf-8",
            )
            result = audit(
                argparse.Namespace(
                    rank=1,
                    baseline_output=classifier,
                    candidate_output=classifier,
                    baseline_taxonomy=taxonomy,
                    candidate_taxonomy=taxonomy,
                    baseline_assignments=baseline_assignments,
                    candidate_assignments=candidate_assignments,
                )
            )
            self.assertTrue(result["classifier_taxonomy_agreement"])
            self.assertFalse(result["complete_agreement"])
            self.assertFalse(result["taxvamb"]["canonical_groups_equal"])


if __name__ == "__main__":
    unittest.main()
