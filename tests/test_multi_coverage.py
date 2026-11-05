from __future__ import annotations

import importlib.machinery
import importlib.util
import sys
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
CLI = PROJECT / "bin" / "eukbinminer"
LOADER = importlib.machinery.SourceFileLoader("eukbinminer_test", str(CLI))
SPEC = importlib.util.spec_from_loader(LOADER.name, LOADER)
assert SPEC is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
LOADER.exec_module(MODULE)

sys.path.insert(0, str(PROJECT / "scripts"))
from split_remag_coverage import SplitError, split_coverage
from filter_prepared_inputs import filter_inputs


class MultiCoverageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.lengths = OrderedDict((identifier, 4) for identifier in ("c1", "c2"))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_depth(self, name: str, text: str) -> Path:
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_generic_multi_column_matrix_is_preserved(self) -> None:
        depth = self.write_depth(
            "matrix.tsv",
            "contig_id\talpha\tbeta\tgamma\n"
            "c1\t1\t2\t3\n"
            "c2\t4\t5\t6\n",
        )
        text, summary = MODULE.normalized_depth_inputs(
            (MODULE.DepthInput(depth),),
            self.lengths,
            "assembly",
        )

        self.assertEqual(
            text.splitlines(),
            [
                "contigName\tcontigLen\ttotalAvgDepth\talpha\talpha-var\t"
                "beta\tbeta-var\tgamma\tgamma-var",
                "c1\t4\t6\t1\t0\t2\t0\t3\t0",
                "c2\t4\t15\t4\t0\t5\t0\t6\t0",
            ],
        )
        self.assertEqual(summary["coverage_columns"], ["alpha", "beta", "gamma"])
        self.assertEqual(summary["coverage_column_count"], 3)
        self.assertEqual(
            summary["total_avg_depth_policy"],
            "sum_of_coverage_columns",
        )

    def test_multi_sample_jgi_pairs_and_variances_are_preserved(self) -> None:
        depth = self.write_depth(
            "jgi.tsv",
            "contigName\tcontigLen\ttotalAvgDepth\tA.bam\tA.bam-var\t"
            "B.bam\tB.bam-var\n"
            "c1\t4\t3\t1\t0.1\t2\t0.2\n"
            "c2\t4\t9\t4\t0.4\t5\t0.5\n",
        )
        text, summary = MODULE.normalized_depth_inputs(
            (MODULE.DepthInput(depth),),
            self.lengths,
            "assembly",
        )

        self.assertEqual(text, depth.read_text(encoding="utf-8"))
        self.assertEqual(summary["coverage_columns"], ["A.bam", "B.bam"])
        self.assertEqual(summary["total_avg_depth_policy"], "provided")

    def test_repeated_named_two_column_files_are_combined(self) -> None:
        alpha = self.write_depth(
            "depth.tsv",
            "contig_id\tdepth\nc1\t1\nc2\t2\n",
        )
        beta = self.write_depth(
            "other.tsv",
            "contig_id\tdepth\nc1\t3\nc2\t4\n",
        )
        text, summary = MODULE.normalized_depth_inputs(
            (
                MODULE.DepthInput(alpha, "alpha"),
                MODULE.DepthInput(beta, "beta"),
            ),
            self.lengths,
            "assembly",
        )

        self.assertEqual(
            text.splitlines(),
            [
                "contigName\tcontigLen\ttotalAvgDepth\talpha\talpha-var\t"
                "beta\tbeta-var",
                "c1\t4\t4\t1\t0\t3\t0",
                "c2\t4\t6\t2\t0\t4\t0",
            ],
        )
        self.assertEqual(summary["depth_source_count"], 2)
        self.assertEqual(summary["coverage_columns"], ["alpha", "beta"])

    def test_invalid_matrix_values_shape_names_and_ids_fail_closed(self) -> None:
        invalid_cases = {
            "negative.tsv": (
                "contig_id\tA\tB\nc1\t1\t-1\nc2\t2\t3\n",
                "finite and non-negative",
            ),
            "nan.tsv": (
                "contig_id\tA\tB\nc1\t1\tNaN\nc2\t2\t3\n",
                "finite and non-negative",
            ),
            "duplicate-column.tsv": (
                "contig_id\tA\tA\nc1\t1\t2\nc2\t2\t3\n",
                "duplicate column",
            ),
            "empty-column.tsv": (
                "contig_id\t\tB\nc1\t1\t2\nc2\t2\t3\n",
                "empty or unsafe",
            ),
            "short-row.tsv": (
                "contig_id\tA\tB\nc1\t1\nc2\t2\t3\n",
                "fields; expected",
            ),
            "duplicate-id.tsv": (
                "contig_id\tA\tB\nc1\t1\t2\nc1\t2\t3\n",
                "duplicate depth ID",
            ),
        }
        for filename, (contents, message) in invalid_cases.items():
            with self.subTest(filename=filename):
                depth = self.write_depth(filename, contents)
                with self.assertRaisesRegex(MODULE.UserError, message):
                    MODULE.normalized_depth_inputs(
                        (MODULE.DepthInput(depth),),
                        self.lengths,
                        "assembly",
                    )

    def test_every_repeated_file_must_cover_every_fasta_contig(self) -> None:
        alpha = self.write_depth(
            "alpha.tsv",
            "contig_id\tdepth\nc1\t1\nc2\t2\n",
        )
        beta = self.write_depth(
            "beta.tsv",
            "contig_id\tdepth\nc1\t3\n",
        )
        with self.assertRaisesRegex(MODULE.UserError, "1 FASTA contigs lack depth"):
            MODULE.normalized_depth_inputs(
                (
                    MODULE.DepthInput(alpha, "alpha"),
                    MODULE.DepthInput(beta, "beta"),
                ),
                self.lengths,
                "assembly",
            )

    def test_remag_adapter_writes_one_headerless_file_per_column(self) -> None:
        matrix = self.write_depth(
            "coverage.tsv",
            "contig_id\talpha\tbeta\nc1\t1\t3\nc2\t2\t4\n",
        )
        out_dir = self.root / "remag"
        metadata = split_coverage(matrix, out_dir)

        self.assertEqual(
            [item["coverage_name"] for item in metadata],
            ["alpha", "beta"],
        )
        self.assertEqual(
            (out_dir / "coverage.0001.tsv").read_text(encoding="utf-8"),
            "c1\t1\nc2\t2\n",
        )
        self.assertEqual(
            (out_dir / "coverage.0002.tsv").read_text(encoding="utf-8"),
            "c1\t3\nc2\t4\n",
        )
        with self.assertRaisesRegex(SplitError, "non-empty output"):
            split_coverage(matrix, out_dir)

    def test_length_filter_preserves_all_coverage_dimensions(self) -> None:
        fasta = self.root / "input.fa"
        fasta.write_text(
            ">c1\nAAAA\n>c2\nCCCC\n",
            encoding="utf-8",
        )
        depth = self.write_depth(
            "jgi.tsv",
            "contigName\tcontigLen\ttotalAvgDepth\talpha\talpha-var\t"
            "beta\tbeta-var\n"
            "c1\t4\t4\t1\t0\t3\t0\n"
            "c2\t4\t6\t2\t0\t4\t0\n",
        )
        output = self.root / "filtered"
        self.assertEqual(filter_inputs(fasta, depth, "assembly", 4, output), "created")
        self.assertEqual(
            (output / "concoct.coverage.tsv").read_text(encoding="utf-8"),
            "contig_id\talpha\tbeta\nc1\t1\t3\nc2\t2\t4\n",
        )
        self.assertEqual(
            (output / "remag.coverage.tsv").read_text(encoding="utf-8"),
            "c1\t4\nc2\t6\n",
        )


if __name__ == "__main__":
    unittest.main()
