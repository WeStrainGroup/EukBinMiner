from __future__ import annotations

import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "materialize_cli_contigs.py"
SPEC = importlib.util.spec_from_file_location("materialize_cli_contigs", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MaterializeCliContigsTests(unittest.TestCase):
    def test_union_and_inclusive_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            euk = root / "euk.fa"
            unknown = root / "unknown.fa.gz"
            euk.write_text(
                ">below\n" + "A" * 1499 + "\n"
                ">boundary description\n" + "C" * 1500 + "\n",
                encoding="utf-8",
            )
            with gzip.open(unknown, "wt", encoding="utf-8") as handle:
                handle.write(">above\n" + "G" * 1501 + "\n")
            depth = root / "depth.tsv"
            depth.write_text(
                "contigName\tdepth\n"
                "below\t1\n"
                "boundary\t2\n"
                "above\t3\n",
                encoding="utf-8",
            )
            metadata = root / "input.json"
            metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": str(euk),
                            "unknown_fasta": str(unknown),
                            "depth_tsv": str(depth),
                        }
                    }
                ),
                encoding="utf-8",
            )
            output = root / "runtime" / "contigs.fa"
            depth_output = root / "runtime" / "depth.tsv"
            summary = root / "runtime" / "summary.json"

            result = MODULE.materialize(
                metadata,
                output,
                summary,
                1500,
                depth_output,
            )

            self.assertEqual(result["total_contigs"], 3)
            self.assertEqual(result["selected_contigs"], 2)
            self.assertEqual(result["selected_bp"], 3001)
            self.assertEqual(
                output.read_text(encoding="utf-8").count(">"),
                3,
            )
            self.assertEqual(
                json.loads(summary.read_text(encoding="utf-8"))[
                    "minimum_contig_bp_inclusive"
                ],
                1500,
            )
            self.assertEqual(
                depth_output.read_text(encoding="utf-8").splitlines(),
                [
                    "contigName\tdepth",
                    "below\t1",
                    "boundary\t2",
                    "above\t3",
                ],
            )
            self.assertEqual(
                result["depth_output"]["duplicate_rows_removed"],
                0,
            )

    def test_duplicate_ids_and_overwrite_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            euk = root / "euk.fa"
            unknown = root / "unknown.fa"
            euk.write_text(">same\nAAAA\n", encoding="utf-8")
            unknown.write_text(">same\nCCCC\n", encoding="utf-8")
            depth = root / "depth.tsv"
            depth.write_text(
                "contigName\tdepth\nsame\t1\nother\t1\n",
                encoding="utf-8",
            )
            metadata = root / "input.json"
            metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": {"path": str(euk)},
                            "unknown_fasta": {"path": str(unknown)},
                            "depth_tsv": {"path": str(depth)},
                        }
                    }
                ),
                encoding="utf-8",
            )
            output = root / "contigs.fa"
            summary = root / "summary.json"
            with self.assertRaisesRegex(ValueError, "duplicate contig ID"):
                MODULE.materialize(metadata, output, summary, 1)

            unknown.write_text(">other\nCCCC\n", encoding="utf-8")
            output.write_text("existing\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                MODULE.materialize(metadata, output, summary, 1)

    def test_contigs_without_depth_are_excluded_and_empty_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            euk = root / "euk.fa"
            unknown = root / "unknown.fa"
            depth = root / "depth.tsv"
            euk.write_text(">euk\nAAAA\n", encoding="utf-8")
            unknown.write_text(">unknown\nCCCC\n", encoding="utf-8")
            depth.write_text(
                "contigName\tdepth\nunrelated\t1\n",
                encoding="utf-8",
            )
            metadata = root / "input.json"
            metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": str(euk),
                            "unknown_fasta": str(unknown),
                            "depth_tsv": str(depth),
                        }
                    }
                ),
                encoding="utf-8",
            )
            output = root / "contigs.fa"
            summary = root / "summary.json"

            result = MODULE.materialize(metadata, output, summary, 1)

            self.assertEqual(result["source_contigs"], 2)
            self.assertEqual(result["dropped_no_depth"], 2)
            self.assertEqual(result["total_contigs"], 0)
            self.assertEqual(output.read_bytes(), b"")

    def test_identical_depth_duplicates_are_canonicalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            euk = root / "euk.fa"
            unknown = root / "unknown.fa"
            euk.write_text(">same\nAAAA\n", encoding="utf-8")
            unknown.write_text(">other\nCCCC\n", encoding="utf-8")
            depth = root / "depth.tsv"
            depth.write_text(
                "contigName\tcontigLen\ttotalAvgDepth\n"
                "same\t4\t1.5\n"
                "same\t4\t1.5\n"
                "other\t4\t2.5\n"
                "unrelated\t5\t3.5\n",
                encoding="utf-8",
            )
            metadata = root / "input.json"
            metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": str(euk),
                            "unknown_fasta": str(unknown),
                            "depth_tsv": str(depth),
                        }
                    }
                ),
                encoding="utf-8",
            )
            output = root / "contigs.fa"
            depth_output = root / "canonical.depth.tsv"
            summary = root / "summary.json"

            result = MODULE.materialize(
                metadata,
                output,
                summary,
                1,
                depth_output,
            )

            self.assertEqual(
                depth_output.read_text(encoding="utf-8").splitlines(),
                [
                    "contigName\tcontigLen\ttotalAvgDepth",
                    "same\t4\t1.5",
                    "other\t4\t2.5",
                    "unrelated\t5\t3.5",
                ],
            )
            self.assertEqual(
                result["depth_source"]["duplicate_rows_removed"],
                1,
            )
            self.assertEqual(
                result["depth_output"]["contig_ids"],
                3,
            )

            conflicting = root / "conflicting.depth.tsv"
            conflicting.write_text(
                "contigName\tdepth\n"
                "same\t1.5\n"
                "same\t9.5\n",
                encoding="utf-8",
            )
            conflicting_metadata = root / "conflicting.input.json"
            conflicting_metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": str(euk),
                            "unknown_fasta": str(unknown),
                            "depth_tsv": str(conflicting),
                        }
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError,
                "conflicting duplicate depth",
            ):
                MODULE.materialize(
                    conflicting_metadata,
                    root / "conflicting.fa",
                    root / "conflicting.summary.json",
                    1,
                    root / "conflicting.canonical.tsv",
                )

    def test_unmatched_malformed_depth_rows_are_audited_and_removed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            euk = root / "euk.fa"
            unknown = root / "unknown.fa"
            euk.write_text(">same\nAAAA\n", encoding="utf-8")
            unknown.write_text(">other\nCCCC\n", encoding="utf-8")
            depth = root / "depth.tsv"
            depth.write_text(
                "contigName\tcontigLen\ttotalAvgDepth\tsample\tsample-var\n"
                "same\t4\t1.5\t1.5\t0.1\n"
                "5.50909\t2.15462\n"
                "other\t4\t2.5\t2.5\t0.2\n",
                encoding="utf-8",
            )
            metadata = root / "input.json"
            metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": str(euk),
                            "unknown_fasta": str(unknown),
                            "depth_tsv": str(depth),
                        }
                    }
                ),
                encoding="utf-8",
            )
            output = root / "contigs.fa"
            depth_output = root / "canonical.depth.tsv"
            summary = root / "summary.json"

            result = MODULE.materialize(
                metadata,
                output,
                summary,
                1,
                depth_output,
            )

            self.assertEqual(
                depth_output.read_text(encoding="utf-8").splitlines(),
                [
                    (
                        "contigName\tcontigLen\ttotalAvgDepth\t"
                        "sample\tsample-var"
                    ),
                    "same\t4\t1.5\t1.5\t0.1",
                    "other\t4\t2.5\t2.5\t0.2",
                ],
            )
            self.assertEqual(
                result["depth_source"]["malformed_rows_removed"],
                1,
            )
            self.assertEqual(
                result["depth_source"]["malformed_row_line_numbers"],
                [3],
            )
            self.assertFalse(
                result["depth_source"][
                    "malformed_row_line_numbers_truncated"
                ]
            )

            selected_malformed = root / "selected-malformed.depth.tsv"
            selected_malformed.write_text(
                "contigName\tcontigLen\ttotalAvgDepth\n"
                "same\t4\n"
                "other\t4\t2.5\n",
                encoding="utf-8",
            )
            selected_metadata = root / "selected-malformed.input.json"
            selected_metadata.write_text(
                json.dumps(
                    {
                        "inputs": {
                            "euk_fasta": str(euk),
                            "unknown_fasta": str(unknown),
                            "depth_tsv": str(selected_malformed),
                        }
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError,
                "malformed depth row for input contig 'same'",
            ):
                MODULE.materialize(
                    selected_metadata,
                    root / "selected-malformed.fa",
                    root / "selected-malformed.summary.json",
                    1,
                    root / "selected-malformed.canonical.tsv",
                )


if __name__ == "__main__":
    unittest.main()
