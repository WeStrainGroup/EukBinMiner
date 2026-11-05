from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

from prepare_inputs import prepare_inputs


class ContigBoundaryTests(unittest.TestCase):
    def test_minimum_contig_length_is_inclusive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "contigs.fa"
            depth = root / "depth.tsv"
            output = root / "prepared"
            fasta.write_text(
                ">c1499\n" + "A" * 1499 + "\n"
                ">c1500\n" + "C" * 1500 + "\n"
                ">c1501\n" + "G" * 1501 + "\n",
                encoding="utf-8",
            )
            depth.write_text(
                "contigName\tcontigLen\ttotalAvgDepth\n"
                "c1499\t1499\t1\n"
                "c1500\t1500\t2\n"
                "c1501\t1501\t3\n",
                encoding="utf-8",
            )

            status = prepare_inputs(
                euk_fasta=fasta,
                unknown_fasta=None,
                depth_tsv=depth,
                sample="boundary",
                out_dir=output,
                minimum_length_inclusive=1500,
            )

            self.assertEqual(status, "created")
            selected = (output / "union.fa").read_text(encoding="utf-8")
            self.assertNotIn(">c1499\n", selected)
            self.assertIn(">c1500\n", selected)
            self.assertIn(">c1501\n", selected)
            metadata = json.loads(
                (output / "preflight.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["filter"]["minimum_contig_bp"], 1500)
            self.assertEqual(metadata["collections"]["union"]["contigs"], 2)


if __name__ == "__main__":
    unittest.main()
