from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
CLI = PROJECT / "bin" / "eukbinminer"
LEGACY_CLI = PROJECT / "bin" / "eukmagminer"


class ReferenceDatasetCliTests(unittest.TestCase):
    def test_legacy_cli_alias_forwards_to_eukbinminer(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(LEGACY_CLI), "--version"],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("deprecated compatibility alias", completed.stderr)

    def test_validate_reports_explicit_qc_reference_datasets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "contigs.fa"
            depth = root / "depth.tsv"
            eukcc_db = root / "eukcc_db"
            busco_db = root / "fungi_odb12"
            fasta.write_text(">c1\nAAAA\n", encoding="utf-8")
            depth.write_text("contig_id\tdepth\nc1\t1\n", encoding="utf-8")
            eukcc_db.mkdir()
            (busco_db / "hmms").mkdir(parents=True)
            (busco_db / "dataset.cfg").write_text(
                "name=fungi\n",
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    str(CLI),
                    "validate",
                    "--sample",
                    "sample_1",
                    "--contigs",
                    str(fasta),
                    "--depth",
                    str(depth),
                    "--binners",
                    "metabat2",
                    "--eukcc-db",
                    str(eukcc_db),
                    "--busco-db",
                    str(busco_db),
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            plan = json.loads(completed.stdout)
            references = plan["qc_references"]
            self.assertEqual(references["eukcc_db"], str(eukcc_db.resolve()))
            self.assertEqual(references["busco_db"], str(busco_db.resolve()))
            self.assertEqual(plan["random_seed"], 6311)

    def test_new_environment_names_take_priority_over_legacy_aliases(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            preferred = root / "preferred"
            legacy = root / "legacy"
            preferred.mkdir()
            legacy.mkdir()
            environment = dict(os.environ)
            environment["EUKBINMINER_EUKCC_DB"] = str(preferred)
            environment["EUKMAGMINER_EUKCC_DB"] = str(legacy)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(CLI),
                    "doctor",
                    "--binners",
                    "metabat2",
                    "--qc-tools",
                    "eukcc",
                ],
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=environment,
            )
            report = json.loads(completed.stdout)
            self.assertEqual(
                report["databases"]["eukcc"]["path"],
                str(preferred.resolve()),
            )


if __name__ == "__main__":
    unittest.main()
