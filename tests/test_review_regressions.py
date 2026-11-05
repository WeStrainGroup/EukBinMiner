from __future__ import annotations
import csv
import fcntl
import json
import shutil
import unittest
from pathlib import Path
from test_multi_coverage import MODULE as cli
import test_mock_end_to_end as mock
from select_within_sample_representatives import fasta_contigs, Candidate, shared_overlap
from merge_busco import parse_summary
from depth_schema import identify_depth_schema


class ReviewRegressionTests(unittest.TestCase):
    setUp = mock.MockEndToEndTests.setUp
    tearDown = mock.MockEndToEndTests.tearDown
    run_cli = mock.MockEndToEndTests.run_cli
    base_run = mock.MockEndToEndTests.base_run

    def arguments(self, output):
        return [*self.base_run(output), "--binners", "metabat2",
                "--qc-tools", "eukcc", "--decision-qc", "eukcc",
                "--eukcc-db", str(self.eukcc_db)]

    def test_reuse_rejects_changes_and_corrupt_mags(self):
        output = self.root / "reuse"
        args = self.arguments(output)
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("already complete", result.stderr)
        result = self.run_cli(*args, "--overlap-threshold", "0.75")
        self.assertEqual(result.returncode, 2)
        self.assertIn("differ", result.stderr)
        original = self.depth.read_text()
        self.depth.write_text(original.replace("c1\t3", "c1\t9"))
        self.assertEqual(self.run_cli(*args).returncode, 2)
        self.depth.write_text(original)
        genome = next((output / "mags").glob("*.fa"))
        genome.write_text(">damaged\nACGT\n")
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 2)
        self.assertIn("incomplete or modified", result.stderr)

    def test_lock_blocks_concurrent_writers(self):
        output = self.root / "locked"
        output.mkdir()
        with (output / ".run.lock").open("a") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            result = self.run_cli(*self.arguments(output))
        self.assertEqual(result.returncode, 2)
        self.assertIn("another process", result.stderr)

    def test_distinct_contig_ids_preserve_identical_sequences(self):
        left = self.root / "left.fa"
        right = self.root / "right.fa"
        left.write_text(">a\nACGT\n>b\nACGT\n")
        right.write_text(">c\nACGT\n")
        lc, rc = fasta_contigs(left), fasta_contigs(right)
        self.assertEqual(len(lc), 2)
        a = Candidate({}, 90, 0, 2, 4, 90, lc, 8)
        b = Candidate({}, 80, 0, 1, 4, 80, rc, 4)
        self.assertEqual(shared_overlap(a, b), (0, 0.0))
        right.write_text(">b\nACGT\n")
        b.contig_bp = fasta_contigs(right)
        self.assertEqual(shared_overlap(a, b), (4, 1.0))

    def test_invalid_quality_cannot_pass(self):
        for completeness, contamination in (("101", "0"), ("90", "-1"), ("nan", "0"), ("90", "inf")):
            passed, score = cli.qc_metrics(completeness, contamination,
                min_completeness=50, max_contamination=10, contamination_penalty=5)
            self.assertEqual((passed, score), ("no", ""))

    def test_jgi_alias_is_applied(self):
        schema = identify_depth_schema(
            ["contigName", "contigLen", "totalAvgDepth", "sample.bam", "sample.bam-var"],
            "assembly", single_coverage_name="mapping_A")
        self.assertEqual(schema.coverage_names, ("mapping_A",))
        self.assertEqual(schema.source_coverage_names, ("sample.bam",))

    def test_invalid_busco_summary_is_rejected(self):
        summary = self.root / "summary.txt"
        for text in ("C:101%[S:100%,D:1%],F:0%,M:0%,n:100",
                     "C:90%[S:80%,D:1%],F:5%,M:5%,n:100",
                     "C:90%[S:89%,D:1%],F:0%,M:0%,n:100"):
            summary.write_text(text)
            self.assertIsNone(parse_summary(summary))

    def test_kunpeng_complete_enodata_is_validated(self):
        output = self.root / "complete-enodata"
        args = self.arguments(output)
        args[args.index("--binners") + 1] = "taxvamb"
        result = self.run_cli(*args, "--kunpeng-db", str(self.kunpeng_db),
                              mock_kunpeng_mode="enodata_complete")
        self.assertEqual(result.returncode, 0, result.stderr)
        provenance = json.loads((output / "taxonomy/kunpeng/provenance.json").read_text())
        self.assertEqual(provenance["enodata_evidence_mode"],
                         "complete_output_report_consistent")

    def test_taxonomy_failure_preserves_other_binners(self):
        output = self.root / "taxonomy-failed"
        args = self.arguments(output)
        args[args.index("--binners") + 1] = "metabat2,vamb,taxvamb"
        result = self.run_cli(*args, "--kunpeng-db", str(self.kunpeng_db),
                              mock_kunpeng_mode="exit")
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads((output / "summary.json").read_text())
        self.assertEqual(summary["completed_binners"], ["metabat2", "vamb"])

    def test_zero_coverage_column_is_logged_and_removed_for_vamb(self):
        self.depth.write_text("contig_id\talpha\tzero\nc1\t3\t0\nc2\t2\t0\nc3\t1\t0\n")
        output = self.root / "zero-column"
        args = self.arguments(output)
        args[args.index("--binners") + 1] = "vamb"
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        metadata = json.loads((output / "inputs/r000001/input.json").read_text())
        self.assertEqual(metadata["coverage_columns"], ["alpha"])
        self.assertEqual(metadata["dropped_all_zero_coverage_columns"], ["zero"])

    def test_missing_vamb_preserves_classic_result(self):
        (self.tools / "vamb").unlink()
        output = self.root / "missing-vamb"
        args = self.arguments(output)
        args[args.index("--binners") + 1] = "metabat2,vamb"
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads((output / "summary.json").read_text())["completed_binners"],
                         ["metabat2"])

    def test_installed_layout_runs_a_classic_binner(self):
        import test_mock_end_to_end
        original_cli = test_mock_end_to_end.CLI
        project = Path(original_cli).parents[1]
        prefix = self.root / "prefix"
        (prefix / "bin").mkdir(parents=True)
        shutil.copy2(original_cli, prefix / "bin/eukbinminer")
        share = prefix / "share/eukbinminer"
        shutil.copytree(project / "scripts", share / "scripts")
        shutil.copytree(project / "envs", share / "envs")
        try:
            test_mock_end_to_end.CLI = prefix / "bin/eukbinminer"
            result = self.run_cli(*self.arguments(self.root / "installed"))
            self.assertEqual(result.returncode, 0, result.stderr)
        finally:
            test_mock_end_to_end.CLI = original_cli

    def test_vamb_invalid_bins_cannot_be_published(self):
        tool = self.tools / "vamb"
        text = tool.read_text()
        original = 'write_fasta(out / "bins" / f"{mode}-bin.fna", records(source))'
        self.assertIn(original, text)
        text = text.replace(original,
                            'write_fasta(out / "bins" / f"{mode}-bin.fna", [("unknown_contig", "ACGT")])')
        tool.write_text(text)
        output = self.root / "invalid-vamb"
        args = self.arguments(output)
        args[args.index("--binners") + 1] = "metabat2,vamb"
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads((output / "summary.json").read_text())
        self.assertEqual(summary["completed_binners"], ["metabat2"])
        self.assertFalse((output / "results/standalone/sample1/results/union/vamb").exists())

    def test_empty_metaeuk_remains_unassessed(self):
        tool = self.tools / "eukcc"
        text = tool.read_text()
        text = text.replace('        bin_path = args[-1]', """
        (out / "run.log").write_text("Metaeuk produced empty fasta.")
        (out / "sample_metaeuk_cleaned.faa").touch()
        (out / "sample_metaeuk.codon.fas").touch()
        raise SystemExit(1)
        bin_path = args[-1]""")
        tool.write_text(text)
        output = self.root / "empty-proteins"
        result = self.run_cli(*self.arguments(output))
        self.assertEqual(result.returncode, 0, result.stderr)
        with (output / "qc/decision_qc.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["eukcc_status"], "metaeuk_empty_proteins")
        self.assertEqual(row["qc_pass"], "no")
        self.assertEqual(row["completeness"], "")

    def test_cached_vamb_abundance_is_checked(self):
        from prepare_inputs import prepare_inputs
        preflight = self.root / "prepared"
        normalized, _ = cli.normalized_depth_inputs(
            (cli.DepthInput(self.depth),), cli.fasta_lengths(self.fasta), "sample1")
        depth = self.root / "normalized.tsv"
        depth.write_text(normalized)
        prepare_inputs(euk_fasta=self.fasta, unknown_fasta=None,
                       depth_tsv=depth, sample="sample1", out_dir=preflight,
                       minimum_length_inclusive=1)
        prepared = self.root / "vamb-input"
        cli.prepare_vamb_input(preflight, prepared, "standalone", "sample1")
        cli.prepare_vamb_input(preflight, prepared, "standalone", "sample1")
        (prepared / "abundance.tsv").write_text("damaged")
        with self.assertRaisesRegex(cli.UserError, "checksum mismatch"):
            cli.prepare_vamb_input(preflight, prepared, "standalone", "sample1")

    def test_selection_preserves_same_binner_and_avoids_transitive_merging(self):
        from select_within_sample_representatives import select
        metadata = self.root / "sample.json"
        metadata.write_text(json.dumps({"rank": 1, "study": "study", "sample": "sample"}))
        rows = []
        for identifier, binner, contigs, completeness in (
            ("a", "remag", ("c1", "c2"), 95),
            ("b", "metabat2", ("c2", "c3"), 90),
            ("c", "concoct", ("c3", "c4"), 85),
        ):
            fasta = self.root / f"{identifier}.fa"
            fasta.write_text("".join(f">{name}\nACGT\n" for name in contigs))
            rows.append({
                "bin_id": identifier, "sample_rank": 1, "study": "study",
                "sample": "sample", "binner": binner, "bin_path": str(fasta),
                "decision_qc": "eukcc", "contamination_source": "eukcc_contamination",
                "qc_pass": "yes", "completeness": completeness, "contamination": 0,
                "contigs": 2, "n50": 4, "bp": 8,
            })
        table = self.root / "qc.tsv"
        for label, candidates, expected in (
            ("chain", rows, ["a", "c"]),
            ("same-binner", [rows[0], {**rows[1], "binner": "remag"}], ["a", "b"]),
        ):
            with table.open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(candidates)
            output = self.root / label
            select(table, metadata, output, 0.5, 5, 50, 10)
            with (output / "representatives.tsv").open() as handle:
                selected = [row["representative_bin_id"] for row in csv.DictReader(handle, delimiter="\t")]
            self.assertEqual(selected, expected)

    def test_failed_current_status_excludes_existing_branch(self):
        results = self.root / "results"
        sample_root = results / "standalone/sample1"
        branch = sample_root / "results/union/metabat2"
        (branch / "bins").mkdir(parents=True)
        (branch / "audit").mkdir()
        (branch / "audit/audit.json").write_text("{}")
        (branch / "fingerprint.json").write_text("{}")
        status = sample_root / "branch_status.tsv"
        status.write_text("input_set\tbinner\tstatus\texit_code\nunion\tmetabat2\tfailed\t1\n")
        self.assertEqual(cli.successful_binners(results, "standalone", "sample1", ("metabat2",)), ())
        status.write_text("input_set\tbinner\tstatus\texit_code\nunion\tmetabat2\tcompleted\t0\n")
        self.assertEqual(cli.successful_binners(results, "standalone", "sample1", ("metabat2",)), ("metabat2",))
