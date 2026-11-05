from __future__ import annotations

import csv
import json
import os
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from mock_k2 import write_options, write_taxonomy


PROJECT = Path(__file__).resolve().parents[1]
CLI = PROJECT / "bin" / "eukbinminer"

MOCK_TOOL = r"""#!/usr/bin/env python3
import csv
import json
import os
import struct
import sys
from collections import deque
from pathlib import Path


name = Path(sys.argv[0]).name
args = sys.argv[1:]


def value(flag):
    return args[args.index(flag) + 1]


def records(path):
    observed = []
    identifier = None
    sequence = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if identifier is not None:
                observed.append((identifier, "".join(sequence)))
            identifier = raw[1:].split()[0]
            sequence = []
        elif raw:
            sequence.append(raw.strip())
    if identifier is not None:
        observed.append((identifier, "".join(sequence)))
    return observed


def write_fasta(path, observed):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in observed),
        encoding="utf-8",
    )


def reverse_complement_bits(value, length):
    output = 0
    for _ in range(length):
        output = (output << 2) | (3 - (value & 3))
        value >>= 2
    return output


def minimizer_event_count(sequence, options_path):
    data = Path(options_path).read_bytes()
    k, l, spaced, toggle = struct.unpack_from("<QQQQ", data)
    mask = (1 << (2 * l)) - 1
    toggle &= mask
    capacity = k - l
    window = deque()
    count = value = valid = events = 0
    for base in sequence.encode("ascii"):
        code = {
            ord("A"): 0,
            ord("C"): 1,
            ord("G"): 2,
            ord("T"): 3,
        }.get(base)
        if code is None:
            window.clear()
            count = value = valid = 0
            continue
        value = ((value << 2) | code) & mask
        valid += 1
        if valid < l:
            continue
        canonical = min(value, reverse_complement_bits(value, l))
        candidate = (canonical & spaced) ^ toggle
        while window and window[-1][0] > candidate:
            window.pop()
        changed = (
            (not window and count >= capacity)
            or count == capacity
        )
        window.append((candidate, count))
        while (
            window
            and count >= capacity
            and window[0][1] < count - capacity
        ):
            window.popleft()
            changed = True
        count += 1
        if count > capacity and changed:
            events += 1
    return events


def binner_records(source, binner):
    observed = records(source)
    if os.environ.get("MOCK_BIN_MODE") != "half":
        return observed
    if binner == "remag":
        return [row for row in observed if row[0] in {"c1", "c2"}]
    if binner == "metabat2":
        return [row for row in observed if row[0] in {"c1", "c3"}]
    return observed


if name == "conda":
    if (
        len(args) < 5
        or args[:2] != ["run", "--no-capture-output"]
        or args[2] not in {"-n", "-p"}
    ):
        raise SystemExit(f"unexpected mock conda invocation: {args!r}")
    os.execvp(args[4], args[4:])
elif name == "seqkit":
    input_list = Path(value("-X"))
    print("file\tnum_seqs\tsum_len\tN50\tmax_len")
    for raw_path in input_list.read_text(encoding="utf-8").splitlines():
        observed = records(raw_path)
        lengths = sorted((len(sequence) for _, sequence in observed), reverse=True)
        total = sum(lengths)
        running = 0
        n50 = 0
        for length in lengths:
            running += length
            if running >= total / 2:
                n50 = length
                break
        print(f"{Path(raw_path).resolve()}\t{len(lengths)}\t{total}\t{n50}\t{max(lengths)}")
elif name == "remag":
    if "--version" in args:
        print("remag 0.4.4")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--random-seed") != expected_seed:
            raise SystemExit("unexpected ReMAG random seed")
        source = args[0]
        out = Path(value("-o"))
        coverage_files = [
            args[index + 1]
            for index, argument in enumerate(args)
            if argument == "-c"
        ]
        expected_coverage_count = os.environ.get(
            "MOCK_EXPECTED_REMAG_COVERAGE_COUNT"
        )
        if (
            expected_coverage_count is not None
            and len(coverage_files) != int(expected_coverage_count)
        ):
            raise SystemExit("unexpected ReMAG coverage file count")
        out.mkdir(parents=True, exist_ok=True)
        (out / "mock_coverage_inputs.json").write_text(
            json.dumps(
                [
                    Path(path).read_text(encoding="utf-8").splitlines()
                    for path in coverage_files
                ]
            ),
            encoding="utf-8",
        )
        write_fasta(out / "bins" / "remag-bin.fa", binner_records(source, "remag"))
elif name == "metabat2":
    if "--help" in args:
        print("MetaBAT 2.18")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--seed") != expected_seed:
            raise SystemExit("unexpected MetaBAT2 random seed")
        source = value("-i")
        prefix = Path(value("-o"))
        write_fasta(
            Path(str(prefix) + ".1.fa"),
            binner_records(source, "metabat2"),
        )
elif name in {
    "concoct",
    "cut_up_fasta.py",
    "merge_cutup_clustering.py",
    "extract_fasta_bins.py",
}:
    if name == "concoct" and "--version" in args:
        print("CONCOCT 1.1.0")
    elif os.environ.get("MOCK_CONCOCT_SUCCESS") == "1":
        if name == "cut_up_fasta.py":
            source = args[0]
            bed = Path(value("-b"))
            bed.parent.mkdir(parents=True, exist_ok=True)
            bed_rows = []
            for identifier, sequence in records(source):
                bed_rows.append(
                    f"{identifier}\t0\t{len(sequence)}\t{identifier}\n"
                )
                print(f">{identifier}\n{sequence}")
            bed.write_text("".join(bed_rows), encoding="utf-8")
        elif name == "concoct":
            coverage = Path(value("--coverage_file"))
            header = coverage.read_text(encoding="utf-8").splitlines()[0]
            expected_columns = os.environ.get("MOCK_EXPECTED_COVERAGE_COLUMNS")
            if expected_columns is not None:
                observed_columns = header.split("\t")[1:]
                if observed_columns != expected_columns.split(","):
                    raise SystemExit("unexpected CONCOCT coverage columns")
            coverage_count = len(header.split("\t")) - 1
            disabled = (
                "--no_cov_normalization" in args
                and "--no_total_coverage" in args
            )
            if coverage_count == 1 and not disabled:
                raise SystemExit("single-column CONCOCT normalization was not disabled")
            if coverage_count > 1 and disabled:
                raise SystemExit("multi-column CONCOCT normalization was disabled")
            out = Path(value("-b"))
            out.mkdir(parents=True, exist_ok=True)
            threshold = value("--length_threshold")
            identifiers = [
                row.split("\t")[0]
                for row in coverage.read_text(encoding="utf-8").splitlines()[1:]
            ]
            (out / f"clustering_gt{threshold}.csv").write_text(
                "contig_id,cluster_id\n"
                + "".join(f"{identifier},0\n" for identifier in identifiers),
                encoding="utf-8",
            )
        elif name == "merge_cutup_clustering.py":
            print(Path(args[0]).read_text(encoding="utf-8"), end="")
        else:
            source = args[0]
            output_path = Path(value("--output_path"))
            write_fasta(output_path / "0.fa", records(source))
    else:
        print("intentional mock CONCOCT failure", file=sys.stderr)
        raise SystemExit(17)
elif name == "eukcc":
    if "--help" in args:
        print("EukCC 2.1.3")
    else:
        out = Path(value("--out"))
        out.mkdir(parents=True, exist_ok=True)
        bin_path = args[-1]
        if "metabat2" in bin_path:
            completeness, contamination = "91", "1"
        elif "taxvamb" in bin_path:
            completeness, contamination = "92", "2"
        else:
            completeness, contamination = "95", "2"
        (out / "quality.csv").write_text(
            "completeness,contamination,ncbi_lng\n"
            f"{completeness},{contamination},Eukaryota;4751\n",
            encoding="utf-8",
        )
elif name == "busco":
    if "--version" in args:
        print("BUSCO 6")
    else:
        out = Path(value("--out_path"))
        bin_id = value("--out")
        out.mkdir(parents=True, exist_ok=True)
        bin_path = value("--in")
        if "metabat2" in bin_path:
            complete, duplicated = "98.0", "1.0"
        else:
            complete, duplicated = "80.0", "2.0"
        single = str(float(complete) - float(duplicated))
        (out / f"short_summary.specific.{bin_id}.txt").write_text(
            f"C:{complete}%[S:{single}%,D:{duplicated}%],"
            f"F:1.0%,M:{99.0 - float(complete)}%,n:100\n",
            encoding="utf-8",
        )
elif name == "vamb":
    if "--version" in args:
        print("Vamb 5.0.4")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--seed") != expected_seed:
            raise SystemExit("unexpected VAMB random seed")
        mode = args[1]
        out = Path(value("--outdir"))
        source = value("--fasta")
        write_fasta(out / "bins" / f"{mode}-bin.fna", records(source))
        (out / "clusters.tsv").write_text("cluster\tcontig\n", encoding="utf-8")
        (out / "log.txt").write_text("mock VAMB complete\n", encoding="utf-8")
elif name == "kun_peng":
    if "--version" in args:
        print("kun_peng 0.7.11")
    else:
        if "--seed" in args or "--random-seed" in args:
            raise SystemExit("Kun-peng must not receive a random seed")
        if os.environ.get("MOCK_KUNPENG_MODE") == "exit":
            print("intentional Kun-peng failure", file=sys.stderr)
            raise SystemExit(23)
        chunk_dir = Path(value("--chunk-dir"))
        output_dir = Path(value("--output-dir"))
        chunk_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        if os.environ.get("MOCK_KUNPENG_MODE") == "missing":
            raise SystemExit(0)
        catalogue = Path(args[-1])
        observed = records(catalogue)
        report = output_dir / "output_1.kreport2"
        mode = os.environ.get("MOCK_KUNPENG_MODE")
        omitted = 1 if mode in {"zero_hit", "enodata", "truncated"} else 0
        classified = len(observed) - (
            1 if mode in {"zero_hit", "enodata"} else 0
        )
        report_rows = []
        if mode in {"zero_hit", "enodata"}:
            report_rows.append(
                f"{100 / len(observed):.2f}\t1\t1\tU\t0\tunclassified\n"
            )
        report_rows.extend(
            [
                f"100.00\t{classified}\t0\tR\t1\troot\n",
                f"100.00\t{classified}\t{classified}\tD\t2759\t"
                "  Eukaryota\n",
            ]
        )
        report.write_text("".join(report_rows), encoding="utf-8")
        output_rows = [
            f"C\t{identifier}\t2759\t{len(sequence)}\t"
            f"2759:{minimizer_event_count(sequence, Path(value('--db')) / 'opts.k2d')}\n"
            for identifier, sequence in observed
        ]
        if omitted:
            output_rows.pop()
        if os.environ.get("MOCK_KUNPENG_MODE") == "duplicate":
            output_rows.append(output_rows[0])
        output = output_dir / "output_1.txt"
        output.write_text(
            "".join(output_rows),
            encoding="utf-8",
        )
        if mode in {"enodata", "enodata_complete"}:
            print("Error: No data available (os error 61)", file=sys.stderr)
            raise SystemExit(1)
else:
    print(f"unsupported mock tool: {name}", file=sys.stderr)
    raise SystemExit(99)
"""


class MockEndToEndTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.tools = self.root / "tools"
        self.tools.mkdir()
        for name in (
            "conda",
            "seqkit",
            "remag",
            "metabat2",
            "concoct",
            "cut_up_fasta.py",
            "merge_cutup_clustering.py",
            "extract_fasta_bins.py",
            "eukcc",
            "busco",
            "vamb",
            "kun_peng",
        ):
            tool = self.tools / name
            tool.write_text(MOCK_TOOL, encoding="utf-8")
            tool.chmod(0o755)

        self.fasta = self.root / "contigs.fa"
        sequences = {
            "c1": "ACGT" * 25,
            "c2": "A" * 99 + "C",
            "c3": "A" * 98 + "CC",
        }
        self.fasta.write_text(
            "".join(
                f">{identifier}\n{sequence}\n"
                for identifier, sequence in sequences.items()
            ),
            encoding="utf-8",
        )
        self.depth = self.root / "depth.tsv"
        self.depth.write_text(
            "contig_id\tdepth\nc1\t3\nc2\t2\nc3\t1\n",
            encoding="utf-8",
        )
        self.taxonomy = self.root / "taxonomy.tsv"
        self.taxonomy.write_text(
            "contigs\tpredictions\n"
            "c1\td__Eukaryota\n"
            "c2\td__Eukaryota\n"
            "c3\td__Eukaryota\n",
            encoding="utf-8",
        )
        self.eukcc_db = self.root / "eukcc_db"
        self.eukcc_db.mkdir()
        (self.eukcc_db / "db.info").write_text("mock\n", encoding="utf-8")
        self.busco_db = self.root / "fungi_odb12"
        (self.busco_db / "hmms").mkdir(parents=True)
        (self.busco_db / "dataset.cfg").write_text("name=fungi\n", encoding="utf-8")
        self.kunpeng_db = self.root / "kunpeng"
        self.kunpeng_db.mkdir()
        write_options(self.kunpeng_db / "opts.k2d")
        write_taxonomy(self.kunpeng_db / "taxo.k2d")
        (self.kunpeng_db / "hash_config.k2d").write_bytes(
            struct.pack("<6Q", 0, 1, 3, 3, 3, 16)
        )
        (self.kunpeng_db / "hash_1.k2d").write_bytes(
            struct.pack("<2Q", 1, 3) + b"\0" * 12
        )
        self.environment = dict(os.environ)
        self.environment["PATH"] = (
            str(self.tools) + os.pathsep + self.environment.get("PATH", "")
        )
        mock_modules = self.root / "mock_python_modules"
        mock_modules.mkdir()
        (mock_modules / "numpy.py").write_text(
            "class _Random:\n"
            "    def seed(self, value):\n"
            "        self.value = value\n"
            "random = _Random()\n",
            encoding="utf-8",
        )
        (mock_modules / "torch.py").write_text(
            "def manual_seed(value):\n"
            "    return value\n"
            "class _Cuda:\n"
            "    @staticmethod\n"
            "    def is_available():\n"
            "        return False\n"
            "cuda = _Cuda()\n",
            encoding="utf-8",
        )
        self.environment["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                (
                    str(mock_modules),
                    self.environment.get("PYTHONPATH", ""),
                ),
            )
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_cli(
        self,
        *arguments: str,
        mock_bin_mode: str | None = None,
        mock_kunpeng_mode: str | None = None,
        expected_seed: int | None = None,
        expected_coverage_columns: tuple[str, ...] | None = None,
        concoct_success: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        environment = dict(self.environment)
        if mock_bin_mode is not None:
            environment["MOCK_BIN_MODE"] = mock_bin_mode
        if mock_kunpeng_mode is not None:
            environment["MOCK_KUNPENG_MODE"] = mock_kunpeng_mode
        if expected_seed is not None:
            environment["MOCK_EXPECTED_SEED"] = str(expected_seed)
        if expected_coverage_columns is not None:
            environment["MOCK_EXPECTED_COVERAGE_COLUMNS"] = ",".join(
                expected_coverage_columns
            )
            environment["MOCK_EXPECTED_REMAG_COVERAGE_COUNT"] = str(
                len(expected_coverage_columns)
            )
        if concoct_success:
            environment["MOCK_CONCOCT_SUCCESS"] = "1"
        result = subprocess.run(
            [sys.executable, str(CLI), *arguments],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )

        if result.returncode and "--output" in arguments:
            output = Path(arguments[arguments.index("--output") + 1])
            for log in sorted(output.rglob("*.log")):
                result.stderr += f"\n--- {log.relative_to(output)} ---\n"
                result.stderr += log.read_text(errors="replace")[-4000:]
        return result

    def base_run(self, output: Path) -> list[str]:
        return [
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            str(self.depth),
            "--output",
            str(output),
            "--retain-min-bin-bp",
            "1",
            "--qc-min-bin-bp",
            "1",
            "--min-contig-bp",
            "1",
            "--threads",
            "1",
        ]

    def test_doctor_checks_current_environment_and_databases(self) -> None:
        completed = self.run_cli(
            "doctor",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--kunpeng-db",
            str(self.kunpeng_db),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertTrue(report["ok"])
        self.assertEqual(report["tool"], "EukBinMiner")
        self.assertEqual(report["execution_mode"], "current_environment_path")
        self.assertGreater(
            report["databases"]["kunpeng"]["shard_bytes"],
            16,
        )
        self.assertEqual(report["databases"]["kunpeng"]["partitions"], 1)

    def test_configured_split_conda_environments_reach_backends(self) -> None:
        output = self.root / "split-conda-environments"
        configured = {
            "CONDA_EXE": str(self.tools / "conda"),
            "REMAG_ENV": str(self.root / "envs" / "mock-remag"),
            "CONCOCT_ENV": str(self.root / "envs" / "mock-concoct"),
            "METABAT2_ENV": str(self.root / "envs" / "mock-metabat2"),
            "EUKCC_ENV": str(self.root / "envs" / "mock-eukcc"),
            "BUSCO_ENV": str(self.root / "envs" / "mock-busco"),
        }
        for suffix, value in configured.items():
            self.environment[f"EUKBINMINER_{suffix}"] = value
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag",
            "--qc-tools",
            "all",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        backend = (output / "config" / "backend.env").read_text(encoding="utf-8")
        for suffix, value in configured.items():
            self.assertIn(f"{suffix}={value}\n", backend)

    def test_partial_binner_failure_busco_decision_and_thresholds(self) -> None:
        output = self.root / "busco-threshold-050"
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag,concoct,metabat2",
            "--qc-tools",
            "all",
            "--decision-qc",
            "busco",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--overlap-threshold",
            "0.50",
            mock_bin_mode="half",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["tool"], "EukBinMiner")
        self.assertEqual(summary["completed_binners"], ["remag", "metabat2"])
        self.assertEqual(summary["representatives"], 1)
        self.assertEqual(summary["alternatives"], 1)
        with (output / "selection" / "representatives.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            representative = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(representative["decision_qc"], "busco")
        self.assertEqual(
            representative["contamination_source"],
            "busco_duplicated",
        )
        self.assertEqual(representative["qc_pass"], "yes")
        self.assertEqual(representative["contamination"], "1")
        self.assertEqual(representative["quality_score"], "93")
        self.assertNotIn("redundancy", representative)
        self.assertNotIn("decision_quality", representative)
        self.assertNotIn("eukcc_quality", representative)
        self.assertEqual(representative["selected_binner"], "metabat2")
        normalized_depth = (
            output / "work" / "input" / "depth.jgi.tsv"
        ).read_text(encoding="utf-8").splitlines()[0]
        self.assertEqual(
            normalized_depth,
            "contigName\tcontigLen\ttotalAvgDepth\t"
            "sample1.bam\tsample1.bam-var",
        )
        mags = sorted((output / "mags").glob("*.fa"))
        self.assertEqual(len(mags), 1)

        strict_output = self.root / "busco-threshold-051"
        strict = self.run_cli(
            *self.base_run(strict_output),
            "--binners",
            "remag,metabat2",
            "--qc-tools",
            "busco",
            "--decision-qc",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--overlap-threshold",
            "0.51",
            mock_bin_mode="half",
        )
        self.assertEqual(strict.returncode, 0, strict.stderr)
        strict_summary = json.loads(
            (strict_output / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(strict_summary["representatives"], 2)
        self.assertEqual(strict_summary["alternatives"], 0)

    def test_repeated_mapping_depths_reach_all_five_binners(self) -> None:
        alpha = self.root / "alpha.depth.tsv"
        beta = self.root / "beta.depth.tsv"
        alpha.write_text(
            "contig_id\tdepth\nc1\t3\nc2\t2\nc3\t1\n",
            encoding="utf-8",
        )
        beta.write_text(
            "contig_id\tdepth\nc1\t30\nc2\t20\nc3\t10\n",
            encoding="utf-8",
        )
        validated = self.run_cli(
            "validate",
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            f"alpha={alpha}",
            "--depth",
            f"beta={beta}",
            "--binners",
            "metabat2",
            "--qc-tools",
            "busco",
            "--decision-qc",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--min-contig-bp",
            "1",
        )
        self.assertEqual(validated.returncode, 0, validated.stderr)
        validation_report = json.loads(validated.stdout)
        self.assertEqual(validation_report["depth_source_count"], 2)
        self.assertEqual(validation_report["coverage_columns"], ["alpha", "beta"])
        self.assertEqual(validation_report["coverage_column_count"], 2)

        output = self.root / "multi-coverage-all-binners"
        completed = self.run_cli(
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            f"alpha={alpha}",
            "--depth",
            f"beta={beta}",
            "--output",
            str(output),
            "--binners",
            "remag,concoct,metabat2,vamb,taxvamb",
            "--taxonomy",
            str(self.taxonomy),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--retain-min-bin-bp",
            "1",
            "--qc-min-bin-bp",
            "1",
            "--min-contig-bp",
            "1",
            "--threads",
            "1",
            expected_coverage_columns=("alpha", "beta"),
            concoct_success=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        failed_diagnostics = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in (output / "results" / "standalone" / "sample1" / "failed").glob(
                "**/*"
            )
            if path.is_file()
        )
        self.assertEqual(
            summary["completed_binners"],
            ["remag", "concoct", "metabat2", "vamb", "taxvamb"],
            completed.stderr
            + (output / "logs" / "binning.classic.log").read_text(
                encoding="utf-8"
            )
            + failed_diagnostics,
        )
        run_metadata = json.loads(
            (output / "config" / "run.json").read_text(encoding="utf-8")
        )
        self.assertEqual(run_metadata["depth_source_count"], 2)
        self.assertEqual(run_metadata["coverage_columns"], ["alpha", "beta"])
        self.assertEqual(run_metadata["coverage_column_count"], 2)

        normalized = (
            output / "work" / "input" / "depth.jgi.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            normalized[0],
            "contigName\tcontigLen\ttotalAvgDepth\talpha\talpha-var\t"
            "beta\tbeta-var",
        )
        self.assertEqual(normalized[1], "c1\t100\t33\t3\t0\t30\t0")

        preflight = output / "results" / "standalone" / "sample1" / "preflight"
        self.assertEqual(
            (preflight / "union.concoct.coverage.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contig_id\talpha\tbeta",
        )
        preflight_metadata = json.loads(
            (preflight / "preflight.json").read_text(encoding="utf-8")
        )
        self.assertEqual(preflight_metadata["depth"]["coverage_column_count"], 2)

        branch_root = (
            output
            / "results"
            / "standalone"
            / "sample1"
            / "results"
            / "union"
        )
        remag_inputs = json.loads(
            (branch_root / "remag" / "mock_coverage_inputs.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            remag_inputs,
            [
                ["c1\t3", "c2\t2", "c3\t1"],
                ["c1\t30", "c2\t20", "c3\t10"],
            ],
        )
        self.assertEqual(
            (branch_root / "concoct" / "fragment_coverage.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contig_id\talpha\tbeta",
        )
        self.assertEqual(
            (output / "inputs" / "r000001" / "abundance.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contigname\talpha\tbeta",
        )
        vamb_metadata = json.loads(
            (output / "inputs" / "r000001" / "input.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(vamb_metadata["coverage_columns"], ["alpha", "beta"])
        self.assertEqual(vamb_metadata["coverage_column_count"], 2)

    def test_single_depth_concoct_keeps_legacy_normalization_policy(self) -> None:
        output = self.root / "single-depth-concoct"
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "concoct",
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            expected_coverage_columns=("sample1.bam",),
            concoct_success=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["completed_binners"], ["concoct"])

    def test_multi_column_matrix_reaches_all_five_binners(self) -> None:
        matrix = self.root / "mapped-samples.coverage.tsv"
        matrix.write_text(
            "contig_id\tdonor_A\tdonor_B\tdonor_C\n"
            "c1\t3\t30\t300\n"
            "c2\t2\t20\t200\n"
            "c3\t1\t10\t100\n",
            encoding="utf-8",
        )
        output = self.root / "multi-column-all-binners"
        completed = self.run_cli(
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            str(matrix),
            "--output",
            str(output),
            "--binners",
            "remag,concoct,metabat2,vamb,taxvamb",
            "--taxonomy",
            str(self.taxonomy),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--retain-min-bin-bp",
            "1",
            "--qc-min-bin-bp",
            "1",
            "--min-contig-bp",
            "1",
            "--threads",
            "1",
            expected_coverage_columns=("donor_A", "donor_B", "donor_C"),
            concoct_success=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(
            summary["completed_binners"],
            ["remag", "concoct", "metabat2", "vamb", "taxvamb"],
        )
        normalized = (
            output / "work" / "input" / "depth.jgi.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            normalized[0],
            "contigName\tcontigLen\ttotalAvgDepth\tdonor_A\tdonor_A-var\t"
            "donor_B\tdonor_B-var\tdonor_C\tdonor_C-var",
        )
        self.assertEqual(normalized[1], "c1\t100\t333\t3\t0\t30\t0\t300\t0")
        preflight = output / "results" / "standalone" / "sample1" / "preflight"
        self.assertEqual(
            (preflight / "union.concoct.coverage.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contig_id\tdonor_A\tdonor_B\tdonor_C",
        )
        self.assertEqual(
            (output / "inputs" / "r000001" / "abundance.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contigname\tdonor_A\tdonor_B\tdonor_C",
        )

    def test_taxvamb_with_user_taxonomy_and_kunpeng(self) -> None:
        taxonomy_output = self.root / "taxvamb-user-taxonomy"
        supplied = self.run_cli(
            *self.base_run(taxonomy_output),
            "--binners",
            "vamb,taxvamb",
            "--taxonomy",
            str(self.taxonomy),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
        )
        self.assertEqual(supplied.returncode, 0, supplied.stderr)
        supplied_summary = json.loads(
            (taxonomy_output / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(supplied_summary["completed_binners"], ["vamb", "taxvamb"])
        self.assertEqual(supplied_summary["representatives"], 1)

        kunpeng_output = self.root / "taxvamb-kunpeng"
        classified = self.run_cli(
            *self.base_run(kunpeng_output),
            "--binners",
            "taxvamb",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
        )
        self.assertEqual(classified.returncode, 0, classified.stderr)
        taxonomy_path = (
            kunpeng_output
            / "taxonomy"
            / "single"
            / "per_sample"
            / "r000001.taxonomy.tsv"
        )
        self.assertTrue(taxonomy_path.is_file())
        self.assertIn("c1", taxonomy_path.read_text(encoding="utf-8"))
        summary = json.loads(
            (kunpeng_output / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(summary["taxonomy"], "kunpeng")
        self.assertEqual(
            summary["taxonomy_provenance"]["backend"],
            "kunpeng",
        )
        self.assertTrue(
            (
                kunpeng_output
                / "taxonomy"
                / "kunpeng"
                / "kunpeng.output.txt.gz"
            ).is_file()
        )

    def test_kunpeng_failures_are_explicit(self) -> None:
        missing_db = self.root / "missing-kunpeng-db"
        missing_db.mkdir()
        invalid = self.run_cli(
            "validate",
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            str(self.depth),
            "--binners",
            "taxvamb",
            "--kunpeng-db",
            str(missing_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--min-contig-bp",
            "1",
        )
        self.assertEqual(invalid.returncode, 2)
        self.assertIn("Kun-peng database is incomplete", invalid.stderr)

        for mode, expected in (
            ("exit", "Kun-peng taxonomy failed"),
            ("missing", "exactly one sequence output"),
            ("duplicate", "duplicate classification"),
            ("truncated", "classified counts differ"),
        ):
            output = self.root / f"kunpeng-{mode}"
            completed = self.run_cli(
                *self.base_run(output),
                "--binners",
                "taxvamb",
                "--kunpeng-db",
                str(self.kunpeng_db),
                "--qc-tools",
                "eukcc",
                "--decision-qc",
                "eukcc",
                "--eukcc-db",
                str(self.eukcc_db),
                mock_kunpeng_mode=mode,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn(expected, completed.stderr)

        zero_hit_output = self.root / "kunpeng-zero-hit"
        completed = self.run_cli(
            *self.base_run(zero_hit_output),
            "--binners",
            "taxvamb",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            mock_kunpeng_mode="zero_hit",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        provenance = json.loads(
            (
                zero_hit_output
                / "taxonomy"
                / "kunpeng"
                / "provenance.json"
            ).read_text(encoding="utf-8")
        )
        compatibility = provenance["compatibility_resolution"]
        self.assertEqual(
            compatibility["raw_omitted_zero_hit_contigs"],
            1,
        )
        self.assertEqual(compatibility["contigs"], 3)
        self.assertEqual(provenance["classifier_exit_code"], 0)
        self.assertFalse(provenance["accepted_nonzero_enodata"])

        enodata_output = self.root / "kunpeng-enodata"
        completed = self.run_cli(
            *self.base_run(enodata_output),
            "--binners",
            "taxvamb",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            mock_kunpeng_mode="enodata",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        provenance = json.loads(
            (
                enodata_output
                / "taxonomy"
                / "kunpeng"
                / "provenance.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(provenance["classifier_exit_code"], 1)
        self.assertTrue(provenance["accepted_nonzero_enodata"])
        self.assertEqual(
            provenance["compatibility_resolution"][
                "raw_omitted_zero_hit_contigs"
            ],
            1,
        )

    def test_kunpeng_rejects_dirty_chunk_directory_and_mutual_exclusion(
        self,
    ) -> None:
        output = self.root / "kunpeng-dirty"
        dirty = output / "work" / "runtime" / "kunpeng-chunks"
        dirty.mkdir(parents=True)
        (dirty / "stale.bin").write_text("stale\n", encoding="utf-8")
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "taxvamb",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("output lacks matching run metadata", completed.stderr)

        mutually_exclusive = self.run_cli(
            "validate",
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            str(self.depth),
            "--binners",
            "taxvamb",
            "--taxonomy",
            str(self.taxonomy),
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--min-contig-bp",
            "1",
        )
        self.assertEqual(mutually_exclusive.returncode, 2)
        self.assertIn("not allowed with argument", mutually_exclusive.stderr)

    def test_zero_qc_candidates_is_a_successful_empty_result(self) -> None:
        output = self.root / "zero-candidates"
        arguments = self.base_run(output)
        qc_index = arguments.index("--qc-min-bin-bp") + 1
        arguments[qc_index] = "1000000"
        completed = self.run_cli(
            *arguments,
            "--binners",
            "metabat2",
            "--qc-tools",
            "all",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["qc_candidates"], 0)
        self.assertEqual(summary["representatives"], 0)
        self.assertEqual(list((output / "mags").glob("*.fa")), [])
        self.assertTrue((output / "mags" / "manifest.tsv").is_file())

    def test_empty_contig_input_is_a_successful_tool_free_result(self) -> None:
        empty_fasta = self.root / "empty.fa"
        empty_fasta.write_bytes(b"")
        unrelated_depth = self.root / "empty.depth.tsv"
        unrelated_depth.write_text(
            "contigName\tcontigLen\ttotalAvgDepth\n"
            "unrelated\t1000\t1\n",
            encoding="utf-8",
        )
        output = self.root / "empty-input"
        completed = self.run_cli(
            "--sample",
            "empty1",
            "--contigs",
            str(empty_fasta),
            "--depth",
            str(unrelated_depth),
            "--output",
            str(output),
            "--binners",
            "all",
            "--qc-tools",
            "all",
            "--decision-qc",
            "eukcc",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--random-seed",
            "6311",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(
            (output / "summary.json").read_text(encoding="utf-8")
        )
        run_metadata = json.loads(
            (output / "config" / "run.json").read_text(encoding="utf-8")
        )
        self.assertTrue(summary["empty_input"])
        self.assertEqual(summary["completed_binners"], [])
        self.assertEqual(summary["taxonomy"], "not_run_empty_input")
        self.assertEqual(summary["discovered_bins"], 0)
        self.assertEqual(summary["representatives"], 0)
        self.assertEqual(run_metadata["contigs"], 0)
        self.assertFalse((output / "logs").exists())
        branch_rows = (
            output
            / "results"
            / "standalone"
            / "empty1"
            / "branch_status.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(branch_rows), 6)
        self.assertTrue(
            all(
                "not_run_empty_input" in row
                for row in branch_rows[1:]
            )
        )
        self.assertTrue((output / "selection" / "representatives.tsv").is_file())
        self.assertTrue((output / "mags" / "manifest.tsv").is_file())

        validated = self.run_cli(
            "validate",
            "--sample",
            "empty1",
            "--contigs",
            str(empty_fasta),
            "--depth",
            str(unrelated_depth),
            "--binners",
            "all",
            "--qc-tools",
            "all",
            "--decision-qc",
            "eukcc",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
        )
        self.assertEqual(validated.returncode, 0, validated.stderr)
        self.assertEqual(json.loads(validated.stdout)["contigs"], 0)

    def test_contigs_below_minimum_are_a_tool_free_empty_result(self) -> None:
        output = self.root / "below-minimum-input"
        arguments = self.base_run(output)
        minimum_index = arguments.index("--min-contig-bp") + 1
        arguments[minimum_index] = "101"
        completed = self.run_cli(
            *arguments,
            "--binners",
            "all",
            "--qc-tools",
            "all",
            "--decision-qc",
            "eukcc",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--random-seed",
            "6311",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(
            (output / "summary.json").read_text(encoding="utf-8")
        )
        run_metadata = json.loads(
            (output / "config" / "run.json").read_text(encoding="utf-8")
        )
        self.assertTrue(summary["empty_input"])
        self.assertEqual(
            summary["empty_input_reason"],
            "no_input_contigs_at_minimum_bp",
        )
        self.assertEqual(summary["completed_binners"], [])
        self.assertEqual(summary["qc_tools"], ["eukcc", "busco"])
        self.assertEqual(summary["random_seed"], 6311)
        self.assertEqual(summary["discovered_bins"], 0)
        self.assertEqual(summary["representatives"], 0)
        self.assertEqual(run_metadata["contigs"], 3)
        self.assertEqual(run_metadata["eligible_contigs"], 0)
        self.assertEqual(run_metadata["eligible_bp"], 0)
        self.assertFalse((output / "logs").exists())
        branch_rows = (
            output
            / "results"
            / "standalone"
            / "sample1"
            / "branch_status.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(branch_rows), 6)
        self.assertTrue(
            all(
                "not_run_empty_input" in row
                for row in branch_rows[1:]
            )
        )
        self.assertTrue((output / ".complete.json").is_file())
        self.assertTrue((output / "mags" / "manifest.tsv").is_file())

    def test_direct_cli_help_and_removed_options(self) -> None:
        help_result = self.run_cli("-h")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("--min-contig-bp", help_result.stdout)
        self.assertIn("--max-contamination", help_result.stdout)
        self.assertIn("--kunpeng-db", help_result.stdout)
        self.assertIn("--random-seed INT", help_result.stdout)
        self.assertIn("usage: eukbinminer", help_result.stdout)
        self.assertIn("default: 6311", help_result.stdout)
        for removed in (
            "--input-min-contig-exclusive",
            "--high-completeness",
            "--medium-completeness",
            "--seqkit",
            "--vamb",
            "--remag-env",
            "--kraken-db",
            "--kraken-memory-mapping",
        ):
            self.assertNotIn(removed, help_result.stdout)

        old_run = self.run_cli("run")
        self.assertNotEqual(old_run.returncode, 0)
        removed_taxonomy = self.run_cli("taxonomy")
        self.assertNotEqual(removed_taxonomy.returncode, 0)

    def test_random_seed_is_unified_and_recorded(self) -> None:
        output = self.root / "random-seed"
        seed = 7311
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag,metabat2,vamb,taxvamb",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--random-seed",
            str(seed),
            expected_seed=seed,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

        run_metadata = json.loads(
            (output / "config" / "run.json").read_text(encoding="utf-8")
        )
        summary = json.loads(
            (output / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(run_metadata["random_seed"], seed)
        self.assertEqual(summary["random_seed"], seed)
        backend_environment = (
            output / "config" / "backend.env"
        ).read_text(encoding="utf-8")
        self.assertIn(f"RANDOM_SEED={seed}\n", backend_environment)

        branch_root = (
            output
            / "results"
            / "standalone"
            / "sample1"
            / "results"
            / "union"
        )
        for binner in ("remag", "metabat2", "vamb", "taxvamb"):
            fingerprint = json.loads(
                (branch_root / binner / "fingerprint.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                fingerprint["fields"]["random_seed"],
                str(seed),
            )

        conflicting_reuse = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag,metabat2,vamb,taxvamb",
            "--kunpeng-db",
            str(self.kunpeng_db),
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--random-seed",
            "8311",
        )
        self.assertEqual(conflicting_reuse.returncode, 2)
        self.assertIn(
            "completed output random seed differs",
            conflicting_reuse.stderr,
        )

        invalid = self.run_cli(
            "validate",
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--depth",
            str(self.depth),
            "--binners",
            "metabat2",
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--random-seed",
            "0",
        )
        self.assertEqual(invalid.returncode, 2)
        self.assertIn("--random-seed must be positive", invalid.stderr)

    def test_quality_thresholds_are_inclusive(self) -> None:
        accepted_output = self.root / "quality-boundary-accepted"
        accepted = self.run_cli(
            *self.base_run(accepted_output),
            "--binners",
            "metabat2",
            "--qc-tools",
            "busco",
            "--decision-qc",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--min-completeness",
            "98",
            "--max-contamination",
            "1",
        )
        self.assertEqual(accepted.returncode, 0, accepted.stderr)
        accepted_summary = json.loads(
            (accepted_output / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(accepted_summary["representatives"], 1)

        rejected_output = self.root / "quality-boundary-rejected"
        rejected = self.run_cli(
            *self.base_run(rejected_output),
            "--binners",
            "metabat2",
            "--qc-tools",
            "busco",
            "--decision-qc",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--min-completeness",
            "98.1",
            "--max-contamination",
            "1",
        )
        self.assertEqual(rejected.returncode, 0, rejected.stderr)
        rejected_summary = json.loads(
            (rejected_output / "summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(rejected_summary["representatives"], 0)

        eukcc_output = self.root / "eukcc-contamination-boundary"
        eukcc = self.run_cli(
            *self.base_run(eukcc_output),
            "--binners",
            "metabat2",
            "--qc-tools",
            "eukcc",
            "--decision-qc",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--min-completeness",
            "91",
            "--max-contamination",
            "1",
        )
        self.assertEqual(eukcc.returncode, 0, eukcc.stderr)
        with (eukcc_output / "selection" / "representatives.tsv").open(
            encoding="utf-8",
            newline="",
        ) as handle:
            representative = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            representative["contamination_source"],
            "eukcc_contamination",
        )
        self.assertEqual(representative["contamination"], "1")


if __name__ == "__main__":
    unittest.main()
