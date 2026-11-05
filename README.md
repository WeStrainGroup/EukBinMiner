# EukBinMiner

EukBinMiner is a standalone command-line workflow for mining,
quality-controlling, and selecting non-redundant eukaryotic bins from one
sample's contigs plus one or more mapped-sample depth/coverage dimensions. It is distributed as a standalone Python/Bash CLI.

## Workflow

```text
contigs + one or more depth/coverage dimensions
  -> Kun-peng taxonomy when TaxVAMB is selected
  -> one or more binners (default: all)
  -> seqkit bin statistics
  -> EukCC and/or BUSCO
  -> within-sample, cross-binner overlap comparison
  -> best representative
  -> non-redundant representative-bin FASTA files
```

Supported binners are `remag`, `concoct`, `metabat2`, `vamb`, and `taxvamb`.
A failed binner does not stop the workflow when another requested binner
produced a usable branch. Bins emitted by the same binner are never compared
to one another.

Two bins from different binners are considered overlapping when

```text
shared_contig_bp / min(bin_a_bp, bin_b_bp) >= overlap_threshold
```

The default threshold is `0.50`. A bin is eligible by default when
completeness is at least 50% and contamination is at most 10%. Eligible
representatives are ordered by quality score
(`completeness - contamination_penalty * contamination`), completeness,
contamination, contig count, N50, and bin ID.

EukCC supplies its native contamination percentage. When BUSCO controls the
decision, EukBinMiner uses BUSCO duplicated percentage as a contamination
proxy. These measurements are useful for the same filtering role but are not
biologically identical.

Reproducible binning uses one public seed, `--random-seed` (default `6311`).
The parsed value is supplied to ReMAG, CONCOCT, MetaBAT2, VAMB, and TaxVAMB
and is recorded in run metadata and every binner branch fingerprint.
Kun-peng does not expose or receive a random-seed argument.

## Calculation settings

| Stage | Applied behavior |
| --- | --- |
| Input | Retain contigs at least 1,500 bp by default. Every retained contig must have coverage in each supplied mapping dimension. |
| ReMAG | Short-read mode, bacterial filtering skipped, seed 6311. The supplied contigs define the candidate input set. |
| CONCOCT | Cut into 10-kb fragments with the last fragment merged. Each fragment inherits its parent contig coverage. Start with half the eligible fragment count, bounded to 2–400 clusters. Fewer than two eligible fragments produce zero bins. |
| CONCOCT coverage | One dimension disables coverage normalization and total coverage. Multiple dimensions use CONCOCT's differential-coverage defaults. |
| MetaBAT2 | Uses the JGI mean/variance pairs. Generic coverage input supplies zero variance placeholders. |
| VAMB / TaxVAMB | Default training schedule, bin splitting disabled, seed applied before model import. All-zero coverage columns are removed and recorded. |
| Raw bins | Minimum retained size defaults to 1 bp. Bins at least 500,000 bp enter QC by default. |
| EukCC | Native completeness and contamination estimates. Verified no-marker or empty-protein outcomes retain an unassessed status. |
| BUSCO | Genome mode with the supplied lineage, offline operation and archived intermediate files. Duplicated BUSCO percentage supplies the decision proxy when selected. |
| Selection | Completeness ≥50%, contamination/proxy ≤10%, score = completeness − 5 × contamination/proxy. Compare only different binners within the current assembly. |

Selection proceeds in descending quality order and compares each candidate
against the representatives already retained. A candidate is suppressed only
when it directly overlaps a retained representative. Each suppressed candidate
records its directly overlapping representative and any ambiguity.
The contig-ID and normalized-sequence pair defines a shared source contig.

## Installation

Clone the source and put `bin/` on `PATH`:

```bash
git clone https://github.com/WeStrainGroup/EukBinMiner.git
export PATH="$PWD/EukBinMiner/bin:${PATH}"
eukbinminer --version
```

EukBinMiner discovers tools from the current `PATH`. The combined environment
specification is `envs/eukbinminer.yml`. The earlier EukCC 2.1.3 / BUSCO 6
dependency solve found incompatible MetaEuk requirements, so this combined
configuration still requires validation. On systems that keep tools in separate named Conda
environments, set `EUKBINMINER_REMAG_ENV`, `EUKBINMINER_CONCOCT_ENV`,
`EUKBINMINER_METABAT2_ENV`, `EUKBINMINER_EUKCC_ENV`, and
`EUKBINMINER_BUSCO_ENV`. The corresponding backend is then run with
`conda run --no-capture-output -n ENV`; an absolute environment path uses
`conda run -p PREFIX`. VAMB/TaxVAMB, Kun-peng, and seqkit remain resolved from
`PATH`. `EUKBINMINER_CONDA_EXE` may be an absolute Conda executable. Legacy
`EUKMAGMINER_*_ENV` names are also accepted. Use `eukbinminer doctor` after
configuring the runtime to check every selected executable and reference
database.

The repository contains a local-development EukBinMiner recipe plus a
companion Kun-peng v0.7.12 recipe built from the tagged MIT-licensed source.
The source uses CC BY-NC-SA 4.0. [GitHub](https://github.com/WeStrainGroup/EukBinMiner). The Conda recipe and combined environment are development
artifacts. A combined EukCC/BUSCO environment still needs a validated packaging solution.
The per-tool specifications are in `envs/`.

The old `eukmagminer` executable remains as a deprecated compatibility alias
for existing scripts. New commands, documentation, summaries, and package
metadata use `EukBinMiner`/`eukbinminer`.

## Reference databases

Pass the unpacked database directories directly:

```text
--eukcc-db /path/to/eukcc_db
--busco-db /path/to/fungi_odb12
```

The equivalent environment variables are `EUKBINMINER_EUKCC_DB` and
`EUKBINMINER_BUSCO_DB`. `--busco-db` must point to the unpacked BUSCO lineage
directory containing `dataset.cfg` and `hmms/`.

TaxVAMB accepts either a two-column taxonomy table (`--taxonomy`) or a
Kun-peng sharded database (`--kunpeng-db`; environment variable
`EUKBINMINER_KUNPENG_DB`). The two inputs are mutually exclusive.
The corresponding legacy `EUKMAGMINER_*` names remain readable when the new
name is unset.

Kun-peng databases contain `hash_config.k2d`, `opts.k2d`, `taxo.k2d`, and a
contiguous set of `hash_1.k2d` through `hash_N.k2d`. `doctor` and `validate`
check the config, shard numbering, headers, capacities, and exact file sizes.
An existing Kraken2 database can be converted once with Kun-peng
`hashshard`; never run the converter against a shared read-only source
directory.

Kun-peng can load one hash shard per worker. Choose the shard size and CPU
count together when planning memory. Stage reference data and inputs according
to the storage policy of your compute system.

When TaxVAMB is selected without `--taxonomy`, EukBinMiner runs Kun-peng
inside the ordinary command:

```text
filtered contigs -> kun_peng classify
  -> Kraken-compatible run-length/LCA resolver
  -> TaxVAMB taxonomy -> TaxVAMB
```

Kun-peng receives `--minimum-hit-groups 3`. Its chunk directory is
task-private below the EukBinMiner runtime directory. EukBinMiner requires
exactly one classification output and one `kreport2`, restores only
report-proven zero-hit omissions, rejects duplicate or unexplained missing
contigs, and validates the index options and compact taxonomy.
Kun-peng v0.7.12 reports one hit per minimizer-change event, whereas Kraken
weighs every k-mer window covered by that event. EukBinMiner restores those run
lengths from the input FASTA and `opts.k2d`, reapplies the three-hit gate, and
resolves the call against `taxo.k2d` before TaxVAMB. Both the compressed raw
Kun-peng output and compatibility output are retained. Provenance records the
number of corrected calls, executable SHA256, reported version, release
provenance, database layout, and output checksums. User-supplied taxonomy TSVs
remain supported.

Kun-peng v0.7.12 omits sequence-output rows for contigs with no database hit,
although those contigs remain counted as unclassified in `kreport2`.
EukBinMiner restores such rows as unclassified only after the report totals,
classified count, raw output, and FASTA IDs jointly prove that every missing
record is a zero-hit omission. Any inconsistent or truncated output still
fails closed. The number restored is recorded in
`raw_omitted_zero_hit_contigs`.

VAMB 5.0.4 records `--seed`, but its Taxometer predictor initializes the
`VAMB2Label` network before applying that seed. EukBinMiner therefore launches
VAMB through a small audited shim that seeds Python, NumPy, and PyTorch before
the VAMB console script is imported. The shim and VAMB both receive the
resolved `--random-seed` value. The seed, seed policy, and launcher checksum
are part of every VAMB/TaxVAMB branch fingerprint.

## Check an installation

`doctor` checks the selected executables and database layout without running
binning:

```bash
eukbinminer doctor \
  --binners all \
  --qc-tools all \
  --eukcc-db /path/to/eukcc_db \
  --busco-db /path/to/fungi_odb12 \
  --kunpeng-db /path/to/kunpeng_db
```

It prints JSON and returns 0 only when every required check passes.
For Kun-peng, the report includes the parsed config, partition count, logical
capacity, and total shard bytes.

## Validate and run

```bash
eukbinminer validate \
  --sample SAMPLE_001 \
  --contigs contigs.fa \
  --depth depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --busco-db /path/to/fungi_odb12 \
  --kunpeng-db /path/to/kunpeng_db

eukbinminer \
  --sample SAMPLE_001 \
  --contigs contigs.fa \
  --depth depth.tsv \
  --output eukbinminer-output \
  --binners all \
  --qc-tools all \
  --decision-qc eukcc \
  --min-contig-bp 1500 \
  --retain-min-bin-bp 1 \
  --qc-min-bin-bp 500000 \
  --min-completeness 50 \
  --max-contamination 10 \
  --contamination-penalty 5 \
  --overlap-threshold 0.50 \
  --random-seed 6311 \
  --eukcc-db /path/to/eukcc_db \
  --busco-db /path/to/fungi_odb12 \
  --kunpeng-db /path/to/kunpeng_db
```

For differential coverage from reads of several samples mapped back to the
same contig assembly, repeat `--depth`. `NAME=TSV` is recommended for generic
two-column files because it records the intended sample dimension explicitly:

```bash
eukbinminer \
  --sample ASSEMBLY_SAMPLE \
  --contigs contigs.fa \
  --depth sample_A=sample_A.depth.tsv \
  --depth sample_B=sample_B.depth.tsv \
  --depth sample_C=sample_C.depth.tsv \
  --output eukbinminer-output \
  --eukcc-db /path/to/eukcc_db \
  --busco-db /path/to/fungi_odb12 \
  --kunpeng-db /path/to/kunpeng_db
```

The equivalent single-matrix form is:

```text
contig_id  sample_A  sample_B  sample_C
contig_1   2.5       0.0       8.1
contig_2   4.0       1.2       3.3
```

and is supplied once as `--depth coverage.tsv`.

`--qc-tools` accepts `eukcc`, `busco`, `eukcc,busco`, or `all`.
`--decision-qc` selects which executed QC tool controls final eligibility and
representative ordering. Contigs enter binning when their length is greater
than or equal to `--min-contig-bp` (default 1500 bp).
`--random-seed` accepts a positive integer and defaults to `6311`.
If that inclusive filter leaves no eligible contigs, the run writes a
checksum-bearing, tool-free empty result; it records all requested binners as
`not_run_empty_input` rather than claiming that they executed.

Each depth input may be a JGI table, a generic two-column table, or a generic
matrix containing an ID column, optional `contigLen`, and one or more numeric
coverage columns. `--depth` may be repeated for independent mapping outputs;
every supplied file must contain every FASTA contig. A `NAME=TSV` label is
valid for a single-coverage file. Without an explicit label, repeated generic
files use a meaningful coverage header when present and otherwise infer the
dimension name from the filename. Input order and column order are preserved.

The canonical JGI table retains all sample depth/variance pairs. For generic
input, variance is zero and `totalAvgDepth` is the sum of all coverage
dimensions; a single JGI input keeps its supplied `totalAvgDepth`. MetaBAT2
receives the complete JGI table. CONCOCT receives every depth column. VAMB/TaxVAMB receive the nonzero
coverage dimensions, with all-zero columns recorded in input metadata. ReMAG 0.4.4 receives one headerless TSV per coverage dimension through
repeated `-c` options, matching its native multi-sample interface. CONCOCT's
single-dimension normalization workaround remains unchanged, while two or
more dimensions use CONCOCT's default differential-coverage normalization.

FASTA IDs, depth IDs, and coverage names must be unique. Rows must have exactly
the declared number of fields. Every coverage, JGI variance, and
`totalAvgDepth` value must be finite and non-negative, and any supplied
`contigLen` must match the FASTA. Duplicate depth IDs are fatal even when rows
are identical. Production adapters may canonicalize read-only legacy inputs,
but the public CLI always validates the canonical copy strictly. The resolved
coverage names/count, source files, total-depth policy, and checksums are
recorded in validation and run metadata.

## Important output

- `mags/`: copied non-redundant representative-bin FASTA files and
  `manifest.tsv` (the path is retained for output-schema compatibility)
- `selection/representatives.tsv`: selected representatives and their metrics
- `selection/alternatives.tsv`: suppressed cross-binner alternatives
- `qc/`: seqkit, EukCC, BUSCO, and decision tables
- `taxonomy/kunpeng/`: raw/compatibility outputs, report, and provenance
- `results/`: retained successful raw binner branches
- `summary.json`: run-level machine-readable summary
- `config/run.json`: resolved inputs, parameters, checksums, and random seed
- `.complete.json`: checksum-bearing completion marker

EukBinMiner never deletes raw bins. Reusing an output directory is allowed
only when its validated completion state is reusable; stale or partial branch
outputs are rejected rather than overwritten.

Kun-peng v0.7.12 can return ENODATA/OS error 61 after producing valid output.
The compatibility path requires the precise error signature and full agreement
between the FASTA, per-contig output and kreport2. Both complete outputs and
report-proven zero-hit omissions are accepted. Exit status and evidence mode
are recorded in taxonomy provenance.

Public QC result fields are `completeness`, `contamination`,
`contamination_source`, `qc_pass`, and `quality_score`.
`contamination_source` is `eukcc_contamination` or `busco_duplicated`.

## Development

Run the local tests with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

The suite exercises mock external tools, real file conversion, QC parsing,
overlap selection, failure isolation, rerun validation and installed layouts.
Historical real-data runs use frozen versions. The changes in this source
release are validated by software regression tests.

## Review changes

- Shared overlap uses the original contig ID plus normalized sequence identity.
  Distinct source IDs with identical sequences retain their individual lengths.
- Reuse checks the input checksums, configuration, CLI/backend code and final
  output checksums. An output-directory lock prevents concurrent writers.
  Reference database paths must point to immutable datasets.
- Invalid QC metrics and inconsistent BUSCO summaries are rejected.
  EukCC bins with verified empty MetaEuk proteins remain explicitly unassessed.
- VAMB drops coverage columns that are zero across every selected contig and
  records their names. A branch with no informative coverage fails.
- A TaxVAMB classification failure allows independently successful binners to
  finish. Shell-stage failures propagate before a branch can be published.
- Installed distributions include backend environment specifications.

## License

CC BY-NC-SA 4.0. See [LICENSE](LICENSE).
External tools and reference datasets retain their respective licenses.
