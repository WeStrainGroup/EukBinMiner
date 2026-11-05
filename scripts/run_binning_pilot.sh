#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${1:-${PROJECT_DIR}/config/SAMN11394100.env}"

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

note() {
    printf '[%s] %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
}

[[ -s "$CONFIG_PATH" ]] || die "Config does not exist or is empty: $CONFIG_PATH"
# shellcheck source=/dev/null
source "$CONFIG_PATH"

: "${STUDY:?Missing STUDY in config}"
: "${SAMPLE:?Missing SAMPLE in config}"
: "${OUT_ROOT:?Missing OUT_ROOT in config}"
: "${PYTHON_EXE:=python3}"
: "${CONDA_EXE:=conda}"
: "${REMAG_ENV:=}"
: "${CONCOCT_ENV:=}"
: "${METABAT2_ENV:=}"
: "${THREADS:=${SLURM_CPUS_PER_TASK:-16}}"
: "${MIN_BIN_SIZE:=1}"
: "${INPUT_MIN_CONTIG_BP:=1500}"
: "${CONCOCT_CHUNK_SIZE:=10000}"
: "${REMAG_MIN_CONTIG:=${INPUT_MIN_CONTIG_BP}}"
: "${CONCOCT_MIN_CONTIG:=${INPUT_MIN_CONTIG_BP}}"
: "${METABAT2_MIN_CONTIG:=${INPUT_MIN_CONTIG_BP}}"
: "${CONCOCT_MAX_CLUSTERS:=400}"
: "${RANDOM_SEED:=6311}"
: "${INPUT_SETS:=union}"
: "${BINNERS:=remag concoct metabat2}"
: "${REMAG_KEEP_INTERMEDIATE:=0}"
: "${CROSS_ACCOUNT_COORDINATOR_USER:=$(id -un)}"

require_positive_integer() {
    local name="$1"
    local value="${!name}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] ||
        die "${name} must be a positive integer, got: ${value}"
}

for integer_setting in \
    THREADS MIN_BIN_SIZE CONCOCT_CHUNK_SIZE REMAG_MIN_CONTIG \
    CONCOCT_MIN_CONTIG METABAT2_MIN_CONTIG CONCOCT_MAX_CLUSTERS \
    RANDOM_SEED; do
    require_positive_integer "$integer_setting"
done
require_positive_integer INPUT_MIN_CONTIG_BP
[[ "$REMAG_KEEP_INTERMEDIATE" == "0" || "$REMAG_KEEP_INTERMEDIATE" == "1" ]] ||
    die "REMAG_KEEP_INTERMEDIATE must be 0 or 1"

if [[ -z "${CONTIG_FASTA:-}" ]]; then
    : "${TIARA_ROOT:?Missing TIARA_ROOT unless CONTIG_FASTA is provided}"
    TIARA_SAMPLE_DIR="${TIARA_ROOT}/${STUDY}/${SAMPLE}"
fi
if [[ -z "${DEPTH_TSV:-}" ]]; then
    : "${RAW_ROOT:?Missing RAW_ROOT unless DEPTH_TSV is provided}"
    RAW_SAMPLE_DIR="${RAW_ROOT}/${STUDY}/${SAMPLE}"
fi
SAMPLE_OUT="${OUT_ROOT}/${STUDY}/${SAMPLE}"
PREP_DIR="${SAMPLE_OUT}/preflight"
RESULTS_DIR="${SAMPLE_OUT}/results"
FAILED_DIR="${SAMPLE_OUT}/failed"
LOG_DIR="${SAMPLE_OUT}/logs"
SUMMARY_DIR="${SAMPLE_OUT}/summary"

mkdir -p "$SAMPLE_OUT" "$RESULTS_DIR" "$FAILED_DIR" "$LOG_DIR"

share_sample_output() {
    [[ -d "$SAMPLE_OUT" ]] || return 0
    chmod -R g+rwX,o-rwx "$SAMPLE_OUT" 2>/dev/null || true
    if command -v setfacl >/dev/null 2>&1; then
        setfacl -R -m \
            "u:${CROSS_ACCOUNT_COORDINATOR_USER}:rwX" \
            "$SAMPLE_OUT" 2>/dev/null || true
        find "$SAMPLE_OUT" -type d -exec setfacl -m \
            "d:u:${CROSS_ACCOUNT_COORDINATOR_USER}:rwx" {} + \
            2>/dev/null || true
    fi
}
trap share_sample_output EXIT

discover_one() {
    local label="$1"
    shift
    local pattern
    local candidate
    local -a candidates=()
    local -a unique_candidates=()

    [[ -d "$TIARA_SAMPLE_DIR" ]] || die "Tiara sample directory not found: $TIARA_SAMPLE_DIR"
    for pattern in "$@"; do
        while IFS= read -r -d '' candidate; do
            candidates+=("$candidate")
        done < <(find "$TIARA_SAMPLE_DIR" -maxdepth 1 -type f -name "$pattern" -print0)
    done

    if ((${#candidates[@]})); then
        while IFS= read -r candidate; do
            unique_candidates+=("$candidate")
        done < <(printf '%s\n' "${candidates[@]}" | LC_ALL=C sort -u)
    fi
    ((${#unique_candidates[@]} == 1)) ||
        die "Expected exactly one ${label} FASTA in ${TIARA_SAMPLE_DIR}; found ${#unique_candidates[@]}. Set its path explicitly in the config."
    printf '%s\n' "${unique_candidates[0]}"
}

if [[ -z "${CONTIG_FASTA:-}" ]]; then
    if [[ -z "${EUK_FASTA:-}" ]]; then
        EUK_FASTA="$(discover_one eukarya \
            'eukarya_*.fa' 'eukarya_*.fasta' 'eukarya_*.fna' \
            'eukarya_*.fa.gz' 'eukarya_*.fasta.gz' 'eukarya_*.fna.gz')"
    fi
    if [[ -z "${UNKNOWN_FASTA:-}" ]]; then
        UNKNOWN_FASTA="$(discover_one unknown \
            'unknown_*.fa' 'unknown_*.fasta' 'unknown_*.fna' \
            'unknown_*.fa.gz' 'unknown_*.fasta.gz' 'unknown_*.fna.gz')"
    fi
fi
if [[ -z "${DEPTH_TSV:-}" ]]; then
    DEPTH_TSV="${RAW_SAMPLE_DIR}/${SAMPLE}_contig_depths.tsv"
fi

if [[ -n "${CONTIG_FASTA:-}" ]]; then
    [[ -s "$CONTIG_FASTA" ]] ||
        die "Contig FASTA does not exist or is empty: $CONTIG_FASTA"
else
    [[ -s "$EUK_FASTA" ]] ||
        die "Eukarya FASTA does not exist or is empty: $EUK_FASTA"
    [[ -s "$UNKNOWN_FASTA" ]] ||
        die "Unknown FASTA does not exist or is empty: $UNKNOWN_FASTA"
fi
[[ -s "$DEPTH_TSV" ]] || die "Depth table does not exist or is empty: $DEPTH_TSV"
command -v "$PYTHON_EXE" >/dev/null 2>&1 || die "Python executable not found: $PYTHON_EXE"
if [[ -n "$REMAG_ENV" || -n "$CONCOCT_ENV" || -n "$METABAT2_ENV" ]]; then
    command -v "$CONDA_EXE" >/dev/null 2>&1 ||
        die "Conda executable not found: $CONDA_EXE"
fi

note "Preflight: validating selected contigs against the depth table"
prepare_args=(
    --depth-tsv "$DEPTH_TSV"
    --sample "$SAMPLE"
    --min-contig-bp "$INPUT_MIN_CONTIG_BP"
    --out-dir "$PREP_DIR"
)
if [[ -n "${CONTIG_FASTA:-}" ]]; then
    prepare_args+=(--contigs "$CONTIG_FASTA")
else
    prepare_args+=(
        --euk-fasta "$EUK_FASTA"
        --unknown-fasta "$UNKNOWN_FASTA"
    )
fi
"$PYTHON_EXE" "${SCRIPT_DIR}/prepare_inputs.py" "${prepare_args[@]}"

union_contig_count="$(
    "$PYTHON_EXE" -c \
        'import json, sys; print(json.load(open(sys.argv[1]))["collections"]["union"]["contig_count"])' \
        "${PREP_DIR}/preflight.json"
)"
if ((union_contig_count == 0)); then
    note "No Tiara union contigs remain at >=${INPUT_MIN_CONTIG_BP} bp; record a valid zero-bin sample"
    STATUS_TMP="$(mktemp "${SAMPLE_OUT}/.branch_status.XXXXXX")"
    printf 'input_set\tbinner\tstatus\texit_code\tresult_path\n' > "$STATUS_TMP"
    summary_args=()
    IFS=' ' read -r -a INPUT_SET_LIST <<< "$INPUT_SETS"
    IFS=' ' read -r -a BINNER_LIST <<< "$BINNERS"
    for input_set in "${INPUT_SET_LIST[@]}"; do
        for binner in "${BINNER_LIST[@]}"; do
            bins_dir="${RESULTS_DIR}/${input_set}/${binner}/bins"
            mkdir -p "$bins_dir"
            printf '%s\t%s\tskipped_no_eligible_contigs\t0\t%s\n' \
                "$input_set" "$binner" "$(dirname -- "$bins_dir")" \
                >> "$STATUS_TMP"
            summary_args+=(--bins-label "${input_set}.${binner}=${bins_dir}")
        done
    done
    mv "$STATUS_TMP" "${SAMPLE_OUT}/branch_status.tsv"
    "$PYTHON_EXE" "${SCRIPT_DIR}/summarize_bins.py" \
        "${summary_args[@]}" \
        --out-dir "$SUMMARY_DIR"
    note "Pilot finished with zero eligible contigs: ${SAMPLE_OUT}"
    exit 0
fi

run_in_env() {
    local env_name="$1"
    shift
    if [[ -n "$env_name" ]]; then
        local selector=-n
        [[ "$env_name" == /* ]] && selector=-p
        "$CONDA_EXE" run --no-capture-output "$selector" "$env_name" "$@"
    else
        "$@"
    fi
}

select_inputs() {
    local input_set="$1"
    local source_fasta
    local source_metabat_depth
    case "$input_set" in
        euk)
            source_fasta="${PREP_DIR}/euk.fa"
            source_metabat_depth="${PREP_DIR}/euk.metabat.depth.tsv"
            ;;
        union)
            source_fasta="${PREP_DIR}/union.fa"
            source_metabat_depth="${PREP_DIR}/union.metabat.depth.tsv"
            ;;
        *)
            die "Unsupported input set: $input_set"
            ;;
    esac
    FILTER_MANIFEST=""
    if ((INPUT_MIN_CONTIG_BP > 1)); then
        FILTER_DIR="${PREP_DIR}/${input_set}.ge${INPUT_MIN_CONTIG_BP}"
        "$PYTHON_EXE" "${SCRIPT_DIR}/filter_prepared_inputs.py" \
            --input-fasta "$source_fasta" \
            --metabat-depth "$source_metabat_depth" \
            --sample "$SAMPLE" \
            --min-contig-bp "$INPUT_MIN_CONTIG_BP" \
            --out-dir "$FILTER_DIR"
        INPUT_FASTA="${FILTER_DIR}/input.fa"
        METABAT_DEPTH="${FILTER_DIR}/metabat.depth.tsv"
        REMAG_COVERAGE="${FILTER_DIR}/remag.coverage.tsv"
        CONCOCT_COVERAGE="${FILTER_DIR}/concoct.coverage.tsv"
        FILTER_MANIFEST="${FILTER_DIR}/filter.json"
    else
        INPUT_FASTA="$source_fasta"
        METABAT_DEPTH="$source_metabat_depth"
        REMAG_COVERAGE="${PREP_DIR}/${input_set}.remag.coverage.tsv"
        CONCOCT_COVERAGE="${PREP_DIR}/${input_set}.concoct.coverage.tsv"
    fi
}

get_tool_identity() {
    local binner="$1"
    local output
    local identity
    local needle
    case "$binner" in
        remag)
            needle="remag"
            output="$(run_in_env "$REMAG_ENV" remag --version 2>&1 || true)"
            ;;
        concoct)
            needle="concoct"
            output="$(run_in_env "$CONCOCT_ENV" concoct --version 2>&1 || true)"
            ;;
        metabat2)
            needle="metabat"
            output="$(run_in_env "$METABAT2_ENV" metabat2 --help 2>&1 || true)"
            ;;
        *)
            die "Unsupported binner for tool identity: $binner"
            ;;
    esac
    identity="$(
        printf '%s\n' "$output" |
            awk -v needle="$needle" '
                index(tolower($0), needle) {
                    print
                    exit
                }
            '
    )"
    if [[ -z "$identity" ]]; then
        identity="$(printf '%s\n' "$output" | awk 'NF { print; exit }')"
    fi
    [[ -n "$identity" ]] ||
        die "Could not determine tool identity for ${binner}"
    printf '%s\n' "$identity"
}

write_branch_fingerprint() {
    local input_set="$1"
    local binner="$2"
    local output_path="$3"
    local coverage_path
    local coverage_column_count
    local env_name
    local env_spec
    local tool_identity
    local -a fingerprint_args=()

    tool_identity="$(get_tool_identity "$binner")"
    fingerprint_args=(
        --field "study=${STUDY}"
        --field "sample=${SAMPLE}"
        --field "input_set=${input_set}"
        --field "binner=${binner}"
        --field "threads=${THREADS}"
        --field "min_bin_size=${MIN_BIN_SIZE}"
        --field "input_min_contig_bp=${INPUT_MIN_CONTIG_BP}"
        --field "tool_identity=${tool_identity}"
        --file "input_fasta=${INPUT_FASTA}"
        --file "preflight_manifest=${PREP_DIR}/preflight.json"
        --file "runner=${SCRIPT_DIR}/run_binning_pilot.sh"
        --file "auditor=${SCRIPT_DIR}/audit_branch.py"
        --file "fingerprint_helper=${SCRIPT_DIR}/branch_fingerprint.py"
    )
    if [[ -n "$FILTER_MANIFEST" ]]; then
        fingerprint_args+=(
            --file "input_filter=${SCRIPT_DIR}/filter_prepared_inputs.py"
            --file "input_filter_manifest=${FILTER_MANIFEST}"
        )
    fi

    case "$binner" in
        remag)
            coverage_path="$CONCOCT_COVERAGE"
            coverage_column_count="$(
                awk -F '\t' 'NR == 1 { print NF - 1; exit }' \
                    "$CONCOCT_COVERAGE"
            )"
            env_name="$REMAG_ENV"
            env_spec="${PROJECT_DIR}/envs/remag.yml"
            fingerprint_args+=(
                --field "environment=${env_name}"
                --field "mode=short-reads"
                --field "skip_bacterial_filter=true"
                --field "min_contig=${REMAG_MIN_CONTIG}"
                --field "random_seed=${RANDOM_SEED}"
                --field "keep_intermediate=${REMAG_KEEP_INTERMEDIATE}"
                --field "coverage_columns=${coverage_column_count}"
                --field "coverage_adapter=one_headerless_tsv_per_column"
                --file "coverage_splitter=${SCRIPT_DIR}/split_remag_coverage.py"
            )
            ;;
        concoct)
            coverage_path="$CONCOCT_COVERAGE"
            coverage_column_count="$(
                awk -F '\t' 'NR == 1 { print NF - 1; exit }' \
                    "$CONCOCT_COVERAGE"
            )"
            env_name="$CONCOCT_ENV"
            env_spec="${PROJECT_DIR}/envs/concoct.yml"
            fingerprint_args+=(
                --field "environment=${env_name}"
                --field "chunk_size=${CONCOCT_CHUNK_SIZE}"
                --field "min_contig_inclusive=${CONCOCT_MIN_CONTIG}"
                --field "length_threshold=$((CONCOCT_MIN_CONTIG - 1))"
                --field "max_clusters=${CONCOCT_MAX_CLUSTERS}"
                --field "random_seed=${RANDOM_SEED}"
                --field "coverage_columns=${coverage_column_count}"
                --field "coverage_normalization=$(
                    if ((coverage_column_count == 1)); then
                        printf 'disabled_single_dimension'
                    else
                        printf 'concoct_default_multi_dimension'
                    fi
                )"
                --file "bin_selector=${SCRIPT_DIR}/select_bins_by_size.py"
            )
            ;;
        metabat2)
            coverage_path="$METABAT_DEPTH"
            env_name="$METABAT2_ENV"
            env_spec="${PROJECT_DIR}/envs/metabat2.yml"
            fingerprint_args+=(
                --field "environment=${env_name}"
                --field "min_contig=${METABAT2_MIN_CONTIG}"
                --field "random_seed=${RANDOM_SEED}"
            )
            ;;
        *)
            die "Unsupported binner for fingerprint: $binner"
            ;;
    esac

    fingerprint_args+=(
        --file "coverage=${coverage_path}"
        --file "environment_spec=${env_spec}"
        --output "$output_path"
    )
    "$PYTHON_EXE" "${SCRIPT_DIR}/branch_fingerprint.py" \
        "${fingerprint_args[@]}"
}

validate_completed_branch() {
    local input_set="$1"
    local binner="$2"
    local final_dir="$3"
    local expected_fingerprint="$4"
    local audit_check_dir
    local audit_file

    [[ -f "${final_dir}/fingerprint.json" ]] ||
        die "Completed branch lacks fingerprint: $final_dir"
    cmp -s "$expected_fingerprint" "${final_dir}/fingerprint.json" ||
        die "Completed branch fingerprint differs for ${input_set}/${binner}. Use a new OUT_ROOT or move the old result aside; it will not be overwritten."
    [[ -d "${final_dir}/bins" ]] ||
        die "Completed branch lacks bins directory: $final_dir"
    for audit_file in audit.json assignments.tsv unbinned.list unbinned.fa; do
        [[ -f "${final_dir}/audit/${audit_file}" ]] ||
            die "Completed branch lacks audit/${audit_file}: $final_dir"
    done

    audit_check_dir="$(mktemp -d "${final_dir}/.audit_check.XXXXXX")"
    (
        cd "$final_dir"
        "$PYTHON_EXE" "${SCRIPT_DIR}/audit_branch.py" \
            --source-fasta "$INPUT_FASTA" \
            --bins-dir bins \
            --label "${input_set}.${binner}" \
            --out-dir "$audit_check_dir"
    )
    for audit_file in audit.json assignments.tsv unbinned.list unbinned.fa; do
        if ! cmp -s \
            "${audit_check_dir}/${audit_file}" \
            "${final_dir}/audit/${audit_file}"; then
            rm -rf -- "$audit_check_dir"
            die "Completed branch audit differs for ${input_set}/${binner}; refusing to reuse a stale or damaged result."
        fi
    done
    rm -rf -- "$audit_check_dir"
}

run_remag() {
    local stage_dir="$1"
    local contig_count
    local coverage_file
    local remag_coverage_dir
    local remag_manifest
    local -a command
    local -a coverage_args=()
    mkdir -p "${stage_dir}/bins"
    remag_coverage_dir="${stage_dir}/remag_coverage"
    "$PYTHON_EXE" "${SCRIPT_DIR}/split_remag_coverage.py" \
        --coverage-matrix "$CONCOCT_COVERAGE" \
        --out-dir "$remag_coverage_dir"
    remag_manifest="${remag_coverage_dir}/manifest.tsv"
    while IFS=$'\t' read -r _index _coverage_name coverage_file _rows _sha256; do
        [[ "$coverage_file" == "file" ]] && continue
        [[ -n "$coverage_file" ]] || continue
        coverage_args+=(-c "${remag_coverage_dir}/${coverage_file}")
    done < "$remag_manifest"
    ((${#coverage_args[@]} > 0)) || die "ReMAG coverage adapter produced no files"
    contig_count="$(awk '/^>/{n += 1} END{print n + 0}' "$INPUT_FASTA")"
    if ((contig_count < 2)); then
        {
            printf 'metric\tvalue\n'
            printf 'state\tskipped_insufficient_contigs\n'
            printf 'contig_count\t%s\n' "$contig_count"
            printf 'minimum_required_contigs\t2\n'
        } > "${stage_dir}/remag_settings.tsv"
        run_in_env "$REMAG_ENV" remag --version \
            > "${stage_dir}/tool_version.txt"
        return 0
    fi
    command=(remag \
        "$INPUT_FASTA" \
        "${coverage_args[@]}" \
        -o "$stage_dir" \
        -m short-reads \
        --skip-bacterial-filter \
        --min-contig-length "$REMAG_MIN_CONTIG" \
        --min-bin-size "$MIN_BIN_SIZE" \
        --random-seed "$RANDOM_SEED" \
        -t "$THREADS")
    if [[ "$REMAG_KEEP_INTERMEDIATE" == "1" ]]; then
        command+=(-k)
    fi
    run_in_env "$REMAG_ENV" "${command[@]}"
    run_in_env "$REMAG_ENV" remag --version > "${stage_dir}/tool_version.txt"
}

run_metabat2() {
    local stage_dir="$1"
    mkdir -p "${stage_dir}/bins"
    run_in_env "$METABAT2_ENV" metabat2 \
        -i "$INPUT_FASTA" \
        -a "$METABAT_DEPTH" \
        -o "${stage_dir}/bins/${SAMPLE}.metabat2" \
        -m "$METABAT2_MIN_CONTIG" \
        -s "$MIN_BIN_SIZE" \
        --seed "$RANDOM_SEED" \
        -t "$THREADS"
    if ! run_in_env "$METABAT2_ENV" metabat2 --help \
        > "${stage_dir}/tool_version.txt" 2>&1; then
        printf 'MetaBAT2 completed, but its help/version probe returned non-zero.\n' \
            >> "${stage_dir}/tool_version.txt"
    fi
}

run_concoct() {
    local stage_dir="$1"
    local fragment_count
    local eligible_fragment_count
    local cluster_count
    local clustering_file
    local coverage_column_count
    local length_threshold
    local branch_state
    local -a concoct_command

    length_threshold=$((CONCOCT_MIN_CONTIG - 1))
    mkdir -p "${stage_dir}/concoct" "${stage_dir}/raw_bins" "${stage_dir}/bins"
    run_in_env "$CONCOCT_ENV" cut_up_fasta.py \
        "$INPUT_FASTA" \
        -c "$CONCOCT_CHUNK_SIZE" \
        -o 0 \
        --merge_last \
        -b "${stage_dir}/contigs_10K.bed" \
        > "${stage_dir}/contigs_10K.fa"

    awk '
        BEGIN { FS = OFS = "\t" }
        NR == FNR {
            if (FNR == 1) {
                header = $0
                column_count = NF
                if (column_count < 2) {
                    print "Coverage table needs at least two columns" > "/dev/stderr"
                    exit 43
                }
                next
            }
            if (NF != column_count) {
                print "Malformed coverage row for " $1 > "/dev/stderr"
                exit 43
            }
            for (column = 2; column <= column_count; column++) {
                coverage[$1, column] = $column
            }
            next
        }
        FNR == 1 {
            print header
        }
        {
            if (!(($1, 2) in coverage)) {
                print "Missing parent coverage for " $1 > "/dev/stderr"
                missing += 1
                next
            }
            printf "%s", $4
            for (column = 2; column <= column_count; column++) {
                printf "%s%s", OFS, coverage[$1, column]
            }
            printf "%s", ORS
        }
        END {
            if (missing > 0) {
                exit 42
            }
        }
    ' "$CONCOCT_COVERAGE" "${stage_dir}/contigs_10K.bed" \
        > "${stage_dir}/fragment_coverage.tsv"

    fragment_count="$(awk '/^>/{n += 1} END{print n + 0}' "${stage_dir}/contigs_10K.fa")"
    eligible_fragment_count="$(
        awk -v minimum="$CONCOCT_MIN_CONTIG" '
            function finish_record() {
                if (seen && sequence_length >= minimum) {
                    eligible += 1
                }
            }
            /^>/ {
                finish_record()
                seen = 1
                sequence_length = 0
                next
            }
            {
                gsub(/[[:space:]]/, "")
                sequence_length += length($0)
            }
            END {
                finish_record()
                print eligible + 0
            }
        ' "${stage_dir}/contigs_10K.fa"
    )"
    cluster_count=0
    branch_state="skipped_insufficient_eligible_fragments"
    if ((eligible_fragment_count >= 2)); then
        cluster_count=$((eligible_fragment_count / 2))
        ((cluster_count < 2)) && cluster_count=2
        ((cluster_count > CONCOCT_MAX_CLUSTERS)) &&
            cluster_count="$CONCOCT_MAX_CLUSTERS"
        branch_state="run"
    fi
    {
        printf 'metric\tvalue\n'
        printf 'state\t%s\n' "$branch_state"
        printf 'fragment_count\t%s\n' "$fragment_count"
        printf 'eligible_fragment_count\t%s\n' "$eligible_fragment_count"
        printf 'min_contig_inclusive\t%s\n' "$CONCOCT_MIN_CONTIG"
        printf 'effective_length_threshold\t%s\n' "$length_threshold"
        printf 'cluster_count\t%s\n' "$cluster_count"
    } > "${stage_dir}/concoct_settings.tsv"
    if ((eligible_fragment_count < 2)); then
        "$PYTHON_EXE" "${SCRIPT_DIR}/select_bins_by_size.py" \
            --input-dir "${stage_dir}/raw_bins" \
            --output-dir "${stage_dir}/bins" \
            --manifest "${stage_dir}/bin_selection.tsv" \
            --min-bp "$MIN_BIN_SIZE"
        run_in_env "$CONCOCT_ENV" concoct --version \
            > "${stage_dir}/tool_version.txt"
        return 0
    fi

    coverage_column_count="$(
        awk -F '\t' 'NR == 1 { print NF - 1; exit }' \
            "${stage_dir}/fragment_coverage.tsv"
    )"
    concoct_command=(concoct \
        --composition_file "${stage_dir}/contigs_10K.fa" \
        --coverage_file "${stage_dir}/fragment_coverage.tsv")
    if ((coverage_column_count == 1)); then
        # Row normalization erases a single dimension.  Keep the historical
        # single-depth behavior, but use CONCOCT defaults for multi-sample
        # differential coverage.
        concoct_command+=(--no_cov_normalization --no_total_coverage)
    fi
    concoct_command+=( \
        --clusters "$cluster_count" \
        --length_threshold "$length_threshold" \
        --threads "$THREADS" \
        --seed "$RANDOM_SEED" \
        -b "${stage_dir}/concoct/")
    run_in_env "$CONCOCT_ENV" "${concoct_command[@]}"

    clustering_file="${stage_dir}/concoct/clustering_gt${length_threshold}.csv"
    [[ -s "$clustering_file" ]] || return 44
    run_in_env "$CONCOCT_ENV" merge_cutup_clustering.py \
        "$clustering_file" \
        > "${stage_dir}/clustering_merged.csv"
    run_in_env "$CONCOCT_ENV" extract_fasta_bins.py \
        "$INPUT_FASTA" \
        "${stage_dir}/clustering_merged.csv" \
        --output_path "${stage_dir}/raw_bins"
    "$PYTHON_EXE" "${SCRIPT_DIR}/select_bins_by_size.py" \
        --input-dir "${stage_dir}/raw_bins" \
        --output-dir "${stage_dir}/bins" \
        --manifest "${stage_dir}/bin_selection.tsv" \
        --min-bp "$MIN_BIN_SIZE"
    run_in_env "$CONCOCT_ENV" concoct --version > "${stage_dir}/tool_version.txt"
}

run_branch() {
    local input_set="$1"
    local binner="$2"
    local parent_dir="${RESULTS_DIR}/${input_set}"
    local final_dir="${parent_dir}/${binner}"
    local stage_dir
    local failed_dir
    local log_file="${LOG_DIR}/${input_set}.${binner}.log"
    local fingerprint_tmp_dir
    local expected_fingerprint
    local rc

    mkdir -p "$parent_dir"
    fingerprint_tmp_dir="$(mktemp -d "${SAMPLE_OUT}/.fingerprint.XXXXXX")"
    expected_fingerprint="${fingerprint_tmp_dir}/fingerprint.json"
    if [[ -f "${final_dir}/.done" ]]; then
        write_branch_fingerprint \
            "$input_set" "$binner" "$expected_fingerprint"
        validate_completed_branch \
            "$input_set" "$binner" "$final_dir" "$expected_fingerprint"
        rm -rf -- "$fingerprint_tmp_dir"
        note "Skip completed branch: ${input_set}/${binner}"
        printf '%s\t%s\tcompleted\t0\t%s\n' "$input_set" "$binner" "$final_dir" >> "$STATUS_TMP"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        return 0
    fi
    if [[ -e "$final_dir" ]]; then
        rm -rf -- "$fingerprint_tmp_dir"
        die "Refusing to overwrite incomplete result directory: $final_dir"
    fi

    stage_dir="$(mktemp -d "${parent_dir}/.${binner}.tmp.XXXXXX")"
    rmdir "$fingerprint_tmp_dir"
    note "Run branch: ${input_set}/${binner}"
    set +e
    (
        set -Eeuo pipefail
        printf 'sample\tinput_set\tbinner\tinput_fasta\tdepth_source\n' > "${stage_dir}/provenance.tsv"
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "$SAMPLE" "$input_set" "$binner" "$INPUT_FASTA" "$DEPTH_TSV" \
            >> "${stage_dir}/provenance.tsv"
        write_branch_fingerprint \
            "$input_set" "$binner" "${stage_dir}/fingerprint.json"
        case "$binner" in
            remag) run_remag "$stage_dir" ;;
            concoct) run_concoct "$stage_dir" ;;
            metabat2) run_metabat2 "$stage_dir" ;;
            *) return 64 ;;
        esac
        (
            cd "$stage_dir"
            "$PYTHON_EXE" "${SCRIPT_DIR}/audit_branch.py" \
                --source-fasta "$INPUT_FASTA" \
                --bins-dir bins \
                --label "${input_set}.${binner}" \
                --out-dir audit
        )
    ) > "$log_file" 2>&1
    rc=$?
    set -e

    if ((rc == 0)); then
        touch "${stage_dir}/.done"
        mv "$stage_dir" "$final_dir"
        printf '%s\t%s\tsuccess\t0\t%s\n' "$input_set" "$binner" "$final_dir" >> "$STATUS_TMP"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        note "Completed branch: ${input_set}/${binner}"
    else
        failed_dir="${FAILED_DIR}/${input_set}.${binner}.$(date -u '+%Y%m%dT%H%M%SZ').$$"
        mv "$stage_dir" "$failed_dir"
        cp "$log_file" "${failed_dir}/branch.log"
        printf '%s\t%s\tfailed\t%s\t%s\n' "$input_set" "$binner" "$rc" "$failed_dir" >> "$STATUS_TMP"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
        note "Branch failed (${rc}); retained diagnostics at ${failed_dir}"
    fi
}

STATUS_TMP="$(mktemp "${SAMPLE_OUT}/.branch_status.XXXXXX")"
printf 'input_set\tbinner\tstatus\texit_code\tresult_path\n' > "$STATUS_TMP"
SUCCESS_COUNT=0
FAILURE_COUNT=0

IFS=' ' read -r -a INPUT_SET_LIST <<< "$INPUT_SETS"
IFS=' ' read -r -a BINNER_LIST <<< "$BINNERS"
for input_set in "${INPUT_SET_LIST[@]}"; do
    select_inputs "$input_set"
    for binner in "${BINNER_LIST[@]}"; do
        run_branch "$input_set" "$binner"
    done
done
mv "$STATUS_TMP" "${SAMPLE_OUT}/branch_status.tsv"

((SUCCESS_COUNT > 0)) || die "All binning branches failed; see ${LOG_DIR} and ${FAILED_DIR}"

summary_args=()
for input_set in "${INPUT_SET_LIST[@]}"; do
    for binner in "${BINNER_LIST[@]}"; do
        bins_dir="${RESULTS_DIR}/${input_set}/${binner}/bins"
        if [[ -d "$bins_dir" ]]; then
            summary_args+=(--bins-label "${input_set}.${binner}=${bins_dir}")
        fi
    done
done

if ((${#summary_args[@]})); then
    "$PYTHON_EXE" "${SCRIPT_DIR}/summarize_bins.py" \
        "${summary_args[@]}" \
        --out-dir "$SUMMARY_DIR"
fi

if ((FAILURE_COUNT > 0)); then
    die "${FAILURE_COUNT} binning branch(es) failed after all requested branches were attempted; successful outputs and their summary were retained. See ${SAMPLE_OUT}/branch_status.tsv."
fi

note "Pilot finished: ${SAMPLE_OUT}"
