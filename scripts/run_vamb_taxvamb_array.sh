#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'
umask 0007

GLOBAL_CONFIG="${1:?Usage: run_vamb_taxvamb_array.sh GLOBAL_CONFIG TAXONOMY_SCOPE}"
REQUESTED_SCOPE="${2:?Usage: run_vamb_taxvamb_array.sh GLOBAL_CONFIG TAXONOMY_SCOPE}"
[[ -s "$GLOBAL_CONFIG" ]] || {
    printf '[ERROR] Config not found: %s\n' "$GLOBAL_CONFIG" >&2
    exit 2
}
# shellcheck source=/dev/null
source "$GLOBAL_CONFIG"

: "${MANIFEST:?Missing MANIFEST}"
: "${PIPELINE_DIR:?Missing PIPELINE_DIR}"
: "${INPUTS_ROOT:?Missing INPUTS_ROOT}"
: "${TAXONOMY_ROOT:?Missing TAXONOMY_ROOT}"
: "${OUT_ROOT:?Missing OUT_ROOT}"
: "${ADDON_ROOT:?Missing ADDON_ROOT}"
: "${ADDON_STATUS_ROOT:?Missing ADDON_STATUS_ROOT}"
: "${PYTHON_EXE:=python3}"
: "${VAMB_EXE:?Missing VAMB_EXE}"
: "${VAMB_TOOL_IDENTITY:?Missing VAMB_TOOL_IDENTITY}"
: "${VAMB_ENV_HISTORY:?Missing VAMB_ENV_HISTORY}"
: "${VAMB_MIN_CONTIG_INCLUSIVE:=1500}"
: "${VAMB_MIN_FASTA_BP:=1}"
: "${RANDOM_SEED:=6311}"
: "${VAMB_TAXVAMB_BINNERS:=vamb taxvamb}"
: "${THREADS:=${SLURM_CPUS_PER_TASK:-8}}"
: "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"

[[ "$RANDOM_SEED" =~ ^[1-9][0-9]*$ ]] || {
    printf '[ERROR] RANDOM_SEED must be a positive integer: %s\n' \
        "$RANDOM_SEED" >&2
    exit 2
}

TAXONOMY_SCOPE="$REQUESTED_SCOPE"
[[ "$TAXONOMY_SCOPE" =~ ^[A-Za-z0-9._-]+$ ]] || {
    printf '[ERROR] Unsafe taxonomy scope: %s\n' "$TAXONOMY_SCOPE" >&2
    exit 2
}

rank="$SLURM_ARRAY_TASK_ID"
row="$(
    awk -F '\t' -v rank="$rank" '
        NR > 1 && $1 == rank {
            print
            found = 1
            exit
        }
        END {
            if (!found) {
                exit 42
            }
        }
    ' "$MANIFEST"
)" || {
    printf '[ERROR] No manifest row for rank %s\n' "$rank" >&2
    exit 2
}
IFS=$'\t' read -r observed_rank study sample _euk _unknown _depth \
    _selection_hash <<< "$row"
[[ "$observed_rank" == "$rank" ]] || {
    printf '[ERROR] Manifest rank mismatch: %s != %s\n' \
        "$observed_rank" "$rank" >&2
    exit 2
}

rank_name="$(printf 'r%06d' "$rank")"
input_dir="${INPUTS_ROOT}/${rank_name}"
input_fasta="${input_dir}/input.fa"
abundance="${input_dir}/abundance.tsv"
input_metadata="${input_dir}/input.json"
taxonomy="${TAXONOMY_ROOT}/${TAXONOMY_SCOPE}/per_sample/${rank_name}.taxonomy.tsv"
required_paths=("$input_fasta" "$abundance" "$input_metadata")
if [[ " ${VAMB_TAXVAMB_BINNERS} " == *" taxvamb "* ]]; then
    required_paths+=("$taxonomy")
fi
for path in "${required_paths[@]}"; do
    [[ -f "$path" ]] || {
        printf '[ERROR] Required input not found: %s\n' "$path" >&2
        exit 2
    }
done

sample_out="${OUT_ROOT}/${study}/${sample}/results/union"
status_dir="${ADDON_STATUS_ROOT}/binner"
staging_root="${ADDON_ROOT}/staging"
mkdir -p "$sample_out" "$status_dir" "$staging_root"

tool_identity="$VAMB_TOOL_IDENTITY"
contig_count="$(
    "$PYTHON_EXE" -c \
        'import json,sys; print(json.load(open(sys.argv[1]))["contigs"])' \
        "$input_metadata"
)"

tighten_owned_permissions() {
    local path="$1"
    local current_user
    current_user="$(id -un)"

    [[ -e "$path" ]] || return 0
    if [[ -d "$path" ]]; then
        if ! find "$path" -xdev -user "$current_user" -type d \
                -exec chmod g+rwX,o-rwx {} + ||
           ! find "$path" -xdev -user "$current_user" -type f \
                -exec chmod g+rw,o-rwx {} +; then
            printf '[WARN] Permission fix-up was incomplete for owned paths under %s\n' \
                "$path" >&2
        fi
        return 0
    fi
    if [[ -O "$path" ]]; then
        if ! chmod g+rw,o-rwx "$path"; then
            printf '[WARN] Permission fix-up failed for owned path %s\n' \
                "$path" >&2
        fi
    else
        printf '[WARN] Permission fix-up skipped for non-owned shared path %s\n' \
            "$path" >&2
    fi
}

write_fingerprint() {
    local binner="$1"
    local output="$2"
    local -a args=(
        --field "study=${study}"
        --field "sample=${sample}"
        --field "input_set=union"
        --field "binner=${binner}"
        --field "tool_identity=${tool_identity}"
        --field "threads=${THREADS}"
        --field "minimum_contig_inclusive=${VAMB_MIN_CONTIG_INCLUSIVE}"
        --field "minimum_fasta_bp=${VAMB_MIN_FASTA_BP}"
        --field "random_seed=${RANDOM_SEED}"
        --field "rng_seed_policy=python_numpy_torch_before_vamb_import"
        --field "binsplitting=disabled"
        --field "default_training_epochs=true"
        --file "input_fasta=${input_fasta}"
        --file "abundance=${abundance}"
        --file "input_metadata=${input_metadata}"
        --file "runner=${PIPELINE_DIR}/scripts/run_vamb_taxvamb_array.sh"
        --file "auditor=${PIPELINE_DIR}/scripts/audit_branch.py"
        --file "environment_spec=${PIPELINE_DIR}/envs/vamb_taxvamb.yml"
        --file "environment_history=${VAMB_ENV_HISTORY}"
        --file "vamb_executable=${VAMB_EXE}"
        --file "vamb_seed_launcher=${PIPELINE_DIR}/scripts/run_seeded_vamb.py"
    )
    if [[ "$binner" == "taxvamb" ]]; then
        args+=(
            --field "taxometer_refinement=enabled"
            --file "taxonomy=${taxonomy}"
        )
    fi
    args+=(--output "$output")
    "$PYTHON_EXE" "${PIPELINE_DIR}/scripts/branch_fingerprint.py" \
        "${args[@]}"
}

validate_existing() {
    local final_dir="$1"
    local expected_fingerprint="$2"
    [[ -s "${final_dir}/fingerprint.json" &&
       -d "${final_dir}/bins" &&
       -s "${final_dir}/audit/audit.json" ]] || return 1
    cmp -s "$expected_fingerprint" "${final_dir}/fingerprint.json" || return 1
    PYTHONPATH="${PIPELINE_DIR}/scripts${PYTHONPATH:+:${PYTHONPATH}}" "$PYTHON_EXE" - "${final_dir}" "$input_fasta" <<'VERIFY'
import json, sys
from pathlib import Path
from audit_branch import audit
import tempfile
final, source = map(Path, sys.argv[1:])
with tempfile.TemporaryDirectory(dir=final.parent, prefix=".reuse-audit.") as work:
    target = Path(work)
    if source.stat().st_size:
        audit(source, final / "bins", "reuse", target)
        actual = json.loads((target / "audit.json").read_text())
        expected = json.loads((final / "audit/audit.json").read_text())
        for key in ("source_sha256", "bins", "assigned_contigs", "assigned_bp"):
            if actual[key] != expected[key]:
                raise SystemExit("VAMB branch audit changed")
VERIFY
}

write_zero_audit() {
    local audit_dir="$1"
    local label="$2"
    mkdir -p "$audit_dir"
    printf 'contig_id\tlength_bp\tstatus\tbin_id\tbin_file\n' \
        > "${audit_dir}/assignments.tsv"
    : > "${audit_dir}/unbinned.list"
    : > "${audit_dir}/unbinned.fa"
    "$PYTHON_EXE" -c '
import hashlib, json, pathlib, sys
source, label, out = map(pathlib.Path, sys.argv[1:4])
payload = {
    "label": str(label),
    "source_fasta": str(source),
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "bins_dir": str(out.parent / "bins"),
    "bin_files": 0,
    "nonempty_bins": 0,
    "source_contigs": 0,
    "source_bp": 0,
    "assigned_contigs": 0,
    "assigned_bp": 0,
    "assigned_contig_fraction": 0.0,
    "assigned_bp_fraction": 0.0,
    "unbinned_contigs": 0,
    "unbinned_bp": 0,
    "bins": [],
}
(out / "audit.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
' "$input_fasta" "$label" "$audit_dir"
}

run_branch() {
    local binner="$1"
    local final_dir="${sample_out}/${binner}"
    local expected_fingerprint
    local stage
    local tool_out
    local failed_log
    local failed_dir

    expected_fingerprint="$(
        mktemp "${staging_root}/.${rank_name}.${binner}.fingerprint.XXXXXX"
    )"
    write_fingerprint "$binner" "$expected_fingerprint"
    if [[ -e "$final_dir" ]]; then
        if validate_existing "$final_dir" "$expected_fingerprint"; then
            rm -f -- "$expected_fingerprint"
            printf '%s already complete: %s\n' "$binner" "$final_dir"
            return 0
        fi
        rm -f -- "$expected_fingerprint"
        printf '[ERROR] Existing branch is stale or partial: %s\n' \
            "$final_dir" >&2
        return 2
    fi

    job_token="${SLURM_JOB_ID:-manual}"
    stage="$(
        mktemp -d \
            "${staging_root}/.${rank_name}.${binner}.${job_token}.XXXXXX"
    )"
    tool_out="${stage}/tool"
    mkdir -p "${stage}/branch/bins"

    if ((contig_count < 2)); then
        {
            printf 'field\tvalue\n'
            printf 'state\tskipped_insufficient_contigs\n'
            printf 'contigs\t%s\n' "$contig_count"
            printf 'minimum_required\t2\n'
        } > "${stage}/branch/run_state.tsv"
        if ((contig_count == 0)); then
            write_zero_audit \
                "${stage}/branch/audit" "union.${binner}"
        else
            "$PYTHON_EXE" "${PIPELINE_DIR}/scripts/audit_branch.py" \
                --source-fasta "$input_fasta" \
                --bins-dir "${stage}/branch/bins" \
                --label "union.${binner}" \
                --out-dir "${stage}/branch/audit"
        fi
    else
        command=(
            "$PYTHON_EXE"
            "${PIPELINE_DIR}/scripts/run_seeded_vamb.py"
            --seed
            "$RANDOM_SEED"
            --
            "$VAMB_EXE"
            bin
        )
        if [[ "$binner" == "vamb" ]]; then
            command+=(default)
        else
            command+=(taxvamb --taxonomy "$taxonomy")
        fi
        command+=(
            --outdir "$tool_out"
            --fasta "$input_fasta"
            --abundance_tsv "$abundance"
            -m "$VAMB_MIN_CONTIG_INCLUSIVE"
            -p "$THREADS"
            --seed "$RANDOM_SEED"
            --minfasta "$VAMB_MIN_FASTA_BP"
            -o
        )
        if ! "${command[@]}"; then
            failed_log="${status_dir}/${rank_name}.${binner}.failed.log"
            if [[ -f "${tool_out}/log.txt" ]]; then
                cp "${tool_out}/log.txt" "$failed_log"
            else
                printf 'VAMB failed before log.txt was created\n' > "$failed_log"
            fi
            failed_dir="${ADDON_ROOT}/failed"
            mkdir -p "$failed_dir"
            mv "$stage" "${failed_dir}/$(basename -- "$stage")"
            rm -f -- "$expected_fingerprint"
            printf '[ERROR] %s failed for %s/%s\n' \
                "$binner" "$study" "$sample" >&2
            return 1
        fi
        [[ -d "${tool_out}/bins" ]] || {
            printf '[ERROR] %s completed without bins directory\n' "$binner" >&2
            failed_dir="${ADDON_ROOT}/failed"
            mkdir -p "$failed_dir"
            mv "$stage" "${failed_dir}/$(basename -- "$stage")"
            rm -f -- "$expected_fingerprint"
            return 1
        }
        rmdir "${stage}/branch/bins"
        mv "${tool_out}/bins" "${stage}/branch/bins"
        while IFS= read -r -d '' path; do
            mv "$path" "${stage}/branch/$(basename -- "$path")"
        done < <(
            find "$tool_out" -maxdepth 1 -type f \
                \( -name '*.tsv' -o -name 'log.txt' \) -print0
        )
        {
            printf 'field\tvalue\n'
            printf 'state\tcompleted\n'
            printf 'contigs\t%s\n' "$contig_count"
        } > "${stage}/branch/run_state.tsv"
        "$PYTHON_EXE" "${PIPELINE_DIR}/scripts/audit_branch.py" \
            --source-fasta "$input_fasta" \
            --bins-dir "${stage}/branch/bins" \
            --label "union.${binner}" \
            --out-dir "${stage}/branch/audit"
    fi

    mv "$expected_fingerprint" "${stage}/branch/fingerprint.json"
    printf '%s\n' "$tool_identity" > "${stage}/branch/tool_version.txt"
    mv "${stage}/branch" "$final_dir"
    # Cross-account waves share this tree through default named-user ACLs.
    # Never recurse with chmod over files created by the other account: POSIX
    # ownership would turn an otherwise valid scientific result into a failed
    # Slurm task. Permission normalization is operational best effort only.
    tighten_owned_permissions "$final_dir"
    rm -rf -- "$stage"
    printf '%s completed: %s\n' "$binner" "$final_dir"
}

vamb_status="not_requested"
taxvamb_status="not_requested"
overall_status=0
IFS=' ' read -r -a requested_binners <<< "$VAMB_TAXVAMB_BINNERS"
((${#requested_binners[@]} > 0)) || {
    printf '[ERROR] VAMB_TAXVAMB_BINNERS is empty\n' >&2
    exit 2
}
for requested_binner in "${requested_binners[@]}"; do
    case "$requested_binner" in
        vamb)
            set +e
            (set -Eeuo pipefail; run_branch vamb)
            vamb_status="$?"
            set -e
            if ((vamb_status != 0)); then
                overall_status=1
            fi
            ;;
        taxvamb)
            set +e
            (set -Eeuo pipefail; run_branch taxvamb)
            taxvamb_status="$?"
            set -e
            if ((taxvamb_status != 0)); then
                overall_status=1
            fi
            ;;
        *)
            printf '[ERROR] Unsupported VAMB_TAXVAMB_BINNERS entry: %s\n' \
                "$requested_binner" >&2
            exit 2
            ;;
    esac
done

status_path="${status_dir}/${rank_name}.status.tsv"
status_stage="$(
    mktemp "${status_dir}/.${rank_name}.status.${SLURM_JOB_ID:-manual}.XXXXXX"
)"
{
    printf 'field\tvalue\n'
    printf 'rank\t%s\n' "$rank"
    printf 'study\t%s\n' "$study"
    printf 'sample\t%s\n' "$sample"
    printf 'taxonomy_scope\t%s\n' "$TAXONOMY_SCOPE"
    printf 'contigs\t%s\n' "$contig_count"
    printf 'vamb_status\t%s\n' "$vamb_status"
    printf 'taxvamb_status\t%s\n' "$taxvamb_status"
    printf 'completed_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
} > "$status_stage"
tighten_owned_permissions "$status_stage"
mv -f -- "$status_stage" "$status_path"

if ((overall_status == 0)); then
    done_path="${status_dir}/${rank_name}.done"
    done_stage="$(
        mktemp "${status_dir}/.${rank_name}.done.${SLURM_JOB_ID:-manual}.XXXXXX"
    )"
    cp "$status_path" "$done_stage"
    tighten_owned_permissions "$done_stage"
    mv -f -- "$done_stage" "$done_path"
    exit 0
fi
exit 1
