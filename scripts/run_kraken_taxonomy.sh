#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

GLOBAL_CONFIG="${1:?Usage: run_kraken_taxonomy.sh GLOBAL_CONFIG SCOPE [MAX_RANK]}"
SCOPE="${2:?Usage: run_kraken_taxonomy.sh GLOBAL_CONFIG SCOPE [MAX_RANK]}"
MAX_RANK="${3:-}"
[[ -s "$GLOBAL_CONFIG" ]] || {
    printf '[ERROR] Config not found: %s\n' "$GLOBAL_CONFIG" >&2
    exit 2
}
# shellcheck source=/dev/null
source "$GLOBAL_CONFIG"

: "${MANIFEST:?Missing MANIFEST}"
: "${PIPELINE_DIR:?Missing PIPELINE_DIR}"
: "${INPUTS_ROOT:?Missing INPUTS_ROOT}"
: "${KRAKEN_ROOT:?Missing KRAKEN_ROOT}"
: "${TAXONOMY_ROOT:?Missing TAXONOMY_ROOT}"
: "${ADDON_STATUS_ROOT:?Missing ADDON_STATUS_ROOT}"
: "${PYTHON_EXE:=python3}"
: "${KRAKEN2_EXE:?Missing KRAKEN2_EXE}"
: "${KRAKEN_DB:?Missing KRAKEN_DB}"
: "${KRAKEN_MIN_HIT_GROUPS:=3}"
: "${KRAKEN_CANARY_USE_MEMORY_MAPPING:=1}"
: "${KRAKEN_FULL_USE_MEMORY_MAPPING:=0}"
: "${THREADS:=${SLURM_CPUS_PER_TASK:-8}}"

[[ "$SCOPE" =~ ^[A-Za-z0-9._-]+$ ]] || {
    printf '[ERROR] Unsafe taxonomy scope: %s\n' "$SCOPE" >&2
    exit 2
}
if [[ -n "$MAX_RANK" && ! "$MAX_RANK" =~ ^[1-9][0-9]*$ ]]; then
    printf '[ERROR] MAX_RANK must be empty or a positive integer: %s\n' \
        "$MAX_RANK" >&2
    exit 2
fi
if [[ "$SCOPE" == canary* ]]; then
    kraken_use_memory_mapping="$KRAKEN_CANARY_USE_MEMORY_MAPPING"
else
    kraken_use_memory_mapping="$KRAKEN_FULL_USE_MEMORY_MAPPING"
fi
[[ "$kraken_use_memory_mapping" == "0" ||
   "$kraken_use_memory_mapping" == "1" ]] || {
    printf '[ERROR] Effective Kraken memory-mapping mode must be 0 or 1: %s\n' \
        "$kraken_use_memory_mapping" >&2
    exit 2
}
for name in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "${KRAKEN_DB}/${name}" ]] || {
        printf '[ERROR] Incomplete Kraken database: %s\n' "${KRAKEN_DB}/${name}" >&2
        exit 2
    }
done

catalogue_dir="${KRAKEN_ROOT}/${SCOPE}.catalogue"
taxonomy_dir="${TAXONOMY_ROOT}/${SCOPE}"
final_run_dir="${KRAKEN_ROOT}/${SCOPE}.run"
mkdir -p "$KRAKEN_ROOT" "$TAXONOMY_ROOT" \
    "${ADDON_STATUS_ROOT}/kraken"

if [[ -s "${taxonomy_dir}/taxonomy.summary.json" &&
      -s "${final_run_dir}/done.tsv" ]]; then
    printf 'Kraken taxonomy already complete: %s\n' "$taxonomy_dir"
    exit 0
fi
[[ ! -e "$taxonomy_dir" ]] || {
    printf '[ERROR] Partial taxonomy output exists: %s\n' "$taxonomy_dir" >&2
    exit 2
}
[[ ! -e "$final_run_dir" ]] || {
    printf '[ERROR] Partial Kraken run output exists: %s\n' "$final_run_dir" >&2
    exit 2
}

catalogue_args=(
    --manifest "$MANIFEST"
    --inputs-root "$INPUTS_ROOT"
    --out-dir "$catalogue_dir"
)
if [[ -n "$MAX_RANK" ]]; then
    catalogue_args+=(--max-rank "$MAX_RANK")
fi
"$PYTHON_EXE" "${PIPELINE_DIR}/scripts/build_kraken_catalogue.py" \
    "${catalogue_args[@]}"
chmod -R g+rwX,o-rwx "$catalogue_dir"

job_token="${SLURM_JOB_ID:-manual}"
stage="$(
    mktemp -d "${KRAKEN_ROOT}/.${SCOPE}.run.${job_token}.XXXXXX"
)"
cleanup_stage() {
    rm -rf -- "$stage"
}
trap cleanup_stage EXIT

kraken_args=(
    --db "$KRAKEN_DB"
    --threads "$THREADS"
    --minimum-hit-groups "$KRAKEN_MIN_HIT_GROUPS"
    --report "${stage}/kraken.report"
    --output "${stage}/kraken.output"
)
if [[ "$kraken_use_memory_mapping" == "1" ]]; then
    kraken_args+=(--memory-mapping)
fi
kraken_args+=("${catalogue_dir}/catalogue.fa")

"$KRAKEN2_EXE" "${kraken_args[@]}"
"$PYTHON_EXE" "${PIPELINE_DIR}/scripts/kraken_to_taxvamb.py" \
    --kraken-output "${stage}/kraken.output" \
    --kraken-report "${stage}/kraken.report" \
    --samples-tsv "${catalogue_dir}/samples.tsv" \
    --out-dir "$taxonomy_dir"

mkdir "${stage}/final"
mv "${stage}/kraken.report" "${stage}/final/kraken.report"
gzip -1 -c "${stage}/kraken.output" \
    > "${stage}/final/kraken.output.gz"
{
    printf 'field\tvalue\n'
    printf 'scope\t%s\n' "$SCOPE"
    printf 'max_rank\t%s\n' "${MAX_RANK:-all}"
    printf 'database\t%s\n' "$KRAKEN_DB"
    printf 'minimum_hit_groups\t%s\n' "$KRAKEN_MIN_HIT_GROUPS"
    printf 'memory_mapping\t%s\n' "$kraken_use_memory_mapping"
    printf 'threads\t%s\n' "$THREADS"
    printf 'kraken_version\t'
    "$KRAKEN2_EXE" --version |
        tr '\n' ' '
    printf '\n'
    for name in hash.k2d opts.k2d taxo.k2d seqid2taxid.map; do
        if [[ -f "${KRAKEN_DB}/${name}" ]]; then
            db_metadata="$(
                "$PYTHON_EXE" -c \
                    'import os, sys; s=os.stat(sys.argv[1]); print(f"{s.st_size}\t{s.st_mtime_ns}")' \
                    "${KRAKEN_DB}/${name}"
            )"
            IFS=$'\t' read -r db_bytes db_mtime_ns <<< "$db_metadata"
            printf '%s_bytes\t%s\n' "$name" "$db_bytes"
            printf '%s_mtime_ns\t%s\n' "$name" "$db_mtime_ns"
        fi
    done
} > "${stage}/final/provenance.tsv"
printf '%s\t%s\t%s\n' \
    "$SCOPE" "${MAX_RANK:-all}" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    > "${stage}/final/done.tsv"
mv "${stage}/final" "$final_run_dir"
chmod -R g+rwX,o-rwx "$taxonomy_dir" "$final_run_dir"

status_path="${ADDON_STATUS_ROOT}/kraken/${SCOPE}.done"
printf '%s\t%s\t%s\n' \
    "$SCOPE" "${MAX_RANK:-all}" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    > "$status_path"
chmod g+rw,o-rwx "$status_path"
