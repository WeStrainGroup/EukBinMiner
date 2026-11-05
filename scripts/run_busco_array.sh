#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
GLOBAL_CONFIG="${1:?Usage: run_busco_array.sh GLOBAL_CONFIG [CANDIDATES] [OUT_DIR]}"
candidate_override="${2:-}"
output_override="${3:-}"
# shellcheck source=/dev/null
source "$GLOBAL_CONFIG"

: "${RUN_ROOT:?Missing RUN_ROOT}"
: "${EUKCC_CANDIDATES:=${RUN_ROOT}/downstream/eukcc_candidates.tsv}"
: "${BUSCO_CANDIDATES:=$EUKCC_CANDIDATES}"
: "${BUSCO_OUT:=${RUN_ROOT}/downstream/busco}"
[[ -z "$candidate_override" ]] || BUSCO_CANDIDATES="$candidate_override"
[[ -z "$output_override" ]] || BUSCO_OUT="$output_override"
: "${BUSCO_ENV:=}"
: "${BUSCO_DB:?Missing BUSCO_DB}"
: "${CONDA_EXE:=conda}"
: "${THREADS:=${SLURM_CPUS_PER_TASK:-4}}"
: "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"

row="$(
    awk -F '\t' -v rank="$SLURM_ARRAY_TASK_ID" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                column[$i] = i
            }
            next
        }
        NR == rank + 1 {
            print $column["bin_id"] "\t" $column["bin_path"]
            found = 1
            exit
        }
        END { if (!found) exit 42 }
    ' "$BUSCO_CANDIDATES"
)" || {
    printf '[ERROR] No BUSCO candidate row %s\n' "$SLURM_ARRAY_TASK_ID" >&2
    exit 2
}
IFS=$'\t' read -r bin_id bin_path <<< "$row"
[[ -n "$bin_id" && -s "$bin_path" ]] || {
    printf '[ERROR] Invalid BUSCO candidate row %s\n' "$SLURM_ARRAY_TASK_ID" >&2
    exit 2
}

mkdir -p "$BUSCO_OUT"
final_dir="${BUSCO_OUT}/${bin_id}"
if [[ -f "${final_dir}/.done" ]]; then
    find "$final_dir" -type f -name 'short_summary*.txt' -size +0c -print -quit |
        grep -q . || {
        printf '[ERROR] Completed BUSCO result lacks summary: %s\n' "$final_dir" >&2
        exit 2
    }
    exit 0
fi
[[ ! -e "$final_dir" ]] || {
    printf '[ERROR] Refusing incomplete BUSCO directory: %s\n' "$final_dir" >&2
    exit 2
}

stage="$(mktemp -d "${BUSCO_OUT}/.${bin_id}.tmp.XXXXXX")"
set +e
command=(busco)
if [[ -n "$BUSCO_ENV" ]]; then
    selector=-n
    [[ "$BUSCO_ENV" == /* ]] && selector=-p
    command=(
        "$CONDA_EXE" run --no-capture-output "$selector" "$BUSCO_ENV"
        "${command[@]}"
    )
fi
"${command[@]}" \
    --in "$bin_path" \
    --out "$bin_id" \
    --out_path "$stage" \
    --mode genome \
    --lineage_dataset "$BUSCO_DB" \
    --cpu "$THREADS" \
    --offline \
    --tar
rc=$?
set -e
if ((rc != 0)); then
    failed="${BUSCO_OUT}/failed.${bin_id}.$(date -u '+%Y%m%dT%H%M%SZ').$$"
    mv "$stage" "$failed"
    printf '[ERROR] BUSCO failed for %s (exit %s); retained %s\n' \
        "$bin_id" "$rc" "$failed" >&2
    exit "$rc"
fi
find "$stage" -type f -name 'short_summary*.txt' -size +0c -print -quit |
    grep -q . || {
    failed="${BUSCO_OUT}/failed.${bin_id}.$(date -u '+%Y%m%dT%H%M%SZ').$$"
    mv "$stage" "$failed"
    printf '[ERROR] BUSCO produced no summary for %s; retained %s\n' \
        "$bin_id" "$failed" >&2
    exit 3
}
touch "${stage}/.done"
mv "$stage" "$final_dir"
