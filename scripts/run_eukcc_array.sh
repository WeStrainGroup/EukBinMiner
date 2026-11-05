#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
GLOBAL_CONFIG="${1:?Usage: run_eukcc_array.sh GLOBAL_CONFIG [CANDIDATES] [OUT_DIR]}"
candidate_override="${2:-}"
output_override="${3:-}"
# shellcheck source=/dev/null
source "$GLOBAL_CONFIG"

: "${RUN_ROOT:?Missing RUN_ROOT}"
: "${BIN_MANIFEST:=${RUN_ROOT}/downstream/bin_manifest.tsv}"
: "${EUKCC_CANDIDATES:=${RUN_ROOT}/downstream/eukcc_candidates.tsv}"
: "${EUKCC_OUT:=${RUN_ROOT}/downstream/eukcc}"
[[ -z "$candidate_override" ]] || EUKCC_CANDIDATES="$candidate_override"
[[ -z "$output_override" ]] || EUKCC_OUT="$output_override"
: "${EUKCC_DB:?Missing EUKCC_DB}"
: "${EUKCC_ENV:=}"
: "${CONDA_EXE:=conda}"
: "${THREADS:=${SLURM_CPUS_PER_TASK:-4}}"
: "${EUKCC_BINS_PER_TASK:=1}"
: "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"

[[ "$EUKCC_BINS_PER_TASK" =~ ^[1-9][0-9]*$ ]] || {
    printf '[ERROR] EUKCC_BINS_PER_TASK must be a positive integer: %s\n' \
        "$EUKCC_BINS_PER_TASK" >&2
    exit 2
}
candidate_count="$(awk 'END { print NR - 1 }' "$EUKCC_CANDIDATES")"
first_rank=$(((SLURM_ARRAY_TASK_ID - 1) * EUKCC_BINS_PER_TASK + 1))
last_rank=$((first_rank + EUKCC_BINS_PER_TASK - 1))
((last_rank > candidate_count)) && last_rank="$candidate_count"
((first_rank <= candidate_count)) || exit 0
mkdir -p "$EUKCC_OUT"

process_candidate() {
    local candidate_rank="$1"
    local row
    local bin_rank bin_id study sample sample_rank binner original_bin_id
    local contigs bp n50 longest sha256 source_path bin_path materialization
    local final_dir stage rc failed

    row="$(
        awk -F '\t' -v rank="$candidate_rank" '
            NR == rank + 1 {
                print
                found = 1
                exit
            }
            END { if (!found) exit 42 }
        ' "$EUKCC_CANDIDATES"
    )" || {
        printf '[ERROR] No EukCC candidate row %s\n' \
            "$candidate_rank" >&2
        return 2
    }
    IFS=$'\t' read -r \
        bin_rank bin_id study sample sample_rank binner original_bin_id \
        contigs bp n50 longest sha256 source_path bin_path materialization \
        <<< "$row"
    [[ "$bin_rank" == "$candidate_rank" && -s "$bin_path" ]] || {
        printf '[ERROR] Invalid EukCC candidate row %s\n' \
            "$candidate_rank" >&2
        return 2
    }

    final_dir="${EUKCC_OUT}/${bin_id}"
    if [[ -f "${final_dir}/.done" ]]; then
        if find "$final_dir" -type f -name '*.csv' -size +0c -print -quit |
            grep -q .; then
            return 0
        fi
        if [[ -s "${final_dir}/status.tsv" ]] &&
            grep -Eq \
                $'^(no_placement_markers|no_suitable_marker_set|placement_alignment_failed|metaeuk_empty_proteins)\t' \
                "${final_dir}/status.tsv"; then
            return 0
        fi
        printf '[ERROR] Completed EukCC result lacks quality CSV or valid status: %s\n' \
            "$final_dir" >&2
        return 2
    fi
    [[ ! -e "$final_dir" ]] || {
        printf '[ERROR] Refusing incomplete EukCC directory: %s\n' \
            "$final_dir" >&2
        return 2
    }

    stage="$(mktemp -d "${EUKCC_OUT}/.${bin_id}.tmp.XXXXXX")"
    set +e
    command=(eukcc single)
    if [[ -n "$EUKCC_ENV" ]]; then
        selector=-n
        [[ "$EUKCC_ENV" == /* ]] && selector=-p
        command=(
            "$CONDA_EXE" run --no-capture-output "$selector" "$EUKCC_ENV"
            "${command[@]}"
        )
    fi
    "${command[@]}" \
        --db "$EUKCC_DB" \
        --out "$stage" \
        --threads "$THREADS" \
        "$bin_path"
    rc=$?
    set -e
    if ((rc != 0)); then
        if grep -FRq 'No placement marker genes found' "$stage"; then
            printf 'eukcc_status\texit_code\tmessage\n' > "${stage}/status.tsv"
            printf 'no_placement_markers\t%s\t%s\n' \
                "$rc" \
                'EukCC found no placement marker genes' >> "${stage}/status.tsv"
            printf '%s\t%s\t%s\n' \
                "$bin_id" "$sha256" "$bin_path" > "${stage}/input.tsv"
            touch "${stage}/.done"
            mv "$stage" "$final_dir"
            printf '[WARN] EukCC found no placement markers for %s; recorded valid status\n' \
                "$bin_id" >&2
            return 0
        fi
        if grep -FRq \
            -e 'No marker gene set could be found' \
            -e 'Could not identify a single suitable marker gene set' \
            "$stage"; then
            printf 'eukcc_status\texit_code\tmessage\n' > "${stage}/status.tsv"
            printf 'no_suitable_marker_set\t%s\t%s\n' \
                "$rc" \
                'EukCC could not identify a suitable marker gene set' \
                >> "${stage}/status.tsv"
            printf '%s\t%s\t%s\n' \
                "$bin_id" "$sha256" "$bin_path" > "${stage}/input.tsv"
            touch "${stage}/.done"
            mv "$stage" "$final_dir"
            printf '[WARN] EukCC found no suitable marker set for %s; recorded valid status\n' \
                "$bin_id" >&2
            return 0
        fi
        # EukCC 2.1.2 can deterministically finish MetaEuk without any
        # predicted proteins.  Such a bin has no defensible EukCC quality
        # estimate; preserve it as explicitly unassessed instead of inventing
        # completeness/contamination values or failing the whole sample.
        # Normalize only the exact log plus zero-length protein/codon outputs,
        # with no quality CSV, so unrelated MetaEuk failures remain fatal.
        if grep -FRq 'Metaeuk produced empty fasta.' "$stage" &&
            find "$stage" -type f -name '*_metaeuk_cleaned.faa' -size 0c \
                -print -quit | grep -q . &&
            find "$stage" -type f -name '*_metaeuk.codon.fas' -size 0c \
                -print -quit | grep -q . &&
            ! find "$stage" -type f -name '*.csv' -size +0c -print -quit |
                grep -q .; then
            printf 'eukcc_status\texit_code\tmessage\n' > "${stage}/status.tsv"
            printf 'metaeuk_empty_proteins\t%s\t%s\n' \
                "$rc" \
                'EukCC MetaEuk produced no proteins; quality is unassessable' \
                >> "${stage}/status.tsv"
            printf '%s\t%s\t%s\n' \
                "$bin_id" "$sha256" "$bin_path" > "${stage}/input.tsv"
            touch "${stage}/.done"
            mv "$stage" "$final_dir"
            printf '[WARN] EukCC MetaEuk produced no proteins for %s; recorded valid unassessed status\n' \
                "$bin_id" >&2
            return 0
        fi
        # EukCC 2.1.2 can exit during hmmalign when exactly one base
        # placement marker is found.  This is deterministic for the affected
        # bin and produces neither an EPA query nor a quality estimate.  Only
        # normalize that narrow filesystem/log signature; later EPA failures,
        # OOMs, and every other non-zero exit remain fatal.
        if grep -FRq \
            'Found 1 marker genes, placing them in the tree using epa-ng' \
            "$stage" &&
            find "$stage" -type f -name 'hmmalign.stdout' -size 0c \
                -print -quit | grep -q . &&
            ! find "$stage" -type f -name 'query.fasta' -print -quit |
                grep -q . &&
            ! find "$stage" -type f -name 'epa-ng.stout' -print -quit |
                grep -q .; then
            printf 'eukcc_status\texit_code\tmessage\n' > "${stage}/status.tsv"
            printf 'placement_alignment_failed\t%s\t%s\n' \
                "$rc" \
                'EukCC single-marker placement alignment failed; quality is unassessable' \
                >> "${stage}/status.tsv"
            printf '%s\t%s\t%s\n' \
                "$bin_id" "$sha256" "$bin_path" > "${stage}/input.tsv"
            touch "${stage}/.done"
            mv "$stage" "$final_dir"
            printf '[WARN] EukCC single-marker placement alignment failed for %s; recorded valid unassessed status\n' \
                "$bin_id" >&2
            return 0
        fi
        failed="${EUKCC_OUT}/failed.${bin_id}.$(date -u '+%Y%m%dT%H%M%SZ').$$"
        mv "$stage" "$failed"
        printf '[ERROR] EukCC failed for %s (exit %s); retained %s\n' \
            "$bin_id" "$rc" "$failed" >&2
        return "$rc"
    fi
    if ! find "$stage" -type f -name '*.csv' -size +0c -print -quit |
        grep -q .; then
        failed="${EUKCC_OUT}/failed.${bin_id}.$(date -u '+%Y%m%dT%H%M%SZ').$$"
        mv "$stage" "$failed"
        printf '[ERROR] EukCC produced no non-empty CSV for %s; retained %s\n' \
            "$bin_id" "$failed" >&2
        return 3
    fi
    printf '%s\t%s\t%s\n' \
        "$bin_id" "$sha256" "$bin_path" > "${stage}/input.tsv"
    touch "${stage}/.done"
    mv "$stage" "$final_dir"
}

failure_count=0
for ((candidate_rank = first_rank; candidate_rank <= last_rank; candidate_rank++)); do
    set +e
    (set -Eeuo pipefail; process_candidate "$candidate_rank")
    candidate_rc="$?"
    set -e
    if ((candidate_rc != 0)); then
        printf '[ERROR] Packed EukCC task %s failed candidate rank %s\n' \
            "$SLURM_ARRAY_TASK_ID" "$candidate_rank" >&2
        failure_count=$((failure_count + 1))
    fi
done
((failure_count == 0)) || {
    printf '[ERROR] Packed EukCC task had %s failed candidate(s)\n' \
        "$failure_count" >&2
    exit 1
}
