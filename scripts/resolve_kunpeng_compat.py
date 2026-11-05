#!/usr/bin/env python3
"""Restore Kraken-compatible k-mer weighting to Kun-peng sequence calls."""

from __future__ import annotations

import argparse
import gzip
import json
import struct
from collections import Counter, deque
from pathlib import Path
from typing import IO, Iterable, Optional, Sequence


DNA_CODES = {
    ord("A"): 0,
    ord("C"): 1,
    ord("G"): 2,
    ord("T"): 3,
    ord("a"): 0,
    ord("c"): 1,
    ord("g"): 2,
    ord("t"): 3,
}
INDEX_OPTIONS_SIZE = 64
TAXONOMY_MAGIC = b"K2TAXDAT"
TAXONOMY_NODE_SIZE = 56
STANDARD_RANKS = (
    ("superkingdom", "d"),
    ("phylum", "p"),
    ("class", "c"),
    ("order", "o"),
    ("family", "f"),
    ("genus", "g"),
    ("species", "s"),
)


def open_text(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def clean_name(value: str) -> str:
    cleaned = " ".join(value.replace(";", " ").replace("\t", " ").split())
    return cleaned or "unnamed"


class K2Options:
    """Portable reader for Kun-peng/Kraken 2 ``opts.k2d``."""

    def __init__(self, path: Path):
        data = path.read_bytes()
        if len(data) != INDEX_OPTIONS_SIZE:
            raise ValueError(
                f"{path}: expected {INDEX_OPTIONS_SIZE} bytes, observed "
                f"{len(data)}"
            )
        (
            self.k,
            self.l,
            self.spaced_seed_mask,
            self.toggle_mask,
        ) = struct.unpack_from("<QQQQ", data, 0)
        self.dna_database = bool(data[32])
        self.minimum_acceptable_hash_value = struct.unpack_from("<Q", data, 40)[0]
        (
            self.reverse_complement_version,
            self.database_version,
            self.database_type,
        ) = struct.unpack_from("<iii", data, 48)
        if (
            not self.dna_database
            or not 1 <= self.l <= 31
            or self.k <= self.l
            or self.k - self.l < 2
            or self.reverse_complement_version != 1
            or self.database_version != 0
            or self.database_type != 0
        ):
            raise ValueError(f"{path}: unsupported index options")
        self.mask = (1 << (2 * self.l)) - 1
        if self.spaced_seed_mask & ~self.mask:
            raise ValueError(f"{path}: minimizer masks exceed l-mer width")
        self.raw_toggle_mask = self.toggle_mask
        self.toggle_mask &= self.mask

    def as_dict(self) -> dict[str, int | bool]:
        return {
            "k": self.k,
            "l": self.l,
            "spaced_seed_mask": self.spaced_seed_mask,
            "toggle_mask": self.raw_toggle_mask,
            "effective_toggle_mask": self.toggle_mask,
            "dna_database": self.dna_database,
            "minimum_acceptable_hash_value": (
                self.minimum_acceptable_hash_value
            ),
            "reverse_complement_version": self.reverse_complement_version,
            "database_version": self.database_version,
            "database_type": self.database_type,
        }


class K2Taxonomy:
    """Read the compact taxonomy used by Kraken 2 and Kun-peng."""

    def __init__(self, path: Path):
        with path.open("rb") as handle:
            if handle.read(len(TAXONOMY_MAGIC)) != TAXONOMY_MAGIC:
                raise ValueError(f"{path}: malformed taxonomy magic")
            header = handle.read(24)
            if len(header) != 24:
                raise ValueError(f"{path}: truncated taxonomy header")
            node_count, name_length, rank_length = struct.unpack("<QQQ", header)
            if node_count < 2:
                raise ValueError(f"{path}: taxonomy has too few nodes")
            nodes = []
            for _ in range(node_count):
                raw_node = handle.read(TAXONOMY_NODE_SIZE)
                if len(raw_node) != TAXONOMY_NODE_SIZE:
                    raise ValueError(f"{path}: truncated taxonomy node array")
                nodes.append(struct.unpack("<7Q", raw_node))
            name_data = handle.read(name_length)
            rank_data = handle.read(rank_length)
            if len(name_data) != name_length or len(rank_data) != rank_length:
                raise ValueError(f"{path}: truncated taxonomy string data")
            if handle.read(1):
                raise ValueError(f"{path}: unexpected trailing taxonomy data")

        self.path = path
        self.node_count = node_count
        self.parents = [int(node[0]) for node in nodes]
        self.name_offsets = [int(node[3]) for node in nodes]
        self.rank_offsets = [int(node[4]) for node in nodes]
        self.external_ids = [int(node[5]) for node in nodes]
        self.name_data = name_data
        self.rank_data = rank_data
        self.external_to_internal = {
            external: internal
            for internal, external in enumerate(self.external_ids)
        }
        if len(self.external_to_internal) != node_count:
            raise ValueError(f"{path}: duplicate external taxonomy IDs")
        for internal, parent in enumerate(self.parents):
            if parent >= node_count:
                raise ValueError(
                    f"{path}: node {internal} has invalid parent {parent}"
                )

    @staticmethod
    def _cstring(data: bytes, offset: int) -> str:
        if offset >= len(data):
            raise ValueError(f"taxonomy string offset is out of bounds: {offset}")
        end = data.find(b"\0", offset)
        if end < 0:
            raise ValueError("taxonomy string is not NUL terminated")
        return data[offset:end].decode("utf-8")

    def name(self, internal: int) -> str:
        return self._cstring(self.name_data, self.name_offsets[internal])

    def rank(self, internal: int) -> str:
        return self._cstring(self.rank_data, self.rank_offsets[internal])

    def lineage_internal(self, internal: int) -> list[int]:
        lineage: list[int] = []
        observed: set[int] = set()
        while internal not in observed:
            lineage.append(internal)
            observed.add(internal)
            parent = self.parents[internal]
            if parent == internal:
                break
            internal = parent
        else:
            raise ValueError("taxonomy contains a parent cycle")
        return lineage

    def lowest_common_ancestor(self, left: int, right: int) -> int:
        if left == 0:
            return right
        if right == 0:
            return left
        right_lineage = set(self.lineage_internal(right))
        return next(
            (
                internal
                for internal in self.lineage_internal(left)
                if internal in right_lineage
            ),
            0,
        )

    def resolve(self, weighted_counts: Counter[int]) -> int:
        """Apply Kraken's ancestor-supported resolve-tree rule."""

        internal_counts: Counter[int] = Counter()
        for external, count in weighted_counts.items():
            try:
                internal = self.external_to_internal[external]
            except KeyError as exc:
                raise ValueError(
                    f"hit taxid {external} is absent from {self.path}"
                ) from exc
            internal_counts[internal] += count

        best_internal = 0
        best_score = 0
        for candidate in sorted(
            internal_counts,
            key=lambda value: self.external_ids[value],
        ):
            score = sum(
                internal_counts.get(ancestor, 0)
                for ancestor in self.lineage_internal(candidate)
            )
            if score > best_score:
                best_internal = candidate
                best_score = score
            elif score == best_score:
                best_internal = self.lowest_common_ancestor(
                    best_internal,
                    candidate,
                )
        return self.external_ids[best_internal]

    def taxvamb_lineage(self, external: int) -> str:
        try:
            internal = self.external_to_internal[external]
        except KeyError as exc:
            raise ValueError(
                f"classification taxid {external} is absent from {self.path}"
            ) from exc
        standard: dict[str, tuple[str, int]] = {}
        for node in reversed(self.lineage_internal(internal)):
            rank = self.rank(node)
            if rank in {item[0] for item in STANDARD_RANKS}:
                standard[rank] = (self.name(node), self.external_ids[node])
        known = [standard.get(rank) for rank, _ in STANDARD_RANKS]
        deepest = max(
            (index for index, value in enumerate(known) if value is not None),
            default=-1,
        )
        if deepest < 0:
            return ""
        output: list[str] = []
        for index in range(deepest + 1):
            value = known[index]
            prefix = STANDARD_RANKS[index][1]
            if value is None:
                next_taxid = next(
                    known[future][1]
                    for future in range(index + 1, deepest + 1)
                    if known[future] is not None
                )
                output.append(f"{prefix}__unclassified_for_{next_taxid}")
            else:
                name, taxid = value
                output.append(
                    f"{prefix}__{clean_name(name)}__taxid_{taxid}"
                )
        return ";".join(output)


def reverse_complement_bits(value: int, length: int) -> int:
    output = 0
    for _ in range(length):
        output = (output << 2) | (3 - (value & 3))
        value >>= 2
    return output


def minimizer_event_weights(sequence: str, options: K2Options) -> list[int]:
    """Reproduce seqkmer 0.1.5 event emission and recover k-mer spans."""

    window: deque[tuple[int, int]] = deque()
    window_count = 0
    lmer_value = 0
    valid_bases = 0
    kmer_index = 0
    previous_event: Optional[int] = None
    weights: list[int] = []
    window_capacity = options.k - options.l

    for base in sequence.encode("ascii"):
        code = DNA_CODES.get(base)
        if code is None:
            window.clear()
            window_count = 0
            lmer_value = 0
            valid_bases = 0
            continue
        lmer_value = ((lmer_value << 2) | code) & options.mask
        valid_bases += 1
        if valid_bases < options.l:
            continue

        canonical = min(
            lmer_value,
            reverse_complement_bits(lmer_value, options.l),
        )
        candidate = (
            (canonical & options.spaced_seed_mask) ^ options.toggle_mask
        )
        while window and window[-1][0] > candidate:
            window.pop()
        changed = (
            (not window and window_count >= window_capacity)
            or window_count == window_capacity
        )
        window.append((candidate, window_count))
        while (
            window
            and window_count >= window_capacity
            and window[0][1] < window_count - window_capacity
        ):
            window.popleft()
            changed = True
        window_count += 1
        if window_count <= window_capacity:
            continue

        if changed:
            if previous_event is not None:
                weights[-1] = kmer_index - previous_event
            weights.append(0)
            previous_event = kmer_index
        kmer_index += 1

    if previous_event is not None:
        weights[-1] = kmer_index - previous_event
    if any(weight <= 0 for weight in weights):
        raise ValueError("failed to reconstruct positive minimizer spans")
    return weights


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: Optional[str] = None
    sequence: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence)
                identifier = line[1:].split(None, 1)[0]
                if not identifier:
                    raise ValueError(f"{path}: empty FASTA ID at line {line_number}")
                if identifier in records:
                    raise ValueError(f"{path}: duplicate FASTA ID {identifier!r}")
                sequence = []
            elif line:
                if identifier is None:
                    raise ValueError(
                        f"{path}: sequence before first FASTA header"
                    )
                sequence.append(line)
    if identifier is not None:
        records[identifier] = "".join(sequence)
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return records


def parse_event_taxids(value: str, *, path: Path, line_number: int) -> list[int]:
    events: list[int] = []
    for token in value.split():
        if token == "|:|":
            raise ValueError(f"{path}: paired hit streams are not supported")
        fields = token.rsplit(":", 1)
        if (
            len(fields) != 2
            or not fields[0].isdigit()
            or not fields[1].isdigit()
            or int(fields[1]) < 1
        ):
            raise ValueError(
                f"{path}: invalid hit token {token!r} at line {line_number}"
            )
        events.extend([int(fields[0])] * int(fields[1]))
    return events


def encode_weighted_hits(
    event_taxids: Sequence[int],
    event_weights: Sequence[int],
) -> tuple[str, Counter[int]]:
    encoded: list[tuple[int, int]] = []
    counts: Counter[int] = Counter()
    for taxid, weight in zip(event_taxids, event_weights):
        if encoded and encoded[-1][0] == taxid:
            encoded[-1] = (taxid, encoded[-1][1] + weight)
        else:
            encoded.append((taxid, weight))
        if taxid:
            counts[taxid] += weight
    return " ".join(f"{taxid}:{weight}" for taxid, weight in encoded), counts


def read_report_counts(path: Path) -> tuple[int, int]:
    """Return raw classified/unclassified sequence counts from kreport2."""

    classified: Optional[int] = None
    unclassified = 0
    observed_unclassified = False
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t", 5)
            if len(fields) != 6:
                raise ValueError(
                    f"{path}: malformed report line {line_number}"
                )
            clade_count = fields[1].strip()
            rank = fields[3].strip()
            taxid = fields[4].strip()
            if not clade_count.isdigit():
                raise ValueError(
                    f"{path}: invalid clade count at line {line_number}"
                )
            count = int(clade_count)
            if rank == "U" and taxid == "0":
                if observed_unclassified:
                    raise ValueError(f"{path}: duplicate unclassified row")
                observed_unclassified = True
                unclassified = count
            elif rank == "R" and taxid == "1":
                if classified is not None:
                    raise ValueError(f"{path}: duplicate root row")
                classified = count
    return (classified or 0), unclassified


def resolve(
    classification_output: Path,
    fasta: Path,
    options_path: Path,
    taxonomy_path: Path,
    output: Path,
    *,
    classification_report: Optional[Path] = None,
    minimum_hit_groups: int = 3,
) -> dict[str, object]:
    if output.exists():
        raise ValueError(f"compatibility output already exists: {output}")
    if minimum_hit_groups < 1:
        raise ValueError("minimum_hit_groups must be positive")
    options = K2Options(options_path)
    taxonomy = K2Taxonomy(taxonomy_path)
    records = read_fasta(fasta)
    observed: set[str] = set()
    corrected_calls = 0
    corrected_status = 0
    classified = 0
    unclassified = 0
    total_events = 0
    total_kmer_windows = 0
    raw_classified = 0
    raw_unclassified = 0
    omitted_zero_hit = 0
    output.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open_text(classification_output) as source, output.open(
            "wt",
            encoding="utf-8",
            newline="\n",
        ) as target:
            for line_number, raw_line in enumerate(source, start=1):
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(
                        f"{classification_output}: expected 5 fields at "
                        f"line {line_number}"
                    )
                raw_status, identifier, raw_taxid, raw_length, raw_hits = fields
                if raw_status not in {"C", "U"}:
                    raise ValueError(
                        f"{classification_output}: invalid status "
                        f"{raw_status!r} at line {line_number}"
                    )
                raw_classified += int(raw_status == "C")
                raw_unclassified += int(raw_status == "U")
                if identifier in observed:
                    raise ValueError(
                        f"duplicate classification for {identifier!r}"
                    )
                observed.add(identifier)
                if identifier not in records:
                    raise ValueError(
                        f"classification contains unknown contig {identifier!r}"
                    )
                sequence = records[identifier]
                if not raw_length.isdigit() or int(raw_length) != len(sequence):
                    raise ValueError(
                        f"{identifier}: classified length {raw_length!r} "
                        f"does not equal FASTA length {len(sequence)}"
                    )
                event_taxids = parse_event_taxids(
                    raw_hits,
                    path=classification_output,
                    line_number=line_number,
                )
                event_weights = minimizer_event_weights(sequence, options)
                if len(event_taxids) != len(event_weights):
                    raise ValueError(
                        f"{identifier}: Kun-peng hit events "
                        f"{len(event_taxids)} != reconstructed events "
                        f"{len(event_weights)}"
                    )
                weighted_hits, weighted_counts = encode_weighted_hits(
                    event_taxids,
                    event_weights,
                )
                hit_groups = sum(taxid != 0 for taxid in event_taxids)
                call = (
                    taxonomy.resolve(weighted_counts)
                    if hit_groups >= minimum_hit_groups and weighted_counts
                    else 0
                )
                status = "C" if call else "U"
                corrected_status += int(status != raw_status)
                corrected_calls += int(
                    status != raw_status or str(call) != raw_taxid
                )
                classified += int(status == "C")
                unclassified += int(status == "U")
                total_events += len(event_taxids)
                total_kmer_windows += sum(event_weights)
                target.write(
                    f"{status}\t{identifier}\t{call}\t{len(sequence)}\t"
                    f"{weighted_hits}\n"
                )
        missing = sorted(set(records) - observed)
        if missing:
            if classification_report is None:
                raise ValueError(
                    "classification is missing FASTA contigs: "
                    + ", ".join(missing[:20])
                )
            report_classified, report_unclassified = read_report_counts(
                classification_report
            )
            if report_classified + report_unclassified != len(records):
                raise ValueError(
                    "Kun-peng report total does not match FASTA contigs: "
                    f"{report_classified + report_unclassified} != "
                    f"{len(records)}"
                )
            if report_classified != raw_classified:
                raise ValueError(
                    "Kun-peng output/report classified counts differ: "
                    f"{raw_classified} != {report_classified}"
                )
            expected_omitted = report_unclassified - raw_unclassified
            if expected_omitted != len(missing):
                raise ValueError(
                    "Kun-peng missing-contig count is inconsistent with "
                    "the report's zero-hit omission count: "
                    f"{len(missing)} != {expected_omitted}"
                )
            with output.open(
                "at",
                encoding="utf-8",
                newline="\n",
            ) as target:
                for identifier in missing:
                    sequence = records[identifier]
                    windows = sum(minimizer_event_weights(sequence, options))
                    hit_stream = f"0:{windows}" if windows else ""
                    target.write(
                        f"U\t{identifier}\t0\t{len(sequence)}\t"
                        f"{hit_stream}\n"
                    )
                    unclassified += 1
                    total_kmer_windows += windows
            omitted_zero_hit = len(missing)
        elif classification_report is not None:
            report_classified, report_unclassified = read_report_counts(
                classification_report
            )
            if (
                report_classified != raw_classified
                or report_unclassified != raw_unclassified
                or report_classified + report_unclassified != len(records)
            ):
                raise ValueError(
                    "Kun-peng output/report sequence counts are inconsistent"
                )
    except Exception:
        output.unlink(missing_ok=True)
        raise

    return {
        "schema_version": 1,
        "algorithm": "kraken2_run_length_compatibility",
        "classification_output": str(classification_output.resolve()),
        "fasta": str(fasta.resolve()),
        "options": str(options_path.resolve()),
        "taxonomy": str(taxonomy_path.resolve()),
        "contigs": len(records),
        "classified": classified,
        "unclassified": unclassified,
        "raw_output_contigs": len(observed),
        "raw_classified": raw_classified,
        "raw_unclassified": raw_unclassified,
        "raw_omitted_zero_hit_contigs": omitted_zero_hit,
        "corrected_calls": corrected_calls,
        "corrected_statuses": corrected_status,
        "minimum_hit_groups": minimum_hit_groups,
        "minimizer_events": total_events,
        "kmer_windows": total_kmer_windows,
        "index_options": options.as_dict(),
        "taxonomy_nodes": taxonomy.node_count,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classification-output", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--options", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--classification-report", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--minimum-hit-groups", type=int, default=3)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    summary = resolve(
        args.classification_output,
        args.fasta,
        args.options,
        args.taxonomy,
        args.output,
        classification_report=args.classification_report,
        minimum_hit_groups=args.minimum_hit_groups,
    )
    encoded = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.summary is not None:
        if args.summary.exists():
            raise ValueError(f"summary already exists: {args.summary}")
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
