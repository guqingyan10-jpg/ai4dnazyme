# Source guide (2026-10-20; documentation comments only).
# Purpose: Search user-supplied RNA reference for recognition-arm matches
# Inputs: Reference FASTA and candidate/mismatch settings
# Generates/returns: Hit records and specificity summaries returned to the design run
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Streaming, substitution-only off-target screen for two RNA-binding arms.

This is a *binding-sequence* screen, not a cleavage or clinical safety predictor.
Both substrate arms must match within a TOTAL mismatch budget. The central two
nucleotides are reported separately and never included in that arm budget.
All sense-strand windows are covered by exact/one-substitution arm seeds: with
at most two substitutions, at least one of the two arms has at most one.

The optional ``pyahocorasick`` package accelerates matching many seed patterns.
Without it, overlapping ``str.find`` searches implement the identical screen.
Ambiguous IUPAC bases remain ambiguous; overlapping ambiguous windows are
skipped, NOT interpreted as matched canonical nucleotides.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import gzip
import hashlib
import io
import re
import time
from typing import Callable, Iterable

from .core import motif_matches, normalize_motif_pattern


_CANONICAL = frozenset("ACGU")
_IUPAC = frozenset("ACGUTNRYSWKMBDHV")
_ACCESSION = re.compile(
    r"(?<![A-Za-z0-9])(?:[NX][MR]_[0-9]+|ENS[A-Z]*T[0-9]+)(?:\.[0-9]+)?"
)
_GENE_ID = re.compile(r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$", re.I)


def _normalize_gene(value: str) -> str:
    value = value.strip().strip("\"'").upper()
    if _GENE_ID.fullmatch(value):
        value = value.split(".", 1)[0]
    return value


def _header_metadata(header: str) -> tuple[str, tuple[str, ...], str, str]:
    """Return accession, normalized gene aliases, display gene, source label."""
    token = header.split(None, 1)[0]
    # Handles ref|NM_...|, Ensembl, and NM_001168.3_verified_from_PS.
    match = _ACCESSION.search(token)
    accession = match.group(0) if match else token
    symbols: list[str] = []
    gene_ids: list[str] = []
    sources: list[str] = []
    for key, value in re.findall(
        r"(?:^|[\s\[])\b(gene_symbol|gene_name|gene)\s*[:=]\s*([^\s\]\[;]+)",
        header,
        flags=re.I,
    ):
        normalized = _normalize_gene(value.rstrip(","))
        if not normalized:
            continue
        (gene_ids if _GENE_ID.fullmatch(normalized) else symbols).append(normalized)
        sources.append(key.lower())
    # NCBI RefSeq descriptive headers conventionally contain (GENE_SYMBOL).
    for value in re.findall(r"\(([A-Za-z][A-Za-z0-9_.-]{0,63})\)", header):
        if not _ACCESSION.fullmatch(value):
            symbols.append(_normalize_gene(value))
            sources.append("parenthetical_symbol")
    # GENCODE: ENST|ENSG|OTTHUMG|OTTHUMT|transcript_name|gene_name|length|...
    parts = token.split("|")
    if len(parts) >= 6 and parts[0].startswith("ENS") and _GENE_ID.fullmatch(parts[1]):
        gene_ids.append(_normalize_gene(parts[1]))
        if parts[5] and parts[5] != "-":
            symbols.append(_normalize_gene(parts[5]))
        sources.append("gencode_pipe")
    aliases = tuple(sorted(set(symbols + gene_ids)))
    display_gene = symbols[0] if symbols else (gene_ids[0] if gene_ids else "")
    return accession, aliases, display_gene, ";".join(sorted(set(sources)))


def _iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    """Read plain/gzip FASTA; reject malformed input rather than silently omit it."""
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rb") as raw, io.TextIOWrapper(raw, encoding="utf-8-sig") as stream:
        header: str | None = None
        sequence: list[str] = []
        for line_number, line in enumerate(stream, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith(">"):
                if header is not None:
                    if not sequence:
                        raise ValueError(f"Empty FASTA record: {header!r}")
                    yield header, "".join(sequence).upper().replace("T", "U")
                header = stripped[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header on line {line_number}")
                sequence = []
                continue
            if header is None:
                raise ValueError(f"Sequence before FASTA header on line {line_number}")
            letters = "".join(stripped.split()).upper()
            invalid = sorted(set(letters) - _IUPAC)
            if invalid:
                raise ValueError(
                    f"Invalid FASTA symbols {invalid!r} on line {line_number}; "
                    "only RNA/DNA IUPAC nucleotide symbols are accepted"
                )
            sequence.append(letters)
        if header is None:
            raise ValueError("FASTA contains no records")
        if not sequence:
            raise ValueError(f"Empty FASTA record: {header!r}")
        yield header, "".join(sequence).upper().replace("T", "U")


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _one_mismatch_seeds(sequence: str, include_substitutions: bool) -> Iterable[str]:
    yield sequence
    if include_substitutions:
        for index, original in enumerate(sequence):
            for base in "ACGU":
                if base != original:
                    yield sequence[:index] + base + sequence[index + 1 :]


def _ambiguous_window_count(positions: list[int], length: int, n_nt: int) -> int:
    """Count the union of invalid start intervals without enumerating every window."""
    last_start = n_nt - length
    if last_start < 0:
        return 0
    total, current_start, current_end = 0, -1, -2
    for position in positions:
        lower, upper = max(0, position - length + 1), min(last_start, position)
        if lower > upper:
            continue
        if lower > current_end + 1:
            total += max(0, current_end - current_start + 1)
            current_start, current_end = lower, upper
        else:
            current_end = max(current_end, upper)
    return total + max(0, current_end - current_start + 1)


def _classification(intended: bool, genes: tuple[str, ...], target_genes: set[str]) -> str:
    if intended:
        return "intended_target"
    if not genes or not target_genes:
        return "unknown_gene"
    return "target_isoform" if set(genes) & target_genes else "other_gene"


def _csv_safe_text(value: str) -> str:
    # Prevent spreadsheet formula execution when users open a hit CSV in Excel.
    return "'" + value if value.lstrip().startswith(("=", "+", "-", "@")) else value


def _compatible_motif(motif: str, pattern: str) -> bool:
    """Declared two-base recognition policy, not proof of biochemical cleavage."""
    return motif_matches(motif, pattern)


def scan_transcriptome(
    candidates: Iterable[dict],
    fasta_path: str | Path,
    target_accession: str,
    target_gene: str = "",
    max_mismatches: int = 2,
    hit_limit: int = 10000,
    progress: Callable[[dict], None] | None = None,
) -> tuple[dict[str, dict], list[dict], dict]:
    """Screen every sense-strand transcript window within the arm mismatch limit.

    ``candidates`` contain ``candidate_id``, ``target_rna``, ``left_length``,
    ``right_length``, ``start_1based`` and ``end_1based``. Arm lengths are >= 1;
    the target sequence length must equal left + 2 + right. Allowed mismatch
    limits are 0, 1, or 2; gaps, bulges, G-U wobble and reverse strands are not
    modeled. A substitution counts as a mismatch even if it permits wobble.

    All counts remain complete even after ``hit_limit`` exported rows are
    reached. Export rows are first hits in FASTA/candidate/start order, not a
    random sample. Exact accession matching is version-sensitive; other
    locations in the intended transcript are counted and labeled separately.
    Gene relations use annotated symbols/IDs, NEVER substring gene guessing.
    Missing annotation, or unresolved target gene, remains ``unknown_gene``.
    When the intended accession occurs in FASTA, its gene aliases supplement
    ``target_gene``; this also enables late, retrospective isoform assignment.
    ``progress`` receives record/nt/elapsed counters every 100 records and at
    completion. No external service, database download, or file writes occur.
    """
    started = time.monotonic()
    if isinstance(max_mismatches, bool) or max_mismatches not in (0, 1, 2):
        raise ValueError("max_mismatches must be 0, 1, or 2")
    if not isinstance(hit_limit, int) or isinstance(hit_limit, bool) or hit_limit < 0:
        raise ValueError("hit_limit must be a nonnegative integer")
    if not str(target_accession).strip():
        raise ValueError("target_accession is required")
    target_accession = _header_metadata(str(target_accession).strip())[0]
    target_genes = {_normalize_gene(target_gene)} if str(target_gene).strip() else set()
    normalized: list[dict] = []
    seen_candidate_ids: set[str] = set()
    for candidate in candidates:
        identifier = str(candidate.get("candidate_id", ""))
        if not identifier or identifier in seen_candidate_ids:
            raise ValueError("Every candidate_id must be nonempty and unique")
        seen_candidate_ids.add(identifier)
        left, right = candidate.get("left_length"), candidate.get("right_length")
        if any(not isinstance(n, int) or isinstance(n, bool) or n < 1 for n in (left, right)):
            raise ValueError(f"Candidate {identifier}: arm lengths must be positive integers")
        rna = str(candidate.get("target_rna", "")).upper().replace("T", "U")
        if len(rna) != left + right + 2 or set(rna) - _CANONICAL:
            raise ValueError(f"Candidate {identifier}: target_rna must be unambiguous and match arm lengths")
        start, end = candidate.get("start_1based"), candidate.get("end_1based")
        if (
            not isinstance(start, int) or isinstance(start, bool) or start < 1
            or not isinstance(end, int) or isinstance(end, bool)
            or end - start + 1 != len(rna)
        ):
            raise ValueError(f"Candidate {identifier}: inconsistent 1-based target coordinates")
        pattern = normalize_motif_pattern(candidate.get("motif_pattern", "AG"))
        normalized.append({**candidate, "candidate_id": identifier, "target_rna": rna, "motif_pattern": pattern})
    if not normalized:
        raise ValueError("At least one candidate is required")
    path = Path(fasta_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    before = path.stat()
    input_hash = _hash_file(path)
    routes: dict[str, set[tuple[int, int]]] = defaultdict(set)
    # At max 1 substitution, one arm MUST be exact; at max 2, one arm <= 1.
    for index, candidate in enumerate(normalized):
        left_length = candidate["left_length"]
        rna = candidate["target_rna"]
        for arm, offset in ((rna[:left_length], 0), (rna[left_length + 2 :], left_length + 2)):
            for seed in _one_mismatch_seeds(arm, max_mismatches == 2):
                routes[seed].add((index, offset))
    ordered_routes = {seed: tuple(sorted(value)) for seed, value in sorted(routes.items())}
    try:
        import ahocorasick  # type: ignore[import-not-found]
    except ImportError:
        automaton = None
        backend = "stdlib_str_find"
    else:
        automaton = ahocorasick.Automaton()
        for seed, values in ordered_routes.items():
            automaton.add_word(seed, (len(seed), values))
        automaton.make_automaton()
        backend = "pyahocorasick"
    bins: list[dict[tuple[bool, tuple[str, ...], str], Counter]] = [defaultdict(Counter) for _ in normalized]
    counts = [Counter() for _ in normalized]
    lengths = sorted({len(c["target_rna"]) for c in normalized})
    window_counts = {length: Counter() for length in lengths}
    exported: list[dict] = []
    n_records = n_nt = n_ambiguous = n_hits = 0
    target_records_found = 0
    duplicate_records = 0
    duplicate_examples: list[str] = []
    accessions_seen: set[str] = set()
    header_genes_for_target: set[str] = set()
    gene_parents: dict[str, str] = {}

    def gene_root(alias: str) -> str:
        gene_parents.setdefault(alias, alias)
        current = alias
        while gene_parents[current] != current:
            current = gene_parents[current]
        while gene_parents[alias] != alias:
            parent = gene_parents[alias]
            gene_parents[alias] = current
            alias = parent
        return current

    def connect_header_aliases(aliases: tuple[str, ...]) -> None:
        # Only link aliases explicitly annotated together in a FASTA header.
        # Prefer symbols over Ensembl IDs for stable human-readable gene keys.
        roots = {gene_root(alias) for alias in aliases}
        if roots:
            canonical = min(roots, key=lambda item: (bool(_GENE_ID.fullmatch(item)), item))
            for root in roots:
                gene_parents[root] = canonical

    def emit_progress(complete: bool = False) -> None:
        if progress is not None:
            progress({
                "records_scanned": n_records, "nucleotides_scanned": n_nt,
                "binding_hits_counted": n_hits, "hit_rows_exported": len(exported),
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "backend": backend, "complete": complete,
            })

    for header, sequence in _iter_fasta(path):
        accession, genes, display_gene, gene_source = _header_metadata(header)
        connect_header_aliases(genes)
        intended = accession == target_accession
        if intended:
            target_records_found += 1
            target_genes.update(genes)
            header_genes_for_target.update(genes)
        if accession in accessions_seen:
            duplicate_records += 1
            if accession not in duplicate_examples and len(duplicate_examples) < 20:
                duplicate_examples.append(accession)
        accessions_seen.add(accession)
        n_records += 1
        n_nt += len(sequence)
        ambiguous = [i for i, base in enumerate(sequence) if base not in _CANONICAL]
        n_ambiguous += len(ambiguous)
        for length in lengths:
            total = max(0, len(sequence) - length + 1)
            skipped = _ambiguous_window_count(ambiguous, length, len(sequence))
            window_counts[length].update(total_windows=total, ambiguous_windows_skipped=skipped, unambiguous_windows=total - skipped)
        starts: list[set[int]] = [set() for _ in normalized]
        if automaton is not None:
            for seed_end, (seed_length, values) in automaton.iter(sequence):
                seed_start = seed_end - seed_length + 1
                for index, offset in values:
                    starts[index].add(seed_start - offset)
        else:
            for seed, values in ordered_routes.items():
                seed_start = sequence.find(seed)
                while seed_start >= 0:
                    for index, offset in values:
                        starts[index].add(seed_start - offset)
                    seed_start = sequence.find(seed, seed_start + 1)
        for index, candidate in enumerate(normalized):
            left_length = candidate["left_length"]
            length = len(candidate["target_rna"])
            target = candidate["target_rna"]
            # Do not allocate one Counter per candidate/transcript with no hits;
            # most transcripts in a full database have no qualifying window.
            bucket = None
            for start in sorted(starts[index]):
                if start < 0 or start + length > len(sequence):
                    continue
                window = sequence[start : start + length]
                if set(window) - _CANONICAL:
                    continue
                left_mm = sum(a != b for a, b in zip(target[:left_length], window[:left_length]))
                right_mm = sum(a != b for a, b in zip(target[left_length + 2 :], window[left_length + 2 :]))
                arm_mm = left_mm + right_mm
                if arm_mm > max_mismatches:
                    continue
                motif = window[left_length : left_length + 2]
                is_ag = motif == "AG"
                compatible = _compatible_motif(motif, candidate["motif_pattern"])
                is_intended_site = intended and start + 1 == candidate["start_1based"]
                if bucket is None:
                    bucket = bins[index][(intended, genes, display_gene)]
                bucket[(motif, arm_mm)] += 1
                counts[index]["binding_hits_total"] += 1
                counts[index]["intended_site_hits" if is_intended_site else "other_sites_on_target_hits"] += int(intended)
                n_hits += 1
                if len(exported) >= hit_limit:
                    continue
                row = {
                    "candidate_id": candidate["candidate_id"],
                    "transcript_accession": accession,
                    "transcript_header": header,
                    "gene": display_gene,
                    "gene_annotation_source": gene_source,
                    "start_1based": start + 1, "end_1based": start + length,
                    "window_rna": window, "target_rna": target,
                    "left_arm_mismatches": left_mm, "right_arm_mismatches": right_mm,
                    "arm_mismatches": arm_mm,
                    "cut_motif": motif, "ag_site": is_ag,
                    "motif_pattern": candidate["motif_pattern"],
                    "compatible_motif": compatible,
                    "motif_compatibility_is_cleavage_proof": False,
                    "is_intended_site": is_intended_site,
                    "_intended": intended, "_genes": genes,
                }
                exported.append(row)
        if n_records == 1 or n_records % 100 == 0:
            emit_progress()
    # Propagate only documented header links, never external or fuzzy alias guesses.
    target_roots = {gene_root(alias) for alias in target_genes}
    target_genes.update(alias for alias in gene_parents if gene_root(alias) in target_roots)
    summaries: dict[str, dict] = {}
    for index, candidate in enumerate(normalized):
        summary = {
            "off_target_status": "screened",
            "other_gene_exact_ag": 0, "other_gene_1mm_ag": 0,
            "other_gene_2mm_ag": 0, "other_gene_non_ag": 0,
            "other_gene_exact_compatible": 0, "other_gene_1mm_compatible": 0,
            "other_gene_2mm_compatible": 0, "other_gene_noncompatible": 0,
            "unknown_gene_hits": 0, "target_isoform_hits": 0,
            "intended_target_hits": 0, "binding_hits_total": counts[index]["binding_hits_total"],
            "intended_site_hits": counts[index]["intended_site_hits"],
            "other_sites_on_target_hits": counts[index]["other_sites_on_target_hits"],
            "unique_other_gene_count": 0,
            "ag_binding_hits_total": 0, "non_ag_binding_hits_total": 0,
            "motif_pattern": candidate["motif_pattern"],
            "compatible_binding_hits_total": 0, "noncompatible_binding_hits_total": 0,
            "unknown_gene_compatible_hits": 0, "target_isoform_compatible_hits": 0,
            "intended_target_compatible_hits": 0,
        }
        other_genes: set[str] = set()
        for (intended, genes, display_gene), bucket in bins[index].items():
            relation = _classification(intended, genes, target_genes)
            total = sum(bucket.values())
            if relation == "other_gene":
                if total:
                    other_genes.add(gene_root(genes[0]))
                for mismatches in (0, 1, 2):
                    label = "exact" if mismatches == 0 else str(mismatches) + "mm"
                    summary[f"other_gene_{label}_ag"] += bucket[("AG", mismatches)]
                    summary[f"other_gene_{label}_compatible"] += sum(n for (motif, mm), n in bucket.items()
                        if mm == mismatches and _compatible_motif(motif, candidate["motif_pattern"]))
                summary["other_gene_non_ag"] += sum(n for (motif, _), n in bucket.items() if motif != "AG")
                summary["other_gene_noncompatible"] += sum(n for (motif, _), n in bucket.items()
                    if not _compatible_motif(motif, candidate["motif_pattern"]))
            else:
                key = {"intended_target": "intended_target_hits", "target_isoform": "target_isoform_hits", "unknown_gene": "unknown_gene_hits"}[relation]
                summary[key] += total
                compatible_key = {"intended_target": "intended_target_compatible_hits",
                                  "target_isoform": "target_isoform_compatible_hits",
                                  "unknown_gene": "unknown_gene_compatible_hits"}[relation]
                summary[compatible_key] += sum(n for (motif, _), n in bucket.items()
                    if _compatible_motif(motif, candidate["motif_pattern"]))
            summary["ag_binding_hits_total"] += sum(n for (motif, _), n in bucket.items() if motif == "AG")
            summary["non_ag_binding_hits_total"] += sum(n for (motif, _), n in bucket.items() if motif != "AG")
            summary["compatible_binding_hits_total"] += sum(n for (motif, _), n in bucket.items()
                if _compatible_motif(motif, candidate["motif_pattern"]))
            summary["noncompatible_binding_hits_total"] += sum(n for (motif, _), n in bucket.items()
                if not _compatible_motif(motif, candidate["motif_pattern"]))
        summary["unique_other_gene_count"] = len(other_genes)
        summary.update(window_counts[len(candidate["target_rna"])])
        summaries[candidate["candidate_id"]] = summary
    for row in exported:
        intended, genes = row.pop("_intended"), row.pop("_genes")
        row["gene_classification"] = _classification(intended, genes, target_genes)
        row["resolved_gene_key"] = gene_root(genes[0]) if genes else ""
        for key, value in row.items():
            if isinstance(value, str):
                row[key] = _csv_safe_text(value)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("Transcriptome FASTA changed during the scan; results are not trustworthy")
    warnings: list[str] = []
    if not target_records_found:
        warnings.append("Exact target accession was not found in the screened FASTA; accession versions are significant")
    if not target_genes:
        warnings.append("Target gene could not be resolved; all non-intended transcript hits are classified unknown_gene")
    if target_gene and header_genes_for_target and _normalize_gene(target_gene) not in header_genes_for_target:
        warnings.append("Provided target_gene is absent from the exact target FASTA header; both provided and header aliases were retained")
    if duplicate_records:
        warnings.append("Duplicate accession records are counted independently; transcript-window counts are not unique-gene counts")
    if n_ambiguous:
        warnings.append("Windows containing ambiguous IUPAC nucleotides were excluded and counted explicitly")
    audit = {
        "analysis_type": "sense-strand substitution-only two-arm binding sequence screen; NOT cleavage or safety validation",
        "fasta_path": str(path), "input_file_sha256": input_hash,
        "input_file_bytes": before.st_size,
        "target_accession": target_accession, "provided_target_gene": str(target_gene),
        "resolved_target_gene_aliases": sorted(target_genes),
        "gene_identity_resolution": "symbols/IDs linked only when annotated together in FASTA headers; otherwise separate annotation keys",
        "exact_target_records_found": target_records_found,
        "records_scanned": n_records, "nucleotides_scanned": n_nt,
        "ambiguous_nucleotides": n_ambiguous,
        "window_counts_by_length": {str(k): dict(v) for k, v in window_counts.items()},
        "duplicate_accession_records": duplicate_records,
        "duplicate_accession_examples": duplicate_examples,
        "candidate_count": len(normalized), "seed_pattern_count": len(ordered_routes),
        "candidate_motif_policies": {candidate["candidate_id"]: candidate["motif_pattern"] for candidate in normalized},
        "motif_compatibility_is_cleavage_proof": False,
        "backend": backend, "max_total_arm_mismatches": max_mismatches,
        "central_motif_counted_in_arm_mismatches": False,
        "reverse_strand_searched": False, "gaps_or_bulges_modeled": False,
        "gu_wobble_modeled": False,
        "binding_hits_counted": n_hits, "hit_rows_exported": len(exported),
        "hit_export_limit": hit_limit, "hit_rows_truncated": n_hits > len(exported),
        "summary_counts_complete": True,
        "hit_export_policy": "first hits in FASTA/candidate/start order; counts are complete even when rows are capped",
        "csv_formula_escape": "string fields starting with =,+,-,@ after whitespace are prefixed with an apostrophe",
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "warnings": warnings,
    }
    emit_progress(complete=True)
    return summaries, exported, audit
