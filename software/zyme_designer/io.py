# Source guide (2026-10-20; documentation comments only).
# Purpose: Read sequences and serialize design outputs
# Inputs: FASTA, result records and destination paths
# Generates/returns: Files requested by CLI/UI callers
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Input/output, provenance and identifier resolution without implicit network access."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

from .core import annotate_header, normalize_rna


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(2**22), b""):
            h.update(block)
    return h.hexdigest()


def fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as f:
        header, parts = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, normalize_rna("".join(parts), allow_ambiguity=True)
                header, parts = line[1:].strip(), []
                if not header:
                    raise ValueError("Empty FASTA header")
            elif header is None:
                raise ValueError("FASTA must begin with a >header")
            else:
                parts.append(line)
        if header is not None:
            yield header, normalize_rna("".join(parts), allow_ambiguity=True)


def resolve_target(path, selector=None):
    found = []
    seen = set()
    for header, seq in fasta_records(path):
        accession, gene = annotate_header(header)
        if selector is None or selector == accession or selector == gene:
            if accession in seen:
                raise ValueError(f"Duplicate target accession in reference: {accession}")
            seen.add(accession)
            found.append({"accession": accession, "gene": gene, "header": header, "sequence": seq})
    if not found:
        raise ValueError(f"No exact accession/gene match for {selector!r}; use a versioned transcript accession")
    if len(found) != 1:
        names = ", ".join(r["accession"] for r in found[:20])
        raise ValueError(f"Ambiguous input ({len(found)} transcripts: {names}); choose a specific versioned transcript")
    return found[0]


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, fields=None):
    fields = fields or sorted({key for row in rows for key in row})
    if not fields:
        fields = ["status"]
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {}
            for key in fields:
                value = row.get(key)
                if isinstance(value, (dict, list, tuple)):
                    value = json.dumps(value, ensure_ascii=False, allow_nan=False)
                # Do not allow arbitrary database headers to become Excel formulae.
                if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
                    value = "'" + value
                clean[key] = value
            writer.writerow(clean)
