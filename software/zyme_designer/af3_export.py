# Source guide (2026-10-20; documentation comments only).
# Purpose: Assemble selected design chains and export AF3 inputs
# Inputs: Candidate, scaffold/trigger/substrate and export settings
# Generates/returns: AF3 input JSON, chain mapping and sequence exports
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Export documented AlphaFold inputs, never submit jobs or invent coordinates."""
from __future__ import annotations

import json
import re
from pathlib import Path

from .core import normalize_rna, sequence_hash
from .io import write_json


def candidate_inputs(candidate: dict, target: dict, scaffold: dict | None = None,
                     flank: int = 0, calcium_count: int = 0) -> tuple[dict, dict, dict]:
    if not isinstance(flank, int) or isinstance(flank, bool) or flank < 0:
        raise ValueError("AlphaFold context flank must be a nonnegative integer")
    if not isinstance(calcium_count, int) or isinstance(calcium_count, bool) or not 0 <= calcium_count <= 10:
        raise ValueError("Explicit modeled Ca ion count must be an integer in 0–10, not a concentration")
    sequence = normalize_rna(target["sequence"], allow_ambiguity=True)
    keys = ("start_1based", "end_1based", "cut_after_1based", "left_length", "right_length")
    if any(not isinstance(candidate.get(key), int) or isinstance(candidate[key], bool) for key in keys):
        raise ValueError("Candidate coordinates and arm lengths must be integer 1-based values")
    first, last, cut = (candidate[key] for key in keys[:3])
    left, right = candidate["left_length"], candidate["right_length"]
    if not (1 <= first <= cut < last <= len(sequence) and left > 0 and right > 0
            and last - first + 1 == left + right + 2 and cut == first + left):
        raise ValueError("Candidate coordinates, excluded central dinucleotide and arm lengths disagree")
    if normalize_rna(candidate["target_rna"]) != sequence[first - 1:last]:
        raise ValueError("Candidate target sequence disagrees with full-transcript coordinates")
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", candidate["candidate_id"]):
        raise ValueError("Candidate ID must be a safe basename")
    for key in ("left_arm_dna", "right_arm_dna", "left_partzyme", "right_partzyme"):
        if candidate.get(key) and (not isinstance(candidate[key], str) or set(candidate[key]) - set("ACGT")):
            raise ValueError("Exported DNA chains must contain only canonical uppercase A/C/G/T")
    start = max(1, candidate["start_1based"] - flank)
    end = min(len(sequence), candidate["end_1based"] + flank)
    substrate = normalize_rna(sequence[start - 1:end])
    complete = bool(candidate.get("left_partzyme") and candidate.get("right_partzyme"))
    entities = [
        ("A", "dna", candidate["left_partzyme"] if complete else candidate["left_arm_dna"],
         "left_partzyme" if complete else "left_substrate_arm_ONLY"),
        ("B", "dna", candidate["right_partzyme"] if complete else candidate["right_arm_dna"],
         "right_partzyme" if complete else "right_substrate_arm_ONLY"),
        ("C", "rna", substrate, "RNA_substrate_context"),
    ]
    if complete and scaffold and scaffold.get("trigger_rna"):
        entities.append(("D", "rna", normalize_rna(scaffold["trigger_rna"]), "trigger"))
    name = candidate["candidate_id"] + ("_scaffold_complex" if complete else "_arms_binding_only")
    local = {"name": name, "modelSeeds": [1], "dialect": "alphafold3", "version": 1,
             "sequences": [{kind: {"id": cid, "sequence": seq}}
                           for cid, kind, seq, role in entities]}
    server_types = {"dna": "dnaSequence", "rna": "rnaSequence"}
    server = {"name": name, "modelSeeds": [], "dialect": "alphafoldserver", "version": 1,
              "sequences": [{server_types[kind]: {"sequence": seq, "count": 1}}
                            for cid, kind, seq, role in entities]}
    if calcium_count:
        local["sequences"].append({"ligand": {"id": [chr(ord("E") + i) for i in range(calcium_count)],
                                               "ccdCodes": ["CA"]}})
        server["sequences"].append({"ion": {"ion": "CA", "count": calcium_count}})
    manifest = {
        "candidate_id": candidate["candidate_id"], "submitted_to_any_service": False,
        "prediction_executed": False, "coordinate_file_generated": False,
        "scope": "fixed_scaffold_complex_candidate" if complete else "binding_arms_without_catalytic_core",
        "target_accession": target["accession"],
        "full_transcript_sequence_sha256_rna": sequence_hash(sequence),
        "RNA_chain": "C", "RNA_chain_transcript_start_1based": start,
        "RNA_chain_transcript_end_1based": end,
        "RNA_chain_site_start_1based": candidate["start_1based"] - start + 1,
        "RNA_chain_site_end_1based": candidate["end_1based"] - start + 1,
        "RNA_chain_cut_after_1based": candidate["cut_after_1based"] - start + 1,
        "chains": [{"id": cid, "type": kind, "role": role, "sequence": seq}
                   for cid, kind, seq, role in entities],
        "modeled_Ca_ion_count": calcium_count,
        "warnings": [
            "These are job requests, not structures or evidence of cleavage.",
            "A short RNA context is not a full-transcript or in-cell accessibility prediction.",
            "No catalytic core is generated from recognition arms alone.",
            "Modeled ion count is NOT a molar salt concentration or validated catalytic coordination.",
            "External predictions/output licenses remain separate from this project's MIT source license.",
            "Low-confidence nucleic-acid predictions do not determine candidate ranking.",
        ],
        "official_format_sources": [
            "https://github.com/google-deepmind/alphafold/blob/main/server/README.md",
            "https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md",
        ],
    }
    return local, server, manifest


def export_candidates(rows: list[dict], target: dict, outdir: Path,
                      scaffold: dict | None = None, top: int = 5,
                      flank: int = 0, calcium_count: int = 0) -> dict:
    if not isinstance(top, int) or isinstance(top, bool) or top < 0:
        raise ValueError("AlphaFold export count must be nonnegative")
    if top == 0:
        return {"status": "disabled", "prediction_executed": False}
    if outdir.exists():
        raise FileExistsError("AlphaFold export directory already exists; will not overwrite")
    outdir.mkdir(parents=True)
    combined, exported = [], []
    for candidate in rows[:top]:
        local, server, manifest = candidate_inputs(candidate, target, scaffold, flank, calcium_count)
        folder = outdir / candidate["candidate_id"]
        folder.mkdir()
        write_json(folder / "local_input.json", local)
        write_json(folder / "server_input.json", [server])
        write_json(folder / "chain_mapping.json", manifest)
        combined.append(server)
        exported.append({"candidate_id": candidate["candidate_id"],
                         "scope": manifest["scope"], "directory": folder.name})
    write_json(outdir / "server_jobs.json", combined)
    readme = ["# AlphaFold 3: Input Export and Structure Viewing", "",
        "This directory contains prediction requests only. The software has not signed in to or submitted jobs to AlphaFold Server, downloaded model weights, or generated predicted coordinates.", "",
        "1. Open https://alphafoldserver.com in your browser, use your own account, and review the service and output terms.",
        "2. Upload a candidate's server_input.json and check each DNA/RNA entity. local_input.json uses the separate local AlphaFold 3 format; the two formats are not interchangeable.",
        "3. After prediction, download the results ZIP. Coordinates are commonly supplied as .cif/mmCIF, not necessarily .pdb.",
        "4. Open viewer/structure.html in the analysis results and load a .pdb or .cif file. Optionally load chain_mapping.json to highlight the RNA cut site.",
        "5. Inspect chain sequences and quality information such as pLDDT/PAE. Low-confidence structures are descriptive and do not establish catalysis, metal coordination, or in-cell accessibility.", "",
        "Arm-only candidates have no catalytic core and are not complete MNAzymes. Complete partzymes and an available trigger are exported only when the user selects a fixed scaffold.",
        "An ion count is not an experimental concentration. No Ca ions are added by default; --af3-calcium-count explicitly sets an exploratory modeled ion count.",
        "A short RNA context does not represent full-length mRNA folding. An external prediction service requires an internet connection, and its output terms are separate from this project's MIT source license.", "",
        "This software does not automatically upload sequences or structures or incur cloud-computing charges."]
    (outdir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    return {"status": "inputs_exported_NOT_predicted", "prediction_executed": False,
            "submitted_to_any_service": False, "candidate_count": len(exported),
            "RNA_context_flank_nt": flank, "modeled_Ca_ion_count": calcium_count,
            "candidates": exported, "directory": outdir.name}
