# Source guide (2026-10-20; documentation comments only).
# Purpose: Parse/map molecular structures and prepare viewer metadata
# Inputs: Coordinate files, candidate and chain/residue mapping
# Generates/returns: Structure diagnostics and viewer/export records
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Optional PDB diagnostic, not a validated 3D accessibility predictor.

The first version intentionally leaves 3D measurements out of candidate ranking.
PDB numbering must be an explicitly mapped transcript or a continuous chain
numbered 1..N. B-factors are interpreted as pLDDT only on an explicit declaration.
"""
from __future__ import annotations

import hashlib
import math
from itertools import product
from pathlib import Path

RNA_RESIDUES = {"A": "A", "C": "C", "G": "G", "U": "U",
                "RA": "A", "RC": "C", "RG": "G", "RU": "U",
                "ADE": "A", "CYT": "C", "GUA": "G", "URA": "U"}


def _parse_pdb(path: Path) -> dict:
    residues = {}
    first_model = None
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("MODEL "):
                if first_model is None:
                    first_model = line[10:14].strip()
                else:
                    break
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[16:17] not in {" ", "A", ""}:
                continue  # Avoid duplicate alternative conformations.
            residue_name = line[17:20].strip().upper()
            if residue_name not in RNA_RESIDUES:
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper()
            if element in {"H", "D"} or (not element and atom_name.lstrip("0123456789").startswith(("H", "D"))):
                continue
            try:
                key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip())
                xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
                bfactor = float(line[60:66]) if line[60:66].strip() else None
            except ValueError as exc:
                raise ValueError("Malformed RNA coordinate/B-factor field in PDB") from exc
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError("Non-finite coordinates in PDB")
            value = residues.setdefault(key, {"base": RNA_RESIDUES[residue_name], "atoms": [], "b_factors": []})
            value["atoms"].append(xyz)
            if bfactor is not None:
                value["b_factors"].append(bfactor)
    if not residues:
        raise ValueError("No canonical RNA residues found; this diagnostic supports PDB, not mmCIF")
    return residues


def pdb_site_metrics(path: str | Path, start_1based: int, end_1based: int,
                     chain: str | None = None, confidence_source: str = "unknown", *,
                     residue_map: dict | None = None,
                     transcript_sequence: str | None = None,
                     neighbor_cutoff_angstrom: float = 8.0) -> dict:
    """Map a recognition region to RNA coordinates and expose limitations.

    residue_map maps (chain, residue_number, insertion_code) tuples to transcript
    1-based indices, or "chain:residue_number:insertion_code" string keys. A map
    is required for gaps, insertion codes, or fragments not starting at residue 1.
    If transcript_sequence is supplied in identity mode, the whole RNA chain
    must exactly match it, preventing fragment and isoform misassignment.
    """
    if start_1based < 1 or end_1based < start_1based:
        raise ValueError("Invalid 1-based target interval")
    if not math.isfinite(neighbor_cutoff_angstrom) or neighbor_cutoff_angstrom <= 0:
        raise ValueError("Neighbor cutoff must be positive and finite")
    source = Path(path)
    residues = _parse_pdb(source)
    chains = sorted({key[0] for key in residues})
    if chain is None:
        if len(chains) != 1:
            raise ValueError("Multiple RNA chains present; explicitly choose a chain")
        chain = chains[0]
    keys = sorted(key for key in residues if key[0] == chain)
    if not keys:
        raise ValueError(f"RNA chain {chain!r} not found")
    warnings = [
        "3D diagnostic is not included in candidate ranking and is not proof of in-cell target exposure.",
        "Neighbor counts are geometric proximity, not solvent accessibility, a binding probability or a rate.",
        "Only first PDB model and canonical RNA heavy atoms are analyzed; proteins, ions and waters are not exposure scores.",
    ]
    if residue_map is None:
        if any(key[2] for key in keys) or [key[1] for key in keys] != list(range(1, len(keys) + 1)):
            raise ValueError("RNA chain is not continuously numbered 1..N; provide an explicit residue_map")
        if transcript_sequence is not None:
            seq = "".join(transcript_sequence.upper().split()).replace("T", "U")
            observed = "".join(residues[key]["base"] for key in keys)
            if observed != seq:
                raise ValueError("PDB chain sequence does not exactly match the complete input transcript")
            mapping_status = "complete_transcript_sequence_verified_identity_1based"
        else:
            mapping_status = "continuous_1based_identity_assumption_not_sequence_verified"
            warnings.append("No transcript sequence supplied: identity numbering is assumed and must be verified by the user.")
        mapping = {key: key[1] for key in keys}
    else:
        mapping = {}
        for key in keys:
            location = residue_map.get(key, residue_map.get(f"{key[0]}:{key[1]}:{key[2]}"))
            if location is not None:
                if not isinstance(location, int) or location < 1:
                    raise ValueError("Explicit transcript indices must be positive integers")
                mapping[key] = location
        if len(set(mapping.values())) != len(mapping):
            raise ValueError("More than one PDB residue maps to the same transcript index")
        if transcript_sequence is not None:
            seq = "".join(transcript_sequence.upper().split()).replace("T", "U")
            for key, position in mapping.items():
                if position > len(seq) or seq[position - 1] != residues[key]["base"]:
                    raise ValueError("Explicit PDB mapping is inconsistent with input transcript sequence")
        mapping_status = "explicit_residue_mapping"
    by_position = {position: key for key, position in mapping.items()}
    is_plddt = confidence_source.lower() in {"plddt", "alphafold_plddt"}
    if not is_plddt:
        warnings.append("B-factors are not interpreted as pLDDT without explicit confidence_source='plddt'.")
    else:
        if any(not math.isfinite(value) or not 0 <= value <= 100 for key in keys for value in residues[key]["b_factors"]):
            raise ValueError("Declared pLDDT contains values outside [0,100]")
        warnings.append("pLDDT is local coordinate confidence, not experimental validation or direct accessibility.")
    rows = []
    cutoff_squared = neighbor_cutoff_angstrom ** 2
    # A cell grid avoids scanning every atom pair for long mRNAs. Each heavy
    # atom can only be within cutoff of its own or one of 26 adjacent cells.
    # It is exact for this distance definition, not an approximate downsample.
    grid = {}
    for other, other_res in residues.items():
        for xyz in other_res["atoms"]:
            cell = tuple(math.floor(value / neighbor_cutoff_angstrom) for value in xyz)
            grid.setdefault(cell, []).append((other, xyz))
    adjacent_offsets = tuple(product((-1, 0, 1), repeat=3))
    for position in range(start_1based, end_1based + 1):
        key = by_position.get(position)
        row = {"transcript_position_1based": position,
               "coordinate_present": key is not None,
               "chain": chain, "pdb_residue_number": key[1] if key else None,
               "base": residues[key]["base"] if key else None,
               "plddt_mean": None, "nonlocal_neighbor_residue_count": None}
        if key is not None:
            residue = residues[key]
            if is_plddt and residue["b_factors"]:
                row["plddt_mean"] = sum(residue["b_factors"]) / len(residue["b_factors"])
            neighbors = set()
            for xyz in residue["atoms"]:
                cell = tuple(math.floor(value / neighbor_cutoff_angstrom) for value in xyz)
                for offset in adjacent_offsets:
                    search_cell = tuple(value + change for value, change in zip(cell, offset))
                    for other, other_xyz in grid.get(search_cell, ()):
                        if other == key or other in neighbors:
                            continue
                        # Exclude own-chain sequence neighbors (within 2 bases).
                        if other in mapping and abs(mapping[other] - position) <= 2:
                            continue
                        if sum((a - b) ** 2 for a, b in zip(xyz, other_xyz)) <= cutoff_squared:
                            neighbors.add(other)
            row["nonlocal_neighbor_residue_count"] = len(neighbors)
        rows.append(row)
    confidence = [row["plddt_mean"] for row in rows if row["plddt_mean"] is not None]
    return {
        "structure_model": "optional_PDB_geometric_diagnostic_no_ranking_effect",
        "pdb_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "mapping_status": mapping_status, "selected_chain": chain,
        "confidence_source": "declared_pLDDT" if is_plddt else "not_interpreted",
        "target_coordinate_coverage_fraction": sum(row["coordinate_present"] for row in rows) / len(rows),
        "target_mean_plddt": sum(confidence) / len(confidence) if confidence else None,
        "neighbor_cutoff_angstrom": neighbor_cutoff_angstrom,
        "residues": rows, "warnings": warnings,
    }
