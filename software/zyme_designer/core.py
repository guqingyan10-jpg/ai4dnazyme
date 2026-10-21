# Source guide (2026-10-20; documentation comments only).
# Purpose: Enumerate windows, build recognition arms and calculate ranking
# Inputs: Target RNA sequence and design/scoring settings
# Generates/returns: Candidate objects and score components returned to CLI/UI
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Deterministic sequence enumeration and interpretable, explicitly heuristic ranking."""
from __future__ import annotations

import hashlib
import math
import re
from collections import Counter

RNA_BASES = frozenset("ACGU")
IUPAC = frozenset("ACGURYSWKMBDHVN")
IUPAC_RNA_SETS = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"), "U": frozenset("U"),
    "R": frozenset("AG"), "Y": frozenset("CU"), "S": frozenset("CG"), "W": frozenset("AU"),
    "K": frozenset("GU"), "M": frozenset("AC"), "B": frozenset("CGU"),
    "D": frozenset("AGU"), "H": frozenset("ACU"), "V": frozenset("ACG"), "N": RNA_BASES,
}
DEFAULT_WEIGHTS = {"accessibility": .35, "binding": .20, "balance": .10,
                   "specificity": .30, "sequence_risk": .05}
DEFAULT_LENGTH_POLICY = {
    "preferred_min": 7, "preferred_max": 9,
    "practical_min": 6, "practical_max": 12,
    "preferred_scale_nt": 3., "practical_scale_nt": 2.,
    "strong_binding_threshold_kcal_mol": -14.,
    "strong_binding_scale_kcal_mol": 6.,
}
# Existing project template; do not silently treat this split core as universal.
BIRC5_SCAFFOLD = {
    "name": "existing_BIRC5_Ca_adapted_8_17_mir21_project_template",
    "left_fixed_prefix": "TCAACATCAGTGAGACGAA",
    "right_fixed_suffix": "TCTCAGCCTGATAAGCTA",
    "trigger_rna": "UAGCUUAUCAGACUGAUGUUGA",
    "left_trigger_length": 11, "right_trigger_length": 11,
    "source": "Author-provided project sequences; core/trigger retained, substrate arms replaced",
    "catalytic_activity_preserved_after_arm_change": "not established",
    "core_motif_pattern": "AG",
    "core_motif_policy": "project_specific_8_17_AG_conservative_policy_not_universal_core",
}


def normalize_rna(sequence, allow_ambiguity=False):
    value = "".join(str(sequence).split()).upper().replace("T", "U")
    valid = IUPAC if allow_ambiguity else RNA_BASES
    if not value or set(value) - valid:
        raise ValueError("Sequence must contain A/C/G/T/U; ambiguous IUPAC accepted only in reference RNA")
    return value


def rc_dna(sequence):
    return normalize_rna(sequence).translate(str.maketrans("ACGU", "TGCA"))[::-1]


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def normalize_motif_pattern(pattern):
    """Normalize a two-nucleotide IUPAC recognition policy; T means RNA U."""
    value = "".join(str(pattern).split()).upper().replace("T", "U")
    if len(value) != 2 or set(value) - IUPAC:
        raise ValueError("motif_pattern must be exactly two RNA IUPAC symbols, e.g. AG, RY, or NN")
    return value


def motif_matches(motif, pattern):
    """Policy compatibility only; not a biochemical cleavage specificity model."""
    pattern = normalize_motif_pattern(pattern)
    motif = str(motif).upper().replace("T", "U")
    return len(motif) == 2 and set(motif) <= RNA_BASES and all(
        base in IUPAC_RNA_SETS[code] for base, code in zip(motif, pattern))


def annotate_header(header):
    accession_match = re.search(r"(?:NM|NR|XM|XR)_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?", header)
    accession = accession_match.group(0) if accession_match else header.split()[0]
    # Ensembl often puts gene:ENSG... before gene_symbol:...; prefer the symbol.
    gene_match = re.search(r"(?:gene_symbol|gene_name)[:=]([\w.-]+)", header)
    if not gene_match:
        gene_match = re.search(r"\bgene[:=]([\w.-]+)", header)
    names = re.findall(r"\(([A-Za-z0-9_.-]+)\)(?=,|$)", header)
    gene = gene_match.group(1) if gene_match else (names[0] if names else "")
    fields = header.split()[0].split("|")
    if not gene and len(fields) >= 6 and fields[0].startswith("ENST"):
        gene = fields[5]  # GENCODE transcript|gene|...|transcript_name|gene_name|...
    if accession.split(".")[0] == "NM_001168" and not gene:
        gene = "BIRC5"
    return accession, gene


def validate_scaffold(scaffold):
    for key in ("left_fixed_prefix", "right_fixed_suffix"):
        if not scaffold.get(key) or set(scaffold[key].upper()) - set("ACGT"):
            raise ValueError(f"Scaffold {key} must be an unmodified DNA sequence")
    result = dict(scaffold)
    result["left_fixed_prefix"] = result["left_fixed_prefix"].upper()
    result["right_fixed_suffix"] = result["right_fixed_suffix"].upper()
    if result.get("trigger_rna"):
        trigger = normalize_rna(result["trigger_rna"])
        la, rb = int(result.get("left_trigger_length", 0)), int(result.get("right_trigger_length", 0))
        if la <= 0 or rb <= 0 or la + rb != len(trigger):
            raise ValueError("Scaffold requires positive trigger arm lengths jointly covering the trigger")
        if result["left_fixed_prefix"][:la] != rc_dna(trigger[-la:]):
            raise ValueError("Left trigger segment is not reverse-complementary to the trigger's right segment")
        if result["right_fixed_suffix"][-rb:] != rc_dna(trigger[:rb]):
            raise ValueError("Right trigger segment is not reverse-complementary to the trigger's left segment")
        result["trigger_rna"] = trigger
        result["trigger_complementarity_checked"] = True
    else:
        result["trigger_complementarity_checked"] = False
    return result


def enumerate_candidates(sequence, left_lengths=(7,), right_lengths=(7,), scaffold=None,
                         *, motif_pattern="AG", cut_positions=None):
    sequence = normalize_rna(sequence, allow_ambiguity=True)
    motif_pattern = normalize_motif_pattern(motif_pattern)
    scaffold = validate_scaffold(scaffold) if scaffold else None
    if scaffold and scaffold.get("core_motif_pattern"):
        core_pattern = normalize_motif_pattern(scaffold["core_motif_pattern"])
        if core_pattern != motif_pattern:
            raise ValueError("Scaffold has a different declared core motif policy; generic recognition does not make this core universal")
    rows, skipped = [], 0
    left_lengths, right_lengths = tuple(left_lengths), tuple(right_lengths)
    if not left_lengths or not right_lengths:
        raise ValueError("At least one arm length per side is required")
    for length in left_lengths + right_lengths:
        if not isinstance(length, int) or isinstance(length, bool) or not 4 <= length <= 20:
            raise ValueError("Arm lengths supported: 4–20 nt; outside this range require a revised model")
    compatible_positions = [i + 1 for i in range(len(sequence) - 1)
                            if motif_matches(sequence[i:i + 2], motif_pattern)]
    if cut_positions is None:
        selected_positions = compatible_positions
    else:
        requested = list(cut_positions)
        if not requested:
            raise ValueError("Explicit cut_positions must not be empty")
        if any(not isinstance(pos, int) or isinstance(pos, bool) or not 1 <= pos < len(sequence) for pos in requested):
            raise ValueError("cut_positions must be 1-based integer positions before the final transcript nucleotide")
        selected_positions = sorted(set(requested))
        if any(pos not in compatible_positions for pos in selected_positions):
            raise ValueError("Explicit cut position does not match the declared canonical two-base motif policy")
    for position in selected_positions:
        center = position - 1
        actual_motif = sequence[center:center + 2]
        if not motif_matches(actual_motif, motif_pattern):
            continue
        for left in sorted(set(left_lengths)):
            for right in sorted(set(right_lengths)):
                start, end = center - left, center + 2 + right
                if start < 0 or end > len(sequence):
                    continue
                target = sequence[start:end]
                if set(target) - RNA_BASES:
                    skipped += 1
                    continue
                left_arm, right_arm = rc_dna(target[:left]), rc_dna(target[-right:])
                row = {"candidate_id": f"{actual_motif}{center + 1}_L{left}_R{right}",
                       "start_1based": start + 1, "end_1based": end,
                       "cut_after_1based": center + 1, "motif": actual_motif,
                       "motif_pattern": motif_pattern,
                       "motif_compatibility_is_cleavage_proof": False,
                       "target_rna": target, "left_length": left, "right_length": right,
                       "left_arm_dna": left_arm, "right_arm_dna": right_arm,
                       "left_partzyme": scaffold["left_fixed_prefix"] + left_arm if scaffold else None,
                       "right_partzyme": right_arm + scaffold["right_fixed_suffix"] if scaffold else None,
                       "scaffold_status": "project_template_retained" if scaffold else "substrate_arms_only",
                       "scaffold_core_motif_policy": scaffold.get("core_motif_pattern", "user_declared_unverified") if scaffold else "no_core_provided",
                       "arm_complementarity_checked": rc_dna(left_arm) == target[:left].replace("U", "T")
                           and rc_dna(right_arm) == target[-right:].replace("U", "T"),
                       "p_unpaired_target": None, "p_unpaired_left_arm": None,
                       "p_unpaired_right_arm": None, "p_unpaired_center": None,
                       "opening_dg_kcal_mol": None, "off_target_status": "not_screened",
                       "structure_score_used": False}
                # Check coordinates, lengths, and exclusion of the central junction.
                assert motif_matches(sequence[row["cut_after_1based"] - 1:row["cut_after_1based"] + 1], motif_pattern)
                assert len(target) == left + right + 2 and row["arm_complementarity_checked"]
                rows.append(row)
    return rows, {"AG_junctions_in_target": sequence.count("AG"),
                  "motif_pattern": motif_pattern, "compatible_junctions_in_target": len(compatible_positions),
                  "cut_positions_policy": "all_compatible" if cut_positions is None else "explicit_user_selected",
                  "selected_cut_positions_1based": selected_positions,
                  "candidate_designs": len(rows), "ambiguous_candidate_windows_skipped": skipped,
                  "cleavage_motif_policy": "AG-only conservative policy for selected 8–17 core, not exhaustive specificity" if motif_pattern == "AG"
                    else "User-declared IUPAC recognition policy; no catalytic core activity or cleavage specificity inferred"}


def _finite(value):
    return value is not None and not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def normalize_length_policy(policy=None):
    """Project-adjustable heuristic preferences, not universal biological limits."""
    policy = {} if policy is None else policy
    if not isinstance(policy, dict) or set(policy) - set(DEFAULT_LENGTH_POLICY):
        raise ValueError("Unknown length-policy keys; use documented preference and penalty parameters")
    result = {**DEFAULT_LENGTH_POLICY, **policy}
    for key in ("preferred_min", "preferred_max", "practical_min", "practical_max"):
        if not isinstance(result[key], int) or isinstance(result[key], bool):
            raise ValueError("Length-policy bounds must be integer nucleotide counts")
    if not 4 <= result["practical_min"] <= result["preferred_min"] <= result["preferred_max"] <= result["practical_max"] <= 20:
        raise ValueError("Require 4 <= practical_min <= preferred_min <= preferred_max <= practical_max <= 20")
    for key in ("preferred_scale_nt", "practical_scale_nt", "strong_binding_scale_kcal_mol"):
        if not _finite(result[key]) or result[key] <= 0:
            raise ValueError("Length/binding penalty scales must be positive finite numbers")
    if not _finite(result["strong_binding_threshold_kcal_mol"]):
        raise ValueError("Strong-binding threshold must be finite")
    return result


def _length_tradeoff(row, policy):
    lengths = [row["left_length"], row["right_length"]]
    if any(not isinstance(length, int) or isinstance(length, bool) or not 4 <= length <= 20 for length in lengths):
        raise ValueError("Ranking requires arm lengths in the supported 4–20 nt interval")
    preferred = [max(policy["preferred_min"] - n, 0, n - policy["preferred_max"]) for n in lengths]
    practical = [max(policy["practical_min"] - n, 0, n - policy["practical_max"]) for n in lengths]
    preferred_cost = sum(distance ** 2 for distance in preferred) / (2 * policy["preferred_scale_nt"] ** 2)
    practical_cost = sum(distance ** 2 for distance in practical) / (2 * policy["practical_scale_nt"] ** 2)
    energies = [row.get("left_binding_dg_kcal_mol"), row.get("right_binding_dg_kcal_mol")]
    binding_known = all(_finite(energy) for energy in energies)
    strong_cost = None if not binding_known else sum(max(0., policy["strong_binding_threshold_kcal_mol"] - energy) ** 2
                       for energy in energies) / (2 * policy["strong_binding_scale_kcal_mol"] ** 2)
    reasons = []
    if any(preferred):
        reasons.append("arm_length_outside_user_preferred_range")
    if any(practical):
        reasons.append("arm_length_outside_user_practical_range")
    if strong_cost is not None and strong_cost > 0:
        reasons.append("very_strong_1M_reference_duplex_stability_penalized_not_a_release_rate")
    if not binding_known:
        reasons.append("strong_binding_penalty_not_assessed_missing_reference_energy")
    if not reasons:
        reasons.append("within_user_preferences_no_excess_reference_binding_penalty")
    total_cost = preferred_cost + practical_cost + (strong_cost if strong_cost is not None else 0.)
    return {"length_preference_penalty_cost": preferred_cost,
            "outside_practical_range_penalty_cost": practical_cost,
            "length_penalty": 1. - math.exp(-preferred_cost - practical_cost),
            "strong_binding_penalty_cost": strong_cost,
            "strong_binding_penalty": None if strong_cost is None else 1. - math.exp(-strong_cost),
            "length_tradeoff_score": math.exp(-total_cost), "length_tradeoff_reasons": reasons,
            "length_tradeoff_is_release_probability": False}


def score_candidates(rows, weights=None, length_policy=None):
    length_policy = normalize_length_policy(length_policy)
    weights = DEFAULT_WEIGHTS if weights is None else weights
    if set(weights) != set(DEFAULT_WEIGHTS) or any(not _finite(v) or v < 0 for v in weights.values()):
        raise ValueError("Provide all five nonnegative finite heuristic score weights")
    if sum(weights.values()) <= 0:
        raise ValueError("At least one ranking weight must be positive")
    for row in rows:
        components = {}
        if _finite(row.get("p_unpaired_target")):
            components["accessibility"] = row["p_unpaired_target"]
        a, b = row.get("left_binding_dg_kcal_mol"), row.get("right_binding_dg_kcal_mol")
        if _finite(a) and _finite(b):
            # Deliberately saturating heuristic, not an inferred occupancy or Kd.
            components["binding"] = min(1., max(0., -(a + b) / 20.))
            components["balance"] = math.exp(-abs(a - b) / 3.)
            if _finite(row.get("opening_dg_kcal_mol")):
                row["binding_plus_opening_proxy_kcal_mol"] = a + b + row["opening_dg_kcal_mol"]
            else:
                row["binding_plus_opening_proxy_kcal_mol"] = None
            row["strong_arm_release_caution"] = min(a, b) < length_policy["strong_binding_threshold_kcal_mol"]
        risk_run = max(row.get("max_self_complementary_run", 0), row.get("max_cross_complementary_run", 0))
        components["sequence_risk"] = math.exp(-max(risk_run - 4, 0) / 3.)
        if row.get("off_target_status") == "screened":
            risk = (10. * row.get("other_gene_exact_compatible", row.get("other_gene_exact_ag", 0))
                    + 2. * row.get("other_gene_1mm_compatible", row.get("other_gene_1mm_ag", 0))
                    + .25 * row.get("other_gene_2mm_compatible", row.get("other_gene_2mm_ag", 0))
                    + .05 * row.get("other_gene_noncompatible", row.get("other_gene_non_ag", 0))
                    + .10 * row.get("unknown_gene_hits", 0)
                    + max(0, row.get("intended_target_hits", 1) - 1))
            components["specificity"] = 1. / (1. + risk)
        coverage = sum(weights[name] for name in components)
        row["ranking_score_before_length_tradeoff"] = 100. * sum(weights[name] * val for name, val in components.items()) / coverage if coverage else 0.
        row.update(_length_tradeoff(row, length_policy))
        row["ranking_score"] = row["ranking_score_before_length_tradeoff"] * row["length_tradeoff_score"]
        row["score_components"] = components
        row["score_weight_coverage"] = coverage / sum(weights.values())
        row["ranking_stage"] = ("screened_with_accessibility" if "specificity" in components and "accessibility" in components
                                else "screened_no_accessibility" if "specificity" in components
                                else "accessibility_only_not_screened" if "accessibility" in components
                                else "sequence_only_not_screened")
        row["ranking_is_cleavage_probability"] = False
        row["pareto_front"] = None
    # Pareto comparison only in the fully scored group; missing is never zero.
    eligible = [r for r in rows if set(r["score_components"]) == set(DEFAULT_WEIGHTS)]
    for row in eligible:
        vector = row["score_components"]
        row["pareto_front"] = not any(
            all(other["score_components"][k] >= v for k, v in vector.items())
            and other["length_tradeoff_score"] >= row["length_tradeoff_score"]
            and (any(other["score_components"][k] > v for k, v in vector.items())
                 or other["length_tradeoff_score"] > row["length_tradeoff_score"])
            for other in eligible if other is not row)
    tiers = {"screened_with_accessibility": 0, "screened_no_accessibility": 1,
             "accessibility_only_not_screened": 2, "sequence_only_not_screened": 3}
    rows.sort(key=lambda row: (tiers[row["ranking_stage"]], -row["ranking_score"], row["cut_after_1based"],
                               row["left_length"], row["right_length"]))
    for rank, row in enumerate(rows, 1):
        row["display_rank"] = rank
    return rows


def ranking_audit(rows, weights=None, length_policy=None):
    length_policy = normalize_length_policy(length_policy)
    return {"method": "deterministic heuristic weighted score + fully-scored Pareto nondominance",
            "weights": dict(DEFAULT_WEIGHTS if weights is None else weights),
            "weights_fitted_to_project_experiments": False,
            "normalization": "100 * sum available weighted components / sum available weights * length_tradeoff_score; stages displayed separately",
            "length_policy": length_policy,
            "length_policy_status": "configurable project heuristic, not universal DNAzyme biology or a measured release model",
            "equations": {"accessibility": "p_unpaired_full_window", "binding": "clip(-(dG_left+dG_right)/20,0,1)",
                          "balance": "exp(-abs(dG_left-dG_right)/3)",
                          "sequence_risk": "exp(-max(max_complementary_run-4,0)/3)",
                          "specificity": "1/(1+10*exact_other_compatible+2*1mm_other_compatible+0.25*2mm_other_compatible+0.05*other_noncompatible+0.10*unknown_gene_hits+extra_target_sites); legacy AG fields only when compatible fields absent",
                          "length_tradeoff_score": "exp[-sum(d_preferred_arm^2)/(2*preferred_scale_nt^2)-sum(d_practical_arm^2)/(2*practical_scale_nt^2)-sum(max(0,threshold-dG_arm)^2)/(2*strong_binding_scale^2)]",
                          "length_distance": "distance below range_min or above range_max, otherwise zero; apply separately to BOTH arms",
                          "strong_binding_missing": "not assessed when energies missing; preserve None and reduce no unknown energy to a measured zero"},
            "stage_counts": dict(Counter(row["ranking_stage"] for row in rows)),
            "scores_are_not_probabilities_or_rates": True,
            "ranking_limit": "best among the evaluated candidate set and declared reference scope, not proven global optimum"}
