# Source guide (2026-10-20; documentation comments only).
# Purpose: Calculate reference RNA/DNA recognition-arm thermodynamics
# Inputs: Arm sequences, temperature and supported ionic settings
# Generates/returns: Thermodynamic score components; no MD trajectory
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Transparent RNA/DNA hybrid reference thermodynamics, NOT enzyme activity.

Nearest-neighbor enthalpies (kcal/mol) and entropies (cal/mol/K) are from
Sugimoto et al. 1995, DOI 10.1021/bi00035a029. Values were cross-checked against
the official Biopython R_DNA_NN1 table on 2026-10-06:
https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py
RNA strand is always written 5'->3'; the paired DNA is its perfect complement.
The original 1 M NaCl reference model is not extrapolated to catalytic Ca2+.
"""
from __future__ import annotations

import hashlib
import json
import math

# Each dinucleotide is on the RNA strand, with U represented explicitly.
# These are measured scientific parameters, not copied software algorithms.
RNA_DNA_NN = {
    "AA": (-7.8, -21.9), "AC": (-5.9, -12.3),
    "AG": (-9.1, -23.5), "AU": (-8.3, -23.9),
    "CA": (-9.0, -26.1), "CC": (-9.3, -23.2),
    "CG": (-16.3, -47.1), "CU": (-7.0, -19.7),
    "GA": (-5.5, -13.5), "GC": (-8.0, -17.1),
    "GG": (-12.8, -31.9), "GU": (-7.8, -21.6),
    "UA": (-7.8, -23.2), "UC": (-8.6, -22.9),
    "UG": (-10.4, -28.4), "UU": (-11.5, -36.4),
}
INIT_H_KCAL_MOL = 1.9
INIT_S_CAL_MOL_K = -3.9
R_CAL_MOL_K = 1.987
MODEL = "Sugimoto1995_RNA_DNA_NN_1M_NaCl_reference_no_Ca_correction"
PARAMETER_SHA256 = hashlib.sha256(json.dumps(
    {"nn": RNA_DNA_NN, "init_h": INIT_H_KCAL_MOL,
     "init_s": INIT_S_CAL_MOL_K}, sort_keys=True,
    separators=(",", ":")).encode("ascii")).hexdigest()


def _normalize(sequence: str, rna: bool = True) -> str:
    seq = "".join(str(sequence).upper().split())
    seq = seq.replace("T", "U") if rna else seq.replace("U", "T")
    alphabet = set("ACGU" if rna else "ACGT")
    if not seq or set(seq) - alphabet:
        raise ValueError("Only non-empty, unmodified A/C/G/U (or T) sequences are supported")
    return seq


def hybrid_reference(sequence: str, temperature_c: float = 37,
                     strand_concentration_uM: float = 0.1) -> dict:
    """Perfect hybrid reference at 1 M NaCl and equimolar strand concentration.

    Tm uses total strand concentration / 4, with two different strand chemistries
    (no self-complementary duplex symmetry factor). No mismatches, dangling ends,
    RNA target folding, cooperative assembly, modifications, or Ca2+ are modeled.
    """
    seq = _normalize(sequence)
    if len(seq) < 2:
        raise ValueError("RNA/DNA nearest-neighbor calculation needs at least two bases")
    if not math.isfinite(float(temperature_c)) or float(temperature_c) <= -273.15:
        raise ValueError("Temperature must be finite and above absolute zero")
    if not math.isfinite(float(strand_concentration_uM)) or strand_concentration_uM <= 0:
        raise ValueError("Reference strand concentration must be positive and finite")
    enthalpy = INIT_H_KCAL_MOL
    entropy = INIT_S_CAL_MOL_K
    for i in range(len(seq) - 1):
        dh, ds = RNA_DNA_NN[seq[i:i + 2]]
        enthalpy += dh
        entropy += ds
    kelvin = float(temperature_c) + 273.15
    dg = enthalpy - kelvin * entropy / 1000.0
    concentration_factor_M = strand_concentration_uM * 1e-6 / 2
    tm = 1000 * enthalpy / (
        entropy + R_CAL_MOL_K * math.log(concentration_factor_M)) - 273.15
    return {
        "rna_sequence_5to3": seq, "enthalpy_kcal_mol": enthalpy,
        "entropy_cal_mol_K": entropy, "binding_dg_kcal_mol": dg,
        "tm_c": tm, "reference_sodium_mM": 1000.0,
        "reference_strand_concentration_uM_each": strand_concentration_uM,
    }


def binding_features(left_rna: str, right_rna: str, temperature_c: float = 37,
                     sodium_mM: float = 300, *,
                     strand_concentration_uM: float = 0.1) -> dict:
    """Reference features for the two independently hybridizing substrate arms.

    Input sodium is recorded, NOT used in an unsupported RNA/DNA or Ca2+ salt
    correction. Energies are not MM/GBSA, binding Kd, cleavage rate, or kcat.
    """
    if not math.isfinite(float(sodium_mM)) or sodium_mM <= 0:
        raise ValueError("Input sodium concentration must be positive and finite")
    left = hybrid_reference(left_rna, temperature_c, strand_concentration_uM)
    right = hybrid_reference(right_rna, temperature_c, strand_concentration_uM)
    warnings = [
        "Perfect RNA/DNA hybrid reference only; no folding, mismatches, dangling ends or modifications.",
        "Catalytic Ca2+ and Tris are not modeled; Tm and deltaG are not predictions in the project buffer.",
        "Two arm deltaG values do not define a cooperative MNAzyme binding free energy, Kd, or cleavage rate.",
        "Tm reference assumes each separate RNA and DNA strand is 0.1 uM unless explicitly overridden; tethered partzyme arms differ.",
    ]
    if not math.isclose(sodium_mM, 1000.0):
        warnings.append(
            "Input Na+ differs from the 1 M reference; no unvalidated salt correction was applied.")
    if min(len(left["rna_sequence_5to3"]), len(right["rna_sequence_5to3"])) < 5 or max(
            len(left["rna_sequence_5to3"]), len(right["rna_sequence_5to3"])) > 13:
        warnings.append("One arm length lies outside the 5-13 nt experimental oligo range in Sugimoto1995.")
    return {
        "left_binding_dg_kcal_mol": left["binding_dg_kcal_mol"],
        "right_binding_dg_kcal_mol": right["binding_dg_kcal_mol"],
        "sum_arm_binding_dg_kcal_mol": left["binding_dg_kcal_mol"] + right["binding_dg_kcal_mol"],
        "arm_binding_imbalance_kcal_mol": abs(left["binding_dg_kcal_mol"] - right["binding_dg_kcal_mol"]),
        "left_tm_c": left["tm_c"], "right_tm_c": right["tm_c"],
        "binding_model": MODEL, "binding_parameter_sha256": PARAMETER_SHA256,
        "binding_reference_sodium_mM": 1000.0,
        "binding_input_sodium_mM": float(sodium_mM),
        "binding_temperature_c": float(temperature_c),
        "tm_reference_strand_concentration_uM_each": strand_concentration_uM,
        "binding_reference_details": {"left": left, "right": right},
        "warnings": warnings,
    }


def _reverse_complement_dna(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _longest_match(a: str, b: str) -> int:
    previous = [0] * (len(b) + 1)
    best = 0
    for ca in a:
        current = [0] * (len(b) + 1)
        for j, cb in enumerate(b, 1):
            if ca == cb:
                current[j] = previous[j - 1] + 1
                best = max(best, current[j])
        previous = current
    return best


def _hairpin_run(seq: str, minimum_loop: int = 3) -> int:
    comp = str.maketrans("ACGT", "TGCA")
    best = 0
    for i in range(len(seq)):
        for j in range(i + minimum_loop + 2, len(seq)):
            run = 0
            while i + run < j - run and j - i - 2 * run - 1 >= minimum_loop:
                if seq[i + run] != seq[j - run].translate(comp):
                    break
                run += 1
                best = max(best, run)
    return best


def complementarity_risk(left_arm_dna: str, right_arm_dna: str,
                         left_partzyme: str | None = None,
                         right_partzyme: str | None = None) -> dict:
    """Sequence-only contiguous Watson-Crick complementarity diagnostic.

    This is not secondary-structure free energy. Complete partzymes may contain
    an intentionally complementary catalytic stem; it is separately disclosed
    and not automatically interpreted as a harmful interaction.
    """
    left, right = _normalize(left_arm_dna, False), _normalize(right_arm_dna, False)
    result = {
        "left_arm_self_dimer_run_nt": _longest_match(left, _reverse_complement_dna(left)),
        "right_arm_self_dimer_run_nt": _longest_match(right, _reverse_complement_dna(right)),
        "left_arm_hairpin_run_nt": _hairpin_run(left),
        "right_arm_hairpin_run_nt": _hairpin_run(right),
        "cross_arm_complementary_run_nt": _longest_match(left, _reverse_complement_dna(right)),
        "complementarity_model": "contiguous_Watson_Crick_sequence_proxy_no_energy",
        "warnings": ["Complementarity runs are flags, not folding probabilities or empirical failure thresholds."],
    }
    result["longest_unintended_complementary_run_nt"] = max(
        value for key, value in result.items() if key.endswith("_run_nt"))
    if left_partzyme is not None and right_partzyme is not None:
        lp, rp = _normalize(left_partzyme, False), _normalize(right_partzyme, False)
        result["full_partzyme_cross_run_nt_including_expected_core"] = _longest_match(
            lp, _reverse_complement_dna(rp))
        result["warnings"].append(
            "Full-partzyme complementarity includes the deliberately paired catalytic core; no automatic rejection was applied.")
    return result
