# Source guide (2026-10-20; documentation comments only).
# Purpose: Call the ViennaRNA local-folding backend
# Inputs: RNA sequence, temperature, window/span and backend
# Generates/returns: Local unpaired probabilities and folding diagnostics
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Pinned ViennaRNA Python local-window backend for the portable desktop app.

This calls the same vrna_probs_window recursions as RNAplfold, not RNAfold MFE
or global probabilities. Callback probabilities are indexed (end_1based, length).

Verified against ViennaRNA 2.7.2 primary source on 2026-10-06:
* src/bin/RNAplfold.c, plfold_callback()/print_up(): callback i is lunp row i.
* src/ViennaRNA/partfunc/pf_window.c: pU[k][startu] corresponds to
  [k-startu+1..k], and return_pU passes k without an index transformation.
* interfaces/pf_window.i: arrays expose a dummy index 0, real values start at 1.

The generic documentation's "5' i"/pfl_fold_up start-position prose is not
used to shift this callback: source and real BIRC5 RNAplfold regression resolve
that ambiguity. Optional CLI-vs-API tests can run when RNAplfold is available.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import math
from pathlib import Path
import threading
import time
from typing import Callable

from .accessibility import parse_lunp

TESTED_VIENNARNA_VERSION = "2.7.2"
PROBABILITY_ROUNDOFF_TOLERANCE = 1e-12
_FOLD_LOCK = threading.Lock()
PRIMARY_SOURCES = [
    "https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/src/bin/RNAplfold.c",
    "https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/src/ViennaRNA/partfunc/pf_window.c",
    "https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/interfaces/pf_window.i",
    "https://viennarna.readthedocs.io/en/latest/partfunc/window.html",
]


def _get_rna_module():
    try:
        rna = importlib.import_module("RNA")
    except ImportError as exc:
        raise RuntimeError(
            "ViennaRNA Python backend is unavailable. Use the bundled desktop "
            "runtime, install ViennaRNA==2.7.2 in a separate environment, or "
            "select the external RNAplfold backend.") from exc
    version = getattr(rna, "__version__", None)
    if version != TESTED_VIENNARNA_VERSION:
        raise RuntimeError(
            f"Python backend was validated for ViennaRNA {TESTED_VIENNARNA_VERSION}; "
            f"detected {version!r}. Use the pinned runtime or external RNAplfold.")
    for name in ("md", "fold_compound", "OPTION_WINDOW", "PROBS_WINDOW_UP", "ANY_LOOP",
                 "params_load_RNA_Turner2004", "parameter_set_rna_turner2004"):
        if not hasattr(rna, name):
            raise RuntimeError(f"ViennaRNA Python runtime lacks required API {name}")
    return rna


def _validated_input(sequence: str, W: int, L: int, u: int,
                     T: float, salt_molar: float | None) -> tuple[str, dict, dict]:
    seq = "".join(str(sequence).upper().split()).replace("T", "U")
    if len(seq) < 4 or set(seq) - set("ACGU"):
        raise ValueError("Local Python folding requires at least 4 unmodified A/C/G/U bases (T is normalized to U)")
    if not all(isinstance(value, int) and not isinstance(value, bool) and value > 0
               for value in (W, L, u)):
        raise ValueError("W, L and u must be positive integers")
    if W < 4 or L > W or u > W:
        raise ValueError("Local folding requires W >= 4, L <= W and u <= W")
    if not math.isfinite(float(T)) or not 0 <= float(T) <= 100:
        raise ValueError("This validated RNA model backend supports finite temperatures from 0 to 100 degrees C")
    if salt_molar is not None and (not math.isfinite(float(salt_molar)) or salt_molar <= 0):
        raise ValueError("RNA-model monovalent salt must be finite and positive (M)")
    requested = {"W": W, "L": L, "u": u, "T": float(T), "salt_molar": salt_molar}
    # This is the short-sequence clamping in RNAplfold.c, explicitly audited.
    effective_W = min(W, len(seq))
    effective = {"W": effective_W, "L": min(L, effective_W),
                 "u": min(u, effective_W), "T": float(T),
                 "salt_molar": float(salt_molar) if salt_molar is not None else 1.021}
    return seq, requested, effective


def _collect_probabilities(rna, sequence: str, effective: dict,
                           progress: Callable[[dict], None] | None = None) -> tuple[dict, dict, dict]:
    """Collect complete all-loop UP callbacks and audit effective model fields.

    Callback exceptions are captured for later propagation: SWIG callbacks can
    otherwise print/swallow Python errors while native recursions continue.
    A module lock serializes native parameter loading and computation.
    """
    probability_rows, errors = {}, []
    numerics = {"roundoff_clipped_probability_count": 0,
                "maximum_roundoff_clip_magnitude": 0.0,
                "roundoff_tolerance": PROBABILITY_ROUNDOFF_TOLERANCE,
                "nonphysical_boundary_values_discarded": 0}
    with _FOLD_LOCK:
        if not rna.params_load_RNA_Turner2004():
            raise RuntimeError("ViennaRNA failed to load built-in RNA Turner2004 parameters")
        md = rna.md()
        defaults = {
            "temperature": effective["T"], "window_size": effective["W"],
            "max_bp_span": effective["L"], "salt": effective["salt_molar"],
            "dangles": 2, "betaScale": 1.0, "pf_smooth": 1,
            "special_hp": 1, "noLP": 0, "noGU": 0, "noGUclosure": 0,
            "logML": 0, "circ": 0, "gquad": 0, "uniq_ML": 0,
            "energy_set": 0, "min_loop_size": 3, "compute_bpp": 1,
            "sfact": 1.07, "saltMLLower": 6, "saltMLUpper": 24,
            "saltDPXInit": 99999, "saltDPXInitFact": -45.324,
            "helical_rise": 2.8, "backbone_length": 6.0,
        }
        for field, value in defaults.items():
            if not hasattr(md, field):
                raise RuntimeError(f"Installed ViennaRNA cannot explicitly set model parameter {field}")
            setattr(md, field, value)
        model_details = {field: getattr(md, field) for field in defaults}
        fc = rna.fold_compound(sequence, md, rna.OPTION_WINDOW)
        if fc is None or not hasattr(fc, "probs_window"):
            raise RuntimeError("ViennaRNA could not construct a local-window fold compound")

        def callback(pr, size, end, maximum, kind, data):
            if errors:
                return
            if not (kind & rna.PROBS_WINDOW_UP):
                return
            if kind & rna.ANY_LOOP != rna.ANY_LOOP:
                errors.append("Unexpected loop-split probability callback; all-loop UP probabilities were requested")
                return
            try:
                if not 1 <= end <= len(sequence) or end in probability_rows:
                    raise ValueError("Callback contains duplicate or out-of-range end position")
                physical_size = min(end, effective["u"])
                # pf_window.c has a documented n < min_loop_size+2 shortcut:
                # it reports full ulength arrays even before the 5' boundary.
                # Preserve only actual [end-length+1,end] transcript segments.
                expected_native_size = effective["u"] if len(sequence) < 5 else physical_size
                if maximum != effective["u"] or size != expected_native_size:
                    raise ValueError("Callback probability dimensions do not match end-indexed UP regions")
                if len(pr) <= size:
                    raise ValueError("Callback probability array omits the required 1-based indices")
                row = []
                numerics["nonphysical_boundary_values_discarded"] += size - physical_size
                for length in range(1, physical_size + 1):
                    if pr[length] is None:
                        raise ValueError("A physically valid unpaired segment has missing probability")
                    value = float(pr[length])
                    if not math.isfinite(value) or value < -PROBABILITY_ROUNDOFF_TOLERANCE or value > 1 + PROBABILITY_ROUNDOFF_TOLERANCE:
                        raise ValueError("ViennaRNA returned a non-finite or out-of-range probability")
                    # Native summation can yield 1+7e-16 for a fully unpaired
                    # poly-A segment. Only bounded floating-point roundoff is
                    # projected onto [0,1]; each occurrence is explicitly audited.
                    if not 0 <= value <= 1:
                        bounded = min(1.0, max(0.0, value))
                        numerics["roundoff_clipped_probability_count"] += 1
                        numerics["maximum_roundoff_clip_magnitude"] = max(
                            numerics["maximum_roundoff_clip_magnitude"], abs(value - bounded))
                        value = bounded
                    row.append(value)
                probability_rows[end] = row
                if progress is not None and (end % 100 == 0 or end == len(sequence)):
                    progress({"phase": "RNA_LOCAL_WINDOW_PF", "positions_completed": len(probability_rows),
                              "positions_total": len(sequence)})
            except Exception as exc:
                errors.append(str(exc))

        succeeded = fc.probs_window(effective["u"], rna.PROBS_WINDOW_UP, callback, None)
    if errors:
        raise RuntimeError("ViennaRNA probability callback failed: " + errors[0])
    if not succeeded:
        raise RuntimeError("ViennaRNA local partition function returned failure")
    if set(probability_rows) != set(range(1, len(sequence) + 1)):
        raise RuntimeError("ViennaRNA did not return probabilities for the complete transcript")
    return probability_rows, model_details, numerics


def run_vienna_python(sequence: str, outdir: str | Path, W: int = 200,
                      L: int = 150, u: int = 16, T: float = 37,
                      salt_molar: float | None = None, *,
                      progress: Callable[[dict], None] | None = None) -> tuple[Path, dict]:
    """Compute and write standard lunp probabilities without a CLI executable.

    Local callbacks are averaged over full windows containing the segment,
    exactly as RNAplfold. No all-loop contexts are dropped, no MFE is used,
    and no dummy/wrapped/unknown segment is invented. No output is overwritten.
    The native backend supports only its regression-validated pinned version.
    """
    seq, requested, effective = _validated_input(sequence, W, L, u, T, salt_molar)
    rna = _get_rna_module()
    out = Path(outdir).resolve()
    seq_hash = hashlib.sha256(seq.encode("ascii")).hexdigest()
    identifier = "target_" + seq_hash[:16]
    lunp = out / (identifier + "_lunp")
    reserved = (lunp, out / "input.fasta", out / "lunp_manifest.json")
    if any(path.exists() for path in reserved):
        raise FileExistsError("Refusing to overwrite an existing Python folding output or manifest; choose a new folder")
    started = time.monotonic()
    rows, model_details, numerics = _collect_probabilities(rna, seq, effective, progress)
    elapsed = time.monotonic() - started
    out.mkdir(parents=True, exist_ok=True)
    lines = ["#unpaired probabilities", "#i\tl=" + "\t".join(str(i) for i in range(1, effective["u"] + 1))]
    for end in range(1, len(seq) + 1):
        values = [format(value, ".17g") for value in rows[end]]
        values += ["NA"] * (effective["u"] - len(values))
        lines.append(str(end) + "\t" + "\t".join(values))
    with lunp.open("x", encoding="ascii", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    with (out / "input.fasta").open("x", encoding="ascii", newline="\n") as handle:
        handle.write(f">{identifier}\n{seq}\n")
    written = parse_lunp(lunp)
    if any(written[(end, length)] != value for end, row in rows.items()
           for length, value in enumerate(row, 1)):
        raise RuntimeError("Written lunp probability round-trip does not reproduce the native calculation")
    parameter_text = rna.parameter_set_rna_turner2004
    if not isinstance(parameter_text, str):
        raise RuntimeError("Pinned RNA Turner2004 parameter text is not available for hashing")
    parameter_sha = hashlib.sha256(parameter_text.encode("utf-8")).hexdigest()
    audit = {
        "software": "ViennaRNA Python API", "backend": "viennarna_python_local_window",
        "version": rna.__version__, "validated_version": TESTED_VIENNARNA_VERSION,
        "native_function": "RNA.fold_compound.probs_window / vrna_probs_window",
        "sequence_sha256_rna": seq_hash, "sequence_length_nt": len(seq),
        "input_fasta": "input.fasta", "lunp_file": lunp.name,
        "lunp_sha256": hashlib.sha256(lunp.read_bytes()).hexdigest(),
        "parameters": effective, "requested_parameters": requested,
        "model_details": model_details, "energy_parameter_set": "RNA_Turner2004",
        "energy_parameters_sha256": parameter_sha,
        "model_configuration_sha256": hashlib.sha256(json.dumps(
            {"version": rna.__version__, "parameters": effective,
             "model_details": model_details, "energy_parameters_sha256": parameter_sha},
            sort_keys=True, separators=(",", ":"), allow_nan=False).encode("ascii")).hexdigest(),
        "salt_status": "explicit_monovalent_RNA_model_only" if salt_molar is not None else "pinned_default_1.021M_not_project_buffer",
        "probability_indexing": "row=end_1based; column=unpaired_length; segment=[end-length+1,end]",
        "callback_context": "PROBS_WINDOW_UP | ANY_LOOP; no UP_SPLIT",
        "probability_print_precision": "17 significant digits; RNAplfold CLI usually prints 7",
        "n_probability_rows": len(rows),
        "n_physical_segments": sum(len(row) for row in rows.values()),
        "numerical_audit": numerics,
        "elapsed_seconds": elapsed,
        "cli_comparison_status": "same native recursion and indexing; see separate regression tests, no CLI invoked in this run",
        "primary_source_urls": PRIMARY_SOURCES,
        "warnings": [
            "RNA-only local secondary-structure ensemble, not MNAzyme binding or cellular exposure.",
            "Salt is a monovalent RNA-model parameter, not Ca2+/Tris or mixed-ion experimental buffer.",
            "Native local probability is not a cleavage efficiency, kinetic barrier, kcat, or off-target cleavage probability.",
            "CLI and API numerical comparisons require the same version, RNA parameters, temperature, window, span and salt.",
        ],
    }
    if any(requested[key] != effective[key] for key in ("W", "L", "u")):
        audit["warnings"].append("Short-sequence W/L/u clamping matched RNAplfold; requested and effective parameters are separately recorded.")
    if numerics["roundoff_clipped_probability_count"]:
        audit["warnings"].append("Native probabilities outside [0,1] by <=1e-12 floating-point roundoff were bounded; count and magnitude are audited.")
    with (out / "lunp_manifest.json").open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    return lunp, audit
