# Source guide (2026-10-20; documentation comments only).
# Purpose: Validate/import/compute local RNA accessibility
# Inputs: Sequence, probability input and folding settings
# Generates/returns: Accessibility values and provenance returned to candidate calculation
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""RNAplfold unpaired probabilities with explicit end-position indexing.

Official format: row i, column length l means bases [i-l+1 .. i], not a
forward window beginning at i. https://viennarna.readthedocs.io/en/latest/man/RNAplfold.html
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

R_KCAL_MOL_K = 0.00198720425864083
PROBABILITY_FLOOR = 1e-300


def parse_lunp(path: str | Path) -> dict[tuple[int, int], float | None]:
    """Parse probability-only text lunp; preserve zero and unknown distinctly.

    RNAplfold -O output is opening energies, not a probability table, and is
    deliberately rejected when its header or numerical values reveal that.
    """
    table = {}
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "opening" in line.lower() and "ener" in line.lower():
                    raise ValueError("Opening-energy output is not an unpaired-probability lunp table")
                continue
            fields = line.split()
            try:
                end = int(fields[0])
            except ValueError as exc:
                raise ValueError(f"Invalid lunp position at line {line_number}") from exc
            if end < 1 or len(fields) < 2:
                raise ValueError(f"Invalid lunp row at line {line_number}")
            for length, token in enumerate(fields[1:], 1):
                key = (end, length)
                if key in table:
                    raise ValueError(f"Duplicate lunp position/length {key}")
                if token.lower() in {"na", "n/a", "nan", "-", "."}:
                    table[key] = None
                    continue
                value = float(token)
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"Invalid probability at line {line_number}, length {length}")
                if length > end:
                    # This region extends before the transcript; it cannot exist.
                    table[key] = None
                else:
                    table[key] = value
    if not table:
        raise ValueError("No unpaired probabilities found in lunp file")
    return table


def _opening(probability: float | None, kelvin: float) -> float | None:
    if probability is None:
        return None
    return -R_KCAL_MOL_K * kelvin * math.log(max(probability, PROBABILITY_FLOOR))


def region_metrics(table: dict, start_1based: int, left_length: int,
                   right_length: int, temperature_c: float = 37) -> dict:
    """Retrieve an unpaired full recognition window, two arms, and central AG.

    Window is left arm + unpaired AG motif + right arm. Opening terms are
    -RT ln averaged unpaired probability, not a kinetic barrier or cellular
    occupancy; the full-window probability is not the product of arm values.
    """
    if not all(isinstance(value, int) and not isinstance(value, bool) and value > 0
               for value in (start_1based, left_length, right_length)):
        raise ValueError("Recognition start and arm lengths must be positive")
    if not math.isfinite(float(temperature_c)) or temperature_c <= -273.15:
        raise ValueError("Temperature must be finite and above absolute zero")
    full_length = left_length + 2 + right_length
    left_end = start_1based + left_length - 1
    motif_end = left_end + 2
    right_end = motif_end + right_length
    locations = {"full_window": (right_end, full_length),
                 "left_arm": (left_end, left_length),
                 "right_arm": (right_end, right_length),
                 "ag": (motif_end, 2)}
    out = {"accessibility_model": "RNAplfold_local_RNA_only_ensemble",
           "accessibility_temperature_c": float(temperature_c),
           "accessibility_window_start_1based": start_1based,
           "accessibility_window_end_1based": right_end,
           "accessibility_window_length_nt": full_length,
           "accessibility_indexing": "(segment_end_1based,segment_length)",
           "warnings": ["RNA-only unpaired probabilities do not model partzyme binding, catalytic Ca2+, RNA-binding proteins or in-cell cleavage."]}
    missing, zero = [], []
    for region, location in locations.items():
        probability = table.get(location)
        if probability is not None and not 0 <= probability <= 1:
            raise ValueError(f"Probability for {region} is outside [0,1]")
        out[f"{region}_p_unpaired"] = probability
        out[f"{region}_opening_dg_kcal_mol"] = _opening(probability, temperature_c + 273.15)
        out[f"{region}_lunp_key"] = list(location)
        if probability is None:
            missing.append(region)
        elif probability == 0:
            zero.append(region)
    if full_length == 16:
        out["full16_p_unpaired"] = out["full_window_p_unpaired"]
    out["missing_probability_regions"] = missing
    out["probability_zero_regions"] = zero
    out["opening_probability_floor"] = PROBABILITY_FLOOR
    if missing:
        out["warnings"].append("Missing/n/a probabilities stay unknown, never interpreted as zero accessibility.")
    if zero:
        out["warnings"].append("Printed zero probability: opening term is capped at -RT ln(1e-300), not a finite measured opening energy.")
    return out


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_rnaplfold(sequence: str, outdir: str | Path, W: int = 200,
                 L: int = 150, u: int = 16, T: float = 37,
                 executable: str = "RNAplfold", salt_molar: float | None = None
                 ) -> tuple[Path, dict]:
    """Run RNAplfold safely with captured version, inputs, logs, output hashes.

    Default salt is explicitly unspecified (the detected program's default),
    not falsely labeled as project buffer. Explicit salt requires version help
    advertising --salt; older tools fail safely instead of silently ignoring it.
    Existing lunp output is never overwritten. No command is shell-evaluated.
    """
    seq = "".join(str(sequence).upper().split()).replace("T", "U")
    if not seq or set(seq) - set("ACGU"):
        raise ValueError("RNAplfold input must contain only unmodified A/C/G/U (or T)")
    if not all(isinstance(value, int) and value > 0 for value in (W, L, u)):
        raise ValueError("W, L and u must be positive integers")
    if L > W or u > W:
        raise ValueError("RNAplfold requires L <= W and requested unpaired length u <= W")
    if not math.isfinite(float(T)) or T <= -273.15:
        raise ValueError("Temperature must be finite and above absolute zero")
    if salt_molar is not None and (not math.isfinite(float(salt_molar)) or salt_molar <= 0):
        raise ValueError("Explicit salt must be positive molar concentration")
    resolved = shutil.which(str(executable))
    if resolved is None:
        raise FileNotFoundError("RNAplfold executable not found; install ViennaRNA or supply a verified lunp file")
    out = Path(outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    seq_hash = hashlib.sha256(seq.encode("ascii")).hexdigest()
    identifier = "target_" + seq_hash[:16]
    lunp = out / (identifier + "_lunp")
    if lunp.exists():
        raise FileExistsError(f"Refusing to overwrite existing RNAplfold output: {lunp}")
    version_run = subprocess.run([resolved, "--version"], text=True, capture_output=True, check=True)
    version = (version_run.stdout + version_run.stderr).strip()
    if salt_molar is not None:
        help_run = subprocess.run([resolved, "--help"], text=True, capture_output=True, check=True)
        if "--salt" not in help_run.stdout + help_run.stderr:
            raise RuntimeError("Installed RNAplfold does not advertise --salt; explicit salt cannot be applied")
    command = [resolved, "-W", str(W), "-L", str(L), "-u", str(u), "-T", str(T)]
    if salt_molar is not None:
        command += ["--salt", str(salt_molar)]
    fasta = f">{identifier}\n{seq}\n"
    (out / "input.fasta").write_text(fasta, encoding="ascii")
    run = subprocess.run(command, input=fasta, text=True, cwd=out,
                         capture_output=True, check=False)
    (out / "RNAplfold.stdout.log").write_text(run.stdout, encoding="utf-8")
    (out / "RNAplfold.stderr.log").write_text(run.stderr, encoding="utf-8")
    if run.returncode:
        raise RuntimeError(f"RNAplfold failed ({run.returncode}); inspect {out / 'RNAplfold.stderr.log'}")
    if not lunp.is_file():
        raise RuntimeError(f"RNAplfold did not create expected lunp output {lunp.name}")
    table = parse_lunp(lunp)
    if max(end for end, _ in table) != len(seq):
        raise ValueError("RNAplfold probability table length does not match input transcript")
    audit = {
        "software": "RNAplfold", "version": version,
        "sequence_sha256_rna": seq_hash, "sequence_length_nt": len(seq),
        "command": command, "input_fasta": "input.fasta", "lunp_file": lunp.name,
        "lunp_sha256": _sha256_file(lunp),
        "parameters": {"W": W, "L": L, "u": u, "T": T,
                       "salt_molar": salt_molar},
        "salt_status": "explicit_RNA_model_only" if salt_molar is not None else "program_default_not_project_buffer",
        "probability_indexing": "row=end_1based; column=unpaired_length",
        "warnings": ["RNAplfold predicts RNA-only ensembles, not catalytic Ca2+ coordination or cellular accessibility."],
    }
    (out / "lunp_manifest.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return lunp, audit
