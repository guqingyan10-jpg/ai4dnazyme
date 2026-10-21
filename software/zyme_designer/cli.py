# Source guide (2026-10-20; documentation comments only).
# Purpose: Run target-site design, optional analyses and exports
# Inputs: Target FASTA; arm/motif settings; optional reference/folding inputs
# Generates/returns: Run directory: candidate CSV/FASTA, JSON, report, manifests and local viewer
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""CLI with auditable stages, optional folding and no implicit model training."""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request

from . import __version__
from .accessibility import parse_lunp, region_metrics, run_rnaplfold
from .core import (BIRC5_SCAFFOLD, DEFAULT_WEIGHTS, enumerate_candidates, ranking_audit,
                   score_candidates, sequence_hash)
from .io import file_hash, read_json, resolve_target, write_csv, write_json
from .offtarget import scan_transcriptome
from .thermo import binding_features, complementarity_risk


def lengths(value):
    try:
        result = tuple(sorted({int(n) for n in value.split(",")}))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use comma-separated integer arm lengths, e.g. 6,7,8") from exc
    if not result or any(n < 4 or n > 20 for n in result):
        raise argparse.ArgumentTypeError("Arm lengths must be within 4–20 nt")
    return result


def _sequence_manifest(manifest, seq_hash):
    recorded = manifest.get("sequence_sha256_rna") or manifest.get("transcript", {}).get("sequence_sha256_rna")
    if recorded != seq_hash:
        raise ValueError("Accessibility manifest's normalized RNA SHA256 does not match the complete input transcript")


def _import_accessibility(args, target, rows, out):
    sequence = target["sequence"]
    audit = {"status": "not_computed", "reason": "No RNAplfold run or verified imported probabilities supplied"}
    by_key = {(r["start_1based"], r["end_1based"], r["left_length"], r["right_length"]): r for r in rows}
    if args.site_table:
        if not args.accessibility_manifest:
            raise ValueError("--site-table requires --accessibility-manifest (RNA hash and parameters)")
        manifest = read_json(args.accessibility_manifest)
        _sequence_manifest(manifest, sequence_hash(sequence))
        if manifest.get("site_table_sha256") != file_hash(args.site_table):
            raise ValueError("Imported site-table SHA256 does not match its accessibility manifest")
        matched = set()
        with args.site_table.open(encoding="utf-8-sig", newline="") as handle:
            for original in csv.DictReader(handle):
                if original.get("accession") and original["accession"] != target["accession"]:
                    continue
                start, end = int(original["start_1based"]), int(original["end_1based"])
                left_raw, right_raw = original.get("left_length"), original.get("right_length")
                if bool(left_raw) != bool(right_raw):
                    raise ValueError("Site-table must specify both left_length and right_length together")
                if left_raw and right_raw:
                    left, right = int(left_raw), int(right_raw)
                elif original.get("arm_length"):
                    left = right = int(original["arm_length"])
                elif original.get("cut_after_1based"):
                    left = int(original["cut_after_1based"]) - start
                    right = end - int(original["cut_after_1based"]) - 1
                elif end - start + 1 == 16:
                    left = right = 7  # Legacy project tables were explicitly 7+AG+7.
                else:
                    raise ValueError("Site-table must specify left/right lengths or cut coordinate for non-legacy windows")
                row = by_key.get((start, end, left, right))
                if row is None:
                    continue
                if original.get("cut_after_1based") and int(original["cut_after_1based"]) != row["cut_after_1based"]:
                    raise ValueError("Imported cut coordinate disagrees with candidate arm lengths")
                if original.get("candidate_id") and original["candidate_id"] != row["candidate_id"]:
                    raise ValueError("Imported candidate ID disagrees with candidate coordinates")
                recorded_sequence = original.get("target_rna") or original.get("site_rna")
                if recorded_sequence != row["target_rna"]:
                    raise ValueError("Imported site-table sequence disagrees with candidate coordinates")
                if row["candidate_id"] in matched:
                    raise ValueError("Duplicate site-table entries for one candidate")
                legacy = left == right == 7
                for key, names in {
                    "p_unpaired_target": ("p_unpaired_target",) + (("p_unpaired_full16",) if legacy else ()),
                    "p_unpaired_left_arm": ("p_unpaired_left_arm",) + (("p_unpaired_left7",) if legacy else ()),
                    "p_unpaired_right_arm": ("p_unpaired_right_arm",) + (("p_unpaired_right7",) if legacy else ()),
                    "p_unpaired_center": ("p_unpaired_center", "p_unpaired_cut_ag", "p_unpaired_center2"),
                }.items():
                    raw = next((original[n] for n in names if original.get(n) not in (None, "")), None)
                    value = None if raw is None else float(raw)
                    if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
                        raise ValueError("Imported probability is outside [0,1]")
                    row[key] = value
                p = row["p_unpaired_target"]
                row["opening_dg_kcal_mol"] = None if p is None else -.00198720425864083 * (args.temperature + 273.15) * math.log(max(p, 1e-300))
                row["accessibility_source"] = "imported_audited_site_table"
                row["opening_probability_floor"] = 1e-300
                row["accessibility_warnings"] = (["Printed zero: capped opening proxy is a numerical limit, not a finite measured energy"] if p == 0 else [])
                matched.add(row["candidate_id"])
        audit = {"status": "imported_site_table", "input_file_sha256": file_hash(args.site_table),
                 "matched_candidates": len(matched),
                 "candidates_with_target_probability": sum(r["p_unpaired_target"] is not None for r in rows),
                 "all_candidate_probabilities_available": all(r["p_unpaired_target"] is not None for r in rows),
                 "sequence_windows_verified": True, "source_manifest": manifest,
                 "warning": "Imported analysis values, not a newly executed RNAplfold run; absent candidates remain unknown"}
    else:
        table = None
        if args.fold:
            if args.fold_backend == "python":
                from .vienna_backend import run_vienna_python
                lunp, audit = run_vienna_python(
                    sequence, out / "folding", W=args.window, L=args.span,
                    u=max(r["left_length"] + r["right_length"] + 2 for r in rows),
                    T=args.temperature, salt_molar=args.rna_salt)
            else:
                lunp, audit = run_rnaplfold(sequence, out / "folding", W=args.window, L=args.span,
                                          u=max(r["left_length"] + r["right_length"] + 2 for r in rows),
                                          T=args.temperature, executable=args.rnaplfold, salt_molar=args.rna_salt)
            table = parse_lunp(lunp)
            audit["status"] = ("new_ViennaRNA_python_run" if args.fold_backend == "python"
                               else "new_RNAplfold_run")
        elif args.lunp:
            if args.accessibility_manifest:
                audit = read_json(args.accessibility_manifest)
                _sequence_manifest(audit, sequence_hash(sequence))
                if audit.get("lunp_sha256") and audit["lunp_sha256"] != file_hash(args.lunp):
                    raise ValueError("lunp SHA256 differs from the supplied manifest")
                if not audit.get("lunp_sha256"):
                    if not args.trust_lunp_sequence:
                        raise ValueError("lunp manifest lacks lunp_sha256; supply file hash or explicit --trust-lunp-sequence")
                    audit = {**audit, "sequence_association": "user_asserted_not_independently_verified",
                             "warnings": [*audit.get("warnings", []), "Manifest lacks lunp file hash; user explicitly asserts sequence association"]}
            elif not args.trust_lunp_sequence:
                raise ValueError("--lunp requires a matching manifest, or explicit --trust-lunp-sequence with a provenance warning")
            else:
                audit = {"sequence_sha256_rna": sequence_hash(sequence),
                         "sequence_association": "user_asserted_not_independently_verified",
                         "warnings": ["User explicitly asserts that lunp came from this exact RNA; no original run manifest supplied"]}
            table = parse_lunp(args.lunp)
            if max(end for end, _ in table) != len(sequence):
                raise ValueError("lunp transcript length differs from the input transcript")
            audit = {**audit, "status": "imported_lunp", "input_file_sha256": file_hash(args.lunp)}
        if table is not None:
            for row in rows:
                metrics = region_metrics(table, row["start_1based"], row["left_length"], row["right_length"], args.temperature)
                row.update({"p_unpaired_target": metrics["full_window_p_unpaired"],
                            "p_unpaired_left_arm": metrics["left_arm_p_unpaired"],
                            "p_unpaired_right_arm": metrics["right_arm_p_unpaired"],
                            "p_unpaired_center": metrics["ag_p_unpaired"],
                            "opening_dg_kcal_mol": metrics["full_window_opening_dg_kcal_mol"],
                            "accessibility_source": audit["status"],
                            "accessibility_missing_regions": metrics["missing_probability_regions"],
                            "opening_probability_floor": metrics["opening_probability_floor"],
                            "accessibility_warnings": metrics["warnings"]})
    return audit


def _validate_import_temperature(args):
    if not args.accessibility_manifest:
        return
    data = read_json(args.accessibility_manifest)
    parameters = data.get("parameters", {})
    temperature = parameters.get("T", parameters.get("T_celsius", parameters.get("temperature_c")))
    if temperature is not None and not math.isclose(float(temperature), args.temperature):
        raise ValueError("Imported RNA probabilities are at a different temperature; do not silently combine conditions")


def design(args):
    started = time.monotonic()
    if args.out.exists():
        raise FileExistsError("Output directory already exists; use a new directory to preserve prior results")
    if args.target and not args.transcriptome:
        raise ValueError("--target accession/gene needs --transcriptome FASTA")
    if args.screen_top < 0 or args.hit_limit < 0:
        raise ValueError("screen-top and hit-limit must be nonnegative")
    if not math.isfinite(args.temperature) or args.temperature <= -273.15:
        raise ValueError("Temperature must be finite and above absolute zero")
    _validate_import_temperature(args)
    weights = read_json(args.weights) if args.weights else DEFAULT_WEIGHTS
    scaffold = (BIRC5_SCAFFOLD if args.scaffold == "builtin-birc5-mir21" else
                read_json(args.scaffold) if args.scaffold else None)
    if args.fetch:
        if not __import__("re").fullmatch(r"(?:NM|NR|XM|XR)_\d+\.\d+", args.fetch):
            raise ValueError("Fetch requires a versioned RefSeq RNA accession; no gene-name guessing")
        query = urllib.parse.urlencode({"db": "nuccore", "id": args.fetch, "rettype": "fasta", "retmode": "text"})
        request = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + query,
                                         headers={"User-Agent": "ZymeSiteDesigner/0.1 (explicit user-requested RNA retrieval)"})
        with urllib.request.urlopen(request, timeout=60) as response:
            text = response.read().decode("utf-8")
        args.out.mkdir(parents=True)
        fetched = args.out / "fetched_target.fasta"
        fetched.write_text(text, encoding="utf-8")
        target = resolve_target(fetched, args.fetch)
        target_source = fetched
    else:
        target_source = args.target_fasta if args.target_fasta else args.transcriptome
        target = resolve_target(target_source, args.target)
    if args.gene:
        if target["gene"] and target["gene"].upper() != args.gene.upper():
            raise ValueError("Provided target gene contradicts the target FASTA annotation")
        target["gene"] = args.gene
    target["sequence_sha256_rna"] = sequence_hash(target["sequence"])
    target["length_nt"] = len(target["sequence"])
    rows, enumeration = enumerate_candidates(target["sequence"], args.left_lengths, args.right_lengths, scaffold,
                                              motif_pattern=args.motif, cut_positions=args.cut_after)
    if not rows:
        raise ValueError("No unambiguous recognition windows for the requested motif, cut positions and arm lengths")
    length_policy = read_json(args.length_policy) if args.length_policy else {}
    if args.preferred_arm_min is not None:
        length_policy["preferred_min"] = args.preferred_arm_min
    if args.preferred_arm_max is not None:
        length_policy["preferred_max"] = args.preferred_arm_max
    args.out.mkdir(parents=True, exist_ok=args.fetch is not None)

    def phase(name, **extra):
        status = {"status": "running", "phase": name, "elapsed_seconds": round(time.monotonic() - started, 3), **extra}
        write_json(args.out / "progress.json", status)
        print(json.dumps(status, ensure_ascii=False), file=sys.stderr, flush=True)

    phase("SEQUENCE_ENUMERATION", candidates=len(rows))
    for row in rows:
        features = binding_features(row["target_rna"][:row["left_length"]], row["target_rna"][-row["right_length"]:],
                                    args.temperature, args.sodium, strand_concentration_uM=args.strand_concentration)
        row.update({k: v for k, v in features.items() if k not in ("warnings", "binding_reference_details")})
        row["binding_warnings"] = features["warnings"]
        risk = complementarity_risk(row["left_arm_dna"], row["right_arm_dna"], row["left_partzyme"], row["right_partzyme"])
        row.update({k: v for k, v in risk.items() if k != "warnings"})
        row["max_self_complementary_run"] = max(risk["left_arm_self_dimer_run_nt"], risk["right_arm_self_dimer_run_nt"],
                                                risk["left_arm_hairpin_run_nt"], risk["right_arm_hairpin_run_nt"])
        row["max_cross_complementary_run"] = risk["cross_arm_complementary_run_nt"]
    phase("RNA_ACCESSIBILITY")
    accessibility_audit = _import_accessibility(args, target, rows, args.out)
    score_candidates(rows, weights, length_policy=length_policy)
    # Sequence/RNA ranking selects a bounded batch for the expensive reference scan.
    reference_audit = {"status": "not_screened", "reason": "No transcriptome FASTA supplied; off-target risk is UNKNOWN"}
    hits = []
    if args.transcriptome:
        selected = list(rows if args.screen_top == 0 else rows[:args.screen_top])
        chosen = {r["candidate_id"] for r in selected}
        for row in rows:
            if row["start_1based"] in args.include_start and row["candidate_id"] not in chosen:
                selected.append(row)
                chosen.add(row["candidate_id"])
        phase("TRANSCRIPTOME_SEQUENCE_SCAN", candidate_count=len(selected), candidate_total=len(rows))
        summaries, hits, reference_audit = scan_transcriptome(
            selected, args.transcriptome, target["accession"], target["gene"], args.mismatches, args.hit_limit,
            progress=lambda data: phase("TRANSCRIPTOME_SEQUENCE_SCAN", **data))
        reference_audit["status"] = "complete_for_selected_candidates"
        reference_audit["candidate_selection"] = {"screen_top_requested": args.screen_top,
            "include_start_1based": args.include_start, "selected_candidate_ids": sorted(chosen),
            "all_candidates_screened": len(selected) == len(rows),
            "selection_before_offtarget": "available RNA accessibility and duplex sequence reference features"}
        for row in rows:
            if row["candidate_id"] in summaries:
                row.update(summaries[row["candidate_id"]])
    score_candidates(rows, weights, length_policy=length_policy)
    structure_results = {}
    if args.pdb:
        from .structure import pdb_site_metrics
        phase("OPTIONAL_3D_DIAGNOSTIC")
        mapping = read_json(args.pdb_map) if args.pdb_map else None
        for row in rows[:args.structure_top]:
            structure_results[row["candidate_id"]] = pdb_site_metrics(
                args.pdb, row["start_1based"], row["end_1based"], chain=args.pdb_chain,
                confidence_source=args.confidence_source, residue_map=mapping,
                transcript_sequence=target["sequence"])
    from .af3_export import export_candidates
    phase("ALPHAFOLD_INPUT_EXPORT")
    af3_audit = export_candidates(rows, target, args.out / "alphafold3_inputs", scaffold,
                                  top=args.af3_top, flank=args.af3_flank,
                                  calcium_count=args.af3_calcium_count)
    metadata = {"software": "Zyme Site Designer", "version": __version__,
                "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "analysis_type": "RNA cleavage-window and DNA substrate-arm prioritization, NOT cleavage-rate prediction",
                "scope_label": args.reference_label or "user-supplied FASTA scope; completeness not independently established",
                "coordinates": "1-based inclusive transcript positions; cut_after is A in A|G",
                "enumeration": enumeration, "accessibility": accessibility_audit,
                "offtarget": reference_audit, "ranking": ranking_audit(rows, weights, length_policy=length_policy),
                "alphafold3": af3_audit,
                "scaffold": scaffold, "structure_influences_ranking": False,
                "input_sha256": {"target_fasta": file_hash(target_source),
                                  **({"transcriptome": reference_audit["input_file_sha256"]} if args.transcriptome else {})},
                "parameters": {"left_lengths": args.left_lengths, "right_lengths": args.right_lengths,
                               "motif_pattern": args.motif, "explicit_cut_after_1based": args.cut_after,
                               "temperature_c": args.temperature, "sodium_mM_requested": args.sodium,
                               "strand_concentration_uM_each_reference": args.strand_concentration,
                               "max_arm_substitutions": args.mismatches, "RNA_window": args.window,
                               "RNA_max_pair_span": args.span},
                "limitations": ["A dinucleotide recognition rule is a user-selected design constraint, not proof that any catalytic core cleaves it. AG is the conservative existing 8–17 project preset.",
                    "A fixed catalytic core plus changed arms is a design strategy, not proof that structure/activity is unchanged.",
                    "Preferred arm lengths and excessive-binding penalties are configurable heuristics, not universal kinetic optima or release-rate estimates.",
                    "Reference RNA/DNA hybrid values are at 1 M NaCl, not the project Ca2+/Tris buffer.",
                    "RNAplfold is a local RNA-only ensemble; cells, RNA-binding proteins and enzyme-induced unfolding are absent.",
                    "Substitution-only off-target screening omits gaps/bulges and treats wobble as a mismatch; no tissue expression weighting.",
                    "Unknown gene annotation and ambiguous RNA windows are not evidence of absence of risk.",
                    "Scores/Pareto are interpretable heuristics, not ML-trained activity, probabilities, absolute rates or guaranteed optima.",
                    "A limited preselected reference scan can miss a superior unselected candidate; use --screen-top 0 for all candidates.",
                    "No new MD or chemical cleavage simulation is performed, and optional 3D diagnostics do not enter ranking."],
                "elapsed_seconds": round(time.monotonic() - started, 3)}
    result = {"metadata": metadata, "target": target, "candidates": rows,
              "offtarget_hits": hits, "structure_diagnostics": structure_results}
    write_json(args.out / "result.json", result)
    write_json(args.out / "run_manifest.json", metadata)
    write_csv(args.out / "candidates.csv", rows)
    write_csv(args.out / "offtarget_hits.csv", hits)
    write_json(args.out / "structure_diagnostics.json", structure_results)
    (args.out / "target.fasta").write_text(f">{target['accession']} gene={target['gene'] or 'unknown'}\n{target['sequence']}\n", encoding="utf-8")
    sequences = []
    for row in rows:
        for key in ("target_rna", "left_arm_dna", "right_arm_dna", "left_partzyme", "right_partzyme"):
            if row.get(key):
                sequences.extend((f">{row['candidate_id']}:{key} start={row['start_1based']} cut_after={row['cut_after_1based']}", row[key]))
    (args.out / "candidate_sequences.fasta").write_text("\n".join(sequences) + "\n", encoding="ascii")
    report = ["# Zyme Site Designer — analysis report", "",
              f"Transcript: {target['accession']} · Gene: {target['gene'] or 'unspecified'} · Length: {target['length_nt']} nt",
              f"Central motif: {args.motif}. Binding arms exclude the two central nucleotides.",
              f"Candidates: {len(rows)} · Reference-screened candidates: {reference_audit.get('candidate_count', 0)}",
              f"Reference: {metadata['scope_label']}", "",
              "## Top candidates", "",
              "| Rank | Target / cut | RNA | Left DNA arm | Right DNA arm | P(unpaired) | Score | Evidence stage |",
              "|---:|---|---|---|---|---:|---:|---|"]
    for row in rows[:20]:
        probability = "Not assessed" if row["p_unpaired_target"] is None else f"{100 * row['p_unpaired_target']:.2f}%"
        report.append(f"| {row['display_rank']} | {row['start_1based']}–{row['end_1based']} / {row['cut_after_1based']}\\|{row['cut_after_1based']+1} | {row['target_rna']} | {row['left_arm_dna']} | {row['right_arm_dna']} | {probability} | {row['ranking_score']:.2f} | {row['ranking_stage']} |")
    report += ["", "## Model scope and limitations", "",
               "Scores prioritize candidate recognition sites; they are not cleavage rates or drug-effect predictions.",
               "RNA/DNA duplex reference ΔG and Tm use a 1 M NaCl model, separate from assembled MNAzyme stability and Ca²⁺ reaction conditions.",
               ""] + ["- " + limitation for limitation in metadata["limitations"]]
    (args.out / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    source_ui = Path(__file__).resolve().parents[1] / "ui"
    if source_ui.is_dir():
        shutil.copytree(source_ui, args.out / "viewer")
        data = json.dumps(result, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
        (args.out / "viewer" / "result-data.js").write_text("window.ZYME_RESULT = " + data + ";\n", encoding="utf-8")
    files = [p for p in args.out.rglob("*") if p.is_file() and p.name != "progress.json"]
    write_csv(args.out / "file_manifest.csv", [{"relative_path": p.relative_to(args.out).as_posix(),
               "size_bytes": p.stat().st_size, "sha256": file_hash(p)} for p in sorted(files)])
    phase("COMPLETE", output_dir=str(args.out.resolve()), status="complete")
    print(str((args.out / "result.json").resolve()))
    return result


def parser():
    p = argparse.ArgumentParser(description="RNA cut-window / DNA substrate-arm design; no kinetic or clinical prediction")
    p.add_argument("--version", action="version", version=__version__)
    commands = p.add_subparsers(dest="command", required=True)
    d = commands.add_parser("design", help="Enumerate, annotate, screen and rank substrate-arm candidates")
    source = d.add_mutually_exclusive_group(required=True)
    source.add_argument("--target-fasta", type=Path, help="One complete transcript FASTA")
    source.add_argument("--target", help="Exact versioned accession or uniquely resolved gene in transcriptome")
    source.add_argument("--fetch", help="Explicit opt-in NCBI fetch of a versioned RNA accession")
    d.add_argument("--gene", help="Optional target gene symbol; conflicts with annotation are rejected")
    d.add_argument("--transcriptome", type=Path, help="Local plain/gzip reference RNA FASTA; required for off-target screening")
    d.add_argument("--reference-label", help="Declared reference species/release/scope, printed in reports")
    d.add_argument("--left-lengths", type=lengths, default=(7,))
    d.add_argument("--right-lengths", type=lengths, default=(7,))
    d.add_argument("--motif", default="AG", help="Two-base RNA IUPAC rule; AG is the existing 8–17 preset, RY/NN are generic user constraints")
    d.add_argument("--cut-after", type=int, action="append", help="Optional explicit 1-based cut position; both central bases remain outside the two arms")
    d.add_argument("--preferred-arm-min", type=int, help="Override configurable heuristic preferred minimum (default 7 nt)")
    d.add_argument("--preferred-arm-max", type=int, help="Override configurable heuristic preferred maximum (default 9 nt)")
    d.add_argument("--length-policy", type=Path, help="JSON with configurable arm-length and excessive-binding tradeoff parameters")
    d.add_argument("--scaffold", help="builtin-birc5-mir21 or JSON with fixed prefix/suffix; default arms only")
    d.add_argument("--temperature", type=float, default=37.)
    d.add_argument("--sodium", type=float, default=300., help="Recorded requested mM; not used for unsupported hybrid salt correction")
    d.add_argument("--strand-concentration", type=float, default=.1, help="Independent RNA/DNA reference strand concentration each, uM")
    fold = d.add_mutually_exclusive_group()
    fold.add_argument("--fold", action="store_true", help="Run installed RNAplfold on the full input transcript")
    fold.add_argument("--lunp", type=Path, help="Import full transcript RNAplfold probability table")
    fold.add_argument("--site-table", type=Path, help="Import audited coordinate/sequence-specific probabilities")
    d.add_argument("--accessibility-manifest", type=Path, help="RNA sequence SHA256 and folding parameters for imported values")
    d.add_argument("--trust-lunp-sequence", action="store_true", help="Explicit user assertion if original lunp sequence manifest is missing")
    d.add_argument("--rnaplfold", default="RNAplfold")
    d.add_argument("--fold-backend", choices=("executable", "python"), default="executable",
                   help="Local-window folding via external RNAplfold or bundled ViennaRNA Python API")
    d.add_argument("--window", type=int, default=200)
    d.add_argument("--span", type=int, default=150)
    d.add_argument("--rna-salt", type=float, help="Optional monovalent RNA-model salt in M; NOT Ca buffer")
    d.add_argument("--screen-top", type=int, default=20, help="Preselect N candidates for reference scan; 0 means ALL")
    d.add_argument("--include-start", type=int, action="append", default=[], help="Additionally scan candidates at this 1-based start")
    d.add_argument("--mismatches", type=int, choices=(0, 1, 2), default=2)
    d.add_argument("--hit-limit", type=int, default=10000, help="Max exported rows; all qualifying hits still COUNTED")
    d.add_argument("--weights", type=Path, help="JSON specifying all five heuristic score weights")
    d.add_argument("--pdb", type=Path, help="Optional mapped complete-transcript RNA PDB diagnostic, no ranking effect")
    d.add_argument("--pdb-chain")
    d.add_argument("--pdb-map", type=Path, help="JSON residue-map: 'chain:residue:insertion' -> transcript position")
    d.add_argument("--confidence-source", choices=("unknown", "plddt"), default="unknown")
    d.add_argument("--structure-top", type=int, default=3)
    d.add_argument("--af3-top", type=int, default=5, help="Export AlphaFold inputs for N leading candidates, never submit (0 disables)")
    d.add_argument("--af3-flank", type=int, default=0, help="RNA context bases added on each side in AF3 input; 0 exports only the recognition window")
    d.add_argument("--af3-calcium-count", type=int, default=0, help="Explicit exploratory modeled Ca ion count (0–10), NOT salt concentration")
    d.add_argument("--out", required=True, type=Path, help="NEW output directory; existing directories rejected")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    output_preexisted = args.out.exists()
    try:
        design(args)
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as exc:
        if not output_preexisted and args.out.is_dir() and (args.out / "progress.json").exists():
            existing = read_json(args.out / "progress.json")
            # Never alter an existing completed result when its reuse was rejected.
            if existing.get("status") != "complete":
                write_json(args.out / "progress.json", {**existing, "status": "failed", "error": str(exc)})
        p.exit(2, f"STOP: {exc}\n")


if __name__ == "__main__":
    main()
