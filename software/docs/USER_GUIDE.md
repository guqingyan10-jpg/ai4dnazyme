# Zyme Site Designer — Source User Guide

Unless stated otherwise, run the commands below from the extracted `software/` directory.

## 1. Set up Python

Python 3.10 or later is required. An isolated environment is recommended.

Windows:

```text
python -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
```

macOS or Linux, using a compatible POSIX shell:

```text
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

The core design code uses the Python standard library. ViennaRNA provides local RNA accessibility calculations. `pyahocorasick` accelerates reference searching; the application also has a standard-library fallback. Tcl/Tk is required by the desktop interface and is not installed by `requirements.txt`.

The install command uses `-e ./src`, so run it from `software/`. Alternatively, enter `software/src/` and install with `python -m pip install -e ".[desktop]"`.

## 2. Launch the desktop interface

```text
python src/desktop_entry.py
```

Follow **Transcript → Design → Results**:

1. Load one complete target FASTA or paste a sequence. Retain the accession version when using a reference transcript.
2. Set the recognition-arm lengths, junction rule and desired analysis options.
3. Select a results directory and run the analysis.
4. Inspect candidates, their evidence status, reports and sequence exports.

For a new transcript, use local ViennaRNA folding or supply probabilities that match the sequence and folding parameters. To assess reference matches, supply a versioned RNA FASTA and a reference-scope label. Without a reference, specificity is unassessed—not zero risk.

## 3. Check source startup and output generation

This example performs sequence design without folding or reference screening:

```text
python -m zyme_designer design --target-fasta src/examples/BIRC5_NM_001168_3.fasta --gene BIRC5 --left-lengths 7 --right-lengths 7 --motif AG --af3-top 0 --out results/BIRC5_demo
```

The exact supplied sequence and settings generate 194 candidates. This is a software check, not a complete research analysis. Accessibility and specificity remain unassessed. The output directory must be new; an existing directory is rejected to protect earlier results.

## 4. Calculate local RNA accessibility

```text
python -m zyme_designer design --target-fasta src/examples/BIRC5_NM_001168_3.fasta --gene BIRC5 --left-lengths 7,8,9 --right-lengths 7,8,9 --motif AG --fold --fold-backend python --temperature 37 --window 200 --span 150 --rna-salt 1.021 --out results/BIRC5_fold
```

Replace the FASTA path to design for another transcript. Imported accessibility must match the input sequence, coordinates and parameter manifest; a cache for one transcript is not a folding result for another.

## 5. Add reference screening

Add these arguments to a design command:

```text
--transcriptome <path-to-reference-RNA-FASTA> --reference-label <species-release-and-scope>
```

`--screen-top 0` scans all candidates. A large reference may substantially increase runtime. A small example reference is not a full human transcriptome. A compatible sequence hit is a screening result, not measured cleavage of that RNA.

For all options:

```text
python -m zyme_designer design --help
```

## 6. Inspect the output directory

| File or directory | Contents |
|---|---|
| `candidates.csv` | Coordinates, recognition arms, accessibility, scores and evidence status |
| `candidate_sequences.fasta` | Candidate recognition-arm or construct sequences |
| `target.fasta` | Target sequence used in the run |
| `result.json` | Structured results |
| `report.md` | Human-readable analysis report |
| `run_manifest.json` | Input provenance and run settings |
| `file_manifest.csv` | Output inventory and integrity information |
| `progress.json` | Job status |
| `viewer/` | Local result and structure viewer |

Reference-hit, construct, structure and AF3 files depend on enabled options. Preserve the input FASTA, result files and manifests when sharing a run.

## 7. Export complete constructs and AF3 inputs

The default recognition-arm output does not contain a catalytic core. Use `--scaffold builtin-birc5-mir21` only when retaining the project's fixed core and trigger modules is intended, or supply an appropriate custom scaffold.

AF3 exports contain job requests and chain mappings, not predicted coordinates. Submit predictions separately through the chosen external workflow. The application does not upload sequences, run AF3 or download model weights. Each AF3 export includes an English `README.md`.

An imported PDB or mmCIF can be inspected in the structure viewer. Match the chain mapping to the selected candidate before interpreting a cut-site highlight. A short RNA context does not represent the folding of the complete transcript.

## 8. Run tests

From `software/src/`:

```text
python -B -m unittest discover -s tests -v
node tests/test_ui.cjs
node tests/test_structure_viewer.cjs
node tests/test_structure_more.cjs
```

Python tests use `unittest`. The three JavaScript checks require Node.js. Tests that need a usable Tk installation or the optional real ViennaRNA backend may skip when those components are absent. The tests actually run for this distribution are summarized in `VALIDATION_REPORT.md`.

## Troubleshooting

| Symptom | Action |
|---|---|
| `RNA` cannot be imported | Install the requirements using the same Python interpreter used to run the application. |
| Tk cannot start | Configure a working Tcl/Tk installation for that interpreter. |
| Output directory already exists | Choose a new `--out` directory; do not delete previous results just to reuse its name. |
| Imported accessibility is rejected | Check the target sequence hash, coordinates and folding manifest. |
| Reference screening is unassessed | Supply a reference RNA FASTA and its scope label. |
| No candidate windows are generated | Check the junction rule, input sequence and available flanking lengths. |
| A run stops | Read the worker log and preserve its input and manifest for diagnosis. |
