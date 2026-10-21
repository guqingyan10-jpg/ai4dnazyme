# Architecture

The desktop and command-line entry points use the same analysis pipeline. The desktop is a thin interaction layer: it validates input, builds explicit argument lists, launches a worker and reads its progress and outputs.

```text
Desktop: Transcript → Design → Results
                 │ validated configuration
                 ▼
        Separate analysis worker ◀── Command-line interface
                 │
                 ├── Transcript resolution and sequence checks
                 ├── Motif / arm enumeration
                 ├── Duplex features and RNA accessibility
                 ├── Optional reference-transcriptome scan
                 ├── Evidence-aware ranking
                 └── CSV / FASTA / JSON / HTML / AF3 input export
```

## Components

| Component | Responsibility |
|---|---|
| `desktop_entry.py` | Shared source/portable entry point and UTF-8 worker logging |
| `zyme_designer/desktop.py` | Native GUI, configuration, worker lifecycle and result inspection |
| `zyme_designer/cli.py` | Pipeline orchestration, progress, provenance and export |
| `zyme_designer/core.py` | Sequence normalization, candidate enumeration and ranking |
| `zyme_designer/vienna_backend.py` | ViennaRNA local-window calculation and model audit |
| `zyme_designer/accessibility.py` | Probability import and region indexing |
| `zyme_designer/thermo.py` | RNA/DNA reference thermodynamics and complementarity diagnostics |
| `zyme_designer/offtarget.py` | Streaming reference scan and hit classification |
| `zyme_designer/af3_export.py` | AlphaFold input generation and chain mapping |
| `zyme_designer/structure.py` | Optional mapped-coordinate diagnostics |
| `ui/` | Local HTML result and structure viewers |
| `tests/` | Backend, interface and numerical regression tests |

## Data and process handling

The GUI uses an independent worker process to keep calculation off the interface event loop. It polls `progress.json` and reads the worker log. Commands are supplied as argument lists rather than interpolated shell commands.

Each desktop run receives a unique result directory; CLI runs reject an existing directory. Inputs are read locally. Cancellation affects the current worker only and retains existing files.

The report consumes `result.json`, or the packaged `result-data.js` representation. Text imported from user files is treated as data. The report does not require an external CDN.

## Extending the project

Add new numerical features in a backend module, include their units and provenance in the run metadata, and test missing-data behaviour. A new ranking feature should document its equation and scope rather than changing an existing field's meaning. Add deterministic fixtures before connecting a feature to the GUI or report.

New transcripts and reference libraries are data inputs, not source-code changes. New catalytic scaffolds require an explicit configuration and their own biological validation.
