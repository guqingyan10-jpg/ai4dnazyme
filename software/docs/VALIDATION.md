# Validation and reproduction

Software validation covers sequence generation, numerical annotation, screening and interface behaviour. Biological evidence from the original project is separate from verification of the software implementation.

## Historical evidence

The archived 2026-10-06 validation environment was Windows 11 with Python 3.13.9 and ViennaRNA 2.7.2. Its Python suite discovered 111 tests: 110 passed and one was skipped. The skipped comparison required an independent RNAplfold executable that was not installed. It is not recorded as a pass. JavaScript checks covered 27 report assertions, 12 structure-parser assertions and 10 additional structure edge-case assertions.

These counts belong to that dated development snapshot. They must not be quoted as the test count of version 1.0.0. Historical command outcomes are recorded in `qa/validation_summary.json` and the accompanying logs.

The real ViennaRNA regression compared 12 region probabilities against archived BIRC5 data, with absolute error below `5.1e-8`. The target was `NM_001168.3`, 2574 nt. At 37 °C, W200/L150 and the recorded 1.021 M RNA-model salt, the original 401–416 nt window reproduced a probability of approximately 0.1009943374. Its cut is 408|409 and its 7 nt arms are `TTGCAAT` and `TGGTTTC`.

## Verified example scope

| Example | Candidates | Full-window probabilities | Reference scope |
|---|---:|---:|---|
| BIRC5 archived W200 demonstration, 7+7 nt | 194 | 3 candidates | Small demonstration FASTA |
| BIRC5 new local calculation, 7/8/9 nt per arm | 1746 | 1746 candidates | Small demonstration FASTA |
| GAPDH NM_002046.7, 7/8/9 nt per arm | 693 | 693 candidates | No reference supplied |

The GAPDH example used a complete 1285 nt transcript and was repeated with the portable version 1.0.0 worker on 2026-10-07. It exercises reuse on a different transcript. It is not a full-human-reference scan or validation of a new catalytic construct.

## Version 1.0.0 acceptance checks

The source checks completed on 2026-10-07 discovered **114 Python tests: 113 passed and one was skipped**. The skipped test was the independent RNAplfold-executable comparison; that executable was unavailable. The **56 JavaScript assertions passed**: 34 report assertions, 12 structure-parser assertions and 10 further structure edge-case assertions.

Three new native GUI tests passed. Geometry checks covered nine window-size/scaling combinations across all four tabs: **36 primary-action visibility cases**. The tested sizes were 900×600, 1024×680 and 1280×800; Tk scaling factors were 4/3, 5/3 and 2. These factors approximate 100%, 125% and 150% point-to-pixel scaling, rather than a hardware test of every monitor. Expanded Advanced settings and English public labels were also checked.

The version 1.0.0 frozen executable passed the same 36 layout cases and its GUI example. Actual worker calculations recovered 1746 BIRC5 designs, 4914 custom RY designs and 693 GAPDH designs; the original BIRC5 window reproduced P(unpaired) = 0.1009943374. The embedded Windows manifest was checked for PerMonitorV2 DPI awareness. These checks used a sanitized executable search path without Python or Conda directories.

The real-browser report checks passed candidate selection, search, sorting, JSON re-import and CSV export. Both 1440 px desktop and 430 px mobile views remained within the viewport. Native desktop screenshots were rendered from this application's own transparent, input-disabled window, not by capturing another application or the desktop.

Source checks and frozen-build checks are recorded separately in `qa/release_v1_0_0/` and `PACKAGED_QA.json`. The Windows archive is extracted to a fresh directory and re-tested before delivery. These are same-machine acceptance checks, not independent clean-machine or all-monitor certification.

## Reproduce the tests

From the source directory:

```text
python -m pip install -e ".[desktop]"
python -B -m unittest discover -s tests -v
node --check ui/app.js
node tests/test_ui.cjs
node --check ui/structure.js
node tests/test_structure_viewer.cjs
node tests/test_structure_more.cjs
```

Retain the Python version, operating system, dependency versions, command, exit status, skipped-test reasons and full logs. A CI configuration is supplied; its presence does not establish that a team's remote CI has run.

## Verify a packaged release

1. Extract the complete archive to a fresh directory.
2. Open the executable and confirm that Run analysis remains visible when the window is resized.
3. Load and run the BIRC5 example; inspect the reference scope and completed output manifest.
4. Run a new-transcript example with local ViennaRNA and compare its input hash, design settings and candidate count.
5. Open the HTML report, inspect the exported CSV/FASTA/JSON and check a local structure input.

A test on one Windows machine is not certification of every Windows version, monitor or independent clean virtual machine. Keep that environment distinction with the acceptance record.
