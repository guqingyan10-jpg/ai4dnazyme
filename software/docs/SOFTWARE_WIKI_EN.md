# Zyme Site Designer

## Overview

Designing an RNA-cleaving MNAzyme begins with a sequence-level question: which region of a complete transcript should the recognition arms target? In our BIRC5 workflow, candidate junctions had to be considered together with RNA folding, arm complementarity and competing sequences. Repeating these operations manually made it difficult to compare designs and retain a consistent record of the assumptions behind each choice.

We developed **Zyme Site Designer** to bring these steps into one reproducible workflow. The software accepts a complete transcript, generates candidate target windows and complementary DNA arms, annotates them with local accessibility and duplex reference features, and optionally screens a defined RNA reference library. Its outputs retain the input scope, parameters and component scores used to prioritise each candidate.

Version 1.0.0 provides an English Windows desktop interface and a command-line pipeline. Users can complete the main workflow in three steps: **Transcript → Design → Results**.

## Design: an explicit, reusable workflow

We separated the recognition-arm design from the catalytic scaffold. The default result therefore contains two DNA binding arms. A user can explicitly select the project's BIRC5 / miR-21 scaffold, or provide a custom scaffold, when fixed core and trigger modules should be retained.

The AG preset selects the conservative junction policy used in our existing 8–17 design. Custom two-base IUPAC rules and explicit cut positions support other sequence-selection policies. Left and right arm lengths are independently adjustable; the desktop default explores 7, 8 and 9 nt on each side rather than imposing one arm length on every site.

The pipeline combines four kinds of evidence:

| Design question | Software output |
|---|---|
| Does the window match the selected recognition policy? | Transcript coordinates, central junction and verified reverse-complementary DNA arms |
| How accessible is the region in the RNA folding model? | Local-window joint unpaired probabilities and model provenance |
| How do the two arms compare? | Independent RNA/DNA reference energies, arm balance and sequence-complementarity diagnostics |
| Which reference sequences may compete for recognition? | Substitution-based hit counts and matched windows within a declared reference FASTA |

We used configurable, visible ranking terms rather than an activity-trained model. Candidate records include both the score and the available evidence stage, so a design without reference screening does not appear to have an assessed zero risk. The full equations and model settings are documented in [Methods](METHODS.md).

## Build: desktop access and interoperable outputs

The native desktop application uses the same analysis pipeline as the command-line interface. An independent worker handles calculation while the interface displays progress and saved logs. Each run produces a separate directory containing candidate CSV, sequence FASTA, result JSON and a parameter/input manifest.

ViennaRNA 2.7.2 supplies the local RNA-window probabilities. A plain or gzip FASTA can be used for reference screening; no full-human-reference download is assumed. The HTML report provides candidate inspection and export without an external CDN. AlphaFold 3 input JSON and chain mappings connect selected designs to a separate structure-prediction workflow, and the local viewer accepts existing PDB or mmCIF coordinates.

The Windows distribution includes its runtime for users who do not have Python or Conda installed. The source remains available for command-line use and further development. The source tree includes tests, documented dependencies and build scripts.

## Test: numerical regression and reuse

The BIRC5 regression reproduces the original **NM_001168.3** target window at **401–416 nt**, with the **408|409** junction. Under the recorded W200/L150/37 °C RNA-model conditions, its joint unpaired probability is approximately **10.10%**, and the recovered 7 nt arms are **TTGCAAT** and **TGGTTTC**. Twelve archived region probabilities agreed with the local ViennaRNA backend to an absolute error below `5.1e-8`.

| Example calculation | Candidate designs | Full-window probability available | Reference scope |
|---|---:|---:|---|
| Archived BIRC5 W200 demonstration, 7+7 nt | 194 | 3 | Small demonstration FASTA |
| BIRC5 new local calculation, 7/8/9 nt per arm | 1746 | 1746 | Small demonstration FASTA |
| GAPDH NM_002046.7 trial, 7/8/9 nt per arm | 693 | 693 | No reference supplied |

The GAPDH trial uses a complete 1285 nt transcript to demonstrate that the workflow is not hard-coded to BIRC5. Its reference-screening state remains unassessed because no reference was supplied.

For version 1.0.0, the source validation discovered **114 Python tests: 113 passed and one was skipped**, and **56 JavaScript assertions passed**. The skip requires an independent RNAplfold executable that was unavailable in that environment. Tests cover sequence and coordinate handling, imported-probability provenance, numerical regression, reference-hit classification, missing-data states, exports and structure parsing. Reproduction commands and dated evidence are provided in [Validation](VALIDATION.md).

## Learn and iterate: from a working pipeline to a usable application

User testing exposed a concrete interface failure: the expandable content area pushed the primary action below the visible window, leaving the user unable to start a design. The initial interface also presented too many controls at once.

We revised the layout to reserve a fixed bottom action bar for **Run analysis**. The main workflow now follows Transcript, Design and Results, while Activity and Advanced settings retain access to logs and detailed parameters. Loading the BIRC5 example is a separate File-menu action; it allows users to inspect the setup before running it. English labels and native Windows DPI-awareness handling address international use and scaled-display rendering.

The new geometry checks passed **36 primary-action visibility cases** across nine window-size/scaling combinations and four tabs. Expanded Advanced settings and English public labels were checked separately. These records turn the reported interface problem into a documented design–build–test–learn iteration, rather than a cosmetic change without verification.

## Use the software

1. Extract the complete Windows ZIP and open `ZymeSiteDesigner.exe`.
2. Load one transcript in Transcript and review the motif, arms and accessibility method in Design.
3. Select Run analysis, inspect Results and use Open report to access the complete HTML report.

For a first example, select **File → Load BIRC5 example**, then **Run analysis**. For a new transcript, load its own FASTA and use local ViennaRNA. A full-reference screen additionally requires a source/version-defined local reference library.

The [README](../README.md) provides source installation and CLI commands; the [user guide](USER_GUIDE.md) covers inputs, results, exports and troubleshooting. The [iGEM submission guide](IGEM_SUBMISSION.md) distinguishes source publication, the software Wiki page and the portable distribution.

## Scope and future development

Zyme Site Designer prioritises recognition designs; it does not infer chemical cleavage rates or certify new catalytic constructs. Accessibility is an RNA-only local ensemble calculation, and hybrid energies are independent duplex reference values rather than whole-MNAzyme MM/GBSA results. Reference screening is confined to the supplied library and selected candidates. The demonstration library is not a complete human transcriptome.

Future development can add experimentally labelled candidate sets and additional recognition models while retaining explicit training/validation separation and input provenance. The present implementation, model assumptions and distribution terms are described in [Limitations](LIMITATIONS.md) and [Third-party notices](../THIRD_PARTY_NOTICES.md).

The project's own source uses the [MIT License](../LICENSE). Bundled dependencies retain their original terms. This local release package does not assert a public repository, DOI or completed official submission; the published Wiki should point to the team's actual iGEM source release and download artifacts after those have been uploaded.
