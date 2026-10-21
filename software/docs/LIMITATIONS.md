# Limitations

## Design and ranking

- Motif rules are design constraints. A motif-compatible junction is not established as cleavable by an arbitrary DNAzyme core.
- Retaining a scaffold while changing recognition arms does not establish unchanged structure or catalytic activity.
- Weights, preferred arm lengths and strong-binding penalties are configurable heuristics. They are not trained probabilities, kinetic optima or product-release rates.
- The program does not calculate `kchem`, `kobs`, `kcat`, catalytic efficiency or dose–response. Candidate prioritisation and fluorescence kinetics answer different questions.

## Physical models

- RNA accessibility describes a local RNA-only secondary-structure ensemble. It omits cellular protein occupancy, translation, modifications and enzyme-induced unfolding.
- RNA-model monovalent salt does not represent explicit Ca²⁺ catalytic chemistry.
- Hybrid reference energies and melting temperatures describe independent perfect RNA/DNA duplexes at the stated reference conditions. They are not whole-MNAzyme binding energies, MM/PBSA or MM/GBSA estimates, dissociation constants or melting temperatures in the reaction buffer.
- Combining the opening-cost and duplex proxies does not establish a complete thermodynamic cycle for the assembled complex.

## Reference screening

- Screening coverage is limited to the supplied FASTA and selected candidates. No complete human reference is bundled.
- The demonstration reference is a small example and must not be described as a whole-transcriptome specificity assessment.
- Matching uses substitutions only: insertions, deletions, bulges, tissue expression and RNA occupancy are not modelled. GU wobble is not treated as a separate binding model.
- A sequence hit describes a potential recognition competitor, not a measured off-target cleavage event. Missing annotation or absent screening does not imply safety.

## Structure functions

- AlphaFold 3 input export does not submit a prediction. A short RNA context or arm-only construct is not a full transcript or full catalytic complex.
- An imported structure needs sequence and residue mapping. B-factors are not automatically pLDDT, and structural diagnostics do not determine the main ranking.

## Distribution and validation

- Windows acceptance records apply to the tested environments, not every computer or monitor. Display-dependent and dependency-dependent tests can be skipped explicitly.
- Project source is MIT-licensed. Bundled third-party components retain their own terms; redistribution and commercial integration require checking those terms and the runtime's permitted distribution conditions.
- A local release package is not evidence of completed public GitLab publication, a DOI deposit or official iGEM submission.
