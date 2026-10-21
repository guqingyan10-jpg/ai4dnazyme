# Methods

Zyme Site Designer separates sequence generation, model annotation, reference screening and ranking. Each run records the inputs, model configuration and evidence coverage used to obtain its result.

## Candidate generation

The target is one transcript with 1-based coordinates. The input alphabet is normalized to RNA. A two-base IUPAC rule identifies eligible junctions; AG is the conservative preset retained from the project's selected 8–17 design. Explicit cuts must match the declared motif.

For a cut after position `c` and left/right arm lengths `l` and `r`:

```text
target window: [c-l, c+1+r]          length: l+r+2
central junction: [c, c+1]
left RNA flank: [c-l, c-1]
right RNA flank: [c+2, c+1+r]
DNA arm: reverse complement of its RNA flank, written 5′→3′
```

Incomplete flanks and ambiguous candidate windows are skipped. Generated coordinates, lengths and complementarity are checked. The supported arm-length interval is 4–20 nt; the desktop default is independently `7,8,9` on each side.

A scaffold is optional. It supplies fixed DNA prefix/suffix segments and, where specified, a trigger with checked segment complementarity. The built-in BIRC5 / miR-21 scaffold declares an AG motif policy. This keeps scaffold selection explicit rather than assigning that core to every generic motif.

## RNA accessibility

The local Python backend uses **ViennaRNA 2.7.2** window partition functions. The starting configuration is W = 200 nt, L = 150 nt and 37 °C. The maximum unpaired segment length is derived from the longest requested target window. Requested and effective settings are both recorded if a short input requires smaller dimensions.

For RNAplfold-style `lunp` output, the value at row `end`, column `length` is the joint unpaired probability of `[end-length+1, end]`. The program extracts probabilities for the complete target window, the two flanks and the central junction. It does not multiply per-base probabilities to approximate a joint window probability.

The opening-cost proxy is:

```text
G_open = -R*T*ln(P_full_window)
```

Temperature is in kelvin and the result is in kcal/mol. Zero probabilities retain their original value; numerical clipping is used only to evaluate a finite logarithm. Missing probabilities remain missing.

Model provenance includes the RNA sequence hash, temperature, window dimensions, monovalent salt setting and energy-parameter hash. The default RNA-model salt of 1.021 M is a model setting, not the project's calcium-containing reaction buffer. Imported tables require matching sequence and provenance checks. See the [ViennaRNA window API](https://viennarna.readthedocs.io/en/latest/partfunc/window.html).

## Duplex and complementarity features

Each arm uses the perfect RNA/DNA nearest-neighbour reference parameters of [Sugimoto et al. (1995)](https://doi.org/10.1021/bi00035a029). Enthalpy and entropy are accumulated over the RNA dinucleotides and initiation term; `G_arm = H - T*S`. Reference melting temperatures use the stated independent strand concentration, defaulting to 0.1 µM per strand, at 1 M NaCl.

The requested experimental sodium concentration is recorded but is not used in an unsupported hybrid/Ca²⁺ correction. Lengths outside the original 5–13 nt oligonucleotide range receive an extrapolation warning.

Sequence-complementarity features report the longest contiguous Watson–Crick matches within and between the DNA arms. Complete-partzyme complementarity is reported separately because an intended core stem may be present.

## Reference screening

Reference input is a local plain or gzip RNA FASTA. The program checks both recognition flanks using a total budget of zero, one or two substitutions. The central motif is classified separately and does not consume the arm mismatch budget.

Exact and single-substitution arm seeds identify windows for full verification. With a two-substitution total budget, at least one arm has no more than one substitution. Optional Aho–Corasick matching accelerates seed lookup; the standard-library fallback retains the same matching policy.

Counts distinguish the intended transcript, other records of the same gene, other annotated genes and unknown annotation. Compatible and non-compatible central motifs are reported separately. By default, the 20 pre-ranked candidates are screened; `--screen-top 0` screens all candidates. Exported hit rows can be limited while total qualifying counts continue to accumulate.

## Ranking

Ranking is deterministic and configurable. Let `a` and `b` be the arm reference free energies and `q` the largest arm complementarity run:

```text
accessibility = P_full_window
binding       = clip(-(a+b)/20, 0, 1)
balance       = exp(-abs(a-b)/3)
sequence_risk = exp(-max(q-4, 0)/3)

R_hits = 10*N_other_exact_compatible
       + 2*N_other_1mm_compatible
       + 0.25*N_other_2mm_compatible
       + 0.05*N_other_noncompatible
       + 0.10*N_unknown_gene_hits
       + max(N_intended_target_hits-1, 0)
specificity = 1/(1+R_hits)
```

Default weights are accessibility 0.35, binding 0.20, balance 0.10, specificity 0.30 and sequence risk 0.05. Available components are renormalized; missing components are not assigned zero:

```text
score_base = 100*sum(w_j*x_j for available j)/sum(w_j for available j)
```

The arm-length trade-off uses distances from preferred and practical intervals, plus a penalty for very strong reference binding:

```text
d(n;low,high) = max(low-n, 0, n-high)
C_preferred = (d(l;7,9)^2 + d(r;7,9)^2)/(2*3^2)
C_practical = (d(l;6,12)^2 + d(r;6,12)^2)/(2*2^2)
C_strong    = (max(0,-14-a)^2 + max(0,-14-b)^2)/(2*6^2)
score_final = score_base*exp(-(C_preferred+C_practical+C_strong))
```

All intervals, coefficients and scales above are heuristics recorded in the manifest. A missing energy leaves the strong-binding penalty unassessed. Evidence stages are ordered before scores: screened with accessibility; screened without accessibility; accessibility only; sequence only. Pareto comparison is restricted to candidates with all five score components and includes the length trade-off.

These scores prioritise designs for further assessment; they are not fitted to fluorescence kinetics. See [Limitations](LIMITATIONS.md) for the scope of the physical proxies.

## Structure export and diagnostics

Leading candidates can be exported as AlphaFold Server and local AlphaFold 3 input JSON, with a chain mapping. The default is five candidates. Exports contain the chosen arms or scaffold and the stated RNA context. External prediction is a separate step. See [AlphaFold Server input documentation](https://github.com/google-deepmind/alphafold/blob/main/server/README.md) and [AlphaFold 3 input documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md).

Optional coordinate diagnostics require sequence/residue mapping. Imported B-factors are interpreted as pLDDT only when that provenance is explicitly supplied. Structure diagnostics do not enter the candidate ranking.
