# Results viewer

The result workspace displays the output of a completed Zyme Site Designer analysis. It works offline in a modern browser. No file is uploaded and no browser extension is required.

## Review a design

1. Open the HTML report from the desktop application. Results in a generated report load automatically. To review a different run, select **Open results** and choose its `result.json`.
2. Use the candidate site map to locate a recognition window within the target transcript. The vertical axis is the predicted joint probability that the full window is unpaired. Grey points indicate unavailable estimates; orange marks the selected design.
3. Search or sort the candidate table, then select a row to inspect its sequences and design features.
4. Select **Export design CSV** to save the complete record for the selected candidate. Expand **Design notes and interpretation** or **Raw design record** for additional detail.
5. Review the reference screen and analysis record before interpreting a ranking. Reference identity, calculation parameters, sequence digests, and warnings are retained with the results.

All displayed sequence coordinates are 1-based. For example, a cut site of `408|409` lies between nucleotides 408 and 409. Arm sequences are shown in the 5′ to 3′ direction. A reported zero is distinct from **Not available** or **Not evaluated**.

## Inspect supplied coordinates

Select **Structure viewer**, then open a PDB or mmCIF file. Drag to rotate and scroll to zoom. Chain visibility, representation, and colour controls are available beside the preview; **Reset view** restores the initial camera and **Save PNG** exports the displayed image at the browser's device-pixel resolution.

An optional `chain_mapping.json` can highlight a recognition region when the supplied chain sequence and residue numbering match the mapping. Otherwise, expand **Highlight a region** and enter the actual structure's chain identifier and residue numbers. Values in this control are not automatically full-transcript coordinates.

Confidence colouring requires explicit confirmation that the file's B-factor field contains pLDDT. Experimental B-factors are not interpreted as pLDDT automatically.

## Interpretation

Candidate scores express relative design priorities. They are not cleavage probabilities. RNA accessibility, isolated-arm hybridization estimates, and reference sequence matches do not independently determine a cleavage rate or biological specificity. Full partzyme sequences require a defined catalytic and trigger scaffold.

The coordinate viewer renders supplied coordinates; it does not predict structures or validate chemical topology. Backbone connections are visual guides. The first model and the preferred alternate location are displayed, and coordinate diagnostics are not used in candidate ranking.

## Display and privacy

The workspace uses native system fonts, scalable SVG site maps, and a device-pixel-aware canvas. Tables scroll within their panels on narrow screens. Detailed notes are collapsed by default to keep the primary review workflow compact. All styles and scripts are bundled locally; there are no remote fonts, CDNs, network requests, or analytics.
