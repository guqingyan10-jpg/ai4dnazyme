# Zyme Site Designer — Local Analysis and Static Wiki Integration

This guide describes the supplied source interfaces and files. It is not a record of a completed upload, public deployment or official submission.

## Separate calculation from result viewing

| Interface | Function | Not provided by that interface |
|---|---|---|
| Native Tk desktop and Python CLI | Sequence input, local folding, ranking and optional reference screening | A hosted calculation service |
| Static `ui/` frontend | Load saved results, compare candidates and inspect existing molecular coordinates | Browser-side RNAplfold or a full-transcriptome calculation backend |

`result.json` is the file interface between calculation and display. `result-data.js` can embed a saved result in a static page. The current implementation does not provide a REST API. Describe the web interface as a result explorer rather than an online full-analysis service.

## Files for the Wiki team

Paths below are relative to `software/src/`.

| Files | Purpose |
|---|---|
| `ui/index.html`, `ui/guide.html`, `ui/structure.html` | Result, help and structure pages |
| `ui/styles.css`, `ui/app.js`, `ui/structure.js` | Active frontend styles and behavior |
| `ui/result-data.js`, `ui/demo-data.js` | Saved-run data and the supplied example; preserve evidence and reference-scope labels |
| `ui/assets/` | Required artwork and local resources |
| `scripts/prepare_wiki_frontend.py` | Stage the frontend and map uploaded artwork URLs |
| `docs/METHODS.md`, `docs/VALIDATION.md`, `docs/LIMITATIONS.md` | Method, validation and scope descriptions |
| `LICENSE`, `THIRD_PARTY_NOTICES.md`, `third_party_licenses/` | Original attribution and license information |

The source-only ZIP does not include a Windows executable or runtime. Installation and local calculation are documented in the outer `README.md` and `USER_GUIDE.md`.

## Stage the frontend

From `software/src/`, the helper can prepare a new staging directory:

```text
python scripts/prepare_wiki_frontend.py --asset-map uploads.json --out wiki_staged
```

The asset-map JSON uses `assets/zyme-mark.svg` and `assets/recognition.svg` as keys, with actual uploaded artwork URLs as values. The helper validates the configured URL form, copies the frontend and adjusts the image content-security policy. It does not upload files, publish a website or change scientific result records. Do not place fabricated upload URLs in a production configuration.

Integrate the staged files into the team's existing Wiki build. Check navigation, relative paths, CSS scope, content-security policy, mobile layout and accessibility. Keep the local analysis source available independently of the result explorer. A suggested tool subroute is `/software/zyme-site-designer/`; it is not a claim that this address has been published.

## Describe the interfaces accurately

Suggested actions are **Explore the BIRC5 example**, **Open existing results** and **Run the full analysis locally**. A file selected from the user's computer is read by the browser; this is not an upload to a calculation server.

Keep the supplied example's actual reference coverage visible. Missing screening must not become a displayed zero-risk result. Recognition-arm exports must not be labeled as complete catalytic constructs unless the relevant scaffold was explicitly selected.

## Publication checks

Before deployment, consult the team's current [Wiki requirements](https://teams.igem.org/go/deliverables/wiki/requirements), [Wiki FAQ](https://teams.igem.org/go/deliverables/wiki/faq) and [software requirements](https://teams.igem.org/go/deliverables/software/requirements). This source documentation does not certify compliance with rules that may change.

Verify the deployed page's resource requests, including interactions after initial load. The supplied frontend uses local assets rather than a third-party script CDN; the deployment must preserve the intended resource and security configuration. Fill in actual team, repository and download addresses only when they exist.

Retain project and dependency licenses. Do not use the project's MIT source license to replace third-party terms. The license and publication requirements for Wiki content should be handled separately from software dependency licensing.
