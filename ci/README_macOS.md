# Zyme Site Designer for macOS

Choose `arm64` for Apple Silicon (M-series) and `x86_64` for Intel Macs.
The packages are built and automatically tested on macOS 15. Other macOS
versions and individual end-user workflows have not been manually tested.

Open the disk image and drag `ZymeSiteDesigner.app` to Applications. The Python,
Tcl/Tk and ViennaRNA runtime is bundled; a separate Python installation is not
required. Results are saved outside the application in the selected output
directory. Keep application-bundle contents together.

## Security and validation

This research-software release uses an ad-hoc signature. It does not have an
Apple Developer ID signature or Apple notarization. macOS security policy may
therefore prevent opening it. Follow your institution's approved installation
procedures; the application and its installer do not disable security controls.

Automatic validation includes source regression tests, packaged GUI startup,
an end-to-end GUI demo, bundled real RNA folding, sequence/AF3-input exports,
signature-integrity checks and running the executable from the mounted disk
image. These tests are not a claim of manual usability testing or Apple approval.

## License and scope

The project's own code is MIT-licensed. Bundled dependencies retain their own
licenses, included inside the application. No complete transcriptome reference,
AlphaFold account credentials or MD trajectories are supplied. Optional reference
screening requires an appropriate user-provided reference.

Software records and downloads: https://zenodo.org/records/23520422
