# ai4dnazyme

## Zyme Site Designer

Source code and reproducible macOS builds for Zyme Site Designer 1.0.0.
The application selects RNA target windows, designs complementary DNA recognition
arms, annotates local accessibility, supports optional reference screening, and
exports sequences and AlphaFold 3 input files.

[Download Windows Installer](https://zenodo.org/records/23520422)

Mac installation packages are released only after the native macOS build and
packaged application tests pass. Inspect the Actions runs and their acceptance
reports for the exact build status. There is no automatic Zenodo publication.

### Repository layout

- `software/`: current English source, interface, examples, tests and licenses.
- `ci/build_macos.py`: native build, packaged GUI/calculation/export tests and DMG checks.
- `.github/workflows/macos.yml`: separate Apple Silicon and Intel cloud builds.
- `ci/README_macOS.md`: installation instructions and validation/signing boundaries.

### Run from source

Use Python 3.10 or later with Tcl/Tk. Restore the losslessly compressed browser
example, install the package, and open the desktop workspace:

```sh
python ci/restore_assets.py
python -m pip install -e './software[desktop]'
python software/desktop_entry.py
```

The compressed `software/ui/demo-data.js.gz` preserves the complete example
byte-for-byte. Restoration verifies its SHA256 before writing `demo-data.js`.
No observations or examples are downsampled for repository size.

The macOS builds use the pinned runtime and libraries declared in the workflow
and `ci/requirements-macos.txt`. The project retains its MIT license; bundled
dependencies keep their own terms. No credentials are required or stored in this
repository. Zenodo upload and final publication approval happen separately.
