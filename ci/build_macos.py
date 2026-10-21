"""Build and validate a native macOS application. Never upload or publish it."""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import plistlib
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "software"
BUILD = ROOT / "build"
QA = BUILD / "qa"
ARTIFACTS = ROOT / "artifacts"


def run(name, command, *, cwd=SOURCE, timeout=300):
    result = subprocess.run([str(x) for x in command], cwd=cwd, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", timeout=timeout)
    QA.mkdir(parents=True, exist_ok=True)
    (QA / f"{name}.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    print(name, "PASS" if result.returncode == 0 else "FAIL", flush=True)
    if result.returncode:
        print((result.stdout + result.stderr)[-10000:], flush=True)
        raise RuntimeError(name + " failed")
    return result.stdout


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            value.update(block)
    return value.hexdigest()


def validate_result(output):
    with (output / "candidates.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 194, len(rows)
    selected = next(r for r in rows if int(r["start_1based"]) == 401 and int(r["end_1based"]) == 416)
    probability = float(selected["p_unpaired_target"])
    assert math.isclose(probability, 0.1009943374, abs_tol=1e-6), probability
    for name in ["result.json", "report.md", "candidate_sequences.fasta", "viewer/index.html"]:
        assert (output / name).is_file(), name
    assert list((output / "alphafold3_inputs").rglob("*.json")), "AF3 inputs absent"
    return {"candidate_count": len(rows), "selected_window": "401-416", "unpaired_probability": probability}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--architecture", choices=["arm64", "x86_64"], required=True)
    args = parser.parse_args()
    if sys.platform != "darwin" or platform.machine() != args.architecture:
        parser.error("Build and test on native macOS for the selected architecture.")
    if BUILD.exists() or ARTIFACTS.exists():
        parser.error("Output directories already exist; start from a clean runner.")
    QA.mkdir(parents=True)
    ARTIFACTS.mkdir()
    import RNA
    import tkinter
    import certifi
    assert RNA.__version__ == "2.7.2"
    os.environ.setdefault("SSL_CERT_FILE", certifi.where())
    sys.path.insert(0, str(SOURCE))
    from zyme_designer import __version__
    version = __version__
    run("source_tests", [sys.executable, "-B", "-m", "unittest", "discover", "-s", "tests", "-v"], timeout=600)
    for name in ["test_ui.cjs", "test_structure_viewer.cjs", "test_structure_more.cjs"]:
        run(name, ["node", "tests/" + name])
    spec = BUILD / "ZymeSiteDesigner.spec"
    data = [(str(SOURCE / folder), folder) for folder in ["ui", "examples", "docs", "third_party_licenses"]]
    data += [(str(SOURCE / name), ".") for name in ["LICENSE", "THIRD_PARTY_NOTICES.md", "README.md", "CITATION.cff"]]
    spec.write_text(
        "from PyInstaller.utils.hooks import collect_all, collect_data_files\n"
        "rna_datas, rna_binaries, rna_hidden = collect_all('RNA')\n"
        f"a = Analysis([{str(SOURCE / 'desktop_entry.py')!r}], pathex=[{str(SOURCE)!r}], "
        f"binaries=rna_binaries, datas={data!r} + rna_datas + collect_data_files('certifi'), "
        "hiddenimports=rna_hidden + ['ahocorasick', 'certifi'], "
        f"runtime_hooks=[{str(ROOT / 'ci/macos_certifi_hook.py')!r}], "
        "hookspath=[], excludes=[], noarchive=False, optimize=0)\n"
        "pyz = PYZ(a.pure)\n"
        "exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name='ZymeSiteDesigner', "
        "debug=False, bootloader_ignore_signals=False, strip=False, upx=False, console=False, "
        f"target_arch={args.architecture!r}, codesign_identity=None, entitlements_file=None)\n"
        "coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name='ZymeSiteDesigner')\n"
        "app = BUNDLE(coll, name='ZymeSiteDesigner.app', "
        f"version={version!r}, bundle_identifier='org.bnuchina2026.zymesitedesigner', "
        "info_plist={'NSHighResolutionCapable': True, 'LSMinimumSystemVersion': '11.0'})\n",
        encoding="utf-8")
    run("pyinstaller", [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
                         "--distpath", BUILD / "dist", "--workpath", BUILD / "work", spec], timeout=1200)
    app = BUILD / "dist/ZymeSiteDesigner.app"
    executable = app / "Contents/MacOS/ZymeSiteDesigner"
    info = plistlib.loads((app / "Contents/Info.plist").read_bytes())
    assert info["CFBundleShortVersionString"] == version
    architecture = run("architecture", ["lipo", "-archs", executable]).strip()
    assert architecture == args.architecture, architecture
    run("codesign", ["codesign", "--verify", "--deep", "--strict", "--verbose=2", app])
    for mode in ["smoke-test", "smoke-demo"]:
        report = QA / f"frozen_{mode}.json"
        run("frozen_" + mode, [executable, "--" + mode, "--smoke-report", report,
                               "--results-dir", BUILD / ("results_" + mode), "--smoke-max-seconds", "180"], timeout=240)
        assert json.loads(report.read_text(encoding="utf-8"))["success"]
    output = BUILD / "frozen_regression"
    run("frozen_real_folding", [executable, "--worker", "design", "--target-fasta",
        SOURCE / "examples/BIRC5_NM_001168_3.fasta", "--left-lengths", "7", "--right-lengths", "7",
        "--fold", "--fold-backend", "python", "--window", "200", "--span", "150", "--rna-salt", "1.021",
        "--scaffold", "builtin-birc5-mir21", "--af3-top", "1", "--out", output], timeout=240)
    regression = validate_result(output)
    staging = BUILD / "disk-image"
    staging.mkdir()
    run("copy_app", ["ditto", app, staging / app.name])
    os.symlink("/Applications", staging / "Applications")
    shutil.copy2(ROOT / "ci/README_macOS.md", staging / "README_macOS.md")
    shutil.copy2(SOURCE / "LICENSE", staging / "LICENSE")
    evidence = {"status": "PASS", "software_version": version,
        "platform": platform.platform(), "architecture": args.architecture,
        "python_version": platform.python_version(), "tk_version": tkinter.TkVersion,
        "dependencies": {name: importlib.metadata.version(name) for name in ["ViennaRNA", "PyInstaller", "pyahocorasick", "certifi"]},
        "source_commit": os.environ.get("GITHUB_SHA", "local"),
        "workflow_run": os.environ.get("GITHUB_RUN_ID", "local"),
        "built_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "frozen_regression": regression, "signature": "ad-hoc, not Developer ID; not notarized",
        "validation": "Native macOS CI, source tests, frozen GUI startup, GUI demo, real folding and exports",
        "manual_end_user_testing": False}
    (staging / "BUILD_INFO.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    dmg = ARTIFACTS / f"ZymeSiteDesigner_v{version}_macOS_{args.architecture}.dmg"
    run("create_dmg", ["hdiutil", "create", "-volname", "Zyme Site Designer", "-srcfolder", staging,
                       "-format", "UDZO", dmg], timeout=300)
    run("verify_dmg", ["hdiutil", "verify", dmg])
    with tempfile.TemporaryDirectory(prefix="zyme-dmg-test-") as mount:
        run("mount_dmg", ["hdiutil", "attach", "-readonly", "-nobrowse", "-mountpoint", mount, dmg])
        try:
            installed = Path(mount) / app.name
            assert (installed / "Contents/Info.plist").read_bytes() == (app / "Contents/Info.plist").read_bytes()
            run("mounted_codesign", ["codesign", "--verify", "--deep", "--strict", installed])
            result = run("mounted_executable", [installed / "Contents/MacOS/ZymeSiteDesigner", "--worker", "--version"])
            assert version in result
        finally:
            run("detach_dmg", ["hdiutil", "detach", mount])
    evidence["dmg_sha256"] = digest(dmg)
    evidence["dmg_bytes"] = dmg.stat().st_size
    evidence["disk_image_mount_and_executable_test"] = "PASS"
    (QA / "ACCEPTANCE.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    (ARTIFACTS / f"ACCEPTANCE_{args.architecture}.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    (ARTIFACTS / f"SHA256SUMS_{args.architecture}.txt").write_text(digest(dmg) + "  " + dmg.name + "\n", encoding="ascii")
    print(json.dumps(evidence, indent=2), flush=True)


if __name__ == "__main__":
    main()
