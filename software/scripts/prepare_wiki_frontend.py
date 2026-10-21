# Source guide (2026-10-20; documentation comments only).
# Purpose: Stage the active software frontend for a Wiki
# Inputs: src/ui files and staging destination
# Generates/returns: Frontend files copied to the selected destination
# Result location: user-selected --out/desktop folder; helpers may return values only.
# Documentation: software/README.md and software/USER_GUIDE.md.
"""Stage the static viewer with explicit, uploaded iGEM artwork URLs.

No upload, repository modification or publishing is performed by this script.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]

def stage(asset_map: dict, out: Path, ui: Path | None = None) -> dict:
    ui = ui or ROOT / "ui"
    required = {"assets/zyme-mark.svg", "assets/recognition.svg"}
    if set(asset_map) != required:
        raise ValueError("Provide exactly the two artwork paths as keys: " + ", ".join(sorted(required)))
    for url in asset_map.values():
        if not isinstance(url, str):
            raise ValueError("Each asset URL must be a string.")
        parts = urlsplit(url)
        if parts.scheme != "https" or parts.netloc != "static.igem.wiki" or parts.username or parts.password:
            raise ValueError("Artwork must use a real https://static.igem.wiki upload URL.")
        if any(char in url for char in ['"', "'", "<", ">", "\n", "\r"]):
            raise ValueError("Invalid characters in an artwork URL.")
    if out.exists():
        raise FileExistsError("Output exists; choose a new directory to preserve earlier versions.")
    shutil.copytree(ui, out)
    for html in out.glob("*.html"):
        content = html.read_text(encoding="utf-8")
        for local, url in asset_map.items():
            content = content.replace('href="' + local + '"', 'href="' + url + '"')
            content = content.replace('src="' + local + '"', 'src="' + url + '"')
        # Only this known asset origin is added; external computation stays blocked.
        content = content.replace("img-src 'self' data:;", "img-src 'self' data: https://static.igem.wiki;")
        html.write_text(content, encoding="utf-8")
    report = {
        "status": "staged_only_not_uploaded", "artwork": asset_map,
        "external_computation": False,
        "check_before_publishing": [
            "Integrate the page and all JS/CSS/data scripts into the team's iGEM-hosted build.",
            "Use real assigned routes and review the parent Wiki's CSP and navigation.",
            "Run the official external-content check, including after loading the example.",
            "Keep the tiny-reference label and existing missing-data markers.",
        ],
    }
    (out / "WIKI_STAGING.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset-map", required=True, type=Path, help="JSON with real URLs returned by iGEM Uploads")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(stage(json.loads(args.asset_map.read_text(encoding="utf-8")), args.out), indent=2))

if __name__ == "__main__":
    main()
