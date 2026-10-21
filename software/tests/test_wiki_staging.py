"""Static deployment adapter tests; never contact or publish to iGEM."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("wiki_stage", ROOT / "scripts/prepare_wiki_frontend.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class WikiStagingTests(unittest.TestCase):
    def test_uploaded_assets_are_rewritten_and_csp_is_narrow(self):
        assets = {"assets/zyme-mark.svg": "https://static.igem.wiki/teams/EXAMPLE_TEST/mark.svg",
                  "assets/recognition.svg": "https://static.igem.wiki/teams/EXAMPLE_TEST/recognition.svg"}
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "staged"
            report = module.stage(assets, out)
            html = (out / "index.html").read_text(encoding="utf-8")
            self.assertIn('src="' + assets["assets/recognition.svg"] + '"', html)
            self.assertIn("connect-src 'none'", html)
            self.assertIn("img-src 'self' data: https://static.igem.wiki;", html)
            self.assertFalse(report["external_computation"])
            self.assertTrue((out / "demo-data.js").is_file())
            with self.assertRaises(FileExistsError):
                module.stage(assets, out)

    def test_external_asset_urls_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                module.stage({"assets/zyme-mark.svg": "https://evil.example/mark.svg",
                              "assets/recognition.svg": "https://static.igem.wiki/example.svg"}, Path(tmp) / "out")

    def test_missing_asset_map_is_not_silently_assumed(self):
        with self.assertRaises(ValueError):
            module.stage({}, Path("never_created_without_asset_map"))

if __name__ == "__main__":
    unittest.main()
