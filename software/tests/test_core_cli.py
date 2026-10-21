"""Regression tests for the designer and protected evidence boundaries."""
from __future__ import annotations

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from zyme_designer.cli import main
from zyme_designer.core import BIRC5_SCAFFOLD, enumerate_candidates, score_candidates, sequence_hash
from zyme_designer.io import resolve_target, write_csv


class CoreTests(unittest.TestCase):
    def test_original_birc5_sequences(self):
        rows, info = enumerate_candidates("AUUGCAAAGGAAACCA", scaffold=BIRC5_SCAFFOLD)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual((row["start_1based"], row["end_1based"], row["cut_after_1based"]), (1, 16, 8))
        self.assertEqual(row["left_arm_dna"], "TTGCAAT")
        self.assertEqual(row["right_arm_dna"], "TGGTTTC")
        self.assertEqual(row["left_partzyme"], "TCAACATCAGTGAGACGAATTGCAAT")
        self.assertEqual(row["right_partzyme"], "TGGTTTCTCTCAGCCTGATAAGCTA")

    def test_changed_substrate_changes_real_arms(self):
        a, _ = enumerate_candidates("AUUGCAAAGGAAACCA")
        b, _ = enumerate_candidates("AACAAUAAGAAGAAAG")
        self.assertNotEqual(a[0]["left_arm_dna"], b[0]["left_arm_dna"])
        self.assertNotEqual(a[0]["right_arm_dna"], b[0]["right_arm_dna"])
        self.assertTrue(a[0]["arm_complementarity_checked"] and b[0]["arm_complementarity_checked"])
        self.assertIsNone(a[0]["left_partzyme"])

    def test_asymmetric_and_ambiguous_windows(self):
        rows, _ = enumerate_candidates("AAAAAGCCCCCC", (4,), (5,))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["target_rna"], "AAAAAGCCCCC")
        rows, audit = enumerate_candidates("AAANAGCCCCC", (4,), (5,))
        self.assertEqual(rows, [])
        self.assertEqual(audit["ambiguous_candidate_windows_skipped"], 1)

    def test_missing_offtarget_not_zero(self):
        rows, _ = enumerate_candidates("AUUGCAAAGGAAACCA")
        score_candidates(rows)
        self.assertEqual(rows[0]["off_target_status"], "not_screened")
        self.assertNotIn("specificity", rows[0]["score_components"])
        self.assertIsNone(rows[0]["pareto_front"])

    def test_zero_accessibility_not_missing(self):
        rows, _ = enumerate_candidates("AUUGCAAAGGAAACCA")
        rows[0]["p_unpaired_target"] = 0.
        score_candidates(rows)
        self.assertEqual(rows[0]["score_components"]["accessibility"], 0.)
        self.assertEqual(rows[0]["ranking_stage"], "accessibility_only_not_screened")

    def test_pareto_and_weight_audit(self):
        rows, _ = enumerate_candidates("AUUGCAAAGGAAACCAUUUUAGCCCCCCC")
        for i, row in enumerate(rows):
            row.update(p_unpaired_target=.9 if i == 0 else .1, off_target_status="screened",
                       left_binding_dg_kcal_mol=-5., right_binding_dg_kcal_mol=-5.,
                       other_gene_exact_ag=0 if i == 0 else 1)
        score_candidates(rows)
        self.assertTrue(rows[0]["pareto_front"])
        self.assertFalse(rows[-1]["pareto_front"])
        self.assertGreater(rows[0]["ranking_score"], rows[-1]["ranking_score"])


class CLITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.target = self.root / "target.fasta"
        self.target.write_text(">NM_001168.3 gene=BIRC5\nAUUGCAAAGGAAACCA\n", encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def call(self, *extra):
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            return main(["design", "--target-fasta", str(self.target), *map(str, extra)])

    def test_cli_without_rnaplfold_marks_unknown(self):
        out = self.root / "result"
        self.call("--out", out)
        result = json.loads((out / "result.json").read_text(encoding="utf-8"))
        self.assertEqual(len(result["candidates"]), 1)
        self.assertIsNone(result["candidates"][0]["p_unpaired_target"])
        self.assertEqual(result["candidates"][0]["off_target_status"], "not_screened")
        self.assertEqual(json.loads((out / "progress.json").read_text())["status"], "complete")

    def test_existing_result_never_modified(self):
        out = self.root / "result"
        self.call("--out", out)
        before = {p.relative_to(out).as_posix(): p.read_bytes() for p in out.rglob("*") if p.is_file()}
        with self.assertRaises(SystemExit):
            self.call("--out", out)
        self.assertEqual(before, {p.relative_to(out).as_posix(): p.read_bytes() for p in out.rglob("*") if p.is_file()})

    def test_screening_distinguishes_central_motif(self):
        reference = self.root / "reference.fasta"
        reference.write_text(self.target.read_text() + ">OTHER gene=OTHER\nAUUGCAACCGAAACCA\n", encoding="utf-8")
        out = self.root / "result"
        self.call("--out", out, "--transcriptome", reference, "--mismatches", 0, "--screen-top", 0)
        data = json.loads((out / "result.json").read_text(encoding="utf-8"))
        row = data["candidates"][0]
        self.assertEqual(row["other_gene_non_ag"], 1)
        self.assertEqual(row["other_gene_exact_ag"], 0)
        self.assertFalse(data["metadata"]["offtarget"]["central_motif_counted_in_arm_mismatches"])

    def test_import_rejects_sequence_hash_mismatch(self):
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"sequence_sha256_rna": "bad"}))
        with self.assertRaises(SystemExit):
            self.call("--out", self.root / "result", "--site-table", self.root / "absent.csv",
                      "--accessibility-manifest", manifest)

    def test_site_probability_is_checked_against_sequence(self):
        from zyme_designer.io import file_hash
        table = self.root / "table.csv"
        write_csv(table, [{"start_1based": 1, "end_1based": 16, "target_rna": "AUUGCAAAGGAAACCA", "p_unpaired_target": .1}])
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"sequence_sha256_rna": sequence_hash("AUUGCAAAGGAAACCA"),
            "site_table_sha256": file_hash(table), "parameters": {"T": 37}}))
        out = self.root / "result"
        self.call("--out", out, "--site-table", table, "--accessibility-manifest", manifest)
        result = json.loads((out / "result.json").read_text(encoding="utf-8"))
        self.assertEqual(result["candidates"][0]["p_unpaired_target"], .1)

    def test_identifier_does_not_select_isoform_silently(self):
        other = self.root / "multi.fasta"
        other.write_text(self.target.read_text() + ">NM_001012271.2 gene=BIRC5\nAUUGCAAAGGAAACCA\n")
        with self.assertRaises(ValueError):
            resolve_target(other, "BIRC5")
        self.assertEqual(resolve_target(other, "NM_001168.3")["accession"], "NM_001168.3")

    def test_csv_formula_escape(self):
        path = self.root / "safe.csv"
        write_csv(path, [{"header": "=BAD()", "energy": -5.1}])
        with path.open(encoding="utf-8-sig", newline="") as f:
            row = next(csv.DictReader(f))
        self.assertEqual(row["header"], "'=BAD()")
        self.assertEqual(float(row["energy"]), -5.1)


if __name__ == "__main__":
    unittest.main()
