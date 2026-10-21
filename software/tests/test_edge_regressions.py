"""Independent regressions for import identity, unknown values and safe exports."""

from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from zyme_designer.cli import main
from zyme_designer.core import annotate_header, score_candidates, sequence_hash
from zyme_designer.io import file_hash, resolve_target, write_csv


class EdgeRegressionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.target = self.root / "target.fasta"

    def tearDown(self):
        self.temporary.cleanup()

    def target_file(self, sequence, header="TARGET.1 gene=TARGETGENE"):
        self.target.write_text(f">{header}\n{sequence}\n", encoding="utf-8")

    def call(self, name, *arguments):
        output = self.root / name
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            main(["design", "--target-fasta", str(self.target), "--out", str(output), *map(str, arguments)])
        return json.loads((output / "result.json").read_text(encoding="utf-8"))

    def import_files(self, sequence, rows):
        table = self.root / "site_table.csv"
        write_csv(table, rows)
        manifest = self.root / "accessibility_manifest.json"
        manifest.write_text(json.dumps({
            "sequence_sha256_rna": sequence_hash(sequence),
            "site_table_sha256": file_hash(table), "parameters": {"T": 37},
        }), encoding="utf-8")
        return ["--site-table", table, "--accessibility-manifest", manifest]

    def lunp_files(self, sequence, include_hash=False):
        lunp = self.root / "target_lunp"
        lines = ["# unpaired probabilities\n"]
        for end in range(1, len(sequence) + 1):
            tokens = ["0.2" if length <= end else "n/a" for length in range(1, 17)]
            lines.append(str(end) + "\t" + "\t".join(tokens) + "\n")
        lunp.write_text("".join(lines), encoding="utf-8")
        manifest = self.root / "lunp_manifest.json"
        data = {"sequence_sha256_rna": sequence_hash(sequence), "parameters": {"T": 37}}
        if include_hash:
            data["lunp_sha256"] = file_hash(lunp)
        manifest.write_text(json.dumps(data), encoding="utf-8")
        return ["--lunp", lunp, "--accessibility-manifest", manifest]

    def test_same_window_two_ag_splits_imports_to_correct_candidates(self):
        sequence = "AAAAAGAGAAAA"
        self.target_file(sequence)
        table_rows = [
            {"candidate_id": "AG5_L4_R6", "start_1based": 1, "end_1based": 12,
             "cut_after_1based": 5, "left_length": 4, "right_length": 6,
             "target_rna": sequence, "p_unpaired_target": .2,
             "p_unpaired_left_arm": .4, "p_unpaired_right_arm": .6,
             "p_unpaired_cut_ag": .21},
            {"candidate_id": "AG7_L6_R4", "start_1based": 1, "end_1based": 12,
             "cut_after_1based": 7, "left_length": 6, "right_length": 4,
             "target_rna": sequence, "p_unpaired_target": .2,
             "p_unpaired_left_arm": .7, "p_unpaired_right_arm": .8,
             "p_unpaired_cut_ag": .22},
        ]
        data = self.call("asymmetric", "--left-lengths", "4,6", "--right-lengths", "4,6",
                         *self.import_files(sequence, table_rows))
        by_id = {row["candidate_id"]: row for row in data["candidates"]}
        self.assertEqual(len(by_id), 4)
        first, second = by_id["AG5_L4_R6"], by_id["AG7_L6_R4"]
        self.assertEqual((first["p_unpaired_left_arm"], first["p_unpaired_right_arm"], first["p_unpaired_center"]), (.4, .6, .21))
        self.assertEqual((second["p_unpaired_left_arm"], second["p_unpaired_right_arm"], second["p_unpaired_center"]), (.7, .8, .22))

    def test_same_window_cut_only_disambiguates_without_candidate_ids(self):
        sequence = "AAAAAGAGAAAA"
        self.target_file(sequence)
        table_rows = [
            {"start_1based": 1, "end_1based": 12, "cut_after_1based": 5,
             "target_rna": sequence, "p_unpaired_target": .2,
             "p_unpaired_left_arm": .4, "p_unpaired_right_arm": .6},
            {"start_1based": 1, "end_1based": 12, "cut_after_1based": 7,
             "target_rna": sequence, "p_unpaired_target": .2,
             "p_unpaired_left_arm": .7, "p_unpaired_right_arm": .8},
        ]
        data = self.call("cut_only", "--left-lengths", "4,6", "--right-lengths", "4,6",
                         *self.import_files(sequence, table_rows))
        by_id = {row["candidate_id"]: row for row in data["candidates"]}
        self.assertEqual(by_id["AG5_L4_R6"]["p_unpaired_left_arm"], .4)
        self.assertEqual(by_id["AG7_L6_R4"]["p_unpaired_left_arm"], .7)

    def test_legacy_seven_base_probabilities_not_applied_to_other_arms(self):
        sequence = "AAAAAGCCCCCC"
        self.target_file(sequence)
        table_rows = [{"start_1based": 1, "end_1based": 12,
                       "cut_after_1based": 5, "left_length": 4, "right_length": 6,
                       "target_rna": sequence, "p_unpaired_target": .25,
                       "p_unpaired_left7": .55, "p_unpaired_right7": .66}]
        data = self.call("legacy_arms", "--left-lengths", "4", "--right-lengths", "6",
                         *self.import_files(sequence, table_rows))
        row = data["candidates"][0]
        self.assertEqual(row["p_unpaired_target"], .25)
        self.assertIsNone(row["p_unpaired_left_arm"])
        self.assertIsNone(row["p_unpaired_right_arm"])

    def test_legacy_full16_probability_not_applied_to_twelve_base_window(self):
        sequence = "AAAAAGCCCCCC"
        self.target_file(sequence)
        table_rows = [{"start_1based": 1, "end_1based": 12,
                       "cut_after_1based": 5, "left_length": 4, "right_length": 6,
                       "target_rna": sequence, "p_unpaired_full16": .55}]
        data = self.call("legacy_full", "--left-lengths", "4", "--right-lengths", "6",
                         *self.import_files(sequence, table_rows))
        self.assertIsNone(data["candidates"][0]["p_unpaired_target"])
        self.assertFalse(data["metadata"]["accessibility"]["all_candidate_probabilities_available"])

    def test_ensembl_symbol_preferred_over_gene_id(self):
        header = "ENST000001.1 cdna gene:ENSG000001.2 gene_symbol:BIRC5"
        self.assertEqual(annotate_header(header), ("ENST000001.1", "BIRC5"))
        self.target_file("AUUGCAAAGGAAACCA", header)
        data = self.call("ensembl", "--gene", "BIRC5")
        self.assertEqual(data["target"]["gene"], "BIRC5")

    def test_gencode_pipe_gene_symbol_is_resolvable(self):
        header = "ENST000001.1|ENSG000001.2|OTTHUMG1|OTTHUMT1|BIRC5-201|BIRC5|16|protein_coding|"
        self.assertEqual(annotate_header(header), ("ENST000001.1", "BIRC5"))
        self.target_file("AUUGCAAAGGAAACCA", header)
        self.assertEqual(resolve_target(self.target, "BIRC5")["accession"], "ENST000001.1")

    def test_lunp_without_output_hash_requires_explicit_trust(self):
        sequence = "AUUGCAAAGGAAACCA"
        self.target_file(sequence)
        with self.assertRaises(SystemExit) as error:
            self.call("unverified_lunp", *self.lunp_files(sequence))
        self.assertEqual(error.exception.code, 2)

    def test_explicit_trust_for_missing_lunp_hash_has_provenance_warning(self):
        sequence = "AUUGCAAAGGAAACCA"
        self.target_file(sequence)
        data = self.call("trusted_lunp", *self.lunp_files(sequence), "--trust-lunp-sequence")
        audit = data["metadata"]["accessibility"]
        self.assertEqual(data["candidates"][0]["p_unpaired_target"], .2)
        evidence = json.dumps(audit).lower()
        self.assertTrue("user_asserted" in evidence or ("warnings" in audit and ("hash" in evidence or "trust" in evidence)))

    def test_blank_target_probability_not_audited_as_available(self):
        sequence = "AUUGCAAAGGAAACCA"
        self.target_file(sequence)
        table_rows = [{"start_1based": 1, "end_1based": 16, "target_rna": sequence,
                       "p_unpaired_target": "", "p_unpaired_left_arm": "", "p_unpaired_right_arm": ""}]
        data = self.call("blank_probabilities", *self.import_files(sequence, table_rows))
        self.assertIsNone(data["candidates"][0]["p_unpaired_target"])
        self.assertFalse(data["metadata"]["accessibility"]["all_candidate_probabilities_available"])

    def test_incomplete_pareto_vector_does_not_crash_or_impute_zero(self):
        common = {"off_target_status": "screened", "left_length": 7, "right_length": 7}
        rows = [
            {**common, "candidate_id": "complete", "cut_after_1based": 1,
             "p_unpaired_target": .2, "left_binding_dg_kcal_mol": -4., "right_binding_dg_kcal_mol": -5.},
            {**common, "candidate_id": "incomplete", "cut_after_1based": 2,
             "p_unpaired_target": .3},
        ]
        ranked = score_candidates(rows)
        self.assertEqual(len(ranked), 2)
        incomplete = next(row for row in ranked if row["candidate_id"] == "incomplete")
        self.assertNotIn("binding", incomplete["score_components"])
        self.assertLess(incomplete["score_weight_coverage"], 1.)

    def test_formula_escape_includes_leading_whitespace(self):
        path = self.root / "safe.csv"
        write_csv(path, [{"header": " \t=UNTRUSTED()", "numeric": -5.1}])
        with path.open(encoding="utf-8-sig", newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertTrue(row["header"].startswith("'"))
        self.assertEqual(float(row["numeric"]), -5.1)


if __name__ == "__main__":
    unittest.main()
