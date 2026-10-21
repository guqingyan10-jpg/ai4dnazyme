"""Generic recognition-policy and bilateral arm trade-off regressions."""

from pathlib import Path
import tempfile
import unittest

from zyme_designer.core import (
    BIRC5_SCAFFOLD, enumerate_candidates, motif_matches,
    normalize_length_policy, normalize_motif_pattern, ranking_audit, score_candidates,
)
from zyme_designer.offtarget import scan_transcriptome


def scored_row(identifier, length=8, dg=-8., right_length=None):
    return {"candidate_id": identifier, "cut_after_1based": 10,
            "left_length": length, "right_length": right_length or length,
            "left_binding_dg_kcal_mol": dg, "right_binding_dg_kcal_mol": dg,
            "p_unpaired_target": .5, "off_target_status": "screened",
            "other_gene_exact_ag": 0, "other_gene_1mm_ag": 0,
            "other_gene_2mm_ag": 0, "other_gene_non_ag": 0}


class GenericRecognitionTests(unittest.TestCase):
    def test_default_ag_sequences_and_ids_remain_compatible(self):
        rows, audit = enumerate_candidates("AUUGCAAAGGAAACCA", scaffold=BIRC5_SCAFFOLD)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["candidate_id"], "AG8_L7_R7")
        self.assertEqual(rows[0]["motif_pattern"], "AG")
        self.assertEqual(rows[0]["left_partzyme"], "TCAACATCAGTGAGACGAATTGCAAT")
        self.assertEqual(audit["selected_cut_positions_1based"], [8])

    def test_nn_enumerates_all_canonical_junctions_with_complete_flanks(self):
        sequence = "AUGCCAUUGACCGUUAA"
        rows, audit = enumerate_candidates(sequence, (4,), (5,), motif_pattern="NN")
        expected = list(range(5, len(sequence) - 5))
        self.assertEqual([row["cut_after_1based"] for row in rows], expected)
        for row in rows:
            start = row["start_1based"] - 1
            self.assertEqual(row["target_rna"], sequence[start:row["end_1based"]])
            self.assertTrue(row["arm_complementarity_checked"])
            self.assertIsNone(row["left_partzyme"])
            self.assertEqual(row["scaffold_core_motif_policy"], "no_core_provided")
            self.assertFalse(row["motif_compatibility_is_cleavage_proof"])
        self.assertEqual(audit["compatible_junctions_in_target"], len(sequence) - 1)

    def test_ry_policy_and_dna_u_normalization(self):
        self.assertTrue(motif_matches("AU", "RY"))
        self.assertTrue(motif_matches("GC", "RY"))
        self.assertFalse(motif_matches("AG", "RY"))
        self.assertFalse(motif_matches("AN", "NN"))
        self.assertEqual(normalize_motif_pattern("RT"), "RU")
        rows, _ = enumerate_candidates("AAAAGTCCCCCC", (4,), (4,), motif_pattern="RY")
        self.assertEqual(rows[0]["motif"], "GU")

    def test_explicit_cut_coordinates_and_invalid_policy_rejected(self):
        sequence = "AAAAAGAGAAAA"
        rows, audit = enumerate_candidates(sequence, (4, 6), (4, 6), motif_pattern="NN", cut_positions=[7, 7])
        self.assertTrue(rows)
        self.assertEqual({r["cut_after_1based"] for r in rows}, {7})
        self.assertEqual(audit["selected_cut_positions_1based"], [7])
        for positions in ([], [0], [len(sequence)], [True], ["5"]):
            with self.subTest(positions=positions), self.assertRaises(ValueError):
                enumerate_candidates(sequence, motif_pattern="NN", cut_positions=positions)
        with self.assertRaises(ValueError):
            enumerate_candidates(sequence, motif_pattern="AG", cut_positions=[6])
        for pattern in ("A", "AGN", "AZ", "", None):
            with self.subTest(pattern=pattern), self.assertRaises(ValueError):
                normalize_motif_pattern(pattern)

    def test_builtin_core_does_not_become_generic_catalyst(self):
        with self.assertRaises(ValueError):
            enumerate_candidates("AUUGCAAAGGAAACCA", scaffold=BIRC5_SCAFFOLD, motif_pattern="NN")
        rows, _ = enumerate_candidates("AUUGCAAAGGAAACCA", motif_pattern="NN")
        self.assertTrue(all(r["left_partzyme"] is None and r["right_partzyme"] is None for r in rows))

    def test_generic_offtarget_compatibility_separate_from_legacy_ag(self):
        candidate = {"candidate_id": "AU5_L4_R4", "target_rna": "AAAAAUCCCC",
                     "left_length": 4, "right_length": 4, "start_1based": 1,
                     "end_1based": 10, "motif_pattern": "RY"}
        text = ">TARGET.1 gene=TARGETGENE\nAAAAAUCCCC\n"
        text += ">OTHER_A gene=A\nAAAAAUCCCC\n>OTHER_B gene=B\nAAAAGCCCCC\n"
        text += ">OTHER_C gene=C\nAAAAAGCCCC\n"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reference.fasta"
            path.write_text(text, encoding="utf-8")
            summaries, hits, audit = scan_transcriptome([candidate], path, "TARGET.1", "TARGETGENE", max_mismatches=0)
        row = summaries[candidate["candidate_id"]]
        self.assertEqual(row["other_gene_exact_ag"], 1)
        self.assertEqual(row["other_gene_non_ag"], 2)
        self.assertEqual(row["other_gene_exact_compatible"], 2)
        self.assertEqual(row["other_gene_noncompatible"], 1)
        self.assertEqual(row["compatible_binding_hits_total"], 3)
        self.assertEqual(row["noncompatible_binding_hits_total"], 1)
        self.assertFalse(next(h for h in hits if h["cut_motif"] == "AG")["compatible_motif"])
        self.assertFalse(audit["motif_compatibility_is_cleavage_proof"])
        score_row = {**scored_row("generic"), **row}
        score_candidates([score_row])
        self.assertAlmostEqual(score_row["score_components"]["specificity"], 1. / (1. + 20. + .05))

    def test_nn_all_central_motifs_are_policy_compatible_not_proven_cleavable(self):
        candidate = {"candidate_id": "NN", "target_rna": "AAAAAUCCCC",
                     "left_length": 4, "right_length": 4, "start_1based": 1,
                     "end_1based": 10, "motif_pattern": "NN"}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reference.fasta"
            path.write_text(">OTHER gene=OTHER\nAAAAAGCCCC\n", encoding="utf-8")
            summaries, hits, _ = scan_transcriptome([candidate], path, "TARGET", "T", max_mismatches=0)
        self.assertEqual(summaries["NN"]["other_gene_exact_compatible"], 1)
        self.assertEqual(summaries["NN"]["other_gene_noncompatible"], 0)
        self.assertFalse(hits[0]["motif_compatibility_is_cleavage_proof"])


class LengthTradeoffTests(unittest.TestCase):
    def test_longer_stronger_is_not_automatically_better(self):
        rows = [scored_row("long", 20, -30.), scored_row("moderate", 8, -8.)]
        score_candidates(rows)
        self.assertEqual(rows[0]["candidate_id"], "moderate")
        self.assertGreater(rows[0]["length_tradeoff_score"], rows[1]["length_tradeoff_score"])
        self.assertGreater(rows[1]["length_penalty"], 0)
        self.assertGreater(rows[1]["strong_binding_penalty"], 0)
        self.assertFalse(rows[1]["length_tradeoff_is_release_probability"])

    def test_both_arms_penalized_not_just_mean_length(self):
        rows = [scored_row("unequal", 4, -8., right_length=12), scored_row("balanced", 8, -8.)]
        score_candidates(rows)
        self.assertEqual(rows[0]["candidate_id"], "balanced")
        self.assertEqual(rows[0]["length_penalty"], 0)
        self.assertGreater(rows[1]["length_penalty"], 0)

    def test_continuous_strong_binding_penalty_increases_with_excess(self):
        rows = [scored_row("at_threshold", 8, -14.), scored_row("slightly_stronger", 8, -15.), scored_row("much_stronger", 8, -20.)]
        score_candidates(rows)
        by_id = {r["candidate_id"]: r for r in rows}
        self.assertEqual(by_id["at_threshold"]["strong_binding_penalty"], 0)
        self.assertGreater(by_id["slightly_stronger"]["strong_binding_penalty"], 0)
        self.assertGreater(by_id["much_stronger"]["strong_binding_penalty"], by_id["slightly_stronger"]["strong_binding_penalty"])

    def test_user_preferences_are_configurable_and_audited(self):
        policy = {"preferred_min": 15, "preferred_max": 17, "practical_min": 14, "practical_max": 18}
        rows = [scored_row("eight", 8, -8.), scored_row("sixteen", 16, -8.)]
        score_candidates(rows, length_policy=policy)
        self.assertEqual(rows[0]["candidate_id"], "sixteen")
        audit = ranking_audit(rows, length_policy=policy)
        self.assertEqual(audit["length_policy"]["preferred_min"], 15)
        self.assertIn("not universal", audit["length_policy_status"])
        self.assertIn("length_tradeoff_score", audit["equations"])

    def test_missing_reference_energy_is_unknown_not_zero_penalty_evidence(self):
        row = scored_row("unknown")
        row.pop("left_binding_dg_kcal_mol")
        score_candidates([row])
        self.assertIsNone(row["strong_binding_penalty"])
        self.assertIsNone(row["strong_binding_penalty_cost"])
        self.assertIn("strong_binding_penalty_not_assessed_missing_reference_energy", row["length_tradeoff_reasons"])

    def test_invalid_policies_are_rejected(self):
        for policy in ({"preferred_min": 3}, {"preferred_scale_nt": 0}, {"practical_min": 8},
                       {"unexpected": 1}, {"strong_binding_scale_kcal_mol": float("nan")}, {"preferred_max": True}):
            with self.subTest(policy=policy), self.assertRaises(ValueError):
                normalize_length_policy(policy)


if __name__ == "__main__":
    unittest.main()
