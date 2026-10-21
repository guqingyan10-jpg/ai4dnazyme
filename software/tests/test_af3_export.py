from pathlib import Path
import json
import tempfile
import unittest

from zyme_designer.af3_export import candidate_inputs, export_candidates
from zyme_designer.core import BIRC5_SCAFFOLD, enumerate_candidates


class AlphaFoldExportTests(unittest.TestCase):
    sequence = "GGGAUUGCAAAGGAAACCAUUU"

    def example(self, scaffold=None):
        rows, _ = enumerate_candidates(self.sequence, scaffold=scaffold)
        target = {"sequence": self.sequence, "accession": "TEST.1"}
        return rows, target

    def test_arms_only_is_not_a_complete_zyme(self):
        rows, target = self.example()
        local, server, manifest = candidate_inputs(rows[0], target)
        self.assertEqual(local["dialect"], "alphafold3")
        self.assertEqual(server["dialect"], "alphafoldserver")
        self.assertEqual(len(local["sequences"]), 3)
        self.assertEqual(manifest["scope"], "binding_arms_without_catalytic_core")
        self.assertFalse(manifest["prediction_executed"])
        self.assertFalse(manifest["submitted_to_any_service"])

    def test_complete_scaffold_and_rna_trigger(self):
        rows, target = self.example(BIRC5_SCAFFOLD)
        local, server, manifest = candidate_inputs(rows[0], target, BIRC5_SCAFFOLD)
        self.assertEqual(len(local["sequences"]), 4)
        self.assertEqual(local["sequences"][0]["dna"]["sequence"], rows[0]["left_partzyme"])
        self.assertEqual(local["sequences"][3]["rna"]["sequence"], BIRC5_SCAFFOLD["trigger_rna"])
        self.assertEqual(manifest["chains"][3]["role"], "trigger")

    def test_flank_clip_and_transcript_mapping(self):
        rows, target = self.example()
        local, _, manifest = candidate_inputs(rows[0], target, flank=100)
        self.assertEqual(local["sequences"][2]["rna"]["sequence"], self.sequence)
        self.assertEqual(manifest["RNA_chain_transcript_start_1based"], 1)
        self.assertEqual(manifest["RNA_chain_cut_after_1based"], rows[0]["cut_after_1based"])

    def test_explicit_calcium_is_not_a_concentration(self):
        rows, target = self.example()
        local, server, manifest = candidate_inputs(rows[0], target, calcium_count=2)
        self.assertEqual(local["sequences"][-1]["ligand"]["ccdCodes"], ["CA"])
        self.assertEqual(local["sequences"][-1]["ligand"]["id"], ["E", "F"])
        self.assertEqual(server["sequences"][-1]["ion"]["count"], 2)
        self.assertEqual(manifest["modeled_Ca_ion_count"], 2)
        with self.assertRaises(ValueError):
            candidate_inputs(rows[0], target, calcium_count=100)

    def test_server_top_level_is_a_list_and_no_overwrite(self):
        rows, target = self.example()
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "export"
            export_candidates(rows, target, out)
            job = json.loads((out / rows[0]["candidate_id"] / "server_input.json").read_text())
            self.assertIsInstance(job, list)
            with self.assertRaises(FileExistsError):
                export_candidates(rows, target, out)


if __name__ == "__main__":
    unittest.main()
