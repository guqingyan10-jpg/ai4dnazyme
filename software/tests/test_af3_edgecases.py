"""Independent, offline QA of export formats and sequence/coordinate boundaries.

Format basis: official AlphaFold Server README and AF3 docs/input.md.
These tests do not submit predictions, install AF3, or validate catalytic activity.
"""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from zyme_designer.af3_export import candidate_inputs, export_candidates
from zyme_designer.core import BIRC5_SCAFFOLD, enumerate_candidates, sequence_hash


class AF3EdgeCases(unittest.TestCase):
    sequence = "CCCGGGAUUGCAAAGGAAACCAUUUGGG"

    def fixture(self, scaffold=None):
        rows, _ = enumerate_candidates(self.sequence, scaffold=scaffold)
        return next(row for row in rows if row["target_rna"] == "AUUGCAAAGGAAACCA"), {
            "sequence": self.sequence, "accession": "OFFLINE_TEST.1"}

    def test_json_schema_distinctions_and_actual_roles(self):
        candidate, target = self.fixture(BIRC5_SCAFFOLD)
        local, server, mapping = candidate_inputs(candidate, target, BIRC5_SCAFFOLD)
        self.assertEqual(local["dialect"], "alphafold3")
        self.assertEqual(local["version"], 1)
        self.assertTrue(local["modelSeeds"])
        self.assertEqual(server["dialect"], "alphafoldserver")
        self.assertEqual(server["version"], 1)
        self.assertEqual(server["modelSeeds"], [])
        self.assertEqual([next(iter(entity)) for entity in local["sequences"]],
                         ["dna", "dna", "rna", "rna"])
        self.assertEqual([next(iter(entity)) for entity in server["sequences"]],
                         ["dnaSequence", "dnaSequence", "rnaSequence", "rnaSequence"])
        self.assertEqual([chain["role"] for chain in mapping["chains"]],
                         ["left_partzyme", "right_partzyme", "RNA_substrate_context", "trigger"])
        for chain in mapping["chains"]:
            self.assertTrue(set(chain["sequence"]) <= set("ACGT" if chain["type"] == "dna" else "ACGU"))
        self.assertFalse(mapping["prediction_executed"])
        self.assertFalse(mapping["coordinate_file_generated"])
        self.assertFalse(mapping["submitted_to_any_service"])

    def test_context_mapping_round_trip_and_full_transcript_hash(self):
        candidate, target = self.fixture()
        local, _, mapping = candidate_inputs(candidate, target, flank=2)
        context = local["sequences"][2]["rna"]["sequence"]
        self.assertEqual(mapping["full_transcript_sequence_sha256_rna"], sequence_hash(self.sequence))
        start = mapping["RNA_chain_site_start_1based"]
        end = mapping["RNA_chain_site_end_1based"]
        self.assertEqual(context[start - 1:end], candidate["target_rna"])
        cut = mapping["RNA_chain_cut_after_1based"]
        self.assertEqual(context[cut - 1:cut + 1], "AG")
        self.assertEqual(mapping["RNA_chain_transcript_start_1based"] + cut - 1,
                         candidate["cut_after_1based"])

    def test_ambiguity_outside_context_is_not_guessed(self):
        candidate, target = self.fixture()
        target["sequence"] = "N" + target["sequence"][1:]
        local, _, mapping = candidate_inputs(candidate, target)
        self.assertNotIn("N", local["sequences"][2]["rna"]["sequence"])
        self.assertEqual(mapping["full_transcript_sequence_sha256_rna"],
                         sequence_hash(target["sequence"]))
        with self.assertRaises(ValueError):
            candidate_inputs(candidate, target, flank=100)

    def test_invalid_counts_are_not_concentrations_or_booleans(self):
        candidate, target = self.fixture()
        for count in (-1, 11, 100, 0.1, "2", True):
            with self.subTest(count=count), self.assertRaises(ValueError):
                candidate_inputs(candidate, target, calcium_count=count)
        for flank in (-1, 0.1, "2", True):
            with self.subTest(flank=flank), self.assertRaises(ValueError):
                candidate_inputs(candidate, target, flank=flank)

    def test_calcium_copies_are_independent_single_ccd_entities(self):
        candidate, target = self.fixture()
        local, server, mapping = candidate_inputs(candidate, target, calcium_count=3)
        calcium = [entity["ligand"] for entity in local["sequences"] if "ligand" in entity]
        identifiers = []
        for entity in calcium:
            self.assertEqual(entity["ccdCodes"], ["CA"],
                             "CCD list describes a multi-component ligand, not copy count")
            ids = entity["id"]
            identifiers.extend(ids if isinstance(ids, list) else [ids])
        self.assertEqual(len(identifiers), 3)
        self.assertEqual(len(set(identifiers)), 3)
        self.assertEqual(server["sequences"][-1], {"ion": {"ion": "CA", "count": 3}})
        self.assertEqual(mapping["modeled_Ca_ion_count"], 3)
        self.assertTrue(any("NOT a molar" in warning for warning in mapping["warnings"]))

    def test_invalid_dna_and_coordinate_claims_rejected(self):
        candidate, target = self.fixture()
        for field, value in (("left_arm_dna", "ATGN"), ("right_arm_dna", "ATGU"),
                             ("start_1based", 0), ("end_1based", len(self.sequence) + 1),
                             ("cut_after_1based", len(self.sequence))):
            altered = deepcopy(candidate)
            altered[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                candidate_inputs(altered, target)
        inconsistent = deepcopy(candidate)
        inconsistent["target_rna"] = "C" * len(candidate["target_rna"])
        with self.assertRaises(ValueError):
            candidate_inputs(inconsistent, target)

    def test_server_file_lists_and_zero_calcium_no_added_entities(self):
        candidate, target = self.fixture()
        local, _, _ = candidate_inputs(candidate, target)
        self.assertEqual(len(local["sequences"]), 3)
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "export"
            export_candidates([candidate], target, out, top=1)
            combined = json.loads((out / "server_jobs.json").read_text(encoding="utf-8"))
            single = json.loads((out / candidate["candidate_id"] / "server_input.json").read_text(encoding="utf-8"))
            self.assertEqual(combined, single)
            self.assertIsInstance(single, list)
            self.assertEqual(len(single), 1)
            with self.assertRaises(FileExistsError):
                export_candidates([candidate], target, out, top=1)


if __name__ == "__main__":
    unittest.main()
