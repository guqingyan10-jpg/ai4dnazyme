import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import subprocess

from zyme_designer.accessibility import parse_lunp, region_metrics, run_rnaplfold
from zyme_designer.structure import pdb_site_metrics
from zyme_designer.thermo import binding_features, complementarity_risk, hybrid_reference, RNA_DNA_NN


def _atom(serial, residue, base="A", chain="A", x=0.0, bfactor=85.0):
    return (f"ATOM  {serial:5d}  C1' {base:>3s} {chain}{residue:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{bfactor:6.2f}           C\n")


class ThermodynamicsTests(unittest.TestCase):
    def test_all_16_directional_hybrid_steps(self):
        self.assertEqual(len(RNA_DNA_NN), 16)
        # AC and CA differ for RNA/DNA; do not use DNA/DNA strand symmetry.
        self.assertNotEqual(hybrid_reference("AC")["binding_dg_kcal_mol"],
                            hybrid_reference("CA")["binding_dg_kcal_mol"])

    def test_birc5_arms_finite_and_reference_disclosed(self):
        result = binding_features("AUUGCAA", "GAAACCA", 37, 300)
        for key in ["left_binding_dg_kcal_mol", "right_binding_dg_kcal_mol", "left_tm_c", "right_tm_c"]:
            self.assertTrue(math.isfinite(result[key]))
        self.assertEqual(result["binding_reference_sodium_mM"], 1000)
        self.assertEqual(result["binding_input_sodium_mM"], 300)
        self.assertTrue(any("salt correction" in warning for warning in result["warnings"]))
        # Sodium input is metadata, not an invented Ca/Na correction.
        alternative = binding_features("AUUGCAA", "GAAACCA", 37, 100)
        self.assertEqual(result["left_binding_dg_kcal_mol"], alternative["left_binding_dg_kcal_mol"])
        json.dumps(result, allow_nan=False)

    def test_manual_nn_and_temperature_units(self):
        result = hybrid_reference("AG")
        self.assertAlmostEqual(result["enthalpy_kcal_mol"], 1.9 - 9.1)
        self.assertAlmostEqual(result["entropy_cal_mol_K"], -3.9 - 23.5)
        self.assertAlmostEqual(result["binding_dg_kcal_mol"], -7.2 + 310.15 * 27.4 / 1000)
        self.assertEqual(hybrid_reference("ATUG") , hybrid_reference("AUUG"))
        with self.assertRaises(ValueError):
            hybrid_reference("ANNN")

    def test_complementarity_diagnostics(self):
        risk = complementarity_risk("AAAACCC", "GGGTTTT")
        self.assertEqual(risk["cross_arm_complementary_run_nt"], 7)
        self.assertEqual(complementarity_risk("ACGT", "TTTT")["left_arm_self_dimer_run_nt"], 4)
        hairpin = complementarity_risk("AAAGGGTTT", "ACACAC")
        self.assertGreaterEqual(hairpin["left_arm_hairpin_run_nt"], 3)


class AccessibilityTests(unittest.TestCase):
    def test_segment_end_convention_birc5(self):
        table = {(416, 16): .1, (407, 7): .3, (416, 7): .4, (409, 2): .5,
                 (401, 16): .9}  # A wrong forward-start implementation would use .9.
        values = region_metrics(table, 401, 7, 7)
        self.assertEqual(values["full_window_p_unpaired"], .1)
        self.assertEqual(values["left_arm_p_unpaired"], .3)
        self.assertEqual(values["ag_p_unpaired"], .5)
        self.assertEqual(values["full_window_lunp_key"], [416, 16])
        self.assertAlmostEqual(values["full_window_opening_dg_kcal_mol"],
                               -0.00198720425864083 * 310.15 * math.log(.1))

    def test_na_zero_missing_not_confused(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample_lunp"
            path.write_text("#unpaired probabilities\n1 0.8 n/a\n2 0 0.25\n3 NA 0.1\n", encoding="utf8")
            table = parse_lunp(path)
        self.assertIsNone(table[(1, 2)])
        self.assertIsNone(table[(3, 1)])
        self.assertEqual(table[(2, 1)], 0)
        values = region_metrics({(4, 4): 0, (1, 1): None, (4, 1): .5}, 1, 1, 1)
        self.assertEqual(values["full_window_p_unpaired"], 0)
        self.assertTrue(math.isfinite(values["full_window_opening_dg_kcal_mol"]))
        self.assertIsNone(values["ag_p_unpaired"])
        self.assertIn("ag", values["missing_probability_regions"])
        self.assertIn("full_window", values["probability_zero_regions"])
        json.dumps(values, allow_nan=False)

    def test_invalid_probabilities_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad_lunp"
            for text in ["1 1.25", "1 -0.1", "# opening energies\n1 0.5", "1 0.5\n1 0.4"]:
                path.write_text(text, encoding="utf8")
                with self.assertRaises(ValueError):
                    parse_lunp(path)

    def test_runner_audits_input_and_uses_no_shell(self):
        def fake_run(command, **kwargs):
            self.assertFalse(kwargs.get("shell", False))
            if command[-1] == "--version":
                return subprocess.CompletedProcess(command, 0, "RNAplfold 2.7.2\n", "")
            fasta = kwargs["input"]
            identifier = fasta.splitlines()[0][1:]
            sequence = fasta.splitlines()[1]
            output = Path(kwargs["cwd"]) / (identifier + "_lunp")
            output.write_text("#unpaired probabilities\n" + "\n".join(
                f"{i} 0.5 n/a" for i in range(1, len(sequence) + 1)), encoding="ascii")
            return subprocess.CompletedProcess(command, 0, "", "")
        with tempfile.TemporaryDirectory() as temporary, patch(
                "zyme_designer.accessibility.shutil.which", return_value="RNAplfold"), patch(
                "zyme_designer.accessibility.subprocess.run", side_effect=fake_run):
            lunp, audit = run_rnaplfold("ACGT", temporary, W=20, L=15, u=2)
            self.assertTrue(lunp.exists())
            self.assertEqual(audit["sequence_length_nt"], 4)
            self.assertEqual(audit["salt_status"], "program_default_not_project_buffer")
            self.assertEqual(len(audit["lunp_sha256"]), 64)
            self.assertTrue((Path(temporary) / "lunp_manifest.json").exists())
            with self.assertRaises(FileExistsError):
                run_rnaplfold("ACGT", temporary, W=20, L=15, u=2)

    def test_salt_unsupported_not_silently_applied(self):
        fake = subprocess.CompletedProcess(["RNAplfold"], 0, "old ViennaRNA without salt\n", "")
        with tempfile.TemporaryDirectory() as temporary, patch(
                "zyme_designer.accessibility.shutil.which", return_value="RNAplfold"), patch(
                "zyme_designer.accessibility.subprocess.run", return_value=fake):
            with self.assertRaises(RuntimeError):
                run_rnaplfold("ACGU", temporary, W=20, L=15, u=2, salt_molar=.3)

    def test_invalid_region_and_sequence_inputs(self):
        for value in [0, 1.5, True]:
            with self.assertRaises(ValueError):
                region_metrics({}, value, 7, 7)
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                run_rnaplfold("ACGU;rm bad", temporary)


class StructureTests(unittest.TestCase):
    def test_bfactor_not_automatically_confidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "model.pdb"
            path.write_text(_atom(1, 1) + _atom(2, 2, "G", x=3) + "END\n", encoding="ascii")
            default = pdb_site_metrics(path, 1, 3)
            explicit = pdb_site_metrics(path, 1, 2, confidence_source="plddt", transcript_sequence="AG")
        self.assertIsNone(default["target_mean_plddt"])
        self.assertEqual(default["target_coordinate_coverage_fraction"], 2 / 3)
        self.assertEqual(explicit["target_mean_plddt"], 85)
        self.assertEqual(explicit["mapping_status"], "complete_transcript_sequence_verified_identity_1based")

    def test_fragment_needs_mapping_and_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fragment.pdb"
            path.write_text(_atom(1, 401) + _atom(2, 402, "G"), encoding="ascii")
            with self.assertRaises(ValueError):
                pdb_site_metrics(path, 401, 402)
            mapped = pdb_site_metrics(path, 401, 402, residue_map={
                ("A", 401, ""): 401, ("A", 402, ""): 402})
            self.assertEqual(mapped["target_coordinate_coverage_fraction"], 1)
            continuous = Path(temporary) / "continuous.pdb"
            continuous.write_text(_atom(1, 1) + _atom(2, 2, "G"), encoding="ascii")
            with self.assertRaises(ValueError):
                pdb_site_metrics(continuous, 1, 2, transcript_sequence="AA")

    def test_exact_nonlocal_neighbor_geometry(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "neighbor.pdb"
            path.write_text("".join(_atom(i, i, x=x) for i, x in enumerate(
                (0, 1, 2, 7.9, 8.1, 100), 1)), encoding="ascii")
            mapped = pdb_site_metrics(path, 1, 1, transcript_sequence="AAAAAA")
        # Residues 2/3 are local sequence neighbors. 4 is within 8 A, 5 is not.
        self.assertEqual(mapped["residues"][0]["nonlocal_neighbor_residue_count"], 1)


if __name__ == "__main__":
    unittest.main()
