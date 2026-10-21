import csv
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from zyme_designer.accessibility import parse_lunp, region_metrics
from zyme_designer.vienna_backend import run_vienna_python


def fake_rna(invalid=False, split=False, missing=False):
    class FakeMd:
        def __init__(self):
            for field in ("temperature", "window_size", "max_bp_span", "salt", "dangles", "betaScale",
                          "pf_smooth", "special_hp", "noLP", "noGU", "noGUclosure", "logML", "circ",
                          "gquad", "uniq_ML", "energy_set", "min_loop_size", "compute_bpp", "sfact",
                          "saltMLLower", "saltMLUpper", "saltDPXInit", "saltDPXInitFact", "helical_rise", "backbone_length"):
                setattr(self, field, 0)

    class FakeFc:
        def __init__(self, sequence, md, flags):
            self.sequence, self.md = sequence, md
            assert flags == 16

        def probs_window(self, u, options, callback, data):
            assert options == 8192
            for end in range(1, len(self.sequence) + 1):
                if missing and end == len(self.sequence):
                    continue
                size = min(u, end)
                values = [None] + [.5 ** length for length in range(1, size + 1)]
                values += [None] * (u - size)
                if invalid and end == 4:
                    values[1] = 1.2
                callback(values, size, end, u, 8192 | (1 if split else 15), data)
            return 1

    return SimpleNamespace(__version__="2.7.2", md=FakeMd, fold_compound=FakeFc,
                           OPTION_WINDOW=16, PROBS_WINDOW_UP=8192, ANY_LOOP=15,
                           params_load_RNA_Turner2004=lambda: 1,
                           parameter_set_rna_turner2004="MOCK PARAMETERS, NOT SCIENTIFIC DATA")


class ViennaBackendMockTests(unittest.TestCase):
    def test_end_indexed_lunp_output_audit_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary, patch(
                "zyme_designer.vienna_backend.importlib.import_module", return_value=fake_rna()):
            lunp, audit = run_vienna_python("ACGUAAGUACGU", temporary, W=10, L=8, u=5, T=37, salt_molar=.3)
            table = parse_lunp(lunp)
            self.assertEqual(table[(5, 3)], .125)
            self.assertIsNone(table[(1, 3)])
            self.assertEqual(audit["parameters"]["salt_molar"], .3)
            self.assertEqual(audit["n_probability_rows"], 12)
            self.assertEqual(audit["model_details"]["noGU"], 0)
            self.assertEqual(audit["model_details"]["dangles"], 2)
            self.assertEqual(len(audit["energy_parameters_sha256"]), 64)
            self.assertTrue((Path(temporary) / "lunp_manifest.json").is_file())
            json.dumps(audit, allow_nan=False)
            with self.assertRaises(FileExistsError):
                run_vienna_python("ACGUAAGUACGU", temporary, W=10, L=8, u=5)

    def test_short_sequence_effective_dimensions_recorded(self):
        with tempfile.TemporaryDirectory() as temporary, patch(
                "zyme_designer.vienna_backend.importlib.import_module", return_value=fake_rna()):
            _, audit = run_vienna_python("AAAAAA", temporary, W=200, L=150, u=16)
        self.assertEqual(audit["parameters"]["W"], 6)
        self.assertEqual(audit["parameters"]["L"], 6)
        self.assertEqual(audit["parameters"]["u"], 6)
        self.assertEqual(audit["requested_parameters"]["W"], 200)
        self.assertEqual(audit["parameters"]["salt_molar"], 1.021)

    def test_unknown_version_and_missing_api_rejected(self):
        for module in [fake_rna(), fake_rna()]:
            if not hasattr(self, "_did_version"):
                module.__version__ = "9.0"
                self._did_version = True
            else:
                del module.PROBS_WINDOW_UP
            with tempfile.TemporaryDirectory() as temporary, patch(
                    "zyme_designer.vienna_backend.importlib.import_module", return_value=module):
                with self.assertRaises(RuntimeError):
                    run_vienna_python("ACGUACGU", temporary, W=8, L=6, u=4)
                self.assertFalse(list(Path(temporary).glob("*_lunp")))

    def test_invalid_callback_never_published_as_real_probability(self):
        for module in [fake_rna(invalid=True), fake_rna(split=True), fake_rna(missing=True)]:
            with tempfile.TemporaryDirectory() as temporary, patch(
                    "zyme_designer.vienna_backend.importlib.import_module", return_value=module):
                with self.assertRaises(RuntimeError):
                    run_vienna_python("ACGUACGU", temporary, W=8, L=6, u=4)
                self.assertFalse(list(Path(temporary).glob("*_lunp")))

    def test_invalid_inputs_rejected_before_native_call(self):
        cases = [dict(sequence="ACGN"), dict(sequence="AC"), dict(W=3), dict(L=30),
                 dict(u=30), dict(T=math.nan), dict(T=-5), dict(T=110),
                 dict(salt_molar=0), dict(salt_molar=math.inf), dict(W=True)]
        with tempfile.TemporaryDirectory() as temporary, patch(
                "zyme_designer.vienna_backend.importlib.import_module") as importer:
            for overrides in cases:
                args = dict(sequence="ACGUACGUACGU", W=10, L=8, u=5, T=37, salt_molar=.3)
                args.update(overrides)
                with self.assertRaises(ValueError):
                    run_vienna_python(outdir=temporary, **args)
            importer.assert_not_called()


@unittest.skipUnless(importlib.util.find_spec("RNA") is not None, "Real ViennaRNA Python package not installed")
class ViennaBackendRealTests(unittest.TestCase):
    def setUp(self):
        import RNA
        self.RNA = RNA
        if RNA.__version__ != "2.7.2":
            self.skipTest("Real regression requires the validated ViennaRNA2.7.2 runtime")

    def test_real_local_native_wrapper_agrees_and_boundary_indexing(self):
        seq = "GGGAAACCCAUUGCAAAGGAAACCACGUGUUUCGCAAGGGAAACCC"
        with tempfile.TemporaryDirectory() as temporary:
            lunp, audit = run_vienna_python(seq, temporary, W=30, L=20, u=5)
            table = parse_lunp(lunp)
        matrix = self.RNA.pfl_fold_up(seq, 5, 30, 20)
        for end in range(1, len(seq) + 1):
            for length in range(1, min(5, end) + 1):
                self.assertAlmostEqual(table[(end, length)], matrix[end][length], places=12)
        self.assertIsNone(table[(1, 2)])
        self.assertNotAlmostEqual(table[(5, 5)], table[(1, 1)])
        self.assertEqual(audit["energy_parameters_sha256"], "2a43345a495850cfd2e0a78c57c6e02085e6df3c53496fe3289dc294d21732ad")

    def test_real_birc5_archived_rnaplfold_numeric_regression(self):
        project = Path(__file__).resolve().parents[1]
        fasta = project / "examples/BIRC5_NM_001168_3.fasta"
        reference = project / "examples/BIRC5_accessibility_W200_selected.csv"
        if not fasta.is_file() or not reference.is_file():
            self.skipTest("BIRC5 example data not present")
        sequence = "".join(fasta.read_text(encoding="utf8").splitlines()[1:])
        with tempfile.TemporaryDirectory() as temporary:
            lunp, audit = run_vienna_python(sequence, temporary, W=200, L=150, u=16, T=37, salt_molar=1.021)
            table = parse_lunp(lunp)
        self.assertEqual(audit["sequence_sha256_rna"], "b9d001c1d678a5597f226cc7db46a2063c20accbd0edd9470f802441a2930988")
        checks = 0
        with reference.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["accession"] != "NM_001168.3":
                    continue
                metrics = region_metrics(table, int(row["start_1based"]), 7, 7, 37)
                for output, archived in (("full_window_p_unpaired", "p_unpaired_full16"),
                                         ("left_arm_p_unpaired", "p_unpaired_left7"),
                                         ("right_arm_p_unpaired", "p_unpaired_right7"),
                                         ("ag_p_unpaired", "p_unpaired_center2")):
                    self.assertLess(abs(metrics[output] - float(row[archived])), 5.1e-8)
                    checks += 1
        self.assertEqual(checks, 12)

    def test_real_temperature_and_salt_are_actually_used(self):
        sequence = "GGGAAACCCAUUGCAAAGGAAACCACGUGUUUCGCAAGGGAAACCC"
        values = []
        for temperature, salt in ((37, 1.021), (37, .3), (25, 1.021)):
            with tempfile.TemporaryDirectory() as temporary:
                lunp, audit = run_vienna_python(sequence, temporary, W=30, L=20, u=5, T=temperature, salt_molar=salt)
                values.append(parse_lunp(lunp)[(10, 5)])
                self.assertAlmostEqual(audit["model_details"]["salt"], salt)
                self.assertAlmostEqual(audit["model_details"]["temperature"], temperature)
        self.assertNotAlmostEqual(values[0], values[1], places=6)
        self.assertNotAlmostEqual(values[0], values[2], places=6)

    def test_real_short_poly_a_and_dna_alphabet_normalization(self):
        for sequence in ("AAAA", "A" * 32):
            with tempfile.TemporaryDirectory() as temporary:
                lunp, audit = run_vienna_python(sequence, temporary)
                table = parse_lunp(lunp)
            self.assertEqual(audit["parameters"]["W"], len(sequence))
            self.assertEqual(audit["n_probability_rows"], len(sequence))
            for value in table.values():
                if value is not None:
                    self.assertAlmostEqual(value, 1, places=12)
            if len(sequence) == 4:
                self.assertEqual(audit["numerical_audit"]["nonphysical_boundary_values_discarded"], 6)
            else:
                self.assertGreater(audit["numerical_audit"]["roundoff_clipped_probability_count"], 0)
                self.assertLess(audit["numerical_audit"]["maximum_roundoff_clip_magnitude"], 1e-12)
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            _, dna = run_vienna_python("GGGAAACCCTTTGGGAAACCC", a)
            _, rna = run_vienna_python("GGGAAACCCUUUGGGAAACCC", b)
        self.assertEqual(dna["sequence_sha256_rna"], rna["sequence_sha256_rna"])
        self.assertEqual(dna["lunp_sha256"], rna["lunp_sha256"])

    @unittest.skipUnless(shutil.which("RNAplfold"), "RNAplfold executable not available for independent CLI comparison")
    def test_optional_real_cli_matches_api_same_version(self):
        executable = shutil.which("RNAplfold")
        version = subprocess.run([executable, "--version"], text=True, capture_output=True, check=True)
        if "2.7.2" not in version.stdout + version.stderr:
            self.skipTest("CLI comparison requires version2.7.2")
        sequence = "GGGAAACCCAUUGCAAAGGAAACCACGUGUUUCGCAAGGGAAACCC"
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            lunp, _ = run_vienna_python(sequence, work / "api", W=30, L=20, u=5, T=37, salt_molar=.3)
            api = parse_lunp(lunp)
            subprocess.run([executable, "-W", "30", "-L", "20", "-u", "5", "-T", "37", "--salt", ".3"],
                           input=f">independent\n{sequence}\n", text=True, cwd=work, capture_output=True, check=True)
            cli = parse_lunp(work / "independent_lunp")
        self.assertEqual(api.keys(), cli.keys())
        for key in api:
            if api[key] is None:
                self.assertIsNone(cli[key])
            else:
                self.assertLess(abs(api[key] - cli[key]), 5.1e-8)


if __name__ == "__main__":
    unittest.main()
