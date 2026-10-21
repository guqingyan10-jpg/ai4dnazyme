"""Desktop configuration/worker regressions, runnable without opening a GUI."""
from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from zyme_designer import __version__
from zyme_designer.desktop import (
    DesktopConfig, FOLD_NONE, FOLD_PYTHON, FOLD_SITE, SCAFFOLD_ARMS,
    SCAFFOLD_BUILTIN, STRATEGY_GENERIC, build_cli_args, demo_config,
    parse_cut_positions, parse_lengths, pasted_fasta_text, resource_root,
    unique_output_path, validate_config, worker_command,
)


class DesktopTests(unittest.TestCase):
    def setUp(self):
        project = Path(__file__).resolve().parents[1]
        self.temp = tempfile.TemporaryDirectory(prefix="zyme_desktop_test_")
        self.root = Path(self.temp.name)
        self.target = self.root / "target with spaces.fasta"
        self.target.write_text(">NM_001168.3 gene=BIRC5\nAUUGCAAAGGAAACCA\n", encoding="utf-8")
        self.base = DesktopConfig(target_fasta=str(self.target), output_parent=str(self.root))

    def tearDown(self):
        self.temp.cleanup()

    def test_default_three_lengths(self):
        self.assertEqual(parse_lengths(self.base.left_lengths), (7, 8, 9))
        self.assertEqual(self.base.scaffold_mode, SCAFFOLD_ARMS)
        self.assertEqual(self.base.fold_mode, FOLD_PYTHON)

    def test_length_duplicates_sorted(self):
        self.assertEqual(parse_lengths("9,7,7,8"), (7, 8, 9))

    def test_invalid_lengths(self):
        for value in ("", "0", "21", "7;rm", "7.5"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_lengths(value)

    def test_raw_paste_normalizes_rna(self):
        self.assertEqual(pasted_fasta_text("ATTGCA AAGGAAACCA", "BIRC5"), ">USER_INPUT gene=BIRC5\nAUUGCAAAGGAAACCA\n")

    def test_fasta_paste_keeps_header(self):
        self.assertTrue(pasted_fasta_text(">NM_001168.3 gene=BIRC5\nATTGCA").startswith(">NM_001168.3 gene=BIRC5\nAUUGCA"))

    def test_paste_multiple_records_rejected(self):
        with self.assertRaises(ValueError):
            pasted_fasta_text(">first\nAG\n>second\nAG")

    def test_unknown_gene_not_invented(self):
        self.assertNotIn("gene=", pasted_fasta_text("AAG"))

    def test_paste_unknown_bases_rejected(self):
        with self.assertRaises(ValueError):
            pasted_fasta_text("AAN")

    def test_missing_input_rejected(self):
        with self.assertRaises(ValueError):
            validate_config(replace(self.base, target_fasta="missing.fasta"))

    def test_gene_header_injection_rejected(self):
        with self.assertRaises(ValueError):
            validate_config(replace(self.base, gene="BIRC5\n>injection"))

    def test_generic_iupac_and_no_scaffold(self):
        config = replace(self.base, strategy=STRATEGY_GENERIC, motif="RY", fold_mode=FOLD_NONE)
        args = build_cli_args(config, self.target, self.root / "new")
        self.assertEqual(args[args.index("--motif") + 1], "RY")
        self.assertNotIn("--scaffold", args)
        self.assertNotIn("--fold", args)

    def test_generic_builtin_guard(self):
        with self.assertRaises(ValueError):
            validate_config(replace(self.base, strategy=STRATEGY_GENERIC, motif="NN", scaffold_mode=SCAFFOLD_BUILTIN))

    def test_817_motif_guard(self):
        with self.assertRaises(ValueError):
            validate_config(replace(self.base, motif="RY"))

    def test_invalid_motif(self):
        for motif in ("A", "AGG", "A;", ""):
            with self.subTest(motif=motif), self.assertRaises(ValueError):
                validate_config(replace(self.base, motif=motif))

    def test_cut_positions(self):
        self.assertEqual(parse_cut_positions("420,408,408"), (408, 420))
        self.assertEqual(parse_cut_positions(""), ())
        with self.assertRaises(ValueError):
            parse_cut_positions("0")

    def test_numerical_guards(self):
        for updates in ({"temperature": float("nan")}, {"rna_salt": 0}, {"span": 201}, {"screen_top": -1}, {"mismatches": 3}, {"preferred_min": 10, "preferred_max": 7}, {"preferred_min": 4}, {"preferred_max": 14}):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                validate_config(replace(self.base, **updates))

    def test_import_requires_manifest(self):
        with self.assertRaises(ValueError):
            validate_config(replace(self.base, fold_mode=FOLD_SITE, accessibility_path=str(self.target)))

    def test_reference_scope_required(self):
        with self.assertRaises(ValueError):
            validate_config(replace(self.base, transcriptome=str(self.target), reference_label=""))

    def test_unique_output_does_not_create(self):
        first, second = unique_output_path(self.root), unique_output_path(self.root)
        self.assertNotEqual(first, second)
        self.assertFalse(first.exists())
        self.assertEqual(first.parent, self.root)

    def test_default_python_folding_arguments(self):
        args = build_cli_args(self.base, self.target, self.root / "new")
        self.assertIn("--fold", args)
        self.assertEqual(args[args.index("--fold-backend") + 1], "python")
        self.assertIn(str(self.target), args)
        self.assertEqual(args[args.index("--preferred-arm-min") + 1], "7")
        self.assertEqual(args[args.index("--preferred-arm-max") + 1], "9")

    def test_worker_uses_argv_not_shell_string(self):
        args = build_cli_args(self.base, self.target, self.root / "new")
        command = worker_command(args, self.root / "log UTF8.log")
        self.assertIsInstance(command, list)
        self.assertIn(str(self.target), command)
        self.assertIn("--worker-log", command)

    def test_frozen_resources_and_worker(self):
        with patch.object(sys, "frozen", True, create=True), patch.object(sys, "_MEIPASS", str(self.root), create=True):
            self.assertEqual(resource_root(), self.root)
            command = worker_command(["--version"], self.root / "worker.log")
            self.assertEqual(command[:2], [sys.executable, "--worker"])
            self.assertNotIn("desktop_entry.py", command)

    def test_demo_explicit_small_library(self):
        config = demo_config(self.root)
        self.assertIn("NOT full human", config.reference_label)
        self.assertEqual(config.left_lengths, "7")
        self.assertEqual(config.fold_mode, FOLD_SITE)
        args = build_cli_args(config, Path(config.target_fasta), self.root / "new")
        self.assertIn("--site-table", args)
        self.assertEqual(args.count("--include-start"), 3)

    def test_worker_version_no_gui_and_utf8_log(self):
        log = self.root / "worker.log"
        command = worker_command(["--version"], log)
        proc = subprocess.run(command, capture_output=True, timeout=20)
        self.assertEqual(proc.returncode, 0, proc.stderr.decode("utf-8", errors="replace"))
        self.assertIn(__version__, log.read_text(encoding="utf-8"))

    def test_worker_invalid_command_nonzero(self):
        log = self.root / "worker.log"
        proc = subprocess.run(worker_command(["not-a-command"], log), capture_output=True, timeout=20)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("invalid choice", log.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
