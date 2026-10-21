"""Independent exhaustive-reference tests for the streaming seed screen."""

from pathlib import Path
import gzip
import random
import tempfile
import types
import unittest
from unittest.mock import patch

from zyme_designer.offtarget import scan_transcriptome


def candidate(sequence="AACAGGUU", identifier="site", left=3, right=3, start=1):
    return {"candidate_id": identifier, "target_rna": sequence,
            "left_length": left, "right_length": right,
            "start_1based": start, "end_1based": start + len(sequence) - 1}


def exhaustive(sequence, candidates, max_mismatches):
    found = set()
    for item in candidates:
        target = item["target_rna"].upper().replace("T", "U")
        left = item["left_length"]
        for start in range(len(sequence) - len(target) + 1):
            window = sequence[start:start + len(target)]
            if set(window) - set("ACGU"):
                continue
            mismatch = sum(a != b for a, b in zip(target[:left], window[:left]))
            mismatch += sum(a != b for a, b in zip(target[left + 2:], window[left + 2:]))
            if mismatch <= max_mismatches:
                found.add((item["candidate_id"], start + 1, mismatch, window[left:left + 2]))
    return found


class OffTargetTests(unittest.TestCase):
    def scan_text(self, text, candidates=None, **kwargs):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reference.fasta"
            path.write_text(text, encoding="utf-8")
            return scan_transcriptome(candidates or [candidate()], path, "TARGET.1", **kwargs)

    def test_randomized_against_exhaustive_reference_all_limits(self):
        rng = random.Random(81017)
        for trial in range(12):
            sequence = "".join(rng.choice("ACGU") for _ in range(130))
            candidates = []
            for index in range(4):
                left, right = rng.randint(1, 5), rng.randint(1, 5)
                target = "".join(rng.choice("ACGU") for _ in range(left)) + "AG"
                target += "".join(rng.choice("ACGU") for _ in range(right))
                candidates.append(candidate(target, str(index), left, right))
            for limit in (0, 1, 2):
                with self.subTest(trial=trial, mismatch_limit=limit):
                    _, rows, audit = self.scan_text(
                        ">TX1 gene=OTHER\n" + sequence, candidates,
                        target_gene="TARGETGENE", max_mismatches=limit, hit_limit=100000,
                    )
                    actual = {(r["candidate_id"], r["start_1based"], r["arm_mismatches"], r["cut_motif"]) for r in rows}
                    expected = exhaustive(sequence, candidates, limit)
                    self.assertEqual(actual, expected)
                    self.assertEqual(audit["binding_hits_counted"], len(expected))

    def test_optional_automaton_route_and_fallback_parity(self):
        """Exercise the optional API route without requiring a compiled extension."""
        class ReferenceAutomaton:
            def __init__(self):
                self.words = {}

            def add_word(self, word, value):
                self.words[word] = value

            def make_automaton(self):
                pass

            def iter(self, sequence):
                events = []
                for word, value in self.words.items():
                    start = sequence.find(word)
                    while start >= 0:
                        events.append((start + len(word) - 1, value))
                        start = sequence.find(word, start + 1)
                return iter(sorted(events, key=lambda item: item[0]))

        text = ">TARGET.1 (T)\nAACAGGUUAACUCGUUACCAGGAU\n"
        candidates = [candidate(), candidate("ACAGGU", "asymmetric", 1, 3)]
        with patch.dict("sys.modules", {"ahocorasick": None}):
            fallback_summaries, fallback_rows, _ = self.scan_text(text, candidates, target_gene="T")
        fake_module = types.SimpleNamespace(Automaton=ReferenceAutomaton)
        with patch.dict("sys.modules", {"ahocorasick": fake_module}):
            aho_summaries, aho_rows, audit = self.scan_text(text, candidates, target_gene="T")
        self.assertEqual(aho_summaries, fallback_summaries)
        self.assertEqual(aho_rows, fallback_rows)
        self.assertEqual(audit["backend"], "pyahocorasick")

    def test_central_motif_never_consumes_arm_budget(self):
        summaries, rows, _ = self.scan_text(
            ">TX1 (OTHER)\nAACUCGUU\n", target_gene="T", max_mismatches=0,
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["arm_mismatches"], 0)
        self.assertEqual(rows[0]["cut_motif"], "UC")
        self.assertEqual(summaries["site"]["other_gene_non_ag"], 1)

    def test_other_gene_all_mismatch_categories_and_gene_count(self):
        text = ">TX1 gene_symbol:OTHER\nAACAGGUU\n>TX2 gene_symbol:OTHER\nACCAGGUU\n"
        text += ">TX3 gene_symbol:SECOND\nACCAGGAU\n"
        summary, _, _ = self.scan_text(text, target_gene="T")
        site = summary["site"]
        self.assertEqual((site["other_gene_exact_ag"], site["other_gene_1mm_ag"], site["other_gene_2mm_ag"]), (1, 1, 1))
        self.assertEqual(site["unique_other_gene_count"], 2)

    def test_late_gene_resolution_same_gene_and_unknown(self):
        text = ">ENST000002.1 gene:ENSG000001.3 gene_symbol:BIRC5\nAACAGGUU\n"
        text += ">TX_UNKNOWN\nAACAGGUU\n>NM_001168.3_verified_from_PS (BIRC5)\nAACAGGUU\n"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rna.fna"
            path.write_text(text, encoding="utf-8")
            summaries, rows, audit = scan_transcriptome([candidate()], path, "NM_001168.3", max_mismatches=0)
        self.assertEqual(summaries["site"]["target_isoform_hits"], 1)
        self.assertEqual(summaries["site"]["intended_target_hits"], 1)
        self.assertEqual(summaries["site"]["unknown_gene_hits"], 1)
        self.assertEqual(rows[0]["gene_classification"], "target_isoform")
        self.assertEqual(rows[2]["transcript_accession"], "NM_001168.3")
        self.assertEqual(audit["exact_target_records_found"], 1)

    def test_gencode_pipe_and_gene_id_alias_matching(self):
        text = ">ENST000001.1|ENSG000001.2|OTTHUMG1|OTTHUMT1|BIRC5-201|BIRC5|8|protein_coding|\nAACAGGUU\n"
        summaries, rows, _ = self.scan_text(text, target_gene="ENSG000001.8", max_mismatches=0)
        self.assertEqual(summaries["site"]["target_isoform_hits"], 1)
        self.assertEqual(rows[0]["gene"], "BIRC5")

    def test_other_gene_identity_links_only_explicit_header_aliases(self):
        text = ">TX1 gene:ENSG000002.1\nAACAGGUU\n"
        text += ">TX2 gene:ENSG000002.8 gene_symbol:OTHER\nAACAGGUU\n"
        text += ">TX3 gene_symbol:OTHER\nAACAGGUU\n"
        summaries, rows, _ = self.scan_text(text, target_gene="T", max_mismatches=0)
        self.assertEqual(summaries["site"]["other_gene_exact_ag"], 3)
        self.assertEqual(summaries["site"]["unique_other_gene_count"], 1)
        self.assertEqual({r["resolved_gene_key"] for r in rows}, {"OTHER"})

    def test_overlapping_hits_export_cap_does_not_cap_counts(self):
        item = candidate("AAAGAA", left=2, right=2)
        sequence = "A" * 30
        summaries, rows, audit = self.scan_text(
            ">TX1 (OTHER)\n" + sequence, [item],
            target_gene="T", max_mismatches=0, hit_limit=2,
        )
        self.assertEqual(summaries["site"]["binding_hits_total"], 25)
        self.assertEqual(summaries["site"]["other_gene_non_ag"], 25)
        self.assertEqual(len(rows), 2)
        self.assertTrue(audit["hit_rows_truncated"])
        self.assertTrue(audit["summary_counts_complete"])
        self.assertEqual(audit["binding_hits_counted"], 25)

    def test_ambiguous_bases_not_filled_and_window_counts(self):
        item = candidate("AAAGAA", left=2, right=2)
        sequence = "AAAAAANAAAAAARAAAAAA"
        summaries, rows, audit = self.scan_text(
            ">TX1 (OTHER)\n" + sequence, [item], target_gene="T", max_mismatches=0,
        )
        expected = exhaustive(sequence, [item], 0)
        self.assertEqual(len(rows), len(expected))
        self.assertEqual(audit["ambiguous_nucleotides"], 2)
        self.assertEqual(summaries["site"]["unambiguous_windows"], len(expected))
        self.assertEqual(summaries["site"]["total_windows"], len(sequence) - 5)
        self.assertEqual(summaries["site"]["ambiguous_windows_skipped"], len(sequence) - 5 - len(expected))

    def test_compressed_file_hash_and_progress(self):
        records = ">TARGET.1 gene=T\nAACAGGUU\n>TX2 gene=OTHER\nAACAGGUU\n"
        events = []
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rna.fna.gz"
            with gzip.open(path, "wt", encoding="utf-8") as stream:
                stream.write(records)
            summaries, rows, audit = scan_transcriptome([candidate()], path, "TARGET.1", max_mismatches=0, progress=events.append)
        self.assertEqual(len(audit["input_file_sha256"]), 64)
        self.assertEqual(audit["records_scanned"], 2)
        self.assertEqual(audit["nucleotides_scanned"], 16)
        self.assertTrue(events[-1]["complete"])
        self.assertEqual(summaries["site"]["intended_site_hits"], 1)

    def test_additional_sites_in_target_and_duplicate_records(self):
        item = candidate("AAAGAA", left=2, right=2)
        summary, rows, audit = self.scan_text(
            ">TARGET.1 (T)\nAAAAAAAA\n>TARGET.1 (T)\nAAAAAAAA\n", [item], max_mismatches=0,
        )
        self.assertEqual(summary["site"]["intended_target_hits"], 6)
        self.assertEqual(summary["site"]["intended_site_hits"], 2)
        self.assertEqual(summary["site"]["other_sites_on_target_hits"], 4)
        self.assertEqual(audit["duplicate_accession_records"], 1)

    def test_no_target_gene_does_not_guess_known_other_gene(self):
        summary, rows, audit = self.scan_text(
            ">TX1 gene=OTHER\nAACAGGUU\n", max_mismatches=0,
        )
        self.assertEqual(summary["site"]["unknown_gene_hits"], 1)
        self.assertEqual(summary["site"]["other_gene_exact_ag"], 0)
        self.assertTrue(audit["warnings"])

    def test_spreadsheet_formula_escape(self):
        summary, rows, _ = self.scan_text(
            ">@unsafe gene=+unsafe\nAACAGGUU\n", target_gene="T", max_mismatches=0,
        )
        self.assertEqual(rows[0]["transcript_accession"], "'@unsafe")
        self.assertEqual(rows[0]["gene"], "'+UNSAFE")

    def test_malformed_fasta_and_candidates_rejected(self):
        for text in ("AACAGGUU", ">\nAACAGGUU", ">TX\n", ">TX\nAAC-AGGUU", ">TX\nAAC1AGGUU", ""):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.scan_text(text)
        bad = candidate()
        bad["target_rna"] = "NACAGGUU"
        with self.assertRaises(ValueError):
            self.scan_text(">TX\nAACAGGUU", [bad])
        with self.assertRaises(ValueError):
            self.scan_text(">TX\nAACAGGUU", max_mismatches=3)
        with self.assertRaises(ValueError):
            self.scan_text(">TX\nAACAGGUU", hit_limit=-1)


if __name__ == "__main__":
    unittest.main()
