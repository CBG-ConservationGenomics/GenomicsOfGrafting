#!/usr/bin/env python3
"""Codon ids have to land on the species tree before HyPhy will run."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import prepare_hyphy_inputs as hyphy  # noqa: E402


TREE = "(Ath:0.1,Sly:0.2,(Nbe:0.3,Cpe:0.4):0.05);\n"
FASTA = """\
>Ath|Arabidopsis_thaliana|NP_1
ATGATG
>Sly|Solanum_lycopersicum|XP_1
ATGATG
>Nbe|Nicotiana_benthamiana|Nb1
ATGATG
>Nbe|Nicotiana_benthamiana|Nb2
ATGATG
>Zz|Nobody|XP_9
ATGATG
"""


class PrepareHyPhyInputsTest(unittest.TestCase):
    def prepare(self, tree: str = TREE, fasta: str = FASTA, foreground: str = "Ath\n"):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        codon = root / "og.codon.fna"
        species = root / "og.treefile"
        fg = root / "foreground_tips.txt"
        codon.write_text(fasta)
        species.write_text(tree)
        fg.write_text(foreground)
        out_fa = root / "hyphy.fna"
        out_tree = root / "fg.nwk"
        status = root / "prep.tsv"
        hyphy.prepare(codon, species, fg, out_fa, out_tree, status)
        return status.read_text(), out_fa, out_tree

    def test_expands_species_tips_and_labels_foreground(self):
        text, fasta, tree = self.prepare()
        status, _note, n_fg, n_seq = text.rstrip("\n").split("\t")
        self.assertEqual(status, "ok")
        self.assertEqual(n_fg, "1")
        self.assertEqual(n_seq, "4")
        self.assertNotIn("Zz|", fasta.read_text())
        parsed = hyphy.parse_newick(tree.read_text())
        tips = {tip.name: tip.foreground for tip in parsed.tips()}
        self.assertEqual(
            tips,
            {
                "Ath|Arabidopsis_thaliana|NP_1": True,
                "Sly|Solanum_lycopersicum|XP_1": False,
                "Nbe|Nicotiana_benthamiana|Nb1": False,
                "Nbe|Nicotiana_benthamiana|Nb2": False,
            },
        )
        written = tree.read_text()
        self.assertIn("Ath|Arabidopsis_thaliana|NP_1{Foreground}:0.1", written)
        self.assertNotIn("Cpe", written)
        self.assertNotIn("(Ath:", written)

    def test_underscore_ids_match_the_species_prefix(self):
        fasta = ">Ath_Arabidopsis_NP\nATG\n>Sly_Solanum_XP\nATG\n>Nbe_Nb1\nATG\n"
        text, _fasta, tree = self.prepare(fasta=fasta)
        self.assertTrue(text.startswith("ok\t"))
        names = {tip.name for tip in hyphy.parse_newick(tree.read_text()).tips()}
        self.assertEqual(names, {"Ath_Arabidopsis_NP", "Sly_Solanum_XP", "Nbe_Nb1"})

    def test_empty_tree_is_skipped(self):
        text, fasta, tree = self.prepare(tree="")
        self.assertTrue(text.startswith("skip\t"))
        self.assertIn("empty tree", text)
        self.assertFalse(fasta.exists())
        self.assertFalse(tree.exists())

    def test_two_sequences_are_skipped(self):
        fasta = ">Ath|a\nATG\n>Sly|b\nATG\n"
        text, _fasta, _tree = self.prepare(fasta=fasta)
        self.assertTrue(text.startswith("skip\t"))
        self.assertIn("n=2", text)

    def test_bootstrap_label_is_not_a_tip(self):
        fasta = ">Ath|a\nATG\n>Sly|b\nATG\n>Nbe|c\nATG\n"
        text, _fasta, tree = self.prepare(
            tree="(Ath:0.1,Sly:0.2,Nbe:0.3)95:0.4;\n",
            fasta=fasta,
        )
        self.assertTrue(text.startswith("ok\t"))
        names = {tip.name for tip in hyphy.parse_newick(tree.read_text()).tips()}
        self.assertNotIn("95", names)
        self.assertEqual(names, {"Ath|a", "Sly|b", "Nbe|c"})

    def test_pvalue_prefers_the_test_result(self):
        text = '{"other": {"p-value": 0.0001}, "test results": {"p-value": 0.042}}'
        self.assertEqual(hyphy.pvalue_from_json(text), "0.042")

    def test_pvalue_uses_the_smallest_site_value_when_there_is_no_test_block(self):
        text = '{"sites": [{"p-value": 0.2}, {"p-value": 0.01}]}'
        self.assertEqual(hyphy.pvalue_from_json(text), "0.01")

    def test_missing_pvalue_is_blank(self):
        self.assertEqual(hyphy.pvalue_from_json('{"error": "no test"}'), "")


if __name__ == "__main__":
    unittest.main()
