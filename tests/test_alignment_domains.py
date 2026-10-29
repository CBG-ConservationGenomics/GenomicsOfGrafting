#!/usr/bin/env python3
"""InterPro domains keep their labels when they are drawn on an alignment."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import write_summary_pages as pages  # noqa: E402


IPS = "\n".join(
    [
        "\t".join(
            [
                "Ath|Arabidopsis_thaliana|NP_195280.1",
                "md5",
                "6",
                "Pfam",
                "PF00046",
                "Homeobox",
                "2",
                "5",
                "0",
                "T",
                "01-01-2020",
                "IPR001356",
                "Homeobox domain",
            ]
        ),
        "\t".join(
            [
                "Ath|Arabidopsis_thaliana|NP_195280.1",
                "md5",
                "6",
                "PANTHER",
                "PTHR1",
                "homeobox-like",
                "3",
                "6",
                "0",
                "T",
                "01-01-2020",
                "IPR001356",
                "Homeobox domain",
            ]
        ),
        "\t".join(
            [
                "Ath|Arabidopsis_thaliana|NP_195280.1",
                "md5",
                "6",
                "Pfam",
                "PF00000",
                "WUS box",
                "1",
                "1",
                "0",
                "T",
                "01-01-2020",
                "-",
                "-",
            ]
        ),
    ]
) + "\n"


class DomainMappingTest(unittest.TestCase):
    def test_overlapping_signatures_of_one_interpro_entry_merge(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interproscan.tsv"
            path.write_text(IPS)
            index = pages.parse_interpro(path)
        domains = index["Ath|Arabidopsis_thaliana|NP_195280.1"]
        homeobox = [item for item in domains if item.accession == "IPR001356"]
        self.assertEqual(len(homeobox), 1)
        self.assertEqual((homeobox[0].start, homeobox[0].end), (2, 6))
        self.assertEqual(homeobox[0].label, "Homeobox domain")
        wus = [item for item in domains if item.label == "WUS box"]
        self.assertEqual(len(wus), 1)

    def test_gaps_do_not_shift_the_domain(self):
        domain = pages.Domain("IPR001356", "Homeobox domain", 1, 4)
        spans = pages.column_spans("MK--TF", [domain], codon=False)
        self.assertEqual(spans, [(0, 5, domain)])

    def test_codon_alignment_uses_three_columns_per_residue(self):
        domain = pages.Domain("IPR001356", "Homeobox domain", 1, 2)
        seq = "ATGAAA------ACCTTC"
        spans = pages.column_spans(seq, [domain], codon=True)
        self.assertEqual(spans, [(0, 5, domain)])

    def test_alignment_html_labels_the_domain_and_matches_a_rewritten_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interproscan.tsv"
            path.write_text(IPS)
            index = pages.parse_interpro(path)
        lookup = lambda name: pages.domains_for(name, index, {})
        html = pages.protein_blocks(
            [("Ath_Arabidopsis_thaliana_NP_195280.1", "MK--TFAA")],
            domains_for=lookup,
        )
        self.assertIn("Homeobox domain", html)
        self.assertIn("WUS box", html)
        self.assertIn('class="dom on"', html)
        self.assertIn("IPR001356", html)

    def test_internal_orthofinder_id_uses_the_sequence_map(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interproscan.tsv"
            path.write_text(IPS)
            index = pages.parse_interpro(path)
        found = pages.domains_for("0_42", index, {"0_42": "Ath|Arabidopsis_thaliana|NP_195280.1"})
        self.assertTrue(any(item.accession == "IPR001356" for item in found))


if __name__ == "__main__":
    unittest.main()
