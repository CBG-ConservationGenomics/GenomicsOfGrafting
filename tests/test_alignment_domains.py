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

    def test_protein_page_can_hide_domains(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "OG0001.fa").write_text(">Ath_Arabidopsis_thaliana_NP_195280.1\nMK--TFAA\n")
            ips = root / "interproscan.tsv"
            ips.write_text(IPS)
            index = pages.parse_interpro(ips)
            html = pages.render_proteins(["OG0001"], root, index, {})
            plain = pages.render_proteins(["OG0001"], root, None, {})
        self.assertIn("InterPro domains", html)
        self.assertIn('data-track="interpro"', html)
        self.assertIn("hide-interpro", html)
        self.assertIn("domain-legend", html)
        self.assertIn('class="domrow kind-interpro"', html)
        self.assertNotIn("data-track=", plain)

    def test_signal_transmembrane_and_glycan_are_separate_tracks(self) -> None:
        signalp = "\n".join(
            [
                "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position",
                "Ath|Arabidopsis_thaliana|NP_195280.1\tSP\t0.01\t0.99\tCS pos: 3-4. Pr: 0.9",
                "Sly|Solanum_lycopersicum|NP_2\tOTHER\t0.99\t0.01\t",
            ]
        ) + "\n"
        topology = "\n".join(
            [
                "id\ttopology",
                "Ath|Arabidopsis_thaliana|NP_195280.1\tSSSMMMMIII",
            ]
        ) + "\n"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "signalp.tsv").write_text(signalp)
            (root / "deeptmhmm.tsv").write_text(topology)
            (root / "OG0001.fa").write_text(">Ath|Arabidopsis_thaliana|NP_195280.1\nMKNGTFAAAA\n")
            signal = pages.parse_signalp(root / "signalp.tsv")
            tm = pages.parse_deeptmhmm(root / "deeptmhmm.tsv")
            glycan = pages.glycans_from_alignments(root, ["OG0001"])
            index = pages.combine_indexes(signal, tm, glycan)
            html = pages.render_proteins(["OG0001"], root, index, {})
        protein = "Ath|Arabidopsis_thaliana|NP_195280.1"
        kinds = {item.kind for item in index[protein]}
        self.assertEqual(kinds, {"signal", "transmembrane", "glycan"})
        signal_bar = next(item for item in index[protein] if item.kind == "signal")
        self.assertEqual((signal_bar.start, signal_bar.end), (1, 3))
        helix = next(item for item in index[protein] if item.kind == "transmembrane")
        self.assertEqual((helix.start, helix.end), (4, 7))
        sequon = next(item for item in index[protein] if item.kind == "glycan")
        self.assertEqual((sequon.start, sequon.end), (3, 5))
        self.assertNotIn("Sly|Solanum_lycopersicum|NP_2", index)
        self.assertIn("Signal peptides", html)
        self.assertIn("Transmembrane helices", html)
        self.assertIn("N-glycosylation", html)
        self.assertIn("kind-signal", html)
        self.assertIn("kind-transmembrane", html)
        self.assertIn("kind-glycan", html)


if __name__ == "__main__":
    unittest.main()
