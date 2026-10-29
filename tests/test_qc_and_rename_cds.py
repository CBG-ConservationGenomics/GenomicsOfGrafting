#!/usr/bin/env python3
"""NCBI cds_from_genomic accessions have to follow the protein id."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import qc_and_rename_cds as qc  # noqa: E402


ORIG_TO_NEW = {
    "NP_195280.1": "Ath|Arabidopsis_thaliana|NP_195280.1",
    "NP_974689.1": "Ath|Arabidopsis_thaliana|NP_974689.1",
    "Solyc04g078650.4.1": "Sly|Solanum_lycopersicum|Solyc04g078650.4.1",
}
GENE_TO_NEW = {
    "NP_195280": "Ath|Arabidopsis_thaliana|NP_195280.1",
    "NP_974689": "Ath|Arabidopsis_thaliana|NP_974689.1",
    "Solyc04g078650": "Sly|Solanum_lycopersicum|Solyc04g078650.4.1",
}


class LookupTest(unittest.TestCase):
    def test_ncbi_header_uses_its_own_protein_id(self):
        header = (
            "lcl|NC_003075.7_cds_NP_195280.1_1 [gene=WOX13] [locus_tag=AT4G35550] "
            "[protein=WUSCHEL related homeobox 13] [protein_id=NP_195280.1] "
            "homolog NP_974689.1"
        )
        found = qc.lookup_new_id(header.split()[0], ORIG_TO_NEW, GENE_TO_NEW, header=header)
        self.assertEqual(found, "Ath|Arabidopsis_thaliana|NP_195280.1")
        self.assertEqual(
            qc.header_accessions(header),
            ["NP_195280.1"],
        )

    def test_accession_embedded_in_the_id_is_enough(self):
        seq_id = "lcl|NC_003075.7_cds_NP_195280.1_1"
        found = qc.lookup_new_id(seq_id, ORIG_TO_NEW, GENE_TO_NEW)
        self.assertEqual(found, "Ath|Arabidopsis_thaliana|NP_195280.1")

    def test_itag_cds_id_still_matches_the_protein(self):
        found = qc.lookup_new_id("Solyc04g078650.4.1", ORIG_TO_NEW, GENE_TO_NEW)
        self.assertEqual(found, "Sly|Solanum_lycopersicum|Solyc04g078650.4.1")


class RenameScriptTest(unittest.TestCase):
    def test_ncbi_and_itag_cds_are_both_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proteins = root / "proteins.faa"
            cds = root / "cds.fna"
            id_map = root / "id_map.tsv"
            out = root / "cleaned.fna"
            report = root / "report.tsv"
            proteins.write_text(
                ">Ath|Arabidopsis_thaliana|NP_195280.1\nMKTF\n"
                ">Sly|Solanum_lycopersicum|Solyc04g078650.4.1\nMKTF\n"
            )
            cds.write_text(
                ">lcl|NC_003075.7_cds_NP_195280.1_1 [gene=WOX13] "
                "[locus_tag=AT4G35550] [protein_id=NP_195280.1]\n"
                "ATGAAAACCTTC\n"
                ">Solyc04g078650.4.1\n"
                "ATGAAAACCTTC\n"
                ">lcl|NC_003075.7_cds_NP_000000.1_1 [protein_id=NP_000000.1]\n"
                "ATGAAAACCTTC\n"
            )
            id_map.write_text(
                "species_id\tgenus_species\toriginal_id\tgene_key\tnew_id\tlength\n"
                "Ath\tArabidopsis_thaliana\tNP_195280.1\tNP_195280\t"
                "Ath|Arabidopsis_thaliana|NP_195280.1\t4\n"
                "Sly\tSolanum_lycopersicum\tSolyc04g078650.4.1\tSolyc04g078650\t"
                "Sly|Solanum_lycopersicum|Solyc04g078650.4.1\t4\n"
            )
            old = sys.argv
            sys.argv = [
                "qc_and_rename_cds.py",
                "--cds",
                str(cds),
                "--proteins",
                str(proteins),
                "--id-map",
                str(id_map),
                "--out-cds",
                str(out),
                "--out-report",
                str(report),
            ]
            try:
                qc.main()
            finally:
                sys.argv = old
            with out.open() as handle:
                kept = {rec.id for rec in SeqIO.parse(handle, "fasta")}
        self.assertEqual(
            kept,
            {
                "Ath|Arabidopsis_thaliana|NP_195280.1",
                "Sly|Solanum_lycopersicum|Solyc04g078650.4.1",
            },
        )


if __name__ == "__main__":
    unittest.main()
