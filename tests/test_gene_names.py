#!/usr/bin/env python3
"""Locus resolution must stay on the gene that owns the identifier."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import gene_names  # noqa: E402


GFF = """\
##gff-version 3
chr4	RefSeq	gene	1	100	.	-	.	ID=gene-AT4G35550;Dbxref=TAIR:AT4G35550,GeneID:829707;Name=WOX13;gene=WOX13;gene_synonym=F8D20.60;locus_tag=AT4G35550
chr4	RefSeq	mRNA	1	100	.	-	.	ID=rna-NM_119720.4;Parent=gene-AT4G35550;gene=WOX13;locus_tag=AT4G35550
chr4	RefSeq	CDS	1	100	.	-	0	ID=cds-NP_195280.1;Parent=rna-NM_119720.4;Name=NP_195280.1;gene=WOX13;locus_tag=AT4G35550;protein_id=NP_195280.1;product=WUSCHEL related homeobox 13
chr4	RefSeq	gene	200	900	.	+	.	ID=gene-AT4G35500;Name=AT4G35500;gene_synonym=F8D20.10;locus_tag=AT4G35500
chr4	RefSeq	mRNA	200	900	.	+	.	ID=rna-NM_202960.1;Parent=gene-AT4G35500;locus_tag=AT4G35500
chr4	RefSeq	CDS	200	900	.	+	0	ID=cds-NP_974689.1;Parent=rna-NM_202960.1;Name=NP_974689.1;locus_tag=AT4G35500;protein_id=NP_974689.1;Note=Protein kinase superfamily protein
chr1	RefSeq	gene	1	100	.	-	.	ID=gene-AT1G71380;Name=CEL3;gene=CEL3;gene_synonym=ATGH9B3,ATCEL3,F3I17.16;locus_tag=AT1G71380
chr1	RefSeq	CDS	1	100	.	-	0	ID=cds-NP_177294.1;Parent=gene-AT1G71380;Name=NP_177294.1;gene=CEL3;locus_tag=AT1G71380;protein_id=NP_177294.1
chr1	RefSeq	gene	200	5000	.	+	.	ID=gene-AT1G71220;Name=EBS1;gene=EBS1;gene_synonym=F3I17.3;locus_tag=AT1G71220
chr1	RefSeq	CDS	200	5000	.	+	0	ID=cds-NP_001323059.1;Parent=gene-AT1G71220;Name=NP_001323059.1;gene=EBS1;locus_tag=AT1G71220;protein_id=NP_001323059.1
chr1	RefSeq	gene	6000	7000	.	+	.	ID=gene-AT1G22880;Name=CEL5;gene=CEL5;locus_tag=AT1G22880
chr1	RefSeq	CDS	6000	7000	.	+	0	ID=cds-NP_001031082.1;Parent=gene-AT1G22880;Name=NP_001031082.1;gene=CEL5;locus_tag=AT1G22880;protein_id=NP_001031082.1;Note=BEST Arabidopsis thaliana protein match is: cellulase 3 (TAIR:AT1G71380.1)
chr2	ITAG	gene	1	10	.	+	.	ID=gene:Solyc04g078650;Name=Solyc04g078650
chr2	ITAG	mRNA	1	10	.	+	.	ID=mRNA:Solyc04g078650.4.1;Parent=gene:Solyc04g078650;Name=Solyc04g078650.4.1;Note=wuschel related homobox 4
"""


def protein(original_id: str, length: int, species_id: str = "Ath") -> dict:
    genus = "Arabidopsis_thaliana" if species_id == "Ath" else "Solanum_lycopersicum"
    return {
        "species_id": species_id,
        "original_id": original_id,
        "gene_key": original_id.split(".")[0],
        "new_id": f"{species_id}|{genus}|{original_id}",
        "length": str(length),
    }


class IdFormsTest(unittest.TestCase):
    def test_bac_gene_number_is_not_an_isoform_suffix(self):
        self.assertEqual(gene_names.id_forms("F8D20.60"), {"f8d20.60"})
        self.assertNotIn("f8d20", gene_names.id_forms("F8D20.60"))
        self.assertNotIn("f3i17", gene_names.id_forms("F3I17.16"))

    def test_recognized_locus_suffixes_still_collapse(self):
        self.assertIn("at4g35550", gene_names.id_forms("AT4G35550.1"))
        self.assertIn("np195280", gene_names.id_forms("NP_195280.1"))
        self.assertIn("solyc04g078650", gene_names.id_forms("Solyc04g078650.4.1"))
        self.assertIn("np195280", gene_names.id_forms("cds-NP_195280.1"))


class ResolveLocusTest(unittest.TestCase):
    def match(self, gene: str, species_id: str, views):
        pieces, like = gene_names.split_pieces(gene)
        query = gene_names.NameQuery(gene, species_id, gene, 1, pieces, like)
        return gene_names.match_in_species(query, species_id, views)

    def test_bac_neighbors_and_note_citations_do_not_steal_the_locus(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.gff"
            path.write_text(GFF)
            index = gene_names.index_gffs([str(path)])
        rows = [
            protein("NP_195280.1", 268),
            protein("NP_974689.1", 440),
            protein("NP_177294.1", 485),
            protein("NP_001323059.1", 1652),
            protein("NP_001031082.1", 600),
            protein("Solyc04g078650.4.1", 320, "Sly"),
        ]
        views = gene_names.attach_proteins(rows, index)

        method, view, _detail = self.match("AT4G35550", "Ath", views)
        self.assertEqual(method, "locus")
        self.assertEqual(view.row["original_id"], "NP_195280.1")

        method, view, _detail = self.match("WOX13", "Ath", views)
        self.assertEqual(method, "symbol")
        self.assertEqual(view.row["original_id"], "NP_195280.1")

        method, view, _detail = self.match("AT1G71380", "Ath", views)
        self.assertEqual(method, "locus")
        self.assertEqual(view.row["original_id"], "NP_177294.1")

        method, view, _detail = self.match("CEL3", "Ath", views)
        self.assertEqual(method, "symbol")
        self.assertEqual(view.row["original_id"], "NP_177294.1")

        method, view, _detail = self.match("GH9B3", "Ath", views)
        self.assertEqual(method, "symbol")
        self.assertEqual(view.row["original_id"], "NP_177294.1")

        method, view, _detail = self.match("WOX4", "Sly", views)
        self.assertEqual(method, "description")
        self.assertEqual(view.row["original_id"], "Solyc04g078650.4.1")


if __name__ == "__main__":
    unittest.main()
