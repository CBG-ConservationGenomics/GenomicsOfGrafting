#!/usr/bin/env python3
"""Checks for the characterized-CDS confirmation step."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import check_candidate_cds as check  # noqa: E402


FIELDS = [
    "query_gene",
    "query_taxon",
    "species_id",
    "match_method",
    "original_id",
    "gene_key",
    "new_id",
    "annotation_text",
]


def write_crosswalk(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            full = {field: "" for field in FIELDS}
            full.update(row)
            writer.writerow(full)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def row(**kwargs) -> check.CrosswalkRow:
    data = {
        "query_gene": "WOX4",
        "query_taxon": "Solanum lycopersicum",
        "species_id": "Sly",
        "match_method": "description",
        "original_id": "Solyc04g078650.4.1",
        "gene_key": "Solyc04g078650",
        "new_id": "Sly|Solanum_lycopersicum|Solyc04g078650.4.1",
        "annotation_text": "",
    }
    data.update(kwargs)
    return check.CrosswalkRow(**data)


class CheckCandidateCdsTest(unittest.TestCase):
    def test_exact_match_ignores_terminal_stop_and_uracil(self):
        candidate = row()
        protein = "MKTF"
        cds = "ATGAAAACCTTC"
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: protein},
            {candidate.new_id: cds},
            [("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "augaaaaccuuctaa")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "match")
        self.assertEqual(rows[0]["protein_identity"], "100.00")
        self.assertEqual(rows[0]["cds_identity"], "100.00")
        self.assertEqual(failures, [])

    def test_synonymous_cds_is_reported_separately_from_the_protein(self):
        candidate = row()
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKTF"},
            {candidate.new_id: "ATGAAAACCTTC"},
            [("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "ATGAAAACCTTTTAA")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "cds_mismatch")
        self.assertEqual(rows[0]["protein_identity"], "100.00")
        self.assertLess(float(rows[0]["cds_identity"]), 98)
        self.assertEqual(len(failures), 1)

    def test_wrong_protein_fails(self):
        candidate = row()
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKKK"},
            {candidate.new_id: "ATGAAAAAAAAA"},
            [("WOX4", "", "WOX4", "ATGAAAACCTTC")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "mismatch")
        self.assertTrue(failures)

    def test_threshold_boundary(self):
        candidate = row()
        protein = "M" + ("A" * 49)
        same_cds = "ATG" + ("GCT" * 49)
        changed_cds = "ATG" + ("GCT" * 48) + "TGT"
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: protein},
            {candidate.new_id: same_cds},
            [("WOX4", "Sly", "WOX4 Sly", changed_cds)],
            98,
            98,
            True,
        )
        # 49/50 identities is exactly 98%.
        self.assertEqual(rows[0]["protein_identity"], "98.00")
        self.assertEqual(rows[0]["status"], "match")
        self.assertEqual(failures, [])

    def test_species_disambiguates_a_shared_symbol(self):
        sly = row()
        ath = row(
            query_taxon="Arabidopsis thaliana",
            species_id="Ath",
            original_id="AT1G46480.1",
            gene_key="AT1G46480",
            new_id="Ath|Arabidopsis_thaliana|AT1G46480.1",
        )
        rows, failures = check.run_check(
            [sly, ath],
            {sly.new_id: "MKTF", ath.new_id: "MKKK"},
            {sly.new_id: "ATGAAAACCTTC", ath.new_id: "ATGAAAAAAAAA"},
            [("WOX4", "Arabidopsis thaliana", "WOX4|Arabidopsis_thaliana", "ATGAAAAAAAAATAA")],
            98,
            98,
            True,
        )
        by_species = {item["species_id"]: item["status"] for item in rows}
        self.assertEqual(by_species["Ath"], "match")
        self.assertEqual(by_species["Sly"], "no_reference")
        self.assertTrue(any("no_reference" in line for line in failures))

    def test_missing_species_on_a_shared_symbol_is_ambiguous(self):
        sly = row()
        ath = row(
            query_taxon="Arabidopsis thaliana",
            species_id="Ath",
            original_id="AT1G46480.1",
            gene_key="AT1G46480",
            new_id="Ath|Arabidopsis_thaliana|AT1G46480.1",
        )
        rows, _failures = check.run_check(
            [sly, ath],
            {sly.new_id: "MKTF", ath.new_id: "MKTF"},
            {sly.new_id: "ATGAAAACCTTC", ath.new_id: "ATGAAAACCTTC"},
            [("WOX4", "", "WOX4", "ATGAAAACCTTC")],
            98,
            98,
            False,
        )
        self.assertEqual(rows[-1]["status"], "ambiguous_header")
        self.assertTrue(all(item["status"] == "no_reference" for item in rows[:-1]))

    def test_locus_header_matches_the_annotation_id(self):
        candidate = row()
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKTF"},
            {candidate.new_id.replace("|", "_"): "ATGAAAACCTTC"},
            [("Solyc04g078650", "", "Solyc04g078650", "ATGAAAACCTTC")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "match")
        self.assertEqual(failures, [])

    def test_embedded_taxon_in_the_header_token(self):
        candidate = row()
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKTF"},
            {candidate.new_id: "ATGAAAACCTTC"},
            [("WOX4_Solanum_lycopersicum", "", "WOX4_Solanum_lycopersicum", "ATGAAAACCTTC")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "match")
        self.assertEqual(failures, [])

    def test_require_all_can_be_relaxed(self):
        candidate = row()
        other = row(query_gene="WOX13", original_id="Solyc02g082670.1", gene_key="Solyc02g082670", new_id="Sly|Solanum_lycopersicum|Solyc02g082670.1")
        _rows, failures = check.run_check(
            [candidate, other],
            {candidate.new_id: "MKTF", other.new_id: "MKKK"},
            {candidate.new_id: "ATGAAAACCTTC", other.new_id: "ATGAAAAAAAAA"},
            [("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "ATGAAAACCTTC")],
            98,
            98,
            False,
        )
        self.assertEqual(failures, [])

    def test_unmatched_header_and_bad_cds(self):
        candidate = row()
        rows, failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKTF"},
            {candidate.new_id: "ATGAAAACCTTC"},
            [
                ("NOT_A_GENE", "Solanum lycopersicum", "NOT_A_GENE Solanum lycopersicum", "ATGAAAACCTTC"),
                ("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "ATGAAAACCT"),
            ],
            98,
            98,
            True,
        )
        statuses = {item["status"] for item in rows}
        self.assertIn("unmatched_header", statuses)
        self.assertIn("bad_cds", statuses)
        self.assertGreaterEqual(len(failures), 2)

    def test_internal_stop_and_duplicate_header(self):
        candidate = row()
        rows, _failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKTF"},
            {candidate.new_id: "ATGAAAACCTTC"},
            [
                ("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "ATGTAAACCTTC"),
                ("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum copy", "ATGAAAACCTTC"),
            ],
            98,
            98,
            True,
        )
        statuses = [item["status"] for item in rows]
        self.assertIn("bad_cds", statuses)
        self.assertIn("duplicate_header", statuses)

    def test_missing_annotation_sequences(self):
        candidate = row()
        rows, failures = check.run_check(
            [candidate],
            {},
            {},
            [("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "ATGAAAACCTTC")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "missing_annotation_protein")
        self.assertTrue(failures)

        rows, _failures = check.run_check(
            [candidate],
            {candidate.new_id: "MKTF"},
            {},
            [("WOX4", "Solanum lycopersicum", "WOX4 Solanum lycopersicum", "ATGAAAACCTTC")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "missing_annotation_cds")

    def test_unresolved_name_does_not_fail_without_a_reference(self):
        missing = row(query_gene="GH9B3", match_method="unmatched", original_id="", gene_key="", new_id="")
        _rows, failures = check.run_check([missing], {}, {}, [], 98, 98, True)
        self.assertEqual(failures, [])

    def test_reference_for_an_unresolved_name_fails(self):
        missing = row(
            query_gene="GH9B3",
            match_method="unmatched",
            original_id="",
            gene_key="",
            new_id="",
        )
        rows, failures = check.run_check(
            [missing],
            {},
            {},
            [("GH9B3", "Solanum lycopersicum", "GH9B3 Solanum lycopersicum", "ATGAAAACCTTC")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "unresolved")
        self.assertTrue(failures)

    def test_cli_writes_a_passing_summary_and_a_strict_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = row()
            crosswalk = root / "crosswalk.tsv"
            write_crosswalk(
                crosswalk,
                [
                    {
                        "query_gene": candidate.query_gene,
                        "query_taxon": candidate.query_taxon,
                        "species_id": candidate.species_id,
                        "match_method": candidate.match_method,
                        "original_id": candidate.original_id,
                        "gene_key": candidate.gene_key,
                        "new_id": candidate.new_id,
                    }
                ],
            )
            proteins = root / "proteins.faa"
            cds = root / "cds.fna"
            reference = root / "reference.fna"
            write_fasta(proteins, [(candidate.new_id, "MKTF")])
            write_fasta(cds, [(candidate.new_id, "ATGAAAACCTTC")])
            write_fasta(reference, [("WOX4 Solanum lycopersicum", "ATGAAAACCTTCTAA")])
            report = root / "report.tsv"
            summary = root / "summary.txt"
            code = check.main(
                [
                    "--crosswalk",
                    str(crosswalk),
                    "--proteins",
                    str(proteins),
                    "--cds",
                    str(cds),
                    "--reference-cds",
                    str(reference),
                    "--require-all",
                    "--out",
                    str(report),
                    "--summary",
                    str(summary),
                ]
            )
            self.assertEqual(code, 0)
            self.assertTrue(summary.read_text().startswith("pass"))
            with report.open() as fh:
                body = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(body[0]["status"], "match")

            write_fasta(reference, [("WOX4 Solanum lycopersicum", "ATGAAAAAAAAA")])
            code = check.main(
                [
                    "--crosswalk",
                    str(crosswalk),
                    "--proteins",
                    str(proteins),
                    "--cds",
                    str(cds),
                    "--reference-cds",
                    str(reference),
                    "--require-all",
                    "--strict",
                    "--out",
                    str(report),
                    "--summary",
                    str(summary),
                ]
            )
            self.assertEqual(code, 1)
            self.assertTrue(summary.read_text().startswith("fail"))
            self.assertEqual(report.read_text().splitlines()[1].split("\t")[12], "mismatch")

    def test_header_parser(self):
        self.assertEqual(check.parse_header("WOX4 Solanum lycopersicum"), ("WOX4", "Solanum lycopersicum"))
        self.assertEqual(check.parse_header("WOX4|Solanum_lycopersicum"), ("WOX4", "Solanum lycopersicum"))
        self.assertEqual(
            check.parse_header("lcl|WOX4 Solanum lycopersicum"),
            ("WOX4", "Solanum lycopersicum"),
        )
        self.assertEqual(
            check.parse_header("Sly|Solanum_lycopersicum|Solyc04g078650.4.1"),
            ("Solyc04g078650.4.1", "Solanum lycopersicum"),
        )
        tabbed = "WOX13\tArabidopsis thaliana AT4G35550.1 TAIR10"
        self.assertEqual(check.parse_header(tabbed), ("WOX13", "Arabidopsis thaliana"))
        self.assertEqual(
            check.parse_header_fields(tabbed),
            ("WOX13", "Arabidopsis thaliana", "AT4G35550.1"),
        )
        self.assertEqual(
            check.parse_header_fields("WOX4\tSolanum lycopersicum Solyc04g078650.4.1 ITAG4.0"),
            ("WOX4", "Solanum lycopersicum", "Solyc04g078650.4.1"),
        )

    def test_tabbed_header_matches_symbol_or_locus(self):
        header = "WOX13\tArabidopsis thaliana AT4G35550.1 TAIR10"
        by_symbol = row(
            query_gene="WOX13",
            query_taxon="Arabidopsis thaliana",
            species_id="Ath",
            original_id="AT4G35550.1",
            gene_key="AT4G35550",
            new_id="Ath|Arabidopsis_thaliana|AT4G35550.1",
        )
        rows, failures = check.run_check(
            [by_symbol],
            {by_symbol.new_id: "MKTF"},
            {by_symbol.new_id: "ATGAAAACCTTC"},
            [("ignored", "ignored", header, "ATGAAAACCTTCTAA")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "match")
        self.assertEqual(failures, [])

        by_locus = row(
            query_gene="AT4G35550",
            query_taxon="Arabidopsis thaliana",
            species_id="Ath",
            match_method="locus",
            original_id="AT4G35550.1",
            gene_key="AT4G35550",
            new_id="Ath|Arabidopsis_thaliana|AT4G35550.1",
        )
        rows, failures = check.run_check(
            [by_locus],
            {by_locus.new_id: "MKTF"},
            {by_locus.new_id: "ATGAAAACCTTC"},
            [("ignored", "ignored", header, "ATGAAAACCTTCTAA")],
            98,
            98,
            True,
        )
        self.assertEqual(rows[0]["status"], "match")
        self.assertEqual(failures, [])

    def test_fasta_reader_keeps_a_tabbed_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reference.fna"
            path.write_text(">WOX13\tArabidopsis thaliana AT4G35550.1 TAIR10\nATGGCTTCTTCTTCTTCTTCTTAA\n")
            gene, taxon, header, _seq = check.load_reference(str(path))[0]
        self.assertEqual(gene, "WOX13")
        self.assertEqual(taxon, "Arabidopsis thaliana")
        self.assertIn("AT4G35550.1", header)
        self.assertNotIn("TAIR10", taxon)


if __name__ == "__main__":
    unittest.main()
