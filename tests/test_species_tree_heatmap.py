#!/usr/bin/env python3
"""Identity heatmap order follows the OrthoFinder species tree."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import curate_results as report  # noqa: E402


def matrix() -> dict:
    return {
        "labels": ["Nbe:a", "Sly:para", "Ath:a", "Sly:ref", "Osa:a", "Cpe:b", "Cpe:a"],
        "full": [
            "Nbe|Nicotiana_benthamiana|a",
            "Sly|Solanum_lycopersicum|para",
            "Ath|Arabidopsis_thaliana|a",
            "Sly|Solanum_lycopersicum|ref",
            "Osa|Oryza_sativa|a",
            "Cpe|Cucurbita_pepo|b",
            "Cpe|Cucurbita_pepo|a",
        ],
        "species": ["Nbe", "Sly", "Ath", "Sly", "Osa", "Cpe", "Cpe"],
        "reference": [False, False, False, True, False, False, False],
        "values": [
            [None, 1, 2, 3, 4, 5, 6],
            [1, None, 7, 8, 9, 10, 11],
            [2, 7, None, 12, 13, 14, 15],
            [3, 8, 12, None, 16, 17, 18],
            [4, 9, 13, 16, None, 19, 20],
            [5, 10, 14, 17, 19, None, 21],
            [6, 11, 15, 18, 20, 21, None],
        ],
    }


class SpeciesTreeHeatmapTest(unittest.TestCase):
    def test_tip_names_match_species_ids(self) -> None:
        known = {"Ath", "Sly", "Nbe", "Cpe"}
        self.assertEqual(report.match_tip("Sly.clean", known), "Sly")
        self.assertEqual(report.match_tip("Sly.clean.faa", known), "Sly")
        self.assertEqual(report.match_tip("Ath_clean", known), "Ath")
        self.assertEqual(report.match_tip("'Nbe.clean'", known), "Nbe")
        self.assertEqual(report.match_tip("Cpe", known), "Cpe")
        self.assertIsNone(report.match_tip("ghost.clean", known))

    def test_tree_order_and_reorder(self) -> None:
        known = {"Ath", "Sly", "Nbe", "Cpe", "Osa"}
        newick = "((Sly.clean:0.03,Nbe.clean:0.04)1.00:0.22,(Ath.clean:0.04,Cpe.clean:0.05)0.95:0.25);"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "SpeciesTree_rooted.txt"
            path.write_text(newick + "\n")
            tree = report.load_species_tree(str(path), known)
        self.assertIsNotNone(tree)
        assert tree is not None
        # Smaller clades are not smaller here; ladderize keeps the original
        # child order when both clades have the same number of tips, and that
        # original order lists the Solanum clade first.
        self.assertEqual(tree["order"], ["Sly", "Nbe", "Ath", "Cpe"])

        def tip_length(node: dict, name: str) -> float | None:
            if not node["children"] and node["name"] == name:
                return node["length"]
            for child in node["children"]:
                found = tip_length(child, name)
                if found is not None:
                    return found
            return None

        self.assertAlmostEqual(tip_length(tree["root"], "Ath"), 0.04)
        ordered = report.reorder_matrix(matrix(), tree["order"])
        self.assertEqual(ordered["species"], ["Sly", "Sly", "Nbe", "Ath", "Cpe", "Cpe", "Osa"])
        self.assertEqual(ordered["reference"], [True, False, False, False, False, False, False])
        self.assertEqual(ordered["labels"][0], "Sly:ref")
        self.assertEqual(ordered["labels"][4], "Cpe:a")
        self.assertEqual(ordered["labels"][5], "Cpe:b")
        # Ath (index 3) vs the Sly reference (index 0) stays 12.
        self.assertEqual(ordered["values"][3][0], 12)
        self.assertIsNone(ordered["values"][1][1])

    def test_unmatched_tip_is_dropped(self) -> None:
        known = {"Ath", "Sly"}
        newick = "((Sly.clean:0.05,ghost.clean:0.06):0.2,Ath.clean:0.3);"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "SpeciesTree_rooted.txt"
            path.write_text(newick + "\n")
            tree = report.load_species_tree(str(path), known)
        self.assertIsNotNone(tree)
        assert tree is not None
        self.assertEqual(tree["order"], ["Ath", "Sly"])
        names = []

        def walk(node: dict) -> None:
            if node["children"]:
                for child in node["children"]:
                    walk(child)
            elif node["name"]:
                names.append(node["name"])

        walk(tree["root"])
        self.assertEqual(names, ["Ath", "Sly"])

    def test_report_embeds_tree_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            families = root / "families.tsv"
            wide = root / "wide.tsv"
            labels = root / "labels.tsv"
            pairs = root / "pairs.tsv"
            tree = root / "SpeciesTree_rooted.txt"
            families.write_text(
                "Orthogroup\tgenes_of_interest\tquery_species_ids\tn_members\t"
                "n_species_with_members\tmean_pct_id_to_ref\tmin_pct_id_to_ref\tmax_pct_id_to_ref\n"
                "OG0001\tWOX4\tSly\t4\t3\t70\t50\t100\n"
            )
            wide.write_text(
                "Orthogroup\tbest_pct_id_Ath\tbest_pct_id_Sly\tbest_pct_id_Nbe\n"
                "OG0001\t80\t100\t90\n"
            )
            labels.write_text(
                "Orthogroup\tlabel\tsequence_id\tspecies_id\tis_reference\n"
                "OG0001\tNbe:a\tNbe|Nicotiana_benthamiana|a\tNbe\t\n"
                "OG0001\tSly:ref\tSly|Solanum_lycopersicum|ref\tSly\tyes\n"
                "OG0001\tAth:a\tAth|Arabidopsis_thaliana|a\tAth\t\n"
            )
            pairs.write_text(
                "Orthogroup\tseq_a\tseq_b\tlabel_a\tlabel_b\tspecies_a\tspecies_b\tpct_identity\treference_id\n"
                "OG0001\tNbe|Nicotiana_benthamiana|a\tSly|Solanum_lycopersicum|ref\tNbe:a\tSly:ref\tNbe\tSly\t70\tSly|Solanum_lycopersicum|ref\n"
                "OG0001\tNbe|Nicotiana_benthamiana|a\tAth|Arabidopsis_thaliana|a\tNbe:a\tAth:a\tNbe\tAth\t60\tSly|Solanum_lycopersicum|ref\n"
                "OG0001\tSly|Solanum_lycopersicum|ref\tAth|Arabidopsis_thaliana|a\tSly:ref\tAth:a\tSly\tAth\t55\tSly|Solanum_lycopersicum|ref\n"
            )
            tree.write_text("((Ath.clean:0.1,Sly.clean:0.2):0.3,Nbe.clean:0.4);\n")
            out = root / "report"
            sys.argv = [
                "curate_results.py",
                "--families",
                str(families),
                "--wide-identity",
                str(wide),
                "--identity-matrix-pairs",
                str(pairs),
                "--identity-matrix-labels",
                str(labels),
                "--species-tree",
                str(tree),
                "--outdir",
                str(out),
            ]
            report.main()
            html = (out / "index.html").read_text()
            payload = html.split('<script id="payload" type="application/json">', 1)[1]
            payload = payload.split("</script>", 1)[0]
            data = json.loads(payload)
        # Nbe is a single-tip clade, so ladderize places it before Ath+Sly.
        self.assertEqual(data["species_tree"]["order"], ["Nbe", "Ath", "Sly"])
        self.assertEqual(
            data["families"][0]["matrix"]["species"],
            ["Nbe", "Ath", "Sly"],
        )
        self.assertIn("OrthoFinder species tree", html)
        self.assertIn("alignTrees", html)


if __name__ == "__main__":
    unittest.main()
