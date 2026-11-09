#!/usr/bin/env python3
"""The structure page lists each family's domain ranges for Mol*."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import write_summary_pages as pages  # noqa: E402


class StructurePageTest(unittest.TestCase):
    def test_structure_page_carries_domains_and_the_viewer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alignment = root / "OG0000001.fa"
            alignment.write_text(">Ath|Arabidopsis_thaliana|NP_1\nMKNGTF\n>Sly|Solanum_lycopersicum|XP_1\nMKNGTF\n")
            families = root / "families.tsv"
            families.write_text("Orthogroup\tgenes_of_interest\nOG0000001\tWOX4\n")
            interpro = root / "interpro.tsv"
            interpro.write_text(
                "Ath|Arabidopsis_thaliana|NP_1\tmd5\t6\tPfam\tPF00046\tHomeobox\t1\t4\t0\tT\t01-01-2020\tIPR001356\tHomeobox domain\n"
            )
            index = pages.parse_interpro(interpro)
            html = pages.render_structures(
                ["OG0000001"],
                root,
                index,
                {},
                pages.orthogroup_titles(str(families)),
            )
        self.assertIn("structures.html", html)
        self.assertIn("molstar@5.13.1", html)
        self.assertIn("OG0000001", html)
        self.assertIn("WOX4", html)
        self.assertIn("Homeobox domain", html)
        start = html.index('id="structure-catalog">') + len('id="structure-catalog">')
        end = html.index("</script>", start)
        catalog = json.loads(html[start:end])
        self.assertEqual(catalog["families"][0]["members"][0]["domains"][0]["start"], 1)
        self.assertEqual(catalog["families"][0]["members"][0]["domains"][0]["end"], 4)
        self.assertIn("pLDDT confidence", html)


if __name__ == "__main__":
    unittest.main()
