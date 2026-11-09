#!/usr/bin/env python3
"""The graft table uses Rootstock, Scion, and TRUE/FALSE, matched to species in the run."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import derive_graft_traits as traits  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]


class GraftOutcomeFormatTest(unittest.TestCase):
    def test_dated_table_uses_scientific_names_and_true_false(self) -> None:
        pairs = traits.load_matrix(str(ROOT / "data" / "graft_outcomes_11-09-2026.tsv"))
        self.assertGreater(len(pairs), 3000)
        first = pairs[0]
        self.assertEqual(first["rootstock"], "Crossandra infundibuliformis")
        self.assertEqual(first["scion"], "Petunia axillaris × P. inflata")
        self.assertEqual(first["outcome"], "compatible")
        self.assertIn("incompatible", {pair["outcome"] for pair in pairs})

    def test_hybrid_does_not_match_a_parent_species(self) -> None:
        catalog = traits.load_species_catalog(str(ROOT / "conf" / "species.yaml"))
        ids = {row["id"] for row in catalog}
        binomials = traits.binomial_index(catalog)
        self.assertEqual(
            traits.pipeline_species_id("Solanum lycopersicum", ids, binomials),
            "Sly",
        )
        self.assertEqual(
            traits.pipeline_species_id("Sly.clean", ids, binomials),
            "Sly",
        )
        hybrid = "Solanum lycopersicum × S. pimpinellifolium"
        self.assertEqual(traits.pipeline_species_id(hybrid, ids, binomials), hybrid)

    def test_coverage_log_marks_species_missing_from_the_table(self) -> None:
        pairs = [
            {"scion": "Sly", "rootstock": "Nbe", "outcome": "compatible"},
            {"scion": "Ath", "rootstock": "Sly", "outcome": "incompatible"},
            {"scion": "Petunia axillaris × P. inflata", "rootstock": "Nbe", "outcome": "compatible"},
        ]
        analyzed = [
            {"species_id": "Sly", "scientific_name": "Solanum lycopersicum"},
            {"species_id": "Nbe", "scientific_name": "Nicotiana benthamiana"},
            {"species_id": "Cpe", "scientific_name": "Cucurbita pepo"},
        ]
        rows = {row["species_id"]: row for row in traits.coverage_rows(analyzed, pairs)}
        self.assertEqual(rows["Sly"]["in_graft_outcomes"], "yes")
        self.assertEqual(rows["Sly"]["n_as_scion"], 1)
        self.assertEqual(rows["Sly"]["n_as_rootstock"], 1)
        self.assertEqual(rows["Sly"]["n_success"], 1)
        self.assertEqual(rows["Sly"]["n_failure"], 1)
        self.assertEqual(rows["Cpe"]["in_graft_outcomes"], "no")
        self.assertEqual(rows["Cpe"]["n_as_scion"], 0)

    def test_a_species_with_no_scion_grafts_is_not_a_broad_scion(self) -> None:
        self.assertEqual(traits.is_broad([], 0.5), 0)
        self.assertEqual(traits.is_broad([True, True, True, True], 0.5), 1)
        self.assertEqual(traits.is_broad([True, False, False], 0.5), 0)

    def test_example_file_uses_the_same_columns(self) -> None:
        pairs = traits.load_matrix(str(ROOT / "data" / "graft_outcomes.example.tsv"))
        self.assertEqual(pairs[0]["rootstock"], "Solanum lycopersicum")
        self.assertEqual(pairs[0]["scion"], "Nicotiana benthamiana")
        self.assertEqual(pairs[0]["outcome"], "compatible")
        self.assertEqual(pairs[1]["outcome"], "incompatible")


if __name__ == "__main__":
    unittest.main()
