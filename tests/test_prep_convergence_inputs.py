#!/usr/bin/env python3
"""PCOC scenario strings come from numbered leaf names, not the alignment."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import prep_convergence_inputs as prep  # noqa: E402


class ScenarioTest(unittest.TestCase):
    def test_one_event_per_foreground_tip(self) -> None:
        numbered = "((Sly_0:0.1,Ath_1:0.2)3:0.05,Nbe_2:0.3);"
        self.assertEqual(prep.scenario_from_numbered(numbered, ["Sly", "Nbe"]), "0/2")

    def test_species_absent_from_the_tree_are_left_out(self) -> None:
        numbered = "(Sly_4:0.2,Ath_5:0.2);"
        self.assertEqual(prep.scenario_from_numbered(numbered, ["Sly", "Cpe"]), "4")

    def test_foreground_keeps_file_order_and_tree_membership(self) -> None:
        self.assertEqual(prep.on_tree(["Cpe", "Sly", "Ath", "Sly"], ["Ath", "Sly"]), ["Sly", "Ath"])


if __name__ == "__main__":
    unittest.main()
