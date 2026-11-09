#!/usr/bin/env python3
"""Prepare the species tree and foreground list for PCOC and CSUBST.

PCOC numbers every node in postorder and, with --noplot, renames each leaf to
``name_<node id>``. A convergent scenario is those node ids, one independent
transition per foreground species, separated by ``/``.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_NUMBERED_LEAF = re.compile(r"([A-Za-z0-9.|-]+)_(\d+)(?=:|,|\))")


def fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            ids.append(line[1:].split()[0])
    return ids


def read_names(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text().splitlines():
        name = line.strip()
        if name and not name.startswith("#"):
            names.append(name.split()[0])
    return names


def on_tree(names: list[str], tips: list[str]) -> list[str]:
    present = set(tips)
    seen: set[str] = set()
    kept: list[str] = []
    for name in names:
        if name in present and name not in seen:
            seen.add(name)
            kept.append(name)
    return kept


def scenario_from_numbered(text: str, foreground: list[str]) -> str:
    """Map PCOC's ``leaf_<node id>`` names onto one event per foreground tip."""
    wanted = set(foreground)
    seen: set[str] = set()
    node_ids: list[str] = []
    for match in _NUMBERED_LEAF.finditer(text):
        leaf, node_id = match.group(1), match.group(2)
        if leaf in wanted and node_id not in seen:
            seen.add(node_id)
            node_ids.append(node_id)
    return "/".join(node_ids)


def _write_status(path: Path, rows: dict[str, str]) -> None:
    path.write_text("".join(f"{key}\t{value}\n" for key, value in rows.items()))


def root_and_prune(tree_path: Path, keep: list[str], out_path: Path) -> list[str]:
    from Bio import Phylo

    trees = list(Phylo.parse(tree_path, "newick"))
    if not trees or trees[0].count_terminals() == 0:
        raise ValueError(f"no tips in {tree_path}")
    tree = trees[0]
    keep_ids = set(keep)
    for tip in list(tree.get_terminals()):
        if tip.name not in keep_ids:
            tree.prune(tip)
    tips = [tip.name for tip in tree.get_terminals() if tip.name]
    if tree.count_terminals() >= 3:
        tree.root_at_midpoint()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, out_path, "newick")
    return tips


def prepare(args: argparse.Namespace) -> int:
    status = args.out_status
    try:
        protein_ids = fasta_ids(args.protein_fasta)
        tips = root_and_prune(args.tree, protein_ids, args.out_tree)
        foreground = on_tree(read_names(args.foreground), tips)
        args.out_foreground.write_text("".join(f"{name}\n" for name in foreground))
        n_codon = 0
        n_codon_fg = 0
        if args.codon_fasta and args.codon_fasta.is_file() and args.codon_fasta.stat().st_size:
            codon_ids = fasta_ids(args.codon_fasta)
            codon_tips = root_and_prune(args.tree, codon_ids, args.out_codon_tree)
            codon_fg = on_tree(foreground, codon_tips)
            args.out_codon_foreground.write_text("".join(f"{name}\n" for name in codon_fg))
            n_codon = len(codon_tips)
            n_codon_fg = len(codon_fg)
        _write_status(
            status,
            {
                "status": "ok",
                "note": "",
                "n_tips": str(len(tips)),
                "n_foreground": str(len(foreground)),
                "n_codon": str(n_codon),
                "n_codon_foreground": str(n_codon_fg),
            },
        )
    except Exception as exc:
        _write_status(
            status,
            {
                "status": "failed",
                "note": str(exc).replace("\t", " ")[:180],
                "n_tips": "0",
                "n_foreground": "0",
                "n_codon": "0",
                "n_codon_foreground": "0",
            },
        )
    return 0


def write_scenario(args: argparse.Namespace) -> int:
    text = args.numbered.read_text()
    foreground = read_names(args.foreground)
    args.out_scenario.write_text(scenario_from_numbered(text, foreground) + "\n")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    prep = sub.add_parser("prepare")
    prep.add_argument("--tree", type=Path, required=True)
    prep.add_argument("--foreground", type=Path, required=True)
    prep.add_argument("--protein-fasta", type=Path, required=True)
    prep.add_argument("--out-tree", type=Path, required=True)
    prep.add_argument("--out-foreground", type=Path, required=True)
    prep.add_argument("--out-status", type=Path, required=True)
    prep.add_argument("--codon-fasta", type=Path)
    prep.add_argument("--out-codon-tree", type=Path)
    prep.add_argument("--out-codon-foreground", type=Path)
    prep.set_defaults(func=prepare)

    scenario = sub.add_parser("scenario")
    scenario.add_argument("--numbered", type=Path, required=True)
    scenario.add_argument("--foreground", type=Path, required=True)
    scenario.add_argument("--out-scenario", type=Path, required=True)
    scenario.set_defaults(func=write_scenario)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
