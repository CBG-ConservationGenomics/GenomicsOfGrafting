#!/usr/bin/env python3
"""Parse IQ-TREE ancestral reconstruction and map AA changes onto branches.

Flags substitutions on phenotype-shift branches when tip→species mapping allows.

Expects IQ-TREE outputs from: iqtree -s aln -m MFP -asr
  - *.state  (ancestral state file)
  - *.treefile
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import Phylo


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--state", required=True, help="IQ-TREE .state file")
    p.add_argument("--tree", required=True, help="IQ-TREE .treefile with node labels")
    p.add_argument("--orthogroup", required=True)
    p.add_argument("--shift-branches", required=True, help="phenotype_shift_branches.tsv")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--min-prob", type=float, default=0.9)
    return p.parse_args()


def load_shift_tipsets(path: str) -> list[set[str]]:
    sets = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tips = {t for t in row["branch_child_tips"].split(",") if t}
            if tips:
                sets.append(tips)
    return sets


def species_of(seq_id: str) -> str:
    if "|" in seq_id:
        return seq_id.split("|")[0]
    return re.sub(r"\.(clean|faa|fasta)$", "", seq_id, flags=re.I)


def parse_state_file(path: Path) -> dict[str, dict[int, tuple[str, float]]]:
    """node -> {site(1-based): (aa, prob)}"""
    nodes: dict[str, dict[int, tuple[str, float]]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("Node"):
                # header may be: Node Site State p_A p_C ...
                if line.startswith("Node"):
                    continue
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            node, site_s, state = parts[0], parts[1], parts[2]
            try:
                site = int(site_s)
            except ValueError:
                continue
            # probability of called state if present
            prob = 1.0
            if len(parts) > 3:
                # find p_<state>
                for tok in parts[3:]:
                    if tok.upper().startswith(f"P_{state.upper()}") or tok.startswith(f"p_{state}"):
                        try:
                            prob = float(tok.split("=")[-1]) if "=" in tok else float(tok)
                        except ValueError:
                            pass
                # IQ-TREE format: columns p_A p_R p_N ... after State
                # Use max of probability columns if numeric
                nums = []
                for tok in parts[3:]:
                    try:
                        nums.append(float(tok))
                    except ValueError:
                        continue
                if nums:
                    prob = max(nums)
            nodes.setdefault(node, {})[site] = (state, prob)
    return nodes


def clade_tip_species(clade) -> set[str]:
    return {species_of(t.name) for t in clade.get_terminals() if t.name}


def main() -> None:
    args = parse_args()
    shift_sets = load_shift_tipsets(args.shift_branches)
    tree = Phylo.read(args.tree, "newick")
    states = parse_state_file(Path(args.state))

    # Map internal clade (by tip set) to node name if IQ-TREE labeled nodes
    # IQ-TREE -asr uses Node1.. in .state matching labels in treefile
    node_by_tips: dict[frozenset[str], str] = {}
    for clade in tree.find_clades():
        if clade.name and not clade.is_terminal():
            tips = frozenset(clade_tip_species(clade))
            node_by_tips[tips] = clade.name
        elif clade.is_terminal() and clade.name:
            node_by_tips[frozenset([species_of(clade.name)])] = clade.name

    # Parent→child contrasts along tree
    rows = []

    def walk(node, parent):
        if parent is not None and parent.name and node.name:
            p_states = states.get(parent.name) or states.get(str(parent.name))
            c_states = states.get(node.name) or states.get(str(node.name))
            # Also try NodeX labels from state keys if names differ
            if p_states and c_states:
                child_tips = clade_tip_species(node)
                on_shift = "yes" if any(s == child_tips or s.issubset(child_tips) or child_tips.issubset(s) for s in shift_sets) else "no"
                sites = sorted(set(p_states) & set(c_states))
                for site in sites:
                    paa, pprob = p_states[site]
                    caa, cprob = c_states[site]
                    if paa == caa or paa == "-" or caa == "-":
                        continue
                    if pprob < args.min_prob or cprob < args.min_prob:
                        continue
                    rows.append(
                        {
                            "Orthogroup": args.orthogroup,
                            "parent_node": parent.name,
                            "child_node": node.name,
                            "child_tips": ",".join(sorted(child_tips)),
                            "site": site,
                            "parent_aa": paa,
                            "child_aa": caa,
                            "parent_prob": f"{pprob:.3f}",
                            "child_prob": f"{cprob:.3f}",
                            "on_phenotype_shift_branch": on_shift,
                            "substitution": f"{paa}{site}{caa}",
                        }
                    )
        for ch in node.clades:
            walk(ch, node)

    walk(tree.root, None)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "Orthogroup",
        "parent_node",
        "child_node",
        "child_tips",
        "site",
        "parent_aa",
        "child_aa",
        "parent_prob",
        "child_prob",
        "on_phenotype_shift_branch",
        "substitution",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_shift = sum(1 for r in rows if r["on_phenotype_shift_branch"] == "yes")
    print(f"[parse_iqtree_asr] changes={len(rows)} on_shift_branches={n_shift} → {args.out_tsv}")


if __name__ == "__main__":
    main()
