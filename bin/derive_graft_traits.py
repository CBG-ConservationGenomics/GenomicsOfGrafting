#!/usr/bin/env python3
"""Derive lineage- and pair-level graft traits from a phenotype matrix.

Inputs
------
- Long TSV: scion, rootstock, outcome
  or square matrix TSV (rows=scion, cols=rootstock)
- OrthoFinder species tree (Newick); tip labels may be proteome filenames
  (e.g. Sly.clean) — we map to species ids.

Outputs
-------
- pair_outcomes.tsv          every directed scion→rootstock pair
- lineage_traits.tsv        per-species rates / binary broad-grafting flags
- phenotype_shift_branches.tsv   branches whose child clade differs in trait
- species_tree.annotated.nwk     tree with tip labels normalized to species ids
- foreground_tips.txt            tips with broad_scion (or chosen trait) = 1
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import Phylo


OUTCOME_MAP = {
    "compatible": "compatible",
    "compat": "compatible",
    "success": "compatible",
    "ok": "compatible",
    "incompatible": "incompatible",
    "incompat": "incompatible",
    "fail": "incompatible",
    "failure": "incompatible",
    "delayed_failure": "delayed_failure",
    "delayed": "delayed_failure",
    "delayed-failure": "delayed_failure",
    "late_failure": "delayed_failure",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True, help="Graft outcome matrix (long or square TSV)")
    p.add_argument("--species-tree", required=True, help="OrthoFinder SpeciesTree Newick")
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--broad-threshold",
        type=float,
        default=0.5,
        help="Min success rate (compatible only) to call broad_scion / broad_rootstock",
    )
    p.add_argument(
        "--success-outcomes",
        default="compatible",
        help="Comma-separated outcomes counted as success (default: compatible)",
    )
    p.add_argument(
        "--foreground-trait",
        default="broad_scion",
        help="Lineage trait used for foreground tip list (broad_scion|broad_rootstock|broad_either)",
    )
    return p.parse_args()


def normalize_outcome(x: str) -> str:
    key = re.sub(r"\s+", "_", x.strip().lower())
    if key not in OUTCOME_MAP:
        raise ValueError(f"Unknown graft outcome: {x!r}")
    return OUTCOME_MAP[key]


def normalize_species_id(label: str, known: set[str] | None = None) -> str:
    """Map OrthoFinder tip / matrix labels to short species ids when possible."""
    s = label.strip().strip("'\"")
    s = re.sub(r"\.(clean|faa|fasta|proteins)$", "", s, flags=re.I)
    if "|" in s:
        s = s.split("|")[0]
    # Genus_species → try first letters + species start, but prefer known ids
    if known:
        if s in known:
            return s
        for k in known:
            if s.lower() == k.lower():
                return k
            if k.lower() in s.lower() or s.lower() in k.lower():
                return k
    return s


def load_matrix(path: str) -> list[dict]:
    """Return list of {scion, rootstock, outcome}."""
    text = Path(path).read_text()
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        raise SystemExit(f"Empty phenotype matrix: {path}")

    # Detect delimiter
    dialect_delim = "\t" if "\t" in lines[0] else ","
    rows = list(csv.reader(lines, delimiter=dialect_delim))
    header = [h.strip() for h in rows[0]]

    long_cols = {c.lower() for c in header}
    if {"scion", "rootstock", "outcome"}.issubset(long_cols):
        # long format
        idx = {h.lower(): i for i, h in enumerate(header)}
        out = []
        for r in rows[1:]:
            if len(r) < 3:
                continue
            out.append(
                {
                    "scion": r[idx["scion"]].strip(),
                    "rootstock": r[idx["rootstock"]].strip(),
                    "outcome": normalize_outcome(r[idx["outcome"]]),
                }
            )
        return out

    # square matrix: first column = scion ids, header[1:] = rootstock ids
    rootstocks = header[1:]
    out = []
    for r in rows[1:]:
        if not r:
            continue
        scion = r[0].strip()
        for j, rs in enumerate(rootstocks):
            if j + 1 >= len(r):
                continue
            val = r[j + 1].strip()
            if not val or val in {".", "NA", "na", "-"}:
                continue
            out.append(
                {
                    "scion": scion,
                    "rootstock": rs.strip(),
                    "outcome": normalize_outcome(val),
                }
            )
    return out


def tip_species_map(tree: Phylo.BaseTree.Tree) -> dict[str, str]:
    """Original tip name → cleaned species id."""
    tips = [c.name for c in tree.get_terminals() if c.name]
    known = set()
    # first pass: strip extensions
    cleaned = {}
    for t in tips:
        cleaned[t] = normalize_species_id(t)
        known.add(cleaned[t])
    return cleaned


def clade_trait_majority(clade, tip_traits: dict[str, int]) -> float | None:
    tips = [t.name for t in clade.get_terminals() if t.name]
    vals = [tip_traits[t] for t in tips if t in tip_traits]
    if not vals:
        return None
    return sum(vals) / len(vals)


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    success = {x.strip() for x in args.success_outcomes.split(",") if x.strip()}

    pairs_raw = load_matrix(args.matrix)
    species_in_matrix = sorted(
        {p["scion"] for p in pairs_raw} | {p["rootstock"] for p in pairs_raw}
    )
    known = set(species_in_matrix)

    # Normalize matrix species ids
    pairs = []
    for p in pairs_raw:
        pairs.append(
            {
                "scion": normalize_species_id(p["scion"], known),
                "rootstock": normalize_species_id(p["rootstock"], known),
                "outcome": p["outcome"],
            }
        )

    tree = Phylo.read(args.species_tree, "newick")
    tip_map = tip_species_map(tree)
    # Relabel tips to species ids
    for tip in tree.get_terminals():
        if tip.name:
            tip.name = normalize_species_id(tip_map.get(tip.name, tip.name), known)

    # Pair table + reciprocity
    pair_index = {(p["scion"], p["rootstock"]): p["outcome"] for p in pairs}
    pair_rows = []
    for p in pairs:
        sc, rs, oc = p["scion"], p["rootstock"], p["outcome"]
        recip = pair_index.get((rs, sc), "")
        reciprocal_same = ""
        if recip:
            reciprocal_same = "yes" if recip == oc else "no"
        pair_rows.append(
            {
                "scion": sc,
                "rootstock": rs,
                "outcome": oc,
                "success": "yes" if oc in success else "no",
                "reciprocal_outcome": recip,
                "reciprocal_same": reciprocal_same,
                "is_self": "yes" if sc == rs else "no",
            }
        )

    with open(out / "pair_outcomes.tsv", "w", newline="") as fh:
        fields = [
            "scion",
            "rootstock",
            "outcome",
            "success",
            "reciprocal_outcome",
            "reciprocal_same",
            "is_self",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    # Lineage traits
    as_scion = defaultdict(list)
    as_stock = defaultdict(list)
    for p in pair_rows:
        if p["is_self"] == "yes":
            continue
        as_scion[p["scion"]].append(p["success"] == "yes")
        as_stock[p["rootstock"]].append(p["success"] == "yes")

    species = sorted(set(as_scion) | set(as_stock) | {t.name for t in tree.get_terminals() if t.name})
    lineage_rows = []
    tip_foreground = {}
    for sp in species:
        sc_vals = as_scion.get(sp, [])
        rs_vals = as_stock.get(sp, [])
        sc_rate = sum(sc_vals) / len(sc_vals) if sc_vals else float("nan")
        rs_rate = sum(rs_vals) / len(rs_vals) if rs_vals else float("nan")
        broad_sc = int(sc_vals and sc_rate >= args.broad_threshold)
        broad_rs = int(rs_vals and rs_rate >= args.broad_threshold)
        broad_either = int(broad_sc or broad_rs)
        lineage_rows.append(
            {
                "species_id": sp,
                "n_as_scion": len(sc_vals),
                "n_as_rootstock": len(rs_vals),
                "scion_success_rate": f"{sc_rate:.4f}" if sc_vals else "",
                "rootstock_success_rate": f"{rs_rate:.4f}" if rs_vals else "",
                "broad_scion": broad_sc,
                "broad_rootstock": broad_rs,
                "broad_either": broad_either,
                "broad_threshold": args.broad_threshold,
            }
        )
        if args.foreground_trait == "broad_scion":
            tip_foreground[sp] = broad_sc
        elif args.foreground_trait == "broad_rootstock":
            tip_foreground[sp] = broad_rs
        else:
            tip_foreground[sp] = broad_either

    with open(out / "lineage_traits.tsv", "w", newline="") as fh:
        fields = list(lineage_rows[0].keys()) if lineage_rows else ["species_id"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(lineage_rows)

    fg_tips = [sp for sp, v in tip_foreground.items() if v == 1]
    (out / "foreground_tips.txt").write_text("\n".join(fg_tips) + ("\n" if fg_tips else ""))

    # Phenotype-shift branches: parent clade mean vs child clade mean of binary trait
    shift_rows = []
    # assign clade names
    for i, clade in enumerate(tree.find_clades(order="level")):
        if clade.is_terminal():
            continue
        parent = None
        # Bio.Phylo doesn't give parent easily — compute child vs sibling contrast via path
    # Walk all non-root clades
    def all_clades_with_parent(tree):
        for path in tree.get_nonterminals(order="preorder"):
            pass
        # Use root traversal
        stack = [(tree.root, None)]
        while stack:
            node, parent = stack.pop()
            yield node, parent
            for ch in node.clades:
                stack.append((ch, node))

    for node, parent in all_clades_with_parent(tree):
        if parent is None or node.is_terminal():
            continue
        child_mean = clade_trait_majority(node, tip_foreground)
        parent_mean = clade_trait_majority(parent, tip_foreground)
        if child_mean is None or parent_mean is None:
            continue
        # shift if crossing 0.5 boundary or absolute delta large
        child_call = 1 if child_mean >= 0.5 else 0
        parent_call = 1 if parent_mean >= 0.5 else 0
        if child_call != parent_call:
            tips = ",".join(sorted(t.name for t in node.get_terminals() if t.name))
            shift_rows.append(
                {
                    "branch_child_tips": tips,
                    "n_tips": len(list(node.get_terminals())),
                    "parent_trait_mean": f"{parent_mean:.3f}",
                    "child_trait_mean": f"{child_mean:.3f}",
                    "parent_call": parent_call,
                    "child_call": child_call,
                    "trait": args.foreground_trait,
                    "direction": f"{parent_call}->{child_call}",
                }
            )

    with open(out / "phenotype_shift_branches.tsv", "w", newline="") as fh:
        fields = [
            "branch_child_tips",
            "n_tips",
            "parent_trait_mean",
            "child_trait_mean",
            "parent_call",
            "child_call",
            "trait",
            "direction",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(shift_rows)

    Phylo.write(tree, out / "species_tree.annotated.nwk", "newick")

    # Also write a simple species distance matrix (patristic) for Mantel covariate
    terminals = [t for t in tree.get_terminals() if t.name]
    names = [t.name for t in terminals]
    with open(out / "species_distance.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species"] + names)
        for a in terminals:
            row = [a.name]
            for b in terminals:
                if a is b:
                    row.append("0")
                else:
                    try:
                        d = tree.distance(a, b)
                    except Exception:
                        d = ""
                    row.append(f"{d:.6f}" if d != "" else "")
            w.writerow(row)

    print(
        f"[derive_graft_traits] pairs={len(pair_rows)} species={len(lineage_rows)} "
        f"shift_branches={len(shift_rows)} foreground_tips={len(fg_tips)} → {out}"
    )


if __name__ == "__main__":
    main()
