#!/usr/bin/env python3
"""Derive lineage- and pair-level graft traits from a phenotype matrix.

Inputs
------
- CSV or TSV: Rootstock, Scion, success (TRUE/FALSE), as in
  data/graft_outcomes_11-09-2026.tsv
- Or the older long table: scion, rootstock, outcome
- Or a square matrix (rows=scion, columns=rootstock)
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
import sys
from collections import defaultdict
from pathlib import Path


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
    "true": "compatible",
    "yes": "compatible",
    "y": "compatible",
    "t": "compatible",
    "1": "compatible",
    "false": "incompatible",
    "no": "incompatible",
    "n": "incompatible",
    "f": "incompatible",
    "0": "incompatible",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True, help="Graft outcome table (Rootstock, Scion, success)")
    p.add_argument("--species-tree", required=True, help="OrthoFinder SpeciesTree Newick")
    p.add_argument(
        "--species-yaml",
        default=None,
        help="conf/species.yaml, used to match Genus species names to the ids in this run",
    )
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
    outcome_col = "outcome" if "outcome" in long_cols else "success" if "success" in long_cols else ""
    if {"scion", "rootstock"}.issubset(long_cols) and outcome_col:
        idx = {h.lower(): i for i, h in enumerate(header)}
        out = []
        for r in rows[1:]:
            if len(r) <= max(idx["scion"], idx["rootstock"], idx[outcome_col]):
                continue
            raw = r[idx[outcome_col]].strip()
            if not raw:
                continue
            out.append(
                {
                    "scion": r[idx["scion"]].strip(),
                    "rootstock": r[idx["rootstock"]].strip(),
                    "outcome": normalize_outcome(raw),
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


def load_species_catalog(path: str | None) -> list[dict]:
    """id, genus, and species from conf/species.yaml."""
    if not path:
        return []
    file = Path(path)
    if not file.is_file() or file.stat().st_size == 0:
        return []
    text = file.read_text()
    try:
        import yaml
    except ImportError:
        yaml = None
    if yaml is not None:
        data = yaml.safe_load(text) or {}
        rows = []
        for item in data.get("species") or []:
            if item.get("id") and item.get("genus") and item.get("species"):
                rows.append(
                    {
                        "id": str(item["id"]).strip(),
                        "genus": str(item["genus"]).strip(),
                        "species": str(item["species"]).strip(),
                    }
                )
        if rows:
            return rows
    rows = []
    current: dict | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        match = re.match(r"-\s+id:\s*(\S+)", stripped)
        if match:
            if current and current.get("id"):
                rows.append(current)
            current = {"id": match.group(1)}
            continue
        if current is None:
            continue
        field = re.match(r"(genus|species):\s*[\"']?([^\"'#]+)", stripped)
        if field:
            current[field.group(1)] = field.group(2).strip()
    if current and current.get("id"):
        rows.append(current)
    return [row for row in rows if row.get("genus") and row.get("species")]


def binomial_index(catalog: list[dict]) -> dict[str, str]:
    """'genus species' in lower case → pipeline species id. Exact names only."""
    index = {}
    for row in catalog:
        index[f"{row['genus']} {row['species']}".lower()] = row["id"]
    return index


def pipeline_species_id(label: str, ids: set[str], binomials: dict[str, str]) -> str:
    """Map a tree tip or a graft-table name onto a pipeline species id when it is one."""
    name = label.strip().strip("'\"")
    name = re.sub(r"\.(clean|faa|fasta|fa|fas|proteins|pep|aa|txt)$", "", name, flags=re.I)
    if "|" in name:
        name = name.split("|", 1)[0]
    for sid in ids:
        if name.lower() == sid.lower():
            return sid
    return binomials.get(" ".join(name.split()).lower(), name)


def coverage_rows(analyzed: list[dict], pairs: list[dict]) -> list[dict]:
    """One row per species in this run, and whether the graft table records it."""
    rows = []
    for species in analyzed:
        sid = species["species_id"]
        hits = [pair for pair in pairs if pair["scion"] == sid or pair["rootstock"] == sid]
        rows.append(
            {
                "species_id": sid,
                "scientific_name": species.get("scientific_name") or "",
                "in_graft_outcomes": "yes" if hits else "no",
                "n_as_scion": sum(pair["scion"] == sid for pair in hits),
                "n_as_rootstock": sum(pair["rootstock"] == sid for pair in hits),
                "n_success": sum(pair["outcome"] == "compatible" for pair in hits),
                "n_failure": sum(pair["outcome"] != "compatible" for pair in hits),
            }
        )
    return rows


def is_broad(values: list[bool], threshold: float) -> int:
    """1 when this species has grafts in the role and the success rate meets the threshold."""
    if not values:
        return 0
    rate = sum(values) / len(values)
    return int(rate >= threshold)


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

    from Bio import Phylo

    pairs_raw = load_matrix(args.matrix)
    catalog = load_species_catalog(args.species_yaml)
    catalog_ids = {row["id"] for row in catalog}
    binomials = binomial_index(catalog)
    names = {row["id"]: f"{row['genus']} {row['species']}" for row in catalog}

    def as_id(label: str) -> str:
        return pipeline_species_id(label, catalog_ids, binomials)

    pairs_all = [
        {"scion": as_id(pair["scion"]), "rootstock": as_id(pair["rootstock"]), "outcome": pair["outcome"]}
        for pair in pairs_raw
    ]

    tree = Phylo.read(args.species_tree, "newick")
    analyzed_ids = []
    for tip in tree.get_terminals():
        if not tip.name:
            continue
        sid = as_id(tip.name)
        tip.name = sid
        if sid not in analyzed_ids:
            analyzed_ids.append(sid)
    analyzed = [
        {"species_id": sid, "scientific_name": names.get(sid, "")}
        for sid in analyzed_ids
    ]
    analyzed_set = set(analyzed_ids)
    coverage = coverage_rows(analyzed, pairs_all)
    with open(out / "analyzed_species.tsv", "w", newline="") as fh:
        fields = [
            "species_id",
            "scientific_name",
            "in_graft_outcomes",
            "n_as_scion",
            "n_as_rootstock",
            "n_success",
            "n_failure",
        ]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(coverage)
    print(f"[derive_graft_traits] {len(coverage)} species analyzed in this run", file=sys.stderr)
    for row in coverage:
        name = row["scientific_name"] or row["species_id"]
        print(
            f"[derive_graft_traits] {row['species_id']} ({name}): "
            f"in graft outcomes: {row['in_graft_outcomes']}"
            f" (scion {row['n_as_scion']}, rootstock {row['n_as_rootstock']},"
            f" success {row['n_success']}, failure {row['n_failure']})",
            file=sys.stderr,
        )

    # Divergence and the trait table use pairs whose two species were both analyzed.
    pairs = [
        pair for pair in pairs_all if pair["scion"] in analyzed_set and pair["rootstock"] in analyzed_set
    ]
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
        broad_sc = is_broad(sc_vals, args.broad_threshold)
        broad_rs = is_broad(rs_vals, args.broad_threshold)
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
