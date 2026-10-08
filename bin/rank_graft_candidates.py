#!/usr/bin/env python3
"""Merge ASR shift-branch hits, divergence tests, and optional selection/convergence into a ranked table."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--asr-dir", required=True, help="Directory of *.asr_substitutions.tsv")
    p.add_argument("--divergence-tests", required=True)
    p.add_argument("--selection-dir", default=None)
    p.add_argument("--convergence-dir", default=None)
    p.add_argument("--out-tsv", required=True)
    return p.parse_args()


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    args = parse_args()
    asr_dir = Path(args.asr_dir)

    shift_counts = defaultdict(int)
    shift_examples = defaultdict(list)
    for f in asr_dir.rglob("*.asr_substitutions.tsv"):
        for row in read_tsv(f):
            if row.get("on_phenotype_shift_branch") == "yes":
                og = row["Orthogroup"]
                shift_counts[og] += 1
                if len(shift_examples[og]) < 5:
                    shift_examples[og].append(row.get("substitution", ""))

    div = {r["Orthogroup"]: r for r in read_tsv(Path(args.divergence_tests))}

    sel_hit = defaultdict(int)
    if args.selection_dir:
        for f in Path(args.selection_dir).rglob("*.hyphy_summary.tsv"):
            for row in read_tsv(f):
                try:
                    if row.get("pvalue") and float(row["pvalue"]) < 0.05:
                        sel_hit[row["Orthogroup"]] += 1
                except ValueError:
                    continue

    conv_ok = defaultdict(int)
    if args.convergence_dir:
        for f in Path(args.convergence_dir).rglob("*.convergence_summary.tsv"):
            for row in read_tsv(f):
                if row.get("status") == "ok":
                    conv_ok[row["Orthogroup"]] += 1

    ogs = sorted(set(div) | set(shift_counts) | set(sel_hit) | set(conv_ok))
    rows = []
    for og in ogs:
        if og == "POOLED":
            continue
        d = div.get(og, {})
        score = 0
        score += min(shift_counts[og], 10)
        if d.get("candidate_beyond_relatedness") == "yes":
            score += 5
        score += min(sel_hit[og], 4)
        score += min(conv_ok[og], 2)
        rows.append(
            {
                "Orthogroup": og,
                "priority_score": score,
                "n_asr_subs_on_shift_branches": shift_counts[og],
                "example_shift_subs": ";".join(shift_examples[og]),
                "divergence_candidate": d.get("candidate_beyond_relatedness", ""),
                "partial_mantel_p": d.get("partial_mantel_p", ""),
                "mean_div_failure": d.get("mean_div_failure", ""),
                "mean_div_success": d.get("mean_div_success", ""),
                "n_selection_hits": sel_hit[og],
                "n_convergence_ok": conv_ok[og],
            }
        )

    rows.sort(key=lambda r: -int(r["priority_score"]))
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "Orthogroup",
        "priority_score",
        "n_asr_subs_on_shift_branches",
        "example_shift_subs",
        "divergence_candidate",
        "partial_mantel_p",
        "mean_div_failure",
        "mean_div_success",
        "n_selection_hits",
        "n_convergence_ok",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"[rank_graft_candidates] ranked={len(rows)} → {args.out_tsv}")


if __name__ == "__main__":
    main()
