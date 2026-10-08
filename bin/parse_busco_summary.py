#!/usr/bin/env python3
"""Parse a BUSCO short_summary file and apply a completeness gate.

Writes a one-row metrics TSV including pass/fail against --min-complete.
Complete % = Single-copy % + Duplicated % (BUSCO 'C').
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SUMMARY_RE = re.compile(
    r"C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:(\d+)"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", required=True, help="BUSCO short_summary*.txt")
    p.add_argument("--species-id", required=True)
    p.add_argument("--lineage", required=True)
    p.add_argument("--min-complete", type=float, default=90.0)
    p.add_argument("--out-tsv", required=True)
    p.add_argument(
        "--out-status",
        required=True,
        help="Write 'pass' or 'fail' as the only content (for Nextflow branching)",
    )
    return p.parse_args()


def parse_summary(path: Path) -> dict:
    text = path.read_text(errors="replace")
    m = SUMMARY_RE.search(text)
    if not m:
        # Fallback: line-oriented BUSCO v5 tables
        vals = {}
        for line in text.splitlines():
            if "Complete BUSCOs (C)" in line:
                vals["complete_n"] = int(line.strip().split()[0])
            elif "Complete and single-copy BUSCOs (S)" in line:
                vals["single_n"] = int(line.strip().split()[0])
            elif "Complete and duplicated BUSCOs (D)" in line:
                vals["duplicated_n"] = int(line.strip().split()[0])
            elif "Fragmented BUSCOs (F)" in line:
                vals["fragmented_n"] = int(line.strip().split()[0])
            elif "Missing BUSCOs (M)" in line:
                vals["missing_n"] = int(line.strip().split()[0])
            elif "Total BUSCO groups searched" in line:
                vals["n"] = int(line.strip().split()[0])
        if "n" in vals and vals["n"] > 0 and "complete_n" in vals:
            c = 100.0 * vals["complete_n"] / vals["n"]
            s = 100.0 * vals.get("single_n", 0) / vals["n"]
            d = 100.0 * vals.get("duplicated_n", 0) / vals["n"]
            f = 100.0 * vals.get("fragmented_n", 0) / vals["n"]
            miss = 100.0 * vals.get("missing_n", 0) / vals["n"]
            return {
                "complete_pct": c,
                "single_pct": s,
                "duplicated_pct": d,
                "fragmented_pct": f,
                "missing_pct": miss,
                "n_markers": vals["n"],
            }
        raise ValueError(f"Could not parse BUSCO summary: {path}")

    c, s, d, f, miss, n = m.groups()
    return {
        "complete_pct": float(c),
        "single_pct": float(s),
        "duplicated_pct": float(d),
        "fragmented_pct": float(f),
        "missing_pct": float(miss),
        "n_markers": int(n),
    }


def main() -> None:
    args = parse_args()
    summary = Path(args.summary)
    if not summary.exists():
        # BUSCO may nest the short summary; search nearby
        hits = list(summary.parent.rglob("short_summary*.txt"))
        if not hits:
            raise SystemExit(f"BUSCO summary not found: {args.summary}")
        summary = hits[0]

    stats = parse_summary(summary)
    passed = stats["complete_pct"] >= args.min_complete
    status = "pass" if passed else "fail"

    row = {
        "species_id": args.species_id,
        "lineage": args.lineage,
        "complete_pct": f"{stats['complete_pct']:.2f}",
        "single_pct": f"{stats['single_pct']:.2f}",
        "duplicated_pct": f"{stats['duplicated_pct']:.2f}",
        "fragmented_pct": f"{stats['fragmented_pct']:.2f}",
        "missing_pct": f"{stats['missing_pct']:.2f}",
        "n_markers": stats["n_markers"],
        "min_complete_pct": args.min_complete,
        "pass": status,
        "summary_file": str(summary),
    }

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(row)

    Path(args.out_status).write_text(status + "\n")
    print(
        f"[parse_busco] {args.species_id}: complete={stats['complete_pct']:.2f}% "
        f"threshold={args.min_complete}% → {status}"
    )


if __name__ == "__main__":
    main()
