#!/usr/bin/env python3
"""Extract protein sequences for members of selected orthogroups."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selected-ogs", required=True)
    p.add_argument("--proteomes", nargs="+", required=True)
    p.add_argument("--out-fasta", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    wanted: set[str] = set()
    with open(args.selected_ogs) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for m in (row.get("members") or "").split(","):
                m = m.strip()
                if m:
                    wanted.add(m)

    records = []
    for fa in args.proteomes:
        for rec in SeqIO.parse(fa, "fasta"):
            if rec.id in wanted:
                records.append(rec)

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.out_fasta, "fasta")
    print(f"[extract_selected_proteins] wrote {len(records)} / {len(wanted)} requested sequences")


if __name__ == "__main__":
    main()
