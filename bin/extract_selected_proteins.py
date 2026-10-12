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


def id_keys(seq_id: str) -> set[str]:
    """Exact ID plus a form with '|' turned into '_'. OrthoFinder sometimes rewrites pipes."""
    keys = {seq_id}
    if "|" in seq_id:
        keys.add(seq_id.replace("|", "_"))
    return keys


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

    if not wanted:
        raise SystemExit(
            "[extract_selected_proteins] ERROR: the orthogroup table has no member proteins. "
            "No literature gene was placed in an orthogroup, so there is nothing to annotate. "
            "Check the genes-of-interest file (gene name, then genus and species) and the "
            "select_orthogroups log for unmatched names."
        )

    by_id: dict[str, object] = {}
    n_proteome = 0
    for fa in args.proteomes:
        for rec in SeqIO.parse(fa, "fasta"):
            n_proteome += 1
            for key in id_keys(rec.id):
                by_id.setdefault(key, rec)

    records = []
    seen: set[str] = set()
    missing: list[str] = []
    for wid in sorted(wanted):
        rec = None
        for key in id_keys(wid):
            rec = by_id.get(key)
            if rec is not None:
                break
        if rec is None:
            missing.append(wid)
            continue
        if rec.id not in seen:
            records.append(rec)
            seen.add(rec.id)

    if not records:
        sample_wanted = ", ".join(list(sorted(wanted))[:5])
        sample_prot = ", ".join(list(by_id)[:5])
        raise SystemExit(
            f"[extract_selected_proteins] ERROR: 0 of {len(wanted)} orthogroup proteins "
            f"were found among {n_proteome} proteome sequences.\n"
            f"  Requested, for example: {sample_wanted}\n"
            f"  Proteome IDs, for example: {sample_prot}"
        )

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.out_fasta, "fasta")
    print(f"[extract_selected_proteins] wrote {len(records)} / {len(wanted)} requested sequences")
    if missing:
        print(
            f"[extract_selected_proteins] WARNING: {len(missing)} IDs not in the proteomes, "
            f"for example: {', '.join(missing[:8])}"
        )


if __name__ == "__main__":
    main()
