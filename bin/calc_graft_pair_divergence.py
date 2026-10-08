#!/usr/bin/env python3
"""Per orthogroup, compute scion–rootstock protein divergence for each graft pair.

Uses OrthoFinder protein MSAs. Picks one sequence per species (longest if multiple).
Reports pairwise % identity and Hamming distance over ungapped columns.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import AlignIO


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs", required=True, help="pair_outcomes.tsv")
    p.add_argument("--og-list", required=True, help="selected orthogroup IDs")
    p.add_argument("--alignments-dir", required=True)
    p.add_argument("--out-tsv", required=True)
    return p.parse_args()


def find_alignment(align_dir: Path, og: str) -> Path | None:
    for name in (f"{og}.fa", f"{og}.faa", f"{og}.fasta", f"{og}.aln"):
        c = align_dir / name
        if c.exists():
            return c
    hits = list(align_dir.rglob(f"{og}.*"))
    for h in hits:
        if h.suffix.lower() in {".fa", ".faa", ".fasta", ".aln"}:
            return h
    return None


def species_of(seq_id: str) -> str:
    return seq_id.split("|")[0] if "|" in seq_id else seq_id.split(".")[0]


def pick_per_species(alignment) -> dict[str, str]:
    best: dict[str, tuple[int, str]] = {}
    for rec in alignment:
        sp = species_of(rec.id)
        seq = str(rec.seq).upper()
        ungapped = len(seq.replace("-", ""))
        if sp not in best or ungapped > best[sp][0]:
            best[sp] = (ungapped, seq)
    return {sp: seq for sp, (_n, seq) in best.items()}


def pair_stats(seq_a: str, seq_b: str) -> tuple[float | None, int, int]:
    assert len(seq_a) == len(seq_b)
    match = compared = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a == b:
            match += 1
    if compared == 0:
        return None, 0, 0
    return 100.0 * match / compared, compared - match, compared


def main() -> None:
    args = parse_args()
    align_dir = Path(args.alignments_dir)

    pairs = []
    with open(args.pairs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("is_self") == "yes":
                continue
            pairs.append(row)

    with open(args.og_list) as fh:
        ogs = [ln.strip() for ln in fh if ln.strip()]

    rows = []
    for og in ogs:
        aln_path = find_alignment(align_dir, og)
        if not aln_path:
            continue
        aln = AlignIO.read(str(aln_path), "fasta")
        by_sp = pick_per_species(aln)
        for p in pairs:
            sc, rs = p["scion"], p["rootstock"]
            if sc not in by_sp or rs not in by_sp:
                rows.append(
                    {
                        "Orthogroup": og,
                        "scion": sc,
                        "rootstock": rs,
                        "outcome": p["outcome"],
                        "success": p["success"],
                        "pct_identity": "",
                        "n_diffs": "",
                        "n_compared": "",
                        "status": "missing_species_in_og",
                    }
                )
                continue
            pct, ndiff, ncmp = pair_stats(by_sp[sc], by_sp[rs])
            rows.append(
                {
                    "Orthogroup": og,
                    "scion": sc,
                    "rootstock": rs,
                    "outcome": p["outcome"],
                    "success": p["success"],
                    "pct_identity": f"{pct:.2f}" if pct is not None else "",
                    "n_diffs": ndiff,
                    "n_compared": ncmp,
                    "status": "ok",
                }
            )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "Orthogroup",
        "scion",
        "rootstock",
        "outcome",
        "success",
        "pct_identity",
        "n_diffs",
        "n_compared",
        "status",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"[calc_graft_pair_divergence] rows={len(rows)} → {args.out_tsv}")


if __name__ == "__main__":
    main()
