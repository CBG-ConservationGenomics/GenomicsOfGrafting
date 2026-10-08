#!/usr/bin/env python3
"""Prepare a one-sequence-per-species FASTA for IQ-TREE / convergence tools."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--alignment", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-species-map", required=True)
    return p.parse_args()


def species_of(seq_id: str) -> str:
    return seq_id.split("|")[0] if "|" in seq_id else seq_id


def main() -> None:
    args = parse_args()
    aln = AlignIO.read(args.alignment, "fasta")
    best: dict[str, SeqRecord] = {}
    for rec in aln:
        sp = species_of(rec.id)
        ungapped = len(str(rec.seq).replace("-", ""))
        if sp not in best or ungapped > len(str(best[sp].seq).replace("-", "")):
            best[sp] = SeqRecord(Seq(str(rec.seq)), id=sp, description=rec.id)

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([best[sp] for sp in sorted(best)], args.out_fasta, "fasta")
    with open(args.out_species_map, "w") as fh:
        fh.write("species_id\toriginal_id\n")
        for sp in sorted(best):
            fh.write(f"{sp}\t{best[sp].description}\n")
    print(f"[prep_og_alignment] species={len(best)} → {args.out_fasta}")


if __name__ == "__main__":
    main()
