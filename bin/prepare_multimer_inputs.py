#!/usr/bin/env python3
"""Prepare AlphaFold-Multimer inputs for scion ligand × rootstock receptor pairs.

Reads:
  - pair_outcomes.tsv (graft matrix derived pairs)
  - protein classifications (from extract_sequence_features.py)
  - selected OG member FASTA or clean proteomes
  - optional selected_orthogroups.tsv to restrict to GOI families

Writes multimer/ directory:
  - sequences/*.fasta   two-chain FASTA (ligand then receptor)
  - manifest.tsv
  - README.txt
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs", required=True)
    p.add_argument("--classifications", required=True)
    p.add_argument("--proteins", nargs="+", required=True)
    p.add_argument("--selected-ogs", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--success-only",
        action="store_true",
        help="If set, only emit pairs with success=yes (default: all non-self pairs)",
    )
    p.add_argument("--include-failures", action="store_true", default=True)
    return p.parse_args()


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)


def species_of(pid: str) -> str:
    return pid.split("|")[0] if "|" in pid else pid


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    seq_dir = out / "sequences"
    seq_dir.mkdir(parents=True, exist_ok=True)

    seqs = {}
    for fa in args.proteins:
        for rec in SeqIO.parse(fa, "fasta"):
            seqs[rec.id] = rec

    # protein → classes
    class_of = {}
    with open(args.classifications) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            class_of[row["protein_id"]] = set(
                c for c in (row.get("classes") or "").split(";") if c
            )

    # optional restrict to selected OG members
    allowed = set(seqs)
    if args.selected_ogs and Path(args.selected_ogs).exists():
        allowed = set()
        with open(args.selected_ogs) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                for m in (row.get("members") or "").split(","):
                    m = m.strip()
                    if m:
                        allowed.add(m)

    ligands_by_sp = defaultdict(list)
    receptors_by_sp = defaultdict(list)
    for pid, classes in class_of.items():
        if pid not in allowed or pid not in seqs:
            continue
        sp = species_of(pid)
        if "ligand" in classes:
            ligands_by_sp[sp].append(pid)
        if "receptor" in classes:
            receptors_by_sp[sp].append(pid)

    pairs = []
    with open(args.pairs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("is_self") == "yes":
                continue
            if args.success_only and row.get("success") != "yes":
                continue
            pairs.append(row)

    manifest = []
    n = 0
    for row in pairs:
        sc, rs = row["scion"], row["rootstock"]
        # scion ligand × rootstock receptor
        for lig in ligands_by_sp.get(sc, []):
            for recp in receptors_by_sp.get(rs, []):
                n += 1
                job = f"{sc}_{safe_name(lig)}__x__{rs}_{safe_name(recp)}"
                # ColabFold multimer: chains separated by ':' in one file OR two sequences
                # Use two-record FASTA; many ColabFold builds treat multi-record as multimer
                out_fa = seq_dir / f"{job}.fasta"
                SeqIO.write(
                    [
                        SeqRecord(seqs[lig].seq, id=f"{lig}", description="ligand"),
                        SeqRecord(seqs[recp].seq, id=f"{recp}", description="receptor"),
                    ],
                    out_fa,
                    "fasta",
                )
                manifest.append(
                    {
                        "job_id": job,
                        "scion": sc,
                        "rootstock": rs,
                        "outcome": row.get("outcome", ""),
                        "success": row.get("success", ""),
                        "ligand_id": lig,
                        "receptor_id": recp,
                        "orientation": "scion_ligand__rootstock_receptor",
                        "fasta_file": f"sequences/{out_fa.name}",
                    }
                )

    with open(out / "manifest.tsv", "w", newline="") as fh:
        fields = [
            "job_id",
            "scion",
            "rootstock",
            "outcome",
            "success",
            "ligand_id",
            "receptor_id",
            "orientation",
            "fasta_file",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(manifest)

    (out / "README.txt").write_text(
        f"""AlphaFold-Multimer handoff
=========================

Jobs: {len(manifest)}
Each FASTA has two chains: ligand (scion) then receptor (rootstock).

Run with LocalColabFold multimer mode, e.g.:

  colabfold_batch sequences/ out_multimer/ --model-type alphafold2_multimer_v3

Then compare interface PAE / contacts between success=yes vs success=no jobs.
Treat differences as hypotheses for graft signaling mismatch.
"""
    )

    print(
        f"[prepare_multimer_inputs] ligands={sum(len(v) for v in ligands_by_sp.values())} "
        f"receptors={sum(len(v) for v in receptors_by_sp.values())} "
        f"jobs={len(manifest)} → {out}"
    )


if __name__ == "__main__":
    main()
