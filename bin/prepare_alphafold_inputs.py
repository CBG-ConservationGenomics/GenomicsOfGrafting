#!/usr/bin/env python3
"""Prepare an AlphaFold/ColabFold handoff package from gene-family results.

Writes:
  alphafold_inputs/
    manifest.tsv          — what to fold + why
    sequences/*.fasta     — one sequence per model (ColabFold-friendly)
    pairs.tsv             — optional reference↔ortholog pairs for later comparison
    README.txt            — how to run on the GPU server
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selected-ogs", required=True)
    p.add_argument("--identity-summary", required=True)
    p.add_argument("--identity-pairwise", required=True)
    p.add_argument("--proteins-fasta", required=True, help="selected_og_proteins.faa")
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-pct-id", type=float, default=40.0)
    p.add_argument("--max-models-per-og", type=int, default=6)
    p.add_argument(
        "--mode",
        choices=["goi_plus_orthologs", "goi_only", "all_members"],
        default="goi_plus_orthologs",
        help="Which proteins to stage for folding",
    )
    return p.parse_args()


def safe_name(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", seq_id)


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    seq_dir = out / "sequences"
    seq_dir.mkdir(parents=True, exist_ok=True)

    seqs = {r.id: r for r in SeqIO.parse(args.proteins_fasta, "fasta")}

    id_summary = {}
    with open(args.identity_summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            id_summary[row["Orthogroup"]] = row

    # best ortholog per species from pairwise table (highest %ID to reference)
    best_by_og_sp: dict[str, dict[str, tuple[float, str]]] = {}
    with open(args.identity_pairwise) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            og = row["Orthogroup"]
            sp = row.get("target_species_id") or ""
            try:
                pct = float(row["pct_identity"])
            except ValueError:
                continue
            if pct < args.min_pct_id:
                continue
            cur = best_by_og_sp.setdefault(og, {})
            if sp not in cur or pct > cur[sp][0]:
                cur[sp] = (pct, row["target_id"])

    manifest_rows = []
    pair_rows = []

    with open(args.selected_ogs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            og = row["Orthogroup"]
            summ = id_summary.get(og, {})
            ref = summ.get("reference_id") or ""
            matched = [m for m in (row.get("matched_members") or "").split(";") if m]
            members = [m.strip() for m in (row.get("members") or "").split(",") if m.strip()]

            chosen: list[tuple[str, str, str]] = []  # id, role, note

            if args.mode == "all_members":
                for m in members:
                    chosen.append((m, "member", ""))
            elif args.mode == "goi_only":
                for m in matched:
                    chosen.append((m, "goi", ""))
            else:
                # goi_plus_orthologs
                refs = matched[:]
                if ref and ref not in refs:
                    refs.insert(0, ref)
                for m in refs:
                    chosen.append((m, "goi_or_reference", ""))
                for sp, (pct, tid) in sorted(best_by_og_sp.get(og, {}).items()):
                    if tid in {c[0] for c in chosen}:
                        continue
                    chosen.append((tid, "ortholog", f"best_{sp};pct_id={pct:.2f}"))

            # cap models per OG (keep GOI/reference first)
            chosen = chosen[: args.max_models_per_og]

            for seq_id, role, note in chosen:
                rec = seqs.get(seq_id)
                if rec is None:
                    continue
                fname = f"{og}__{safe_name(seq_id)}.fasta"
                SeqIO.write(rec, seq_dir / fname, "fasta")
                manifest_rows.append(
                    {
                        "Orthogroup": og,
                        "genes_of_interest": row.get("genes_of_interest", ""),
                        "query_species_ids": row.get("query_species_ids", ""),
                        "sequence_id": seq_id,
                        "role": role,
                        "note": note,
                        "length": len(rec.seq),
                        "fasta_file": f"sequences/{fname}",
                        "reference_id": ref,
                        "mean_pct_id_to_ref": summ.get("mean_pct_id_to_ref", ""),
                    }
                )
                if role == "ortholog" and ref:
                    pair_rows.append(
                        {
                            "Orthogroup": og,
                            "reference_id": ref,
                            "ortholog_id": seq_id,
                            "ortholog_species_id": seq_id.split("|")[0] if "|" in seq_id else "",
                            "note": note,
                            "reference_fasta": f"sequences/{og}__{safe_name(ref)}.fasta"
                            if (seq_dir / f"{og}__{safe_name(ref)}.fasta").exists()
                            or ref in seqs
                            else "",
                            "ortholog_fasta": f"sequences/{fname}",
                        }
                    )

    with open(out / "manifest.tsv", "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "genes_of_interest",
            "query_species_ids",
            "sequence_id",
            "role",
            "note",
            "length",
            "fasta_file",
            "reference_id",
            "mean_pct_id_to_ref",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    with open(out / "pairs.tsv", "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "reference_id",
            "ortholog_id",
            "ortholog_species_id",
            "note",
            "reference_fasta",
            "ortholog_fasta",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    (out / "README.txt").write_text(
        f"""AlphaFold / ColabFold handoff package
====================================

Prepared by GraftingGeneFamilies (CPU pipeline).
Copy this directory to your GPU server, then run the GPU workflow:

  nextflow run workflows/alphafold.nf -profile gpu \\
    --alphafold_inputs /path/to/this/directory \\
    --outdir alphafold_results

Or run ColabFold directly on sequences/:

  colabfold_batch sequences/ out_colabfold/

Contents
--------
  manifest.tsv     one row per model to fold
  pairs.tsv        reference↔ortholog pairs for structural comparison
  sequences/       individual FASTA files (one protein each)
  README.txt       this file

Filters used
------------
  mode              = {args.mode}
  min_pct_id        = {args.min_pct_id}
  max_models_per_og = {args.max_models_per_og}
  n_models staged   = {len(manifest_rows)}
  n_pairs           = {len(pair_rows)}
"""
    )

    print(
        f"[prepare_alphafold_inputs] staged {len(manifest_rows)} sequences, "
        f"{len(pair_rows)} pairs → {out}"
    )


if __name__ == "__main__":
    main()
