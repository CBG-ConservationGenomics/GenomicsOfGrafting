#!/usr/bin/env python3
"""Prepare an AlphaFold/ColabFold handoff package from gene-family results.

Optionally ranks / filters orthogroups using graft-phenotype priority
(graft_candidate_ranking.tsv from ASR, divergence tests, HyPhy, convergence)
so GPU folding focuses on the most informative families first.

Writes:
  alphafold_inputs/
    manifest.tsv          — what to fold + why (sorted by priority)
    pairs.tsv             — reference↔ortholog pairs for later comparison
    sequences/*.fasta     — one sequence per model (ColabFold-friendly)
    og_priority.tsv       — OG-level ranking used for this handoff
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
    p.add_argument(
        "--priority-ranking",
        default=None,
        help="graft_candidate_ranking.tsv (from RANK_GRAFT_CANDIDATES)",
    )
    p.add_argument(
        "--min-priority-score",
        type=int,
        default=0,
        help="Drop OGs with priority_score below this when ranking is provided",
    )
    p.add_argument(
        "--max-ogs",
        type=int,
        default=0,
        help="After ranking, keep only the top N OGs (0 = no cap)",
    )
    return p.parse_args()


def safe_name(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", seq_id)


def load_priority(path: str | None) -> dict[str, dict]:
    """Return Orthogroup → ranking row. Empty/placeholder files → {}."""
    if not path or not Path(path).exists():
        return {}
    text = Path(path).read_text()
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines or "Orthogroup" not in lines[0]:
        return {}
    out = {}
    for row in csv.DictReader(lines, delimiter="\t"):
        og = row.get("Orthogroup")
        if not og or og == "POOLED":
            continue
        try:
            row["_score"] = int(float(row.get("priority_score") or 0))
        except ValueError:
            row["_score"] = 0
        out[og] = row
    return out


def empty_pri(score: int = 0) -> dict:
    return {
        "priority_score": str(score) if score else "0",
        "n_asr_subs_on_shift_branches": "",
        "divergence_candidate": "",
        "n_selection_hits": "",
        "example_shift_subs": "",
        "_score": score,
    }


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

    priority = load_priority(args.priority_ranking)
    used_priority = bool(priority)

    selected_rows = []
    with open(args.selected_ogs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            selected_rows.append(row)

    scored: list[tuple[int, dict, dict]] = []
    for row in selected_rows:
        og = row["Orthogroup"]
        if used_priority:
            pri = priority.get(og) or empty_pri(0)
            score = int(pri.get("_score", 0))
            if score < args.min_priority_score:
                continue
        else:
            pri = empty_pri(0)
            pri["priority_score"] = ""
            score = 0
        scored.append((score, row, pri))

    scored.sort(key=lambda x: (-x[0], x[1]["Orthogroup"]))
    if used_priority and args.max_ogs > 0:
        scored = scored[: args.max_ogs]

    og_pri_rows = []
    for rank_i, (score, row, pri) in enumerate(scored, 1):
        og_pri_rows.append(
            {
                "fold_rank": rank_i,
                "Orthogroup": row["Orthogroup"],
                "priority_score": pri.get("priority_score", score),
                "n_asr_subs_on_shift_branches": pri.get("n_asr_subs_on_shift_branches", ""),
                "divergence_candidate": pri.get("divergence_candidate", ""),
                "n_selection_hits": pri.get("n_selection_hits", ""),
                "example_shift_subs": pri.get("example_shift_subs", ""),
                "genes_of_interest": row.get("genes_of_interest", ""),
            }
        )

    manifest_rows = []
    pair_rows = []

    for fold_rank, (score, row, pri) in enumerate(scored, 1):
        og = row["Orthogroup"]
        summ = id_summary.get(og, {})
        ref = summ.get("reference_id") or ""
        matched = [m for m in (row.get("matched_members") or "").split(";") if m]
        members = [m.strip() for m in (row.get("members") or "").split(",") if m.strip()]

        chosen: list[tuple[str, str, str]] = []

        if args.mode == "all_members":
            for m in members:
                chosen.append((m, "member", ""))
        elif args.mode == "goi_only":
            for m in matched:
                chosen.append((m, "goi", ""))
        else:
            refs = matched[:]
            if ref and ref not in refs:
                refs.insert(0, ref)
            for m in refs:
                chosen.append((m, "goi_or_reference", ""))
            orths = []
            for sp, (pct, tid) in best_by_og_sp.get(og, {}).items():
                if tid in {c[0] for c in chosen}:
                    continue
                orths.append((sp, pct, tid))
            # Divergence candidates: fold more divergent partners first (structure diffs)
            if pri.get("divergence_candidate") == "yes":
                orths.sort(key=lambda x: x[1])
            else:
                orths.sort(key=lambda x: -x[1])
            for sp, pct, tid in orths:
                chosen.append((tid, "ortholog", f"best_{sp};pct_id={pct:.2f}"))

        chosen = chosen[: args.max_models_per_og]

        priority_note = ""
        if used_priority:
            bits = [f"priority={score}", f"fold_rank={fold_rank}"]
            if pri.get("n_asr_subs_on_shift_branches"):
                bits.append(f"asr_shift={pri['n_asr_subs_on_shift_branches']}")
            if pri.get("divergence_candidate") == "yes":
                bits.append("div_candidate")
            if pri.get("n_selection_hits"):
                bits.append(f"sel={pri['n_selection_hits']}")
            if pri.get("example_shift_subs"):
                bits.append(f"eg={pri['example_shift_subs']}")
            priority_note = ";".join(bits)

        for seq_id, role, note in chosen:
            rec = seqs.get(seq_id)
            if rec is None:
                continue
            fname = f"{og}__{safe_name(seq_id)}.fasta"
            SeqIO.write(rec, seq_dir / fname, "fasta")
            combined_note = ";".join(x for x in (note, priority_note) if x)
            manifest_rows.append(
                {
                    "fold_rank": fold_rank,
                    "priority_score": pri.get("priority_score", "") if used_priority else "",
                    "Orthogroup": og,
                    "genes_of_interest": row.get("genes_of_interest", ""),
                    "query_species_ids": row.get("query_species_ids", ""),
                    "sequence_id": seq_id,
                    "role": role,
                    "note": combined_note,
                    "length": len(rec.seq),
                    "fasta_file": f"sequences/{fname}",
                    "reference_id": ref,
                    "mean_pct_id_to_ref": summ.get("mean_pct_id_to_ref", ""),
                    "n_asr_subs_on_shift_branches": pri.get("n_asr_subs_on_shift_branches", ""),
                    "divergence_candidate": pri.get("divergence_candidate", ""),
                    "n_selection_hits": pri.get("n_selection_hits", ""),
                }
            )
            if role == "ortholog" and ref:
                pair_rows.append(
                    {
                        "Orthogroup": og,
                        "fold_rank": fold_rank,
                        "priority_score": pri.get("priority_score", "") if used_priority else "",
                        "reference_id": ref,
                        "ortholog_id": seq_id,
                        "ortholog_species_id": seq_id.split("|")[0] if "|" in seq_id else "",
                        "note": combined_note,
                        "reference_fasta": f"sequences/{og}__{safe_name(ref)}.fasta",
                        "ortholog_fasta": f"sequences/{fname}",
                    }
                )

    with open(out / "manifest.tsv", "w", newline="") as fh:
        fields = [
            "fold_rank",
            "priority_score",
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
            "n_asr_subs_on_shift_branches",
            "divergence_candidate",
            "n_selection_hits",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    with open(out / "pairs.tsv", "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "fold_rank",
            "priority_score",
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

    with open(out / "og_priority.tsv", "w", newline="") as fh:
        fields = [
            "fold_rank",
            "Orthogroup",
            "priority_score",
            "n_asr_subs_on_shift_branches",
            "divergence_candidate",
            "n_selection_hits",
            "example_shift_subs",
            "genes_of_interest",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(og_pri_rows)

    (out / "README.txt").write_text(
        f"""AlphaFold / ColabFold handoff package
====================================

Prepared by GraftingGeneFamilies (CPU pipeline) AFTER graft-phenotype
ranking when available (ASR shift-branch hits, divergence-beyond-relatedness,
HyPhy selection, convergence).

Copy this directory to your GPU server, then run:

  nextflow run workflows/alphafold.nf -profile gpu \\
    --alphafold_inputs /path/to/this/directory \\
    --outdir alphafold_results

Fold order follows fold_rank in manifest.tsv / og_priority.tsv
(highest phenotype priority first).

Contents
--------
  manifest.tsv     one row per model (sorted by fold_rank)
  pairs.tsv        reference↔ortholog pairs for structural comparison
  og_priority.tsv  OG-level scores used to prioritize this handoff
  sequences/       individual FASTA files (one protein each)
  README.txt       this file

Filters used
------------
  mode                 = {args.mode}
  min_pct_id           = {args.min_pct_id}
  max_models_per_og    = {args.max_models_per_og}
  priority_ranking     = {"yes" if used_priority else "no (all selected OGs)"}
  min_priority_score   = {args.min_priority_score}
  max_ogs              = {args.max_ogs or "none"}
  n_ogs staged         = {len(og_pri_rows)}
  n_models staged      = {len(manifest_rows)}
  n_pairs              = {len(pair_rows)}
"""
    )

    print(
        f"[prepare_alphafold_inputs] staged {len(manifest_rows)} sequences "
        f"({len(og_pri_rows)} OGs), {len(pair_rows)} pairs; "
        f"priority={'on' if used_priority else 'off'} → {out}"
    )


if __name__ == "__main__":
    main()
