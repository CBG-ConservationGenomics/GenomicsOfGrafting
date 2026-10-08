#!/usr/bin/env python3
"""Summarize ColabFold/AlphaFold outputs into a TSV using the handoff manifest."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True)
    p.add_argument("--model-dirs", nargs="+", required=True)
    p.add_argument("--out-tsv", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    by_stem = {}
    with open(args.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            stem = Path(row["fasta_file"]).stem
            by_stem[stem] = row

    out_rows = []
    for d in args.model_dirs:
        path = Path(d)
        stem = path.name.replace("_out", "")
        meta = by_stem.get(stem, {})
        scores = sorted(path.glob("*_scores_rank_*.json")) + sorted(path.glob("*_scores.json"))
        plddt = pae = ""
        if scores:
            data = json.loads(scores[0].read_text())
            if "plddt" in data and data["plddt"]:
                vals = data["plddt"]
                plddt = f"{(sum(vals) / len(vals)):.2f}"
            if "max_pae" in data:
                pae = str(data["max_pae"])
            elif data.get("pae"):
                flat = [x for row in data["pae"] for x in row]
                pae = f"{(sum(flat) / len(flat)):.2f}" if flat else ""
        pdbs = list(path.glob("*.pdb")) + list(path.glob("*.cif"))
        out_rows.append(
            {
                "Orthogroup": meta.get("Orthogroup", ""),
                "sequence_id": meta.get("sequence_id", stem),
                "role": meta.get("role", ""),
                "mean_plddt": plddt,
                "pae_summary": pae,
                "model_dir": path.name,
                "n_structures": len(pdbs),
            }
        )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "sequence_id",
            "role",
            "mean_plddt",
            "pae_summary",
            "model_dir",
            "n_structures",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"[summarize_alphafold] {len(out_rows)} models → {args.out_tsv}")


if __name__ == "__main__":
    main()
