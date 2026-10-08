#!/usr/bin/env python3
"""Map sequence features onto AlphaFold models and flag substitutions in features.

Joins:
  - sequence_features.tsv (1-based seq coords)
  - AF PDB/CIF (residue numbers; pLDDT in B-factors)
  - optional substitutions.tsv from compare_structures.py

Outputs:
  - structure_features.tsv     feature residues with pLDDT
  - substitutions_in_features.tsv  partner diffs overlapping features
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", required=True)
    p.add_argument("--models-dir", required=True)
    p.add_argument("--substitutions", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-plddt", type=float, default=50.0)
    return p.parse_args()


def safe_name(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", seq_id)


def find_structure(models_dir: Path, protein_id: str) -> Path | None:
    frag = safe_name(protein_id)
    candidates = []
    for pat in (f"*{frag}*.pdb", f"*{frag}*.cif"):
        candidates.extend(models_dir.rglob(pat))
    if not candidates:
        return None

    def rank_key(p: Path):
        name = p.name.lower()
        m = re.search(r"rank[_-]?0*(\d+)", name)
        rank = int(m.group(1)) if m else 999
        return (rank, 0 if "relaxed" in name and "unrelaxed" not in name else 1, str(p))

    return sorted(candidates, key=rank_key)[0]


def load_res_plddt(path: Path) -> dict[int, float]:
    parser = PDBParser(QUIET=True) if path.suffix.lower() == ".pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    out = {}
    model = next(structure.get_models())
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=False) or res.id[0] != " ":
                continue
            num = res.id[1]
            if "CA" in res:
                out[num] = float(res["CA"].get_bfactor())
            else:
                atoms = list(res.get_atoms())
                if atoms:
                    out[num] = float(sum(a.get_bfactor() for a in atoms) / len(atoms))
    return out


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    models_dir = Path(args.models_dir)

    features = []
    with open(args.features) as fh:
        features = list(csv.DictReader(fh, delimiter="\t"))

    # cache structures
    plddt_cache: dict[str, dict[int, float]] = {}
    struct_path: dict[str, str] = {}

    feat_rows = []
    for f in features:
        if f.get("feature_type") == "protein_class":
            continue
        pid = f["protein_id"]
        try:
            start, end = int(f["start"]), int(f["end"])
        except (KeyError, ValueError):
            continue
        if pid not in plddt_cache:
            sp = find_structure(models_dir, pid)
            if sp is None:
                plddt_cache[pid] = {}
                struct_path[pid] = ""
            else:
                plddt_cache[pid] = load_res_plddt(sp)
                struct_path[pid] = str(sp)

        vals = []
        for res in range(start, end + 1):
            if res in plddt_cache[pid]:
                vals.append(plddt_cache[pid][res])
        mean_p = sum(vals) / len(vals) if vals else ""
        feat_rows.append(
            {
                **f,
                "structure_file": Path(struct_path[pid]).name if struct_path.get(pid) else "",
                "n_residues_with_coords": len(vals),
                "mean_plddt": f"{mean_p:.1f}" if vals else "",
                "low_confidence": (
                    "yes"
                    if vals and (sum(vals) / len(vals)) < args.min_plddt
                    else ("no" if vals else "missing_structure")
                ),
            }
        )

    with open(out / "structure_features.tsv", "w", newline="") as fh:
        fields = list(feat_rows[0].keys()) if feat_rows else [
            "protein_id",
            "feature_type",
            "start",
            "end",
            "mean_plddt",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(feat_rows)

    # substitutions × features
    sub_out = []
    if args.substitutions and Path(args.substitutions).exists():
        # index features by protein and residue
        by_prot_res: dict[str, list[dict]] = defaultdict(list)
        for f in features:
            if f.get("feature_type") == "protein_class":
                continue
            try:
                a, b = int(f["start"]), int(f["end"])
            except (KeyError, ValueError):
                continue
            for res in range(a, b + 1):
                by_prot_res[f["protein_id"]].append((res, f))

        with open(args.substitutions) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                # compare_structures columns: reference_id, ortholog_id, ref_seq_pos, orth_seq_pos, ...
                ref = row.get("reference_id") or ""
                orth = row.get("ortholog_id") or ""
                try:
                    ref_pos = int(row.get("ref_seq_pos") or 0)
                    orth_pos = int(row.get("orth_seq_pos") or 0)
                except ValueError:
                    continue
                ref_feats = [f for res, f in by_prot_res.get(ref, []) if res == ref_pos]
                orth_feats = [f for res, f in by_prot_res.get(orth, []) if res == orth_pos]
                if not ref_feats and not orth_feats:
                    continue
                labels = sorted(
                    {
                        x.get("feature_type", "")
                        for x in ref_feats + orth_feats
                        if x.get("feature_type")
                    }
                )
                sub_out.append(
                    {
                        **row,
                        "feature_types": ";".join(labels),
                        "ref_feature_names": ";".join(
                            sorted({x.get("feature_name", "") for x in ref_feats if x.get("feature_name")})
                        )[:500],
                        "orth_feature_names": ";".join(
                            sorted({x.get("feature_name", "") for x in orth_feats if x.get("feature_name")})
                        )[:500],
                        "in_functional_feature": "yes",
                    }
                )

    with open(out / "substitutions_in_features.tsv", "w", newline="") as fh:
        if sub_out:
            w = csv.DictWriter(fh, fieldnames=list(sub_out[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(sub_out)
        else:
            fh.write(
                "Orthogroup\treference_id\tortholog_id\tsubstitution\tfeature_types\tin_functional_feature\n"
            )

    print(
        f"[map_features_to_structure] structure_features={len(feat_rows)} "
        f"subs_in_features={len(sub_out)} → {out}"
    )


if __name__ == "__main__":
    main()
