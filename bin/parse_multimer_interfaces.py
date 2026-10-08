#!/usr/bin/env python3
"""Parse AlphaFold/ColabFold multimer models: interface contacts + cross-chain PAE.

For each multimer job directory (or PDB + scores JSON), report:
  - n_interface_contacts (min heavy-atom distance <= cutoff)
  - mean/min interface PAE (from scores JSON pae matrix)
  - mean interface pLDDT
  - iptm/ptm if present

Join to multimer manifest.tsv to compare compatible vs incompatible grafts.

ColabFold scores JSON typically contains:
  { "plddt": [...], "pae": [[...]], "iptm": ..., "ptm": ... }
PAE is indexed by absolute residue index across concatenated chains (0-based).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models-dir", required=True, help="Root dir of ColabFold multimer outputs")
    p.add_argument("--manifest", required=True, help="multimer/manifest.tsv")
    p.add_argument("--outdir", required=True)
    p.add_argument("--distance-cutoff", type=float, default=5.0, help="Å heavy-atom contact cutoff")
    p.add_argument("--pae-cutoff", type=float, default=15.0, help="Max PAE (Å) to count confident contact")
    p.add_argument("--min-plddt", type=float, default=50.0)
    return p.parse_args()


def load_structure(path: Path):
    parser = PDBParser(QUIET=True) if path.suffix.lower() == ".pdb" else MMCIFParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def chain_residues(structure):
    """Return ordered list of (chain_id, resnum, res, ca_plddt) for polymer AAs."""
    model = next(structure.get_models())
    rows = []
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=False) or res.id[0] != " ":
                continue
            plddt = None
            if "CA" in res:
                plddt = float(res["CA"].get_bfactor())
            else:
                atoms = list(res.get_atoms())
                if atoms:
                    plddt = float(sum(a.get_bfactor() for a in atoms) / len(atoms))
            rows.append((chain.id, res.id[1], res, plddt))
    return rows


def min_heavy_distance(res_a, res_b) -> float:
    best = math.inf
    for a in res_a.get_atoms():
        if a.element == "H":
            continue
        for b in res_b.get_atoms():
            if b.element == "H":
                continue
            d = a - b
            if d < best:
                best = d
    return best


def find_pae_matrix(scores: dict):
    if "pae" in scores:
        return scores["pae"]
    if "predicted_aligned_error" in scores:
        return scores["predicted_aligned_error"]
    return None


def load_scores_json(path: Path) -> dict:
    data = json.loads(path.read_text())
    # some AF-DB style wrap in a list
    if isinstance(data, list) and data:
        data = data[0]
    return data


def find_best_model(job_dir: Path) -> tuple[Path | None, Path | None]:
    pdbs = list(job_dir.rglob("*.pdb")) + list(job_dir.rglob("*.cif"))
    if not pdbs:
        # maybe job_dir itself contains files
        pdbs = list(job_dir.glob("*.pdb")) + list(job_dir.glob("*.cif"))
    if not pdbs:
        return None, None

    def rank_key(p: Path):
        m = re.search(r"rank[_-]?0*(\d+)", p.name.lower())
        rank = int(m.group(1)) if m else 999
        return (rank, p.name)

    pdb = sorted(pdbs, key=rank_key)[0]
    # scores json near pdb
    scores = None
    stem = pdb.name
    candidates = list(pdb.parent.glob("*_scores*.json")) + list(pdb.parent.glob("*pae*.json"))
    candidates += list(pdb.parent.glob("*_predicted_aligned_error*.json"))
    # prefer scores with pae matching rank
    rank_m = re.search(r"rank[_-]?0*(\d+)", pdb.name.lower())
    if rank_m:
        r = rank_m.group(1)
        ranked = [c for c in candidates if f"rank_{int(r):03d}" in c.name or f"rank_{r}" in c.name]
        if ranked:
            candidates = ranked + candidates
    if candidates:
        scores = candidates[0]
    return pdb, scores


def analyze_complex(
    pdb_path: Path,
    scores_path: Path | None,
    distance_cutoff: float,
    pae_cutoff: float,
    min_plddt: float,
) -> dict:
    structure = load_structure(pdb_path)
    residues = chain_residues(structure)
    chains = sorted({c for c, _, _, _ in residues})
    if len(chains) < 2:
        return {
            "n_chains": len(chains),
            "n_contacts": 0,
            "n_confident_contacts": 0,
            "status": "need_ge_2_chains",
        }

    # Absolute index in PAE matrix = order of residues as in structure (AF usually A then B)
    # Build list in structure order
    abs_index = {(c, n): i for i, (c, n, _, _) in enumerate(residues)}

    pae = None
    plddt_json = None
    iptm = ptm = ""
    if scores_path and scores_path.exists():
        scores = load_scores_json(scores_path)
        pae = find_pae_matrix(scores)
        plddt_json = scores.get("plddt")
        iptm = scores.get("iptm", scores.get("ipTM", ""))
        ptm = scores.get("ptm", scores.get("pTM", ""))

    # contacts between different chains
    contacts = []
    by_chain = {ch: [] for ch in chains}
    for c, n, res, pld in residues:
        by_chain[c].append((n, res, pld))

    # only first two chains for ligand-receptor
    c1, c2 = chains[0], chains[1]
    pae_vals = []
    plddt_vals = []
    n_conf = 0

    for n1, r1, p1 in by_chain[c1]:
        for n2, r2, p2 in by_chain[c2]:
            d = min_heavy_distance(r1, r2)
            if d > distance_cutoff:
                continue
            i = abs_index[(c1, n1)]
            j = abs_index[(c2, n2)]
            pae_ij = None
            if pae is not None:
                try:
                    # symmetric mean of PAE(i,j) and PAE(j,i)
                    pae_ij = 0.5 * (float(pae[i][j]) + float(pae[j][i]))
                except Exception:
                    pae_ij = None
            pld = None
            if p1 is not None and p2 is not None:
                pld = 0.5 * (p1 + p2)
            elif plddt_json is not None:
                try:
                    pld = 0.5 * (float(plddt_json[i]) + float(plddt_json[j]))
                except Exception:
                    pld = None

            contacts.append((n1, n2, d, pae_ij, pld))
            if pae_ij is not None:
                pae_vals.append(pae_ij)
            if pld is not None:
                plddt_vals.append(pld)
            if (
                pae_ij is not None
                and pae_ij <= pae_cutoff
                and (pld is None or pld >= min_plddt)
            ):
                n_conf += 1

    def mean(xs):
        return sum(xs) / len(xs) if xs else ""

    return {
        "n_chains": len(chains),
        "chain_ids": f"{c1},{c2}",
        "n_contacts": len(contacts),
        "n_confident_contacts": n_conf,
        "mean_contact_distance": f"{mean([c[2] for c in contacts]):.2f}" if contacts else "",
        "mean_interface_pae": f"{mean(pae_vals):.2f}" if pae_vals else "",
        "min_interface_pae": f"{min(pae_vals):.2f}" if pae_vals else "",
        "mean_interface_plddt": f"{mean(plddt_vals):.1f}" if plddt_vals else "",
        "iptm": iptm,
        "ptm": ptm,
        "pdb_file": pdb_path.name,
        "scores_file": scores_path.name if scores_path else "",
        "status": "ok",
    }


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    models = Path(args.models_dir)

    manifest = list(csv.DictReader(open(args.manifest), delimiter="\t"))
    rows = []
    contact_detail = []

    for job in manifest:
        job_id = job["job_id"]
        # find directory matching job_id
        candidates = [models / job_id, models / f"{job_id}_out", models / "models" / job_id]
        candidates += list(models.rglob(f"*{job_id}*"))
        job_dir = None
        for c in candidates:
            if c.is_dir():
                job_dir = c
                break
            if c.is_file() and c.suffix.lower() in {".pdb", ".cif"}:
                job_dir = c.parent
                break
        if job_dir is None:
            rows.append({**job, "status": "model_not_found", "n_contacts": "", "mean_interface_pae": ""})
            continue

        pdb, scores = find_best_model(job_dir)
        if pdb is None:
            rows.append({**job, "status": "pdb_missing", "n_contacts": "", "mean_interface_pae": ""})
            continue

        stats = analyze_complex(
            pdb, scores, args.distance_cutoff, args.pae_cutoff, args.min_plddt
        )
        rows.append({**job, **stats})

    fields = sorted({k for r in rows for k in r.keys()})
    # prefer a readable column order
    preferred = [
        "job_id",
        "scion",
        "rootstock",
        "outcome",
        "success",
        "ligand_id",
        "receptor_id",
        "n_contacts",
        "n_confident_contacts",
        "mean_interface_pae",
        "min_interface_pae",
        "mean_interface_plddt",
        "mean_contact_distance",
        "iptm",
        "ptm",
        "status",
        "pdb_file",
        "scores_file",
    ]
    fields = [c for c in preferred if c in fields] + [c for c in fields if c not in preferred]

    with open(out / "multimer_interface_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # Group comparison: success vs failure
    def collect(flag: str, key: str):
        vals = []
        for r in rows:
            if r.get("success") == flag and r.get("status") == "ok" and r.get(key) not in ("", None):
                try:
                    vals.append(float(r[key]))
                except ValueError:
                    pass
        return vals

    summary_lines = ["metric\tsuccess_mean\tfailure_mean\tn_success\tn_failure"]
    for key in ("mean_interface_pae", "n_confident_contacts", "n_contacts", "mean_interface_plddt"):
        s = collect("yes", key)
        f = collect("no", key)
        sm = f"{sum(s)/len(s):.3f}" if s else ""
        fm = f"{sum(f)/len(f):.3f}" if f else ""
        summary_lines.append(f"{key}\t{sm}\t{fm}\t{len(s)}\t{len(f)}")
    (out / "compatible_vs_incompatible_interface.tsv").write_text(
        "\n".join(summary_lines) + "\n"
    )

    print(
        f"[parse_multimer_interfaces] jobs={len(manifest)} "
        f"ok={sum(1 for r in rows if r.get('status')=='ok')} → {out}"
    )


if __name__ == "__main__":
    main()
