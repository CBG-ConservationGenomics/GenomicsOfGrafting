#!/usr/bin/env python3
"""Compute per-residue solvent accessibility on AlphaFold models (defense focus).

Order of methods:
  1. Bio.PDB.SASA.ShrakeRupley (Biopython ≥1.80) — preferred
  2. FreeSASA Python API if installed
  3. CA/CB neighbor exposure proxy

Outputs:
  residue_surface.tsv
  substitutions_on_surface.tsv  (optional)
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa

# Empirical max ASA (Å²) for RSA normalization (Tien et al. 2013)
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167,
    "Q": 225, "E": 223, "G": 104, "H": 224, "I": 197,
    "L": 201, "K": 236, "M": 224, "F": 240, "P": 159,
    "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174, "X": 200,
}
AA3 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models-dir", required=True)
    p.add_argument("--proteins", nargs="*", default=None)
    p.add_argument("--classifications", default=None, help="Restrict to defense class if set")
    p.add_argument("--substitutions", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument("--rsa-surface-cutoff", type=float, default=0.25)
    return p.parse_args()


def safe_name(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", seq_id)


def find_structure(models_dir: Path, protein_id: str) -> Path | None:
    frag = safe_name(protein_id)
    hits = list(models_dir.rglob(f"*{frag}*.pdb")) + list(models_dir.rglob(f"*{frag}*.cif"))
    if not hits:
        # try trailing locus token
        locus = protein_id.split("|")[-1]
        hits = list(models_dir.rglob(f"*{locus}*.pdb")) + list(models_dir.rglob(f"*{locus}*.cif"))
    if not hits:
        return None

    def rank_key(p: Path):
        m = re.search(r"rank[_-]?0*(\d+)", p.name.lower())
        return (int(m.group(1)) if m else 999, p.name)

    return sorted(hits, key=rank_key)[0]


def load_structure(path: Path):
    parser = PDBParser(QUIET=True) if path.suffix.lower() == ".pdb" else MMCIFParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def as_pdb_path(path: Path) -> Path:
    if path.suffix.lower() == ".pdb":
        return path
    from Bio.PDB import PDBIO

    struct = load_structure(path)
    tmp = path.with_suffix(".tmp_sasa.pdb")
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(tmp))
    return tmp


def residues_from_struct(struct) -> dict[int, tuple[str, object]]:
    out = {}
    model = next(struct.get_models())
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=False) or res.id[0] != " ":
                continue
            aa = AA3.get(res.get_resname().upper(), "X")
            out[res.id[1]] = (aa, res)
    return out


def rsa_biopython(path: Path) -> dict[int, tuple[str, float]]:
    from Bio.PDB.SASA import ShrakeRupley

    work = as_pdb_path(path)
    try:
        struct = load_structure(work)
        sr = ShrakeRupley()
        sr.compute(struct, level="R")
        out = {}
        for num, (aa, res) in residues_from_struct(struct).items():
            sasa = float(getattr(res, "sasa", 0.0) or 0.0)
            out[num] = (aa, sasa / MAX_ASA.get(aa, 200))
        return out
    finally:
        if work != path and work.exists():
            work.unlink(missing_ok=True)


def rsa_freesasa(path: Path) -> dict[int, tuple[str, float]]:
    import freesasa  # type: ignore

    work = as_pdb_path(path)
    try:
        structure = freesasa.Structure(str(work))
        result = freesasa.calc(structure)
        # Aggregate atom SASA by residue number (chain A assumed / first model)
        res_sasa: dict[int, float] = {}
        res_aa: dict[int, str] = {}
        n = structure.nAtoms()
        for i in range(n):
            resnum = structure.residueNumber(i)
            try:
                resnum = int(str(resnum).strip())
            except ValueError:
                continue
            res_sasa[resnum] = res_sasa.get(resnum, 0.0) + result.atomArea(i)
            name = structure.residueName(i).strip().upper()
            res_aa[resnum] = AA3.get(name, "X")
        return {
            num: (res_aa.get(num, "X"), sasa / MAX_ASA.get(res_aa.get(num, "X"), 200))
            for num, sasa in res_sasa.items()
        }
    finally:
        if work != path and work.exists():
            work.unlink(missing_ok=True)


def rsa_neighbor(path: Path) -> dict[int, tuple[str, float]]:
    struct = load_structure(path)
    residues = []
    for num, (aa, res) in residues_from_struct(struct).items():
        atom = res["CB"] if "CB" in res else (res["CA"] if "CA" in res else None)
        if atom is None:
            continue
        residues.append((num, aa, atom.coord))
    out = {}
    for i, (num, aa, coord) in enumerate(residues):
        nneigh = 0
        for j, (_n2, _a2, c2) in enumerate(residues):
            if i == j:
                continue
            dx, dy, dz = coord[0] - c2[0], coord[1] - c2[1], coord[2] - c2[2]
            if dx * dx + dy * dy + dz * dz <= 144.0:
                nneigh += 1
        out[num] = (aa, max(0.0, 1.0 - (nneigh / 30.0)))
    return out


def compute_surface(path: Path) -> tuple[dict[int, tuple[str, float]], str]:
    try:
        return rsa_biopython(path), "biopython_sasa"
    except Exception:
        pass
    try:
        return rsa_freesasa(path), "freesasa"
    except Exception:
        pass
    return rsa_neighbor(path), "neighbor_fallback"


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    models = Path(args.models_dir)

    protein_ids: set[str] = set(args.proteins or [])
    if args.classifications and Path(args.classifications).exists():
        with open(args.classifications) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        defense = {
            row["protein_id"]
            for row in rows
            if "defense" in (row.get("classes") or "").split(";")
        }
        if defense:
            protein_ids = defense if not protein_ids else (protein_ids & defense)
        elif not protein_ids:
            protein_ids = {row["protein_id"] for row in rows}
    if not protein_ids:
        raise SystemExit("No protein IDs — pass --classifications and/or --proteins")

    rows = []
    method = ""
    for pid in sorted(protein_ids):
        sp = find_structure(models, pid)
        if sp is None:
            continue
        surf, method = compute_surface(sp)
        for resnum, (aa, score) in sorted(surf.items()):
            rows.append(
                {
                    "protein_id": pid,
                    "resnum": resnum,
                    "aa": aa,
                    "score_type": "rsa_or_exposure",
                    "score": f"{score:.3f}",
                    "is_surface": "yes" if score >= args.rsa_surface_cutoff else "no",
                    "method": method,
                    "structure_file": sp.name,
                }
            )

    with open(out / "residue_surface.tsv", "w", newline="") as fh:
        fields = [
            "protein_id",
            "resnum",
            "aa",
            "score_type",
            "score",
            "is_surface",
            "method",
            "structure_file",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    sub_out = []
    if args.substitutions and Path(args.substitutions).exists():
        surface_idx = {
            (r["protein_id"], int(r["resnum"])): r for r in rows if r["is_surface"] == "yes"
        }
        with open(args.substitutions) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ref = row.get("reference_id") or ""
                orth = row.get("ortholog_id") or ""
                try:
                    rp = int(row.get("ref_seq_pos") or 0)
                    op = int(row.get("orth_seq_pos") or 0)
                except ValueError:
                    continue
                ref_s = surface_idx.get((ref, rp))
                orth_s = surface_idx.get((orth, op))
                if not ref_s and not orth_s:
                    continue
                sub_out.append(
                    {
                        **row,
                        "ref_is_surface": "yes" if ref_s else "no",
                        "orth_is_surface": "yes" if orth_s else "no",
                        "ref_surface_score": ref_s["score"] if ref_s else "",
                        "orth_surface_score": orth_s["score"] if orth_s else "",
                    }
                )

    with open(out / "substitutions_on_surface.tsv", "w", newline="") as fh:
        if sub_out:
            w = csv.DictWriter(fh, fieldnames=list(sub_out[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(sub_out)
        else:
            fh.write("reference_id\tortholog_id\tsubstitution\tref_is_surface\n")

    print(
        f"[compute_surface_rsa] residues={len(rows)} surface_subs={len(sub_out)} "
        f"method={method or 'none'} → {out}"
    )


if __name__ == "__main__":
    main()
