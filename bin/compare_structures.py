#!/usr/bin/env python3
"""Compare AlphaFold models for reference↔ortholog pairs (CPU stage).

Inputs
------
- alphafold_inputs/pairs.tsv (+ manifest.tsv optional)
- alphafold model directory tree (ColabFold/AF outputs)
- OrthoFinder MultipleSequenceAlignments/ (for substitution mapping)

Outputs
-------
- structure_comparison_summary.tsv  per pair: RMSD, n_subs, optional TM-score
- substitutions.tsv                 per AA difference with structure residue #s + pLDDT
- superimposed/                     ortholog PDBs fitted onto reference (optional)
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path

from Bio import AlignIO
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Superimposer
from Bio.PDB.Polypeptide import is_aa


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs", required=True, help="alphafold_inputs/pairs.tsv")
    p.add_argument("--models-dir", required=True, help="Directory with ColabFold/AF model outputs")
    p.add_argument(
        "--alignments-dir",
        required=True,
        help="OrthoFinder MultipleSequenceAlignments directory",
    )
    p.add_argument("--outdir", required=True)
    p.add_argument("--manifest", default=None, help="Optional alphafold_inputs/manifest.tsv")
    p.add_argument("--min-plddt", type=float, default=50.0, help="Flag subs below this pLDDT")
    p.add_argument("--write-superimposed", action="store_true", default=True)
    p.add_argument("--no-write-superimposed", action="store_false", dest="write_superimposed")
    p.add_argument(
        "--usalign-bin",
        default="USalign",
        help="Optional USalign/TMalign binary on PATH for TM-score (empty to skip)",
    )
    return p.parse_args()


def safe_name(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", seq_id)


def find_alignment(align_dir: Path, og: str) -> Path | None:
    for name in (f"{og}.fa", f"{og}.faa", f"{og}.fasta", f"{og}.aln"):
        p = align_dir / name
        if p.exists():
            return p
    hits = list(align_dir.rglob(f"{og}.*"))
    for h in hits:
        if h.suffix.lower() in {".fa", ".faa", ".fasta", ".aln"}:
            return h
    return None


def find_structure(models_dir: Path, og: str, seq_id: str) -> Path | None:
    """Locate best-ranked PDB/CIF for a staged sequence."""
    stem = f"{og}__{safe_name(seq_id)}"
    candidates: list[Path] = []

    # Preferred layout from workflows/alphafold.nf
    for base in [
        models_dir / f"{stem}_out",
        models_dir / "models" / f"{stem}_out",
        models_dir / stem,
        models_dir / "models" / stem,
    ]:
        if base.is_dir():
            candidates.extend(base.rglob("*.pdb"))
            candidates.extend(base.rglob("*.cif"))

    # Fallback: search by safe name fragment
    if not candidates:
        frag = safe_name(seq_id)
        candidates.extend(models_dir.rglob(f"*{frag}*.pdb"))
        candidates.extend(models_dir.rglob(f"*{frag}*.cif"))

    if not candidates:
        return None

    def rank_key(p: Path) -> tuple:
        name = p.name.lower()
        # prefer rank_001 / rank1, then unrelaxed/relaxed
        m = re.search(r"rank[_-]?0*(\d+)", name)
        rank = int(m.group(1)) if m else 999
        relaxed = 0 if "relaxed" in name and "unrelaxed" not in name else 1
        return (rank, relaxed, len(name), str(p))

    return sorted(candidates, key=rank_key)[0]


def load_structure(path: Path):
    parser = PDBParser(QUIET=True) if path.suffix.lower() == ".pdb" else MMCIFParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def protein_residues(structure):
    """Ordered polymer residues (AA only) from first model, all chains concatenated."""
    model = next(structure.get_models())
    residues = []
    for chain in model:
        for res in chain:
            if is_aa(res, standard=False) and res.id[0] == " ":
                residues.append(res)
    return residues


def res_plddt(res) -> float | None:
    """AlphaFold stores pLDDT in the B-factor column of CA (or atoms)."""
    if "CA" in res:
        return float(res["CA"].get_bfactor())
    atoms = list(res.get_atoms())
    if not atoms:
        return None
    return float(sum(a.get_bfactor() for a in atoms) / len(atoms))


def msa_sequences(aln_path: Path) -> dict[str, str]:
    aln = AlignIO.read(str(aln_path), "fasta")
    return {r.id: str(r.seq).upper() for r in aln}


def ungapped_index_map(aligned: str) -> dict[int, int]:
    """Map 0-based ungapped sequence index → 0-based alignment column."""
    out = {}
    seq_i = 0
    for col_i, aa in enumerate(aligned):
        if aa != "-":
            out[seq_i] = col_i
            seq_i += 1
    return out


def column_to_seq_index(aligned: str) -> dict[int, int]:
    """Map alignment column → 0-based ungapped index (gaps omitted)."""
    out = {}
    seq_i = 0
    for col_i, aa in enumerate(aligned):
        if aa != "-":
            out[col_i] = seq_i
            seq_i += 1
    return out


def compare_pair(
    og: str,
    ref_id: str,
    orth_id: str,
    ref_struct_path: Path,
    orth_struct_path: Path,
    aln_path: Path,
    min_plddt: float,
) -> tuple[dict, list[dict], object | None]:
    """Return summary row, substitution rows, and superimposed orth structure (or None)."""
    seqs = msa_sequences(aln_path)
    if ref_id not in seqs or orth_id not in seqs:
        # try matching by trailing original id
        def resolve(wanted: str) -> str | None:
            if wanted in seqs:
                return wanted
            tail = wanted.split("|")[-1]
            for k in seqs:
                if k == wanted or k.endswith(tail) or k.split("|")[-1] == tail:
                    return k
            return None

        rkey, okey = resolve(ref_id), resolve(orth_id)
        if not rkey or not okey:
            raise ValueError(f"MSA missing {ref_id} and/or {orth_id} in {aln_path}")
        ref_id, orth_id = rkey, okey

    ref_aln = seqs[ref_id]
    orth_aln = seqs[orth_id]
    if len(ref_aln) != len(orth_aln):
        raise ValueError(f"MSA length mismatch for {og}")

    ref_struct = load_structure(ref_struct_path)
    orth_struct = load_structure(orth_struct_path)
    ref_res = protein_residues(ref_struct)
    orth_res = protein_residues(orth_struct)

    ref_col = column_to_seq_index(ref_aln)
    orth_col = column_to_seq_index(orth_aln)

    # Build paired CA atoms for columns with AA in both sequences
    ref_atoms = []
    orth_atoms = []
    subs = []
    n_aligned = 0
    n_ident = 0

    for col in range(len(ref_aln)):
        a = ref_aln[col]
        b = orth_aln[col]
        if a == "-" or b == "-":
            continue
        if col not in ref_col or col not in orth_col:
            continue
        ri = ref_col[col]
        oi = orth_col[col]
        if ri >= len(ref_res) or oi >= len(orth_res):
            # Structure/sequence length mismatch (truncated model etc.)
            continue
        r_res = ref_res[ri]
        o_res = orth_res[oi]
        n_aligned += 1
        if a == b:
            n_ident += 1
        else:
            r_pld = res_plddt(r_res)
            o_pld = res_plddt(o_res)
            low = False
            if r_pld is not None and r_pld < min_plddt:
                low = True
            if o_pld is not None and o_pld < min_plddt:
                low = True
            subs.append(
                {
                    "Orthogroup": og,
                    "reference_id": ref_id,
                    "ortholog_id": orth_id,
                    "msa_column": col + 1,
                    "ref_seq_pos": ri + 1,
                    "orth_seq_pos": oi + 1,
                    "ref_resnum": r_res.id[1],
                    "orth_resnum": o_res.id[1],
                    "ref_aa": a,
                    "orth_aa": b,
                    "substitution": f"{a}{ri + 1}{b}",
                    "ref_plddt": f"{r_pld:.1f}" if r_pld is not None else "",
                    "orth_plddt": f"{o_pld:.1f}" if o_pld is not None else "",
                    "low_confidence": "yes" if low else "no",
                }
            )

        if "CA" in r_res and "CA" in o_res:
            ref_atoms.append(r_res["CA"])
            orth_atoms.append(o_res["CA"])

    rmsd = ""
    n_atoms = len(ref_atoms)
    superimposed = None
    if n_atoms >= 3:
        sup = Superimposer()
        sup.set_atoms(ref_atoms, orth_atoms)
        rmsd = f"{sup.rms:.3f}"
        # Apply the same rotation/translation to all ortholog atoms for writing
        orth_copy = load_structure(orth_struct_path)
        rot, tran = sup.rotran
        for atom in orth_copy.get_atoms():
            atom.transform(rot, tran)
        superimposed = orth_copy

    pct_id = f"{(100.0 * n_ident / n_aligned):.2f}" if n_aligned else ""

    summary = {
        "Orthogroup": og,
        "reference_id": ref_id,
        "ortholog_id": orth_id,
        "reference_structure": str(ref_struct_path),
        "ortholog_structure": str(orth_struct_path),
        "alignment_file": aln_path.name,
        "n_aligned_positions": n_aligned,
        "n_identical": n_ident,
        "msa_pct_identity": pct_id,
        "n_substitutions": len(subs),
        "n_low_confidence_subs": sum(1 for s in subs if s["low_confidence"] == "yes"),
        "ca_rmsd": rmsd,
        "n_ca_atoms_superposed": n_atoms,
        "tm_score_ref": "",
        "tm_score_orth": "",
        "usalign_status": "not_run",
        "status": "ok",
    }
    return summary, subs, superimposed


def run_usalign(usalign_bin: str, ref_pdb: Path, orth_pdb: Path) -> tuple[str, str, str]:
    """Return (tm_score_ref, tm_score_orth, status)."""
    if not usalign_bin:
        return "", "", "skipped"
    exe = shutil.which(usalign_bin) or shutil.which("TMalign")
    if not exe:
        return "", "", "binary_not_found"
    try:
        proc = subprocess.run(
            [exe, str(orth_pdb), str(ref_pdb)],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        out = proc.stdout + "\n" + proc.stderr
        # USalign/TMalign: "TM-score= 0.xxxx (if normalized by length of Chain_1)"
        scores = re.findall(r"TM-score\s*=\s*([0-9.]+)", out)
        if len(scores) >= 2:
            return scores[0], scores[1], "ok"
        if len(scores) == 1:
            return scores[0], "", "ok_partial"
        return "", "", "parse_failed"
    except Exception as exc:  # noqa: BLE001
        return "", "", f"error:{exc}"


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    super_dir = out / "superimposed"
    if args.write_superimposed:
        super_dir.mkdir(parents=True, exist_ok=True)

    models_dir = Path(args.models_dir)
    align_dir = Path(args.alignments_dir)

    summaries = []
    all_subs = []

    with open(args.pairs) as fh:
        pairs = list(csv.DictReader(fh, delimiter="\t"))

    if not pairs:
        print("[compare_structures] WARNING: no pairs in pairs.tsv")

    io = PDBIO()

    for pair in pairs:
        og = pair["Orthogroup"]
        ref_id = pair["reference_id"]
        orth_id = pair["ortholog_id"]

        aln_path = find_alignment(align_dir, og)
        ref_pdb = find_structure(models_dir, og, ref_id)
        orth_pdb = find_structure(models_dir, og, orth_id)

        if not aln_path or not ref_pdb or not orth_pdb:
            summaries.append(
                {
                    "Orthogroup": og,
                    "reference_id": ref_id,
                    "ortholog_id": orth_id,
                    "reference_structure": str(ref_pdb or ""),
                    "ortholog_structure": str(orth_pdb or ""),
                    "alignment_file": str(aln_path.name if aln_path else ""),
                    "n_aligned_positions": "",
                    "n_identical": "",
                    "msa_pct_identity": "",
                    "n_substitutions": "",
                    "n_low_confidence_subs": "",
                    "ca_rmsd": "",
                    "n_ca_atoms_superposed": "",
                    "tm_score_ref": "",
                    "tm_score_orth": "",
                    "usalign_status": "",
                    "status": "missing_"
                    + "+".join(
                        [
                            x
                            for x, ok in [
                                ("alignment", aln_path),
                                ("ref_structure", ref_pdb),
                                ("orth_structure", orth_pdb),
                            ]
                            if not ok
                        ]
                    ),
                }
            )
            continue

        try:
            summary, subs, superimposed = compare_pair(
                og, ref_id, orth_id, ref_pdb, orth_pdb, aln_path, args.min_plddt
            )
        except Exception as exc:  # noqa: BLE001
            summaries.append(
                {
                    "Orthogroup": og,
                    "reference_id": ref_id,
                    "ortholog_id": orth_id,
                    "reference_structure": str(ref_pdb),
                    "ortholog_structure": str(orth_pdb),
                    "alignment_file": aln_path.name,
                    "n_aligned_positions": "",
                    "n_identical": "",
                    "msa_pct_identity": "",
                    "n_substitutions": "",
                    "n_low_confidence_subs": "",
                    "ca_rmsd": "",
                    "n_ca_atoms_superposed": "",
                    "tm_score_ref": "",
                    "tm_score_orth": "",
                    "usalign_status": "",
                    "status": f"error:{exc}",
                }
            )
            continue

        tm1, tm2, tm_status = run_usalign(args.usalign_bin, ref_pdb, orth_pdb)
        summary["tm_score_ref"] = tm1
        summary["tm_score_orth"] = tm2
        summary["usalign_status"] = tm_status

        if args.write_superimposed and superimposed is not None:
            out_pdb = super_dir / f"{og}__{safe_name(orth_id)}__on__{safe_name(ref_id)}.pdb"
            io.set_structure(superimposed)
            io.save(str(out_pdb))
            summary["superimposed_pdb"] = str(out_pdb.name)
        else:
            summary["superimposed_pdb"] = ""

        summaries.append(summary)
        all_subs.extend(subs)

    sum_fields = [
        "Orthogroup",
        "reference_id",
        "ortholog_id",
        "reference_structure",
        "ortholog_structure",
        "alignment_file",
        "n_aligned_positions",
        "n_identical",
        "msa_pct_identity",
        "n_substitutions",
        "n_low_confidence_subs",
        "ca_rmsd",
        "n_ca_atoms_superposed",
        "tm_score_ref",
        "tm_score_orth",
        "usalign_status",
        "superimposed_pdb",
        "status",
    ]
    with open(out / "structure_comparison_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=sum_fields, delimiter="\t", extrasaction="ignore", restval=""
        )
        w.writeheader()
        w.writerows(summaries)

    sub_fields = [
        "Orthogroup",
        "reference_id",
        "ortholog_id",
        "msa_column",
        "ref_seq_pos",
        "orth_seq_pos",
        "ref_resnum",
        "orth_resnum",
        "ref_aa",
        "orth_aa",
        "substitution",
        "ref_plddt",
        "orth_plddt",
        "low_confidence",
    ]
    with open(out / "substitutions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sub_fields, delimiter="\t")
        w.writeheader()
        w.writerows(all_subs)

    n_ok = sum(1 for s in summaries if s.get("status") == "ok")
    print(
        f"[compare_structures] pairs={len(pairs)} ok={n_ok} "
        f"substitutions={len(all_subs)} → {out}"
    )


if __name__ == "__main__":
    main()
