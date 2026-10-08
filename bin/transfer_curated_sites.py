#!/usr/bin/env python3
"""Transfer curated residue numbers from a reference protein to orthologs via MSA.

Curated sites are defined on ONE reference sequence (e.g. tomato GH9 catalytic
residues). This script finds the OrthoFinder protein MSA for the OG, maps each
reference residue to an alignment column, then to the ungapped position in every
other sequence — never copying raw residue numbers across species.

Input YAML (conf/functional_feature_rules.yaml curated_sites section), e.g.:

  curated_sites:
    - orthogroup: OG0001234          # optional but recommended
      reference_id: Sly|Solanum_lycopersicum|Solyc03g123456.1.1
      feature: active_site
      residues: [198, 201, 310]      # 1-based on the REFERENCE ungapped protein
      note: "GH9 catalytic triad (example)"

Or dict form keyed by reference_id (legacy):

  curated_sites:
    "Sly|...|Solyc...":
      - feature: active_site
        residues: [198, 201]
        orthogroup: OG0001234
        note: "..."

Output: transferred_sites.tsv suitable to merge into sequence_features.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import AlignIO

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rules", required=True, help="YAML with curated_sites")
    p.add_argument("--alignments-dir", required=True, help="OrthoFinder MSA directory")
    p.add_argument("--out-tsv", required=True)
    p.add_argument(
        "--selected-ogs",
        default=None,
        help="Optional Orthogroups.tsv or selected_orthogroups.tsv to locate reference→OG",
    )
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


def load_curated(rules_path: str) -> list[dict]:
    if yaml is None:
        raise SystemExit("PyYAML required")
    data = yaml.safe_load(Path(rules_path).read_text()) or {}
    raw = data.get("curated_sites") or {}
    entries = []
    if isinstance(raw, list):
        entries = raw
    elif isinstance(raw, dict):
        for ref_id, sites in raw.items():
            for s in sites or []:
                e = dict(s)
                e["reference_id"] = ref_id
                entries.append(e)
    return entries


def find_og_for_reference(selected_ogs: str | None, ref_id: str) -> str | None:
    if not selected_ogs or not Path(selected_ogs).exists():
        return None
    with open(selected_ogs) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            members = row.get("members") or ""
            if ref_id in members or ref_id in (row.get("matched_members") or ""):
                return row.get("Orthogroup")
    return None


def seq_pos_to_col(aligned: str, seq_pos_1based: int) -> int | None:
    """Map 1-based ungapped residue → 0-based alignment column."""
    n = 0
    for col, aa in enumerate(aligned):
        if aa != "-":
            n += 1
            if n == seq_pos_1based:
                return col
    return None


def col_to_seq_pos(aligned: str, col: int) -> int | None:
    """Map 0-based alignment column → 1-based ungapped residue (None if gap)."""
    if col < 0 or col >= len(aligned) or aligned[col] == "-":
        return None
    return sum(1 for aa in aligned[: col + 1] if aa != "-")


def resolve_seq_id(wanted: str, ids: list[str]) -> str | None:
    if wanted in ids:
        return wanted
    for i in ids:
        if i.endswith(wanted) or wanted.endswith(i) or wanted in i or i in wanted:
            return i
    # match by trailing locus
    wtail = wanted.split("|")[-1]
    for i in ids:
        if i.split("|")[-1] == wtail or i.split("|")[-1].startswith(wtail.split(".")[0]):
            return i
    return None


def main() -> None:
    args = parse_args()
    entries = load_curated(args.rules)
    align_dir = Path(args.alignments_dir)
    out_rows = []

    for e in entries:
        ref_id = e.get("reference_id")
        residues = e.get("residues") or []
        feature = e.get("feature") or "curated_site"
        note = e.get("note") or ""
        og = e.get("orthogroup") or find_og_for_reference(args.selected_ogs, ref_id)
        if not ref_id or not residues:
            continue
        if not og:
            out_rows.append(
                {
                    "reference_id": ref_id,
                    "target_id": "",
                    "orthogroup": "",
                    "feature_type": feature,
                    "ref_resnum": "",
                    "target_resnum": "",
                    "ref_aa": "",
                    "target_aa": "",
                    "status": "orthogroup_unknown",
                    "note": note,
                }
            )
            continue

        aln_path = find_alignment(align_dir, og)
        if not aln_path:
            for r in residues:
                out_rows.append(
                    {
                        "reference_id": ref_id,
                        "target_id": "",
                        "orthogroup": og,
                        "feature_type": feature,
                        "ref_resnum": r,
                        "target_resnum": "",
                        "ref_aa": "",
                        "target_aa": "",
                        "status": "alignment_missing",
                        "note": note,
                    }
                )
            continue

        aln = AlignIO.read(str(aln_path), "fasta")
        seqs = {rec.id: str(rec.seq).upper() for rec in aln}
        ref_key = resolve_seq_id(ref_id, list(seqs))
        if not ref_key:
            for r in residues:
                out_rows.append(
                    {
                        "reference_id": ref_id,
                        "target_id": "",
                        "orthogroup": og,
                        "feature_type": feature,
                        "ref_resnum": r,
                        "target_resnum": "",
                        "ref_aa": "",
                        "target_aa": "",
                        "status": "reference_not_in_msa",
                        "note": note,
                    }
                )
            continue

        ref_aln = seqs[ref_key]
        for res in residues:
            col = seq_pos_to_col(ref_aln, int(res))
            if col is None:
                out_rows.append(
                    {
                        "reference_id": ref_key,
                        "target_id": "",
                        "orthogroup": og,
                        "feature_type": feature,
                        "ref_resnum": res,
                        "target_resnum": "",
                        "ref_aa": "",
                        "target_aa": "",
                        "status": "ref_residue_out_of_range",
                        "note": note,
                    }
                )
                continue
            ref_aa = ref_aln[col]
            for tid, taln in seqs.items():
                tpos = col_to_seq_pos(taln, col)
                if tpos is None:
                    status = "gap_in_target"
                    taa = "-"
                    tpos_out = ""
                else:
                    status = "ok"
                    taa = taln[col]
                    tpos_out = tpos
                out_rows.append(
                    {
                        "reference_id": ref_key,
                        "target_id": tid,
                        "orthogroup": og,
                        "feature_type": feature,
                        "ref_resnum": res,
                        "target_resnum": tpos_out,
                        "ref_aa": ref_aa,
                        "target_aa": taa,
                        "msa_column": col + 1,
                        "status": status,
                        "note": note,
                    }
                )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "orthogroup",
        "reference_id",
        "target_id",
        "feature_type",
        "ref_resnum",
        "target_resnum",
        "ref_aa",
        "target_aa",
        "msa_column",
        "status",
        "note",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    n_ok = sum(1 for r in out_rows if r["status"] == "ok")
    print(f"[transfer_curated_sites] rows={len(out_rows)} ok={n_ok} → {args.out_tsv}")


if __name__ == "__main__":
    main()
