#!/usr/bin/env python3
"""Calculate pairwise protein % identity from OrthoFinder OG alignments."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--og-list", required=True, help="File with Orthogroup IDs")
    p.add_argument(
        "--alignments-dir",
        required=True,
        help="Directory with OG*.fa alignments (OrthoFinder MultipleSequenceAlignments)",
    )
    p.add_argument("--selected-ogs-tsv", required=True, help="Output of select_orthogroups.py")
    p.add_argument(
        "--reference-species-id",
        default="auto",
        help="Species id for %%ID reference, or 'auto' to use the matched GOI member",
    )
    p.add_argument("--out-pairwise", required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


def percent_identity(seq_a: str, seq_b: str) -> float | None:
    """Identity over positions where neither sequence has a gap."""
    assert len(seq_a) == len(seq_b)
    matches = 0
    compared = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a.upper() == b.upper():
            matches += 1
    if compared == 0:
        return None
    return 100.0 * matches / compared


def find_alignment(align_dir: Path, og: str) -> Path | None:
    candidates = [
        align_dir / f"{og}.fa",
        align_dir / f"{og}.faa",
        align_dir / f"{og}.fasta",
        align_dir / f"{og}.aln",
    ]
    for c in candidates:
        if c.exists():
            return c
    # recursive fallback
    hits = list(align_dir.rglob(f"{og}.*"))
    for h in hits:
        if h.suffix.lower() in {".fa", ".faa", ".fasta", ".aln"}:
            return h
    return None


def pick_reference(
    ids: list[str],
    ref_species: str,
    matched_members: list[str],
    query_species_ids: list[str],
) -> str | None:
    """Choose %ID reference: matched GOI member (auto) or a fixed species."""
    id_set = set(ids)

    # 1) Always prefer an actual matched GOI member present in the alignment
    for mid in matched_members:
        if mid in id_set:
            if ref_species in ("", "auto", None):
                return mid
            if mid.startswith(f"{ref_species}|"):
                return mid

    # 2) Fixed species fallback
    if ref_species and ref_species != "auto":
        for mid in ids:
            if mid.startswith(f"{ref_species}|"):
                return mid

    # 3) Species inferred from the GOI query
    for sp in query_species_ids:
        if not sp:
            continue
        for mid in matched_members:
            if mid in id_set and mid.startswith(f"{sp}|"):
                return mid
        for mid in ids:
            if mid.startswith(f"{sp}|"):
                return mid

    # 4) Any matched member, else first sequence
    for mid in matched_members:
        if mid in id_set:
            return mid
    return ids[0] if ids else None


def main() -> None:
    args = parse_args()
    align_dir = Path(args.alignments_dir)

    og_meta = {}
    with open(args.selected_ogs_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            og_meta[row["Orthogroup"]] = row

    with open(args.og_list) as fh:
        ogs = [ln.strip() for ln in fh if ln.strip()]

    pairwise_rows = []
    summary_rows = []

    for og in ogs:
        aln_path = find_alignment(align_dir, og)
        meta = og_meta.get(og, {})
        matched = [m for m in (meta.get("matched_members") or "").split(";") if m]
        query_species = [s for s in (meta.get("query_species_ids") or "").split(";") if s]

        if aln_path is None:
            summary_rows.append(
                {
                    "Orthogroup": og,
                    "genes_of_interest": meta.get("genes_of_interest", ""),
                    "query_species_ids": meta.get("query_species_ids", ""),
                    "reference_id": "",
                    "n_sequences": 0,
                    "mean_pct_id_to_ref": "",
                    "min_pct_id_to_ref": "",
                    "max_pct_id_to_ref": "",
                    "alignment_file": "",
                    "status": "alignment_missing",
                }
            )
            continue

        alignment: MultipleSeqAlignment = AlignIO.read(str(aln_path), "fasta")
        ids = [r.id for r in alignment]
        seqs = {r.id: str(r.seq) for r in alignment}
        ref = pick_reference(ids, args.reference_species_id, matched, query_species)
        if ref is None:
            continue

        ids_to_ref = []
        for other in ids:
            if other == ref:
                continue
            pid = percent_identity(seqs[ref], seqs[other])
            if pid is None:
                continue
            ids_to_ref.append(pid)
            pairwise_rows.append(
                {
                    "Orthogroup": og,
                    "genes_of_interest": meta.get("genes_of_interest", ""),
                    "query_species_ids": meta.get("query_species_ids", ""),
                    "query_id": ref,
                    "target_id": other,
                    "target_species_id": other.split("|")[0] if "|" in other else "",
                    "pct_identity": f"{pid:.2f}",
                    "alignment_file": str(aln_path.name),
                }
            )

        summary_rows.append(
            {
                "Orthogroup": og,
                "genes_of_interest": meta.get("genes_of_interest", ""),
                "query_species_ids": meta.get("query_species_ids", ""),
                "reference_id": ref,
                "n_sequences": len(ids),
                "mean_pct_id_to_ref": f"{(sum(ids_to_ref) / len(ids_to_ref)):.2f}" if ids_to_ref else "",
                "min_pct_id_to_ref": f"{min(ids_to_ref):.2f}" if ids_to_ref else "",
                "max_pct_id_to_ref": f"{max(ids_to_ref):.2f}" if ids_to_ref else "",
                "alignment_file": aln_path.name,
                "status": "ok",
            }
        )

    Path(args.out_pairwise).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_pairwise, "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "genes_of_interest",
            "query_species_ids",
            "query_id",
            "target_id",
            "target_species_id",
            "pct_identity",
            "alignment_file",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(pairwise_rows)

    with open(args.out_summary, "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "genes_of_interest",
            "query_species_ids",
            "reference_id",
            "n_sequences",
            "mean_pct_id_to_ref",
            "min_pct_id_to_ref",
            "max_pct_id_to_ref",
            "alignment_file",
            "status",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    print(f"[calc_percent_identity] pairwise={len(pairwise_rows)} summary={len(summary_rows)}")


if __name__ == "__main__":
    main()
