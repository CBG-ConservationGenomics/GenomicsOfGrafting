#!/usr/bin/env python3
"""Calculate pairwise protein % identity from OrthoFinder OG alignments."""

from __future__ import annotations

import argparse
import csv
import re
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
    p.add_argument("--out-matrix-pairs", required=True, help="All-by-all identity, one pair per row")
    p.add_argument("--out-matrix-labels", required=True, help="Short label for each protein in a matrix")
    p.add_argument("--out-matrix-dir", required=True, help="Directory of square identity matrices")
    return p.parse_args()


def short_label(seq_id: str) -> str:
    """Species plus locus, without the isoform suffix when the locus pattern is known."""
    parts = seq_id.split("|")
    species = parts[0] if len(parts) >= 2 else ""
    tail = parts[-1]
    match = re.match(
        r"^(Solyc\d+g\d+|AT[1-5MC]G\d+|Niben\d+Scf\w+g\d+|Niben\w+g\d+|LOC\d+|(?:XP|NP|WP|YP)_\d+)",
        tail,
        re.IGNORECASE,
    )
    locus = match.group(1) if match else re.sub(r"\.\d+$", "", tail)
    return f"{species}:{locus}" if species else locus


def unique_labels(seq_ids: list[str]) -> list[str]:
    """Keep paralogs distinct when they share a locus label."""
    bases = [short_label(seq_id) for seq_id in seq_ids]
    counts: dict[str, int] = {}
    for base in bases:
        counts[base] = counts.get(base, 0) + 1
    labels = []
    for seq_id, base in zip(seq_ids, bases):
        if counts[base] == 1:
            labels.append(base)
            continue
        tail = seq_id.split("|")[-1]
        species = seq_id.split("|")[0] if "|" in seq_id else ""
        labels.append(f"{species}:{tail}" if species else tail)
    return labels


def order_sequences(seq_ids: list[str], reference: str) -> list[str]:
    """Characterized protein first, then its same-species homologs, then the other species."""
    ref_species = reference.split("|")[0] if "|" in reference else ""

    def sort_key(seq_id: str) -> tuple:
        species = seq_id.split("|")[0] if "|" in seq_id else seq_id
        return (
            0 if seq_id == reference else 1,
            0 if species == ref_species else 1,
            species,
            seq_id,
        )

    return sorted(seq_ids, key=sort_key)


def matrix_rows(
    og: str,
    seq_ids: list[str],
    seqs: dict[str, str],
    reference: str,
) -> tuple[list[dict], list[dict], list[str]]:
    """All-by-all identity, a label table, and the ordered sequence ids."""
    ordered = order_sequences(seq_ids, reference)
    labels = unique_labels(ordered)
    pairs = []
    label_rows = []
    for seq_id, label in zip(ordered, labels):
        species = seq_id.split("|")[0] if "|" in seq_id else ""
        label_rows.append(
            {
                "Orthogroup": og,
                "label": label,
                "sequence_id": seq_id,
                "species_id": species,
                "is_reference": "yes" if seq_id == reference else "",
            }
        )
    for i, left in enumerate(ordered):
        for right in ordered[i + 1 :]:
            pid = percent_identity(seqs[left], seqs[right])
            if pid is None:
                continue
            left_i = ordered.index(left)
            right_i = ordered.index(right)
            pairs.append(
                {
                    "Orthogroup": og,
                    "seq_a": left,
                    "seq_b": right,
                    "label_a": labels[left_i],
                    "label_b": labels[right_i],
                    "species_a": left.split("|")[0] if "|" in left else "",
                    "species_b": right.split("|")[0] if "|" in right else "",
                    "pct_identity": f"{pid:.2f}",
                    "reference_id": reference,
                }
            )
    return pairs, label_rows, ordered


def write_square_matrix(
    path: Path,
    og: str,
    ordered: list[str],
    labels: list[str],
    seqs: dict[str, str],
    reference: str,
) -> None:
    lookup = {seq_id: label for seq_id, label in zip(ordered, labels)}
    with path.open("w", newline="") as fh:
        fh.write(f"# Orthogroup\t{og}\n")
        fh.write(f"# reference\t{lookup.get(reference, reference)}\t{reference}\n")
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sequence", *labels])
        for left in ordered:
            row = [lookup[left]]
            for right in ordered:
                if left == right:
                    row.append("")
                    continue
                pid = percent_identity(seqs[left], seqs[right])
                row.append("" if pid is None else f"{pid:.1f}")
            writer.writerow(row)


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


def align_id_for(target: str, ids: list[str]) -> str | None:
    """Find an alignment sequence that is the same protein as `target`."""
    if not target:
        return None
    id_set = set(ids)
    if target in id_set:
        return target
    flat = target.replace("|", "_")
    if flat in id_set:
        return flat
    tail = target.split("|")[-1]
    for seq_id in ids:
        if seq_id == tail or seq_id.endswith("|" + tail) or seq_id.replace("|", "_") == flat:
            return seq_id
    return None


def pick_reference(
    ids: list[str],
    ref_species: str,
    matched_members: list[str],
    query_species_ids: list[str],
    reference_member: str = "",
) -> str | None:
    """Prefer the literature gene, then a matched GOI, then a fixed species."""
    id_set = set(ids)

    hit = align_id_for(reference_member, ids)
    if hit:
        return hit

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
    matrix_pairs: list[dict] = []
    matrix_labels: list[dict] = []
    matrix_dir = Path(args.out_matrix_dir)
    matrix_dir.mkdir(parents=True, exist_ok=True)
    (matrix_dir / "README.txt").write_text(
        "One square percent-identity matrix per orthogroup.\n"
        "Each row and column is one protein. Paralogs from the same species are separate columns.\n"
        "The characterized protein is listed first. The diagonal is blank.\n"
        "Labels are species:locus. See identity_matrix_labels.tsv for the full sequence id.\n"
    )

    for og in ogs:
        aln_path = find_alignment(align_dir, og)
        meta = og_meta.get(og, {})
        matched = [m for m in (meta.get("matched_members") or "").split(";") if m]
        query_species = [s for s in (meta.get("query_species_ids") or "").split(";") if s]
        reference_member = meta.get("reference_member") or ""

        if aln_path is None:
            summary_rows.append(
                {
                    "Orthogroup": og,
                    "genes_of_interest": meta.get("genes_of_interest", ""),
                    "query_species_ids": meta.get("query_species_ids", ""),
                    "reference_id": reference_member,
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
        ref = pick_reference(
            ids, args.reference_species_id, matched, query_species, reference_member
        )
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

        pairs, label_rows, ordered = matrix_rows(og, ids, seqs, ref)
        matrix_pairs.extend(pairs)
        matrix_labels.extend(label_rows)
        write_square_matrix(
            matrix_dir / f"{og}.identity_matrix.tsv",
            og,
            ordered,
            [row["label"] for row in label_rows],
            seqs,
            ref,
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

    pair_fields = [
        "Orthogroup",
        "seq_a",
        "seq_b",
        "label_a",
        "label_b",
        "species_a",
        "species_b",
        "pct_identity",
        "reference_id",
    ]
    with open(args.out_matrix_pairs, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pair_fields, delimiter="\t")
        w.writeheader()
        w.writerows(matrix_pairs)
    with open(args.out_matrix_labels, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["Orthogroup", "label", "sequence_id", "species_id", "is_reference"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(matrix_labels)

    print(
        f"[calc_percent_identity] pairwise={len(pairwise_rows)} summary={len(summary_rows)} "
        f"matrix_pairs={len(matrix_pairs)}"
    )


if __name__ == "__main__":
    main()
