#!/usr/bin/env python3
"""QC/rename CDS to match cleaned proteome IDs from qc_and_rename_proteome.py.

Keeps CDS that map to retained proteins (via --id-map) and assigns the same
new IDs (Sly|Genus_species|...). Trims transcripts/cDNA to the ORF that
translates to the matched protein when UTRs are present.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cds", required=True)
    p.add_argument("--proteins", required=True, help="Cleaned/renamed protein FASTA")
    p.add_argument("--id-map", required=True, help="Proteome id_map.tsv")
    p.add_argument("--out-cds", required=True)
    p.add_argument("--out-report", required=True)
    p.add_argument("--trim-to-protein", action="store_true", default=True)
    p.add_argument("--no-trim-to-protein", action="store_false", dest="trim_to_protein")
    return p.parse_args()


def open_text_auto(path: str):
    p = Path(path)
    with p.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(p, "rt")
    return p.open("rt")


def find_orf_matching_protein(cds: str, protein: str) -> str | None:
    """Return nt substring whose translation equals protein."""
    cds = cds.upper().replace("U", "T")
    protein = protein.upper().rstrip("*")
    if not cds or not protein:
        return None
    expected = len(protein) * 3

    # Exact CDS length + frame 0
    if len(cds) >= expected:
        chunk = cds[:expected]
        aa = str(Seq(chunk).translate(to_stop=False)).rstrip("*")
        if aa == protein:
            return chunk

    for frame in (0, 1, 2):
        seq = cds[frame:]
        seq = seq[: len(seq) // 3 * 3]
        if len(seq) < expected:
            continue
        aa = str(Seq(seq).translate(to_stop=False))
        idx = aa.find(protein)
        if idx >= 0:
            return seq[idx * 3 : idx * 3 + expected]
    return None


def lookup_new_id(seq_id: str, orig_to_new: dict[str, str], gene_to_new: dict[str, str]) -> str | None:
    sid = seq_id.split()[0]
    if sid in orig_to_new:
        return orig_to_new[sid]
    if seq_id in orig_to_new:
        return orig_to_new[seq_id]
    for orig, new in orig_to_new.items():
        if sid == orig or sid.startswith(orig + ".") or orig.startswith(sid):
            return new
    for gk, new in gene_to_new.items():
        if sid == gk or sid.startswith(gk):
            return new
    return None


def main() -> None:
    args = parse_args()

    orig_to_new: dict[str, str] = {}
    gene_to_new: dict[str, str] = {}
    with open(args.id_map) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            orig_to_new[row["original_id"]] = row["new_id"]
            gene_to_new[row["gene_key"]] = row["new_id"]

    proteins = {
        r.id: str(r.seq).upper().rstrip("*") for r in SeqIO.parse(args.proteins, "fasta")
    }

    with open_text_auto(args.cds) as handle:
        cds_records = list(SeqIO.parse(handle, "fasta"))

    kept = []
    seen = set()
    n_ok = n_trim = n_fail = n_dup = n_nomatch = 0

    for rec in cds_records:
        new_id = lookup_new_id(rec.id, orig_to_new, gene_to_new)
        if not new_id or new_id not in proteins:
            n_nomatch += 1
            continue
        if new_id in seen:
            n_dup += 1
            continue

        cds_seq = str(rec.seq).upper().replace("U", "T")
        prot = proteins[new_id]

        if args.trim_to_protein:
            final = find_orf_matching_protein(cds_seq, prot)
            if final is None:
                n_fail += 1
                continue
            if final != cds_seq and len(final) != len(cds_seq):
                n_trim += 1
            else:
                n_ok += 1
        else:
            if len(cds_seq) != len(prot) * 3:
                n_fail += 1
                continue
            final = cds_seq
            n_ok += 1

        rec.id = new_id
        rec.name = new_id
        rec.description = ""
        rec.seq = Seq(final)
        kept.append(rec)
        seen.add(new_id)

    Path(args.out_cds).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(kept, args.out_cds, "fasta")

    missing = sorted(set(proteins) - seen)
    with open(args.out_report, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "n_cds_input",
                "n_proteins",
                "n_cds_kept",
                "n_exact_or_ok",
                "n_trimmed_orf",
                "n_orf_fail",
                "n_no_protein_match",
                "n_duplicate",
                "n_proteins_missing_cds",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerow(
            {
                "n_cds_input": len(cds_records),
                "n_proteins": len(proteins),
                "n_cds_kept": len(kept),
                "n_exact_or_ok": n_ok,
                "n_trimmed_orf": n_trim,
                "n_orf_fail": n_fail,
                "n_no_protein_match": n_nomatch,
                "n_duplicate": n_dup,
                "n_proteins_missing_cds": len(missing),
            }
        )

    print(
        f"[qc_and_rename_cds] kept={len(kept)}/{len(proteins)} proteins "
        f"(trim={n_trim}, orf_fail={n_fail}, no_match={n_nomatch})"
    )


if __name__ == "__main__":
    main()
