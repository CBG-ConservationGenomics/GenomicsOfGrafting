#!/usr/bin/env python3
"""Build codon MSAs from OrthoFinder protein alignments + matching CDS (pal2nal-style).

For each selected orthogroup protein MSA, project gaps onto CDS codons:
  protein '-' → '---' in codon alignment; otherwise next 3 nt from CDS.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--og-list", required=True)
    p.add_argument("--alignments-dir", required=True, help="OrthoFinder MultipleSequenceAlignments")
    p.add_argument(
        "--cds-fastas",
        nargs="+",
        required=True,
        help="One or more cleaned CDS FASTAs (IDs match protein MSA headers)",
    )
    p.add_argument("--outdir", required=True)
    p.add_argument("--out-report", required=True)
    p.add_argument(
        "--min-sequences",
        type=int,
        default=3,
        help="Skip OGs with fewer successfully projected sequences",
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


def project_codon(protein_aln: str, cds: str) -> str | None:
    """Map protein aligned string + ungapped CDS → codon alignment string."""
    cds = cds.upper().replace("U", "T")
    protein_aln = protein_aln.upper()
    aa_ungapped = protein_aln.replace("-", "").rstrip("*")
    expected = len(aa_ungapped) * 3
    if len(cds) < expected:
        return None
    # Allow longer CDS only if prefix matches translation
    cds_use = cds[:expected]
    try:
        translated = str(Seq(cds_use).translate(to_stop=False)).rstrip("*")
    except Exception:
        return None
    if translated != aa_ungapped:
        # try ignoring trailing partial
        if not aa_ungapped.startswith(translated) and not translated.startswith(aa_ungapped):
            return None
        if len(translated) < len(aa_ungapped):
            return None

    out = []
    i = 0  # codon index into cds_use
    for aa in protein_aln:
        if aa == "-":
            out.append("---")
        elif aa == "*":
            # stop in alignment — consume codon if present
            if i * 3 + 3 <= len(cds_use):
                out.append(cds_use[i * 3 : i * 3 + 3])
                i += 1
            else:
                out.append("---")
        else:
            if i * 3 + 3 > len(cds_use):
                return None
            out.append(cds_use[i * 3 : i * 3 + 3])
            i += 1
    return "".join(out)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    align_dir = Path(args.alignments_dir)

    cds_index = {}
    for fa in args.cds_fastas:
        for rec in SeqIO.parse(fa, "fasta"):
            cds_index[rec.id] = str(rec.seq).upper().replace("U", "T")

    with open(args.og_list) as fh:
        ogs = [ln.strip() for ln in fh if ln.strip()]

    report = []
    n_written = 0

    for og in ogs:
        aln_path = find_alignment(align_dir, og)
        if not aln_path:
            report.append({"Orthogroup": og, "n_in_msa": 0, "n_codon": 0, "status": "alignment_missing"})
            continue

        aln = AlignIO.read(str(aln_path), "fasta")
        codon_recs = []
        n_miss_cds = 0
        n_proj_fail = 0

        for rec in aln:
            cds = cds_index.get(rec.id)
            if not cds:
                n_miss_cds += 1
                continue
            projected = project_codon(str(rec.seq), cds)
            if not projected:
                n_proj_fail += 1
                continue
            codon_recs.append(SeqRecord(Seq(projected), id=rec.id, description=""))

        status = "ok"
        if len(codon_recs) < args.min_sequences:
            status = "too_few_sequences"
        else:
            SeqIO.write(codon_recs, outdir / f"{og}.codon.fna", "fasta")
            n_written += 1

        report.append(
            {
                "Orthogroup": og,
                "n_in_msa": len(aln),
                "n_codon": len(codon_recs),
                "n_missing_cds": n_miss_cds,
                "n_project_fail": n_proj_fail,
                "status": status,
            }
        )

    Path(args.out_report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_report, "w", newline="") as fh:
        fields = [
            "Orthogroup",
            "n_in_msa",
            "n_codon",
            "n_missing_cds",
            "n_project_fail",
            "status",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(report)

    print(f"[build_codon_alignments] wrote {n_written}/{len(ogs)} codon MSAs → {outdir}")


if __name__ == "__main__":
    main()
