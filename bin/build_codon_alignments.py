#!/usr/bin/env python3
"""Build codon MSAs from OrthoFinder protein alignments + matching CDS (pal2nal-style).

For each selected orthogroup protein MSA, project gaps onto CDS codons:
  protein '-' → '---' in codon alignment; otherwise next 3 nt from CDS.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--og-list", required=True)
    p.add_argument("--alignments-dir", required=True, help="OrthoFinder MultipleSequenceAlignments")
    p.add_argument("--cds-fastas", nargs="+", required=True, help="Cleaned CDS FASTAs")
    p.add_argument(
        "--orthofinder-results",
        default="",
        help="OrthoFinder Results directory, used to read SequenceIDs.txt",
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


def id_keys(seq_id: str) -> list[str]:
    """Forms of one sequence id, longest first.

    OrthoFinder often rewrites '|' as '_' in alignment headers, or keeps only
    the original accession. CDS files use the proteome id
    (Sly|Genus_species|Solyc...).
    """
    seq_id = seq_id.strip()
    if not seq_id:
        return []
    keys = [seq_id, seq_id.replace("|", "_")]
    tail = seq_id.split("|")[-1]
    if tail and tail not in keys:
        keys.append(tail)
    parts = seq_id.replace("|", "_").split("_")
    for i in range(1, len(parts)):
        suffix = "_".join(parts[i:])
        if len(suffix) >= 4 and suffix not in keys:
            keys.append(suffix)
    return keys


def load_sequence_id_map(root: Path) -> dict[str, str]:
    """OrthoFinder internal id (0_42) -> original FASTA accession."""
    mapping: dict[str, str] = {}
    if not root.exists():
        return mapping
    paths = [root] if root.is_file() else list(root.glob("**/SequenceIDs.txt"))
    for path in paths:
        if not path.is_file():
            continue
        for line in path.read_text(errors="replace").splitlines():
            if ":" not in line:
                continue
            internal, original = line.split(":", 1)
            original = original.strip().split()[0].rstrip(",")
            if internal.strip() and original:
                mapping[internal.strip()] = original
    return mapping


def index_cds(fastas: list[str]) -> dict[str, str]:
    index: dict[str, str] = {}
    for fa in fastas:
        for rec in SeqIO.parse(fa, "fasta"):
            seq = str(rec.seq).upper().replace("U", "T")
            for key in id_keys(rec.id):
                index.setdefault(key, seq)
    return index


def lookup_cds(seq_id: str, index: dict[str, str], id_map: dict[str, str]) -> str | None:
    candidates = [seq_id]
    if seq_id in id_map:
        candidates.append(id_map[seq_id])
    seen: set[str] = set()
    for candidate in candidates:
        for key in id_keys(candidate):
            if key in seen:
                continue
            seen.add(key)
            hit = index.get(key)
            if hit:
                return hit
    return None


def find_alignment(align_dir: Path, og: str) -> Path | None:
    preferred = [
        align_dir / "MultipleSequenceAlignments" / f"{og}.fa",
        align_dir / f"{og}.fa",
        align_dir / f"{og}.faa",
        align_dir / f"{og}.fasta",
        align_dir / f"{og}.aln",
    ]
    for candidate in preferred:
        if candidate.exists():
            return candidate
    hits = list(align_dir.rglob(f"{og}.*"))
    msa_hits = [hit for hit in hits if "MultipleSequenceAlignments" in hit.parts]
    for hit in msa_hits or hits:
        if hit.suffix.lower() in {".fa", ".faa", ".fasta", ".aln"}:
            return hit
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

    cds_index = index_cds(args.cds_fastas)
    id_map = load_sequence_id_map(Path(args.orthofinder_results)) if args.orthofinder_results else {}
    if not cds_index:
        print("[build_codon_alignments] WARNING: no CDS sequences were read", file=sys.stderr)

    with open(args.og_list) as fh:
        ogs = [ln.strip() for ln in fh if ln.strip()]

    report = []
    n_written = 0

    for og in ogs:
        aln_path = find_alignment(align_dir, og)
        if not aln_path:
            kind = "unassigned_no_alignment" if og.startswith("UNASSIGNED_") else "alignment_missing"
            report.append(
                {
                    "Orthogroup": og,
                    "n_in_msa": 0,
                    "n_codon": 0,
                    "n_missing_cds": 0,
                    "n_project_fail": 0,
                    "status": kind,
                }
            )
            continue

        aln = AlignIO.read(str(aln_path), "fasta")
        codon_recs = []
        n_miss_cds = 0
        n_proj_fail = 0

        for rec in aln:
            cds = lookup_cds(rec.id, cds_index, id_map)
            if not cds:
                n_miss_cds += 1
                if n_miss_cds == 1:
                    example_msa = rec.id
                continue
            projected = project_codon(str(rec.seq), cds)
            if not projected:
                n_proj_fail += 1
                continue
            codon_recs.append(SeqRecord(Seq(projected), id=rec.id, description=""))

        if n_miss_cds and n_miss_cds == len(aln):
            example_cds = next(iter(cds_index), "")
            print(
                f"[build_codon_alignments] {og}: 0/{len(aln)} alignment sequences matched a CDS. "
                f"Example alignment id: {example_msa}. Example CDS id: {example_cds}",
                file=sys.stderr,
            )
        if len(codon_recs) < args.min_sequences:
            status = "too_few_sequences"
        else:
            status = "ok"
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
