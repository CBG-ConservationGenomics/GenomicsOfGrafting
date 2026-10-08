#!/usr/bin/env python3
"""QC plant proteomes: drop isoforms and alt-haplotig proteins, then rename IDs."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


ISOFORM_SUFFIX = re.compile(
    r"(?:\.(\d+)|_iso(\d+)|-T(\d+)|\.t(\d+)|_t(\d+))$", re.IGNORECASE
)
GENE_FROM_PROTEIN = re.compile(
    r"^(?P<gene>.+?)(?:\.\d+|_iso\d+|-T\d+|\.t\d+|_t\d+)?$", re.IGNORECASE
)
LOCUS_RES = [
    re.compile(r"^(Solyc\d+g\d+)", re.IGNORECASE),
    re.compile(r"^(AT[1-5MC]G\d+)", re.IGNORECASE),
    re.compile(r"^(Niben\d+Scf\w+g\d+)", re.IGNORECASE),
    re.compile(r"^(Niben\w+g\d+)", re.IGNORECASE),
    re.compile(r"^(Nb\w+g\d+)", re.IGNORECASE),
    re.compile(r"^(Cp[\w.]+g\d+)", re.IGNORECASE),
    re.compile(r"^(LOC\d+)", re.IGNORECASE),
    re.compile(r"^((?:XP|NP|WP|YP)_\d+)", re.IGNORECASE),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True, help="Input protein FASTA")
    p.add_argument("--species-id", required=True, help="Short species id, e.g. Sly")
    p.add_argument("--genus", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-map", required=True, help="TSV mapping original→new IDs")
    p.add_argument("--gff", default=None, help="Optional GFF3 for haplotig filtering")
    p.add_argument("--keep-longest-isoform", action="store_true", default=True)
    p.add_argument("--no-keep-longest-isoform", action="store_false", dest="keep_longest_isoform")
    p.add_argument("--drop-alt-haplotigs", action="store_true", default=True)
    p.add_argument("--no-drop-alt-haplotigs", action="store_false", dest="drop_alt_haplotigs")
    p.add_argument(
        "--haplotig-patterns",
        default="hap|haplotig|alt_?hap|_H2|unanchored_hap",
        help="Regex (case-insensitive) matching alt-haplotig seqids",
    )
    return p.parse_args()


def gene_key(seq_id: str) -> str:
    """Collapse isoform-like suffixes to a gene locus key."""
    for rx in LOCUS_RES:
        m = rx.match(seq_id)
        if m:
            return m.group(1)
    m = GENE_FROM_PROTEIN.match(seq_id)
    return m.group("gene") if m else seq_id


def load_haplotig_seqids(gff_path: str | None, pattern: re.Pattern[str]) -> set[str]:
    bad: set[str] = set()
    if not gff_path:
        return bad
    path = Path(gff_path)
    if not path.exists():
        return bad
    gene_on_bad_seq: set[str] = set()
    with path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype, _s, _e, _sc, _st, _ph, attrs = parts
            if pattern.search(seqid):
                bad.add(seqid)
                # capture gene/protein names on haplotig scaffolds
                for key in ("ID", "Name", "protein_id", "gene", "locus_tag"):
                    m = re.search(rf"(?:^|;){key}=([^;]+)", attrs)
                    if m:
                        gene_on_bad_seq.add(m.group(1))
    bad.update(gene_on_bad_seq)
    return bad


def is_haplotig_header(header: str, pattern: re.Pattern[str]) -> bool:
    return bool(pattern.search(header))


def open_text_auto(path: str):
    """Open plain text or gzip (detected by magic bytes 1f 8b)."""
    p = Path(path)
    with p.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(p, "rt")
    return p.open("rt")


def main() -> None:
    args = parse_args()
    hap_re = re.compile(args.haplotig_patterns, re.IGNORECASE)
    hap_ids = load_haplotig_seqids(args.gff, hap_re) if args.drop_alt_haplotigs else set()

    with open_text_auto(args.fasta) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    after_hap = []
    for rec in records:
        if args.drop_alt_haplotigs:
            header = f"{rec.id} {rec.description}"
            if rec.id in hap_ids or is_haplotig_header(rec.id, hap_re):
                continue
            # Drop only if a haplotig-like token appears in the ID itself or known GFF set
            if hap_ids and any(tok in hap_ids for tok in header.split()):
                continue
        after_hap.append(rec)

    if args.keep_longest_isoform:
        by_gene: dict[str, list] = defaultdict(list)
        for rec in after_hap:
            by_gene[gene_key(rec.id)].append(rec)
        kept = []
        for _gene, isoforms in by_gene.items():
            isoforms.sort(key=lambda r: (len(r.seq), r.id), reverse=True)
            kept.append(isoforms[0])
    else:
        kept = after_hap

    genus_sp = f"{args.genus}_{args.species}"
    out_records = []
    map_rows = []
    for rec in sorted(kept, key=lambda r: r.id):
        orig = rec.id.replace(" ", "_")
        new_id = f"{args.species_id}|{genus_sp}|{orig}"
        map_rows.append(
            {
                "species_id": args.species_id,
                "genus_species": genus_sp,
                "original_id": orig,
                "gene_key": gene_key(orig),
                "new_id": new_id,
                "length": len(rec.seq),
            }
        )
        rec.id = new_id
        rec.name = new_id
        rec.description = ""
        out_records.append(rec)

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_records, args.out_fasta, "fasta")

    with open(args.out_map, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["species_id", "genus_species", "original_id", "gene_key", "new_id", "length"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(map_rows)

    print(
        f"[qc_and_rename] {args.species_id}: "
        f"input={len(records)} after_haplotig={len(after_hap)} "
        f"output={len(out_records)}"
    )


if __name__ == "__main__":
    main()
