#!/usr/bin/env python3
"""Record where a proteome download came from, and refuse an organism mismatch.

NCBI Datasets writes assembly_data_report.jsonl into the genome package.
That report is the organism, accession, and annotation actually retrieved,
which can differ from the genus and species in conf/species.yaml.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


COLUMNS = [
    "species_id",
    "configured_taxon",
    "source",
    "proteins_source",
    "cds_source",
    "gff_source",
    "proteins_file",
    "cds_file",
    "gff_file",
    "retrieved_utc",
    "tool",
    "remote_organism",
    "remote_taxid",
    "remote_accession",
    "assembly_name",
    "annotation_name",
    "annotation_date",
    "expected_protein_coding",
    "n_proteins",
    "n_cds",
    "proteins_md5",
    "proteins_bytes",
    "sample_header",
    "warning",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--species-id", required=True)
    p.add_argument("--genus", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--source", required=True, choices=("ncbi", "url", "local"))
    p.add_argument("--proteins-source", required=True, help="Accession, URL, or local path requested")
    p.add_argument("--cds-source", default="")
    p.add_argument("--gff-source", default="")
    p.add_argument("--proteins", required=True, help="Staged protein FASTA")
    p.add_argument("--cds", default="")
    p.add_argument("--gff", default="")
    p.add_argument("--proteins-file", default="", help="Path inside the download, or the URL")
    p.add_argument("--cds-file", default="")
    p.add_argument("--gff-file", default="")
    p.add_argument("--assembly-report", default="", help="NCBI assembly_data_report.jsonl")
    p.add_argument("--tool", default="")
    p.add_argument("--out", required=True)
    return p.parse_args()


def first(obj: dict, *paths: tuple[str, ...]):
    for path in paths:
        cur = obj
        ok = True
        for key in path:
            if not isinstance(cur, dict) or key not in cur:
                ok = False
                break
            cur = cur[key]
        if ok and cur not in (None, ""):
            return cur
    return ""


def load_assembly_report(path: str) -> dict:
    file = Path(path) if path else None
    if file is None or not file.is_file() or file.stat().st_size == 0:
        return {}
    for line in file.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            return json.loads(line)
        except json.JSONDecodeError:
            continue
    return {}


def count_fasta(path: str) -> tuple[int, str]:
    file = Path(path) if path else None
    if file is None or not file.is_file() or file.stat().st_size == 0:
        return 0, ""
    n = 0
    header = ""
    with file.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
                if not header:
                    header = line[1:].strip()
    return n, header[:160]


def md5(path: str) -> tuple[str, int]:
    file = Path(path)
    digest = hashlib.md5()
    size = 0
    with file.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            size += len(chunk)
            digest.update(chunk)
    return digest.hexdigest(), size


def organism_matches(genus: str, organism: str) -> bool:
    if not organism:
        return True
    return re.search(rf"\b{re.escape(genus)}\b", organism, re.IGNORECASE) is not None


def main() -> None:
    args = parse_args()
    report = load_assembly_report(args.assembly_report)
    remote_organism = str(
        first(report, ("organism", "organismName"), ("organism", "organism_name"))
    )
    remote_taxid = str(first(report, ("organism", "taxId"), ("organism", "tax_id")))
    remote_accession = str(first(report, ("accession",), ("currentAccession",)))
    assembly_name = str(
        first(report, ("assemblyInfo", "assemblyName"), ("assembly_info", "assembly_name"))
    )
    annotation_name = str(first(report, ("annotationInfo", "name"), ("annotation_info", "name")))
    annotation_date = str(
        first(report, ("annotationInfo", "releaseDate"), ("annotation_info", "release_date"))
    )
    expected = first(
        report,
        ("annotationInfo", "stats", "geneCounts", "proteinCoding"),
        ("annotation_info", "stats", "gene_counts", "protein_coding"),
    )
    expected_text = str(expected) if expected != "" else ""

    n_proteins, sample_header = count_fasta(args.proteins)
    n_cds, _cds_header = count_fasta(args.cds)
    proteins_md5, proteins_bytes = md5(args.proteins)

    warnings: list[str] = []
    if remote_organism and not organism_matches(args.genus, remote_organism):
        warnings.append(
            f"remote organism '{remote_organism}' does not match configured "
            f"{args.genus} {args.species}"
        )
    if expected_text.isdigit() and n_proteins < int(expected_text) * 0.5:
        warnings.append(
            f"protein FASTA has {n_proteins} sequences but the annotation reports "
            f"{expected_text} protein-coding genes"
        )

    row = {
        "species_id": args.species_id,
        "configured_taxon": f"{args.genus} {args.species}",
        "source": args.source,
        "proteins_source": args.proteins_source,
        "cds_source": args.cds_source,
        "gff_source": args.gff_source,
        "proteins_file": args.proteins_file or Path(args.proteins).name,
        "cds_file": args.cds_file,
        "gff_file": args.gff_file,
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tool": args.tool,
        "remote_organism": remote_organism,
        "remote_taxid": remote_taxid,
        "remote_accession": remote_accession or (args.proteins_source if args.source == "ncbi" else ""),
        "assembly_name": assembly_name,
        "annotation_name": annotation_name,
        "annotation_date": annotation_date,
        "expected_protein_coding": expected_text,
        "n_proteins": n_proteins,
        "n_cds": n_cds,
        "proteins_md5": proteins_md5,
        "proteins_bytes": proteins_bytes,
        "sample_header": sample_header,
        "warning": "; ".join(warnings),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    summary = (
        f"[provenance] {args.species_id}: source={args.source} "
        f"requested={args.proteins_source} proteins={n_proteins}"
    )
    if remote_organism:
        summary += f" organism='{remote_organism}'"
        if remote_accession:
            summary += f" accession={remote_accession}"
        if expected_text:
            summary += f" annotation_proteins={expected_text}"
    print(summary)
    if warnings:
        print("[provenance] ERROR: " + "; ".join(warnings), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
