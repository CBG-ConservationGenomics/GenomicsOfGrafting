#!/usr/bin/env python3
"""Confirm each candidate against a characterized CDS from that species.

The gene list is a name and a species. Those names are resolved onto proteins
in the downloaded annotations. This check is optional evidence that the
resolved annotation sequence is the gene that was characterized.

The reference FASTA holds one coding sequence per candidate. A tab separates
the gene name from the species. A locus id and annotation name may follow
the species, separated by spaces:

  >WOX13<TAB>Arabidopsis thaliana AT4G35550.1 TAIR10
  >WOX4 Solanum lycopersicum
  >WOX4|Solanum_lycopersicum

The reference CDS is translated and aligned to the cleaned annotation protein.
The nucleotides are aligned to the cleaned annotation CDS. Both identities
have to reach the given thresholds.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_names import gene_key, id_forms, locus_species, norm  # noqa: E402


STOP_CODONS = {"TAA", "TAG", "TGA"}
DB_PREFIXES = {"lcl", "gb", "ref", "emb", "dbj", "gi", "sp", "tr"}
FIELDS = [
    "query_gene",
    "query_taxon",
    "species_id",
    "new_id",
    "original_id",
    "reference_header",
    "reference_cds_bp",
    "annotation_cds_bp",
    "reference_aa",
    "annotation_aa",
    "protein_identity",
    "cds_identity",
    "status",
    "detail",
]
# A provided reference that cannot be confirmed stops a strict run.
# no_reference is added by the caller when every resolved gene must have one.
FAIL_STATUSES = {
    "cds_mismatch",
    "mismatch",
    "missing_annotation_protein",
    "missing_annotation_cds",
    "unmatched_header",
    "ambiguous_header",
    "duplicate_header",
    "bad_cds",
}


@dataclass
class CrosswalkRow:
    query_gene: str
    query_taxon: str
    species_id: str
    match_method: str
    original_id: str
    gene_key: str
    new_id: str
    annotation_text: str

    @property
    def resolved(self) -> bool:
        return bool(self.new_id) and self.match_method not in {"", "unmatched", "ambiguous"}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crosswalk", required=True, help="gene_id_crosswalk.tsv from select_orthogroups")
    p.add_argument("--proteins", nargs="+", required=True, help="Cleaned protein FASTA files")
    p.add_argument("--cds", nargs="+", required=True, help="Cleaned CDS FASTA files")
    p.add_argument("--reference-cds", required=True, help="Characterized CDS FASTA")
    p.add_argument("--min-protein-identity", type=float, default=98.0)
    p.add_argument("--min-cds-identity", type=float, default=98.0)
    p.add_argument(
        "--require-all",
        action="store_true",
        help="Fail when a resolved candidate has no reference CDS",
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 when a check fails. The report is still written.",
    )
    p.add_argument("--out", required=True, help="Report TSV")
    p.add_argument("--summary", required=True, help="One-line verdict plus detail")
    return p.parse_args(argv)


def taxon_key(text: str) -> str:
    parts = text.replace("|", " ").replace("_", " ").split()
    if len(parts) >= 2:
        return f"{parts[0].lower()}_{parts[1].lower()}"
    if len(parts) == 1:
        return parts[0].lower()
    return ""


def load_crosswalk(path: str) -> list[CrosswalkRow]:
    file = Path(path)
    if not file.is_file():
        raise SystemExit(f"[check_candidate_cds] ERROR: crosswalk not found: {path}")
    with file.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[check_candidate_cds] ERROR: {path} has no candidate rows")
    missing = {"query_gene", "new_id"} - set(rows[0])
    if missing:
        raise SystemExit(
            f"[check_candidate_cds] ERROR: {path} is missing columns: {', '.join(sorted(missing))}"
        )
    loaded = []
    for row in rows:
        loaded.append(
            CrosswalkRow(
                query_gene=(row.get("query_gene") or "").strip(),
                query_taxon=(row.get("query_taxon") or "").strip(),
                species_id=(row.get("species_id") or "").strip(),
                match_method=(row.get("match_method") or "").strip(),
                original_id=(row.get("original_id") or "").strip(),
                gene_key=(row.get("gene_key") or "").strip(),
                new_id=(row.get("new_id") or "").strip(),
                annotation_text=(row.get("annotation_text") or "").strip(),
            )
        )
    return loaded


def load_fasta(paths: list[str]) -> dict[str, str]:
    """Map FASTA ids to sequence text. Empty and non-FASTA files are skipped."""
    seqs: dict[str, str] = {}
    for path in paths:
        file = Path(path)
        if not file.is_file() or file.stat().st_size == 0:
            continue
        with file.open() as fh:
            start = ""
            for line in fh:
                stripped = line.strip()
                if stripped:
                    start = stripped
                    break
            if not start.startswith(">"):
                continue
            fh.seek(0)
            for rec in SeqIO.parse(fh, "fasta"):
                seqs.setdefault(rec.id, str(rec.seq))
    return seqs


def load_reference(path: str) -> list[tuple[str, str, str, str]]:
    """Return (gene, taxon, header, sequence) for each reference record."""
    file = Path(path)
    if not file.is_file() or file.stat().st_size == 0:
        raise SystemExit(f"[check_candidate_cds] ERROR: reference CDS is empty or missing: {path}")
    records = []
    with file.open() as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            header = (rec.description or rec.id).strip()
            gene, taxon, _locus = parse_header_fields(header)
            if not gene:
                raise SystemExit(f"[check_candidate_cds] ERROR: reference CDS has an empty header in {path}")
            records.append((gene, taxon, header, str(rec.seq)))
    if not records:
        raise SystemExit(f"[check_candidate_cds] ERROR: no sequences in {path}")
    return records


def looks_like_locus(token: str) -> bool:
    """True for AT4G35550.1 or Solyc04g078650.4.1, not for an annotation name like TAIR10."""
    token = token.strip().strip(",")
    if not token:
        return False
    if locus_species(token):
        return True
    return gene_key(token) != token


def species_and_locus(rest: str) -> tuple[str, str]:
    """Genus and species, then an optional locus id. A trailing annotation name is ignored."""
    words = rest.replace("_", " ").split()
    if not words:
        return "", ""
    if len(words) == 1:
        return words[0], ""
    # "Ath AT4G35550.1 TAIR10" — the second token is the locus, not a species epithet.
    if looks_like_locus(words[1]):
        taxon = words[0]
        tail = words[1:]
    else:
        taxon = f"{words[0]} {words[1]}"
        tail = words[2:]
    locus = next((token.strip(",") for token in tail if looks_like_locus(token)), "")
    return taxon, locus


def parse_header_fields(text: str) -> tuple[str, str, str]:
    """Return (gene, species, locus).

    The gene and the rest of the header are separated by a tab:

      WOX13<TAB>Arabidopsis thaliana AT4G35550.1 TAIR10

    The species is the genus and species. AT4G35550.1 is a locus id and
    TAIR10 is the annotation name; neither is part of the species. A header
    without a tab still accepts a space or a pipe between the gene and species.
    """
    text = text.strip()
    if not text:
        return "", "", ""
    lowered = text.lower()
    for prefix in DB_PREFIXES:
        token = prefix + "|"
        if lowered.startswith(token):
            text = text.split("|", 1)[1].strip()
            break
    if "\t" not in text and "|" in text:
        parts = [part.strip() for part in text.split("|")]
        if len(parts) == 3 and parts[0] and parts[2]:
            return parts[2], parts[1].replace("_", " "), ""
        if len(parts) == 2 and parts[0] and parts[1]:
            return parts[0], parts[1].replace("_", " "), ""
    if "\t" in text:
        gene, rest = text.split("\t", 1)
    else:
        words = text.split(None, 1)
        gene = words[0]
        rest = words[1] if len(words) > 1 else ""
    gene = gene.strip()
    taxon, locus = species_and_locus(rest.strip())
    return gene, taxon, locus


def parse_header(text: str) -> tuple[str, str]:
    """Split a header into the gene token and the species."""
    gene, taxon, _locus = parse_header_fields(text)
    return gene, taxon


def row_taxon_keys(row: CrosswalkRow) -> set[str]:
    keys = set()
    if row.species_id:
        keys.add(row.species_id.lower())
    if row.query_taxon:
        keys.add(taxon_key(row.query_taxon))
    parts = row.new_id.split("|") if row.new_id else []
    if len(parts) >= 2:
        keys.add(parts[0].lower())
        keys.add(taxon_key(parts[1]))
    return {key for key in keys if key}


def peel_embedded_taxon(gene: str, rows: list[CrosswalkRow]) -> tuple[str, str]:
    """WOX4_Solanum_lycopersicum → (WOX4, Solanum lycopersicum)."""
    lowered = gene.lower()
    keys = []
    for row in rows:
        for key in row_taxon_keys(row):
            if "_" in key:
                keys.append(key)
    for key in sorted(set(keys), key=len, reverse=True):
        suffix = "_" + key
        if lowered.endswith(suffix) and len(gene) > len(suffix):
            return gene[: -len(suffix)], key.replace("_", " ")
    return gene, ""


def gene_matches(row: CrosswalkRow, gene: str) -> bool:
    if norm(gene) == norm(row.query_gene):
        return True
    forms = id_forms(gene)
    if not forms:
        return False
    known = id_forms(row.query_gene) | id_forms(row.original_id) | id_forms(row.gene_key)
    if row.new_id:
        known |= id_forms(row.new_id.split("|")[-1])
    return bool(forms & known)


def taxon_matches(row: CrosswalkRow, taxon: str) -> bool:
    if not taxon:
        return True
    key = taxon_key(taxon)
    if not key:
        return True
    return key in row_taxon_keys(row)


def find_seq(seqs: dict[str, str], new_id: str) -> str:
    if not new_id:
        return ""
    if new_id in seqs:
        return seqs[new_id]
    underscored = new_id.replace("|", "_")
    if underscored in seqs:
        return seqs[underscored]
    last = new_id.split("|")[-1]
    return seqs.get(last, "")


def clean_cds(seq: str) -> tuple[str, str]:
    """Return (coding sequence without a terminal stop, error detail)."""
    text = "".join(seq.split()).upper().replace("U", "T").replace("-", "")
    if not text:
        return "", "sequence is empty"
    if len(text) >= 3 and text[-3:] in STOP_CODONS:
        text = text[:-3]
    odd = [base for base in text if base not in "ACGT"]
    if odd:
        return "", "sequence contains bases other than A, C, G, and T"
    if len(text) % 3 != 0:
        return "", "length is not a multiple of 3 after removing a terminal stop"
    return text, ""


def translate_cds(seq: str) -> str:
    return str(Seq(seq).translate(table=1, to_stop=False)).rstrip("*")


def percent_identity(left: str, right: str, aligner: PairwiseAligner) -> float:
    if not left or not right:
        return 0.0
    if left == right:
        return 100.0
    alignment = aligner.align(left, right)[0]
    length = alignment.length
    if not length:
        return 0.0
    return 100.0 * alignment.counts().identities / length


def blank_row(**values: str) -> dict[str, str]:
    row = {field: "" for field in FIELDS}
    row.update(values)
    return row


def identity_text(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.2f}"


def compare_pair(
    row: CrosswalkRow,
    header: str,
    reference_seq: str,
    proteins: dict[str, str],
    cds_seqs: dict[str, str],
    min_protein: float,
    min_cds: float,
    aligner: PairwiseAligner,
) -> dict[str, str]:
    """Compare one characterized CDS with the annotation sequences for one candidate."""
    base = blank_row(
        query_gene=row.query_gene,
        query_taxon=row.query_taxon,
        species_id=row.species_id,
        new_id=row.new_id,
        original_id=row.original_id,
        reference_header=header,
    )
    if not row.resolved:
        base["status"] = "unresolved"
        base["detail"] = "candidate name was not resolved to an annotation protein"
        return base

    coding, error = clean_cds(reference_seq)
    if error:
        base["status"] = "bad_cds"
        base["detail"] = f"characterized CDS: {error}"
        return base
    protein = translate_cds(coding)
    if "*" in protein:
        base["reference_cds_bp"] = str(len(coding))
        base["reference_aa"] = str(len(protein.replace("*", "")))
        base["status"] = "bad_cds"
        base["detail"] = "characterized CDS has an internal stop codon"
        return base

    annotation_protein = find_seq(proteins, row.new_id).upper().rstrip("*")
    raw_cds = find_seq(cds_seqs, row.new_id)
    annotation_cds = ""
    cds_problem = ""
    if raw_cds:
        annotation_cds, cds_problem = clean_cds(raw_cds)
        if cds_problem:
            annotation_cds = ""

    base["reference_cds_bp"] = str(len(coding))
    base["reference_aa"] = str(len(protein))
    if not annotation_protein:
        base["status"] = "missing_annotation_protein"
        base["detail"] = f"no cleaned protein for {row.new_id}"
        return base

    protein_identity = percent_identity(protein, annotation_protein, aligner)
    base["annotation_aa"] = str(len(annotation_protein))
    base["protein_identity"] = identity_text(protein_identity)
    protein_ok = protein_identity + 1e-9 >= min_protein

    # A wrong protein is a mismatch even when its CDS was dropped at QC.
    # Missing CDS alone, with a protein that already matches, stays
    # missing_annotation_cds.
    if not protein_ok:
        detail = f"annotation protein identity {protein_identity:.2f}% is below {min_protein:g}%"
        if not annotation_cds:
            missing = f"annotation CDS: {cds_problem}" if cds_problem else f"no cleaned CDS for {row.new_id}"
            detail = f"{detail}; {missing}"
        base["status"] = "mismatch"
        base["detail"] = detail
        return base

    if not annotation_cds:
        base["status"] = "missing_annotation_cds"
        base["detail"] = (
            f"annotation CDS: {cds_problem}" if cds_problem else f"no cleaned CDS for {row.new_id}"
        )
        return base

    cds_identity = percent_identity(coding, annotation_cds, aligner)
    base["annotation_cds_bp"] = str(len(annotation_cds))
    base["cds_identity"] = identity_text(cds_identity)

    cds_ok = cds_identity + 1e-9 >= min_cds
    if cds_ok:
        base["status"] = "match"
        base["detail"] = "annotation CDS and protein match the characterized CDS"
    else:
        base["status"] = "cds_mismatch"
        base["detail"] = (
            f"annotation protein matches; CDS identity {cds_identity:.2f}% "
            f"is below {min_cds:g}%"
        )
    return base


def run_check(
    crosswalk: list[CrosswalkRow],
    proteins: dict[str, str],
    cds_seqs: dict[str, str],
    references: list[tuple[str, str, str, str]],
    min_protein: float,
    min_cds: float,
    require_all: bool,
) -> tuple[list[dict[str, str]], list[str]]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    prepared = []
    for gene, taxon, header, seq in references:
        parsed_gene, parsed_taxon, locus = parse_header_fields(header) if header else ("", "", "")
        if parsed_gene:
            gene, taxon = parsed_gene, parsed_taxon
        if not taxon:
            gene, taxon = peel_embedded_taxon(gene, crosswalk)
        symbol_hits = [
            index
            for index, row in enumerate(crosswalk)
            if gene_matches(row, gene) and taxon_matches(row, taxon)
        ]
        locus_hits = []
        if locus and norm(locus) != norm(gene):
            locus_hits = [
                index
                for index, row in enumerate(crosswalk)
                if gene_matches(row, locus) and taxon_matches(row, taxon)
            ]
        if symbol_hits and locus_hits:
            both = [index for index in symbol_hits if index in set(locus_hits)]
            hits = both or symbol_hits + [index for index in locus_hits if index not in symbol_hits]
        else:
            hits = symbol_hits or locus_hits
        prepared.append((gene, taxon, header, seq, hits))

    used: dict[int, str] = {}
    assigned: dict[int, tuple[str, str, str]] = {}
    extras: list[dict[str, str]] = []
    for gene, taxon, header, seq, hits in prepared:
        if not hits:
            extras.append(
                blank_row(
                    query_gene=gene,
                    query_taxon=taxon,
                    reference_header=header,
                    status="unmatched_header",
                    detail="header does not match a gene name or locus id in the candidate list",
                )
            )
            continue
        if len(hits) > 1:
            names = ", ".join(
                f"{crosswalk[index].query_gene} ({crosswalk[index].query_taxon or crosswalk[index].species_id})"
                for index in hits[:6]
            )
            extras.append(
                blank_row(
                    query_gene=gene,
                    query_taxon=taxon,
                    reference_header=header,
                    status="ambiguous_header",
                    detail=f"header matches more than one candidate: {names}",
                )
            )
            continue
        index = hits[0]
        if index in used:
            extras.append(
                blank_row(
                    query_gene=crosswalk[index].query_gene,
                    query_taxon=crosswalk[index].query_taxon,
                    species_id=crosswalk[index].species_id,
                    new_id=crosswalk[index].new_id,
                    original_id=crosswalk[index].original_id,
                    reference_header=header,
                    status="duplicate_header",
                    detail=f"another characterized CDS already matched this candidate ({used[index]})",
                )
            )
            continue
        used[index] = header
        assigned[index] = (header, seq, taxon)

    rows: list[dict[str, str]] = []
    for index, row in enumerate(crosswalk):
        if index not in assigned:
            rows.append(
                blank_row(
                    query_gene=row.query_gene,
                    query_taxon=row.query_taxon,
                    species_id=row.species_id,
                    new_id=row.new_id,
                    original_id=row.original_id,
                    status="no_reference" if row.resolved else "unresolved",
                    detail=(
                        "no characterized CDS was provided for this candidate"
                        if row.resolved
                        else "candidate name was not resolved to an annotation protein"
                    ),
                )
            )
            continue
        header, seq, _taxon = assigned[index]
        rows.append(compare_pair(row, header, seq, proteins, cds_seqs, min_protein, min_cds, aligner))
    rows.extend(extras)

    failures = []
    for row in rows:
        status = row["status"]
        failed = status in FAIL_STATUSES
        if status == "no_reference" and require_all:
            failed = True
        if status == "unresolved" and row["reference_header"]:
            failed = True
        if not failed:
            continue
        label = row["query_gene"] or row["reference_header"] or "(unnamed)"
        taxon = row["query_taxon"] or row["species_id"]
        where = f"{label} ({taxon})" if taxon else label
        numbers = []
        if row["protein_identity"]:
            numbers.append(f"protein {row['protein_identity']}%")
        if row["cds_identity"]:
            numbers.append(f"CDS {row['cds_identity']}%")
        suffix = f" [{', '.join(numbers)}]" if numbers else ""
        failures.append(f"{where}: {status}{suffix}. {row['detail']}")
    return rows, failures


def write_report(path: str, rows: list[dict[str, str]]) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: str, rows: list[dict[str, str]], failures: list[str]) -> None:
    n_match = sum(1 for row in rows if row["status"] == "match")
    if failures:
        lines = [
            "fail",
            f"{len(failures)} of {len(rows)} candidate CDS checks failed ({n_match} matched).",
            *failures[:20],
        ]
    else:
        lines = [
            "pass",
            f"{n_match} of {len(rows)} candidate CDS checks matched the annotation CDS and protein.",
        ]
    Path(path).write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not 0 <= args.min_protein_identity <= 100:
        raise SystemExit("[check_candidate_cds] ERROR: --min-protein-identity must be between 0 and 100")
    if not 0 <= args.min_cds_identity <= 100:
        raise SystemExit("[check_candidate_cds] ERROR: --min-cds-identity must be between 0 and 100")

    crosswalk = load_crosswalk(args.crosswalk)
    proteins = load_fasta(args.proteins)
    cds_seqs = load_fasta(args.cds)
    references = load_reference(args.reference_cds)
    rows, failures = run_check(
        crosswalk,
        proteins,
        cds_seqs,
        references,
        args.min_protein_identity,
        args.min_cds_identity,
        args.require_all,
    )
    write_report(args.out, rows)
    write_summary(args.summary, rows, failures)
    n_match = sum(1 for row in rows if row["status"] == "match")
    print(
        f"[check_candidate_cds] {n_match}/{len(rows)} matched; "
        f"{len(failures)} failed → {args.out}"
    )
    if failures and args.strict:
        print("[check_candidate_cds] ERROR: characterized CDS check failed:", file=sys.stderr)
        for line in failures[:20]:
            print(f"  {line}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
