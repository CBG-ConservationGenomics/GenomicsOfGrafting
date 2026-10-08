#!/usr/bin/env python3
"""Select OrthoFinder orthogroups that contain genes of interest.

Genes-of-interest file formats (mix freely):

  # bare ID — species inferred from ID pattern when possible
  Solyc01g000570
  AT1G01010
  Niben101Scf01786g00001

  # explicit species (recommended for ambiguous IDs like XP_*)
  Ath:AT1G01010
  Nbe:Niben101Scf01786g00001
  Cpe:XP_023545678.1

  # tab-separated: gene_id <tab> species_id
  AT1G01010    Ath
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


# Locus-style prefixes used to collapse isoforms → gene keys
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

ISOFORM_SUFFIX = re.compile(
    r"(?:\.\d+|_iso\d+|-T\d+|\.t\d+|_t\d+)$", re.IGNORECASE
)

# Infer pipeline species id from a bare gene identifier
INFER_SPECIES = [
    (re.compile(r"^Solyc\d+g\d+", re.I), "Sly"),
    (re.compile(r"^AT[1-5MC]G\d+", re.I), "Ath"),
    (re.compile(r"^Niben", re.I), "Nbe"),
    (re.compile(r"^Nb[vVlL\d]", re.I), "Nbe"),
    (re.compile(r"^Cp[\d.]*LG", re.I), "Cpe"),
    (re.compile(r"^Cpe\|", re.I), "Cpe"),
    (re.compile(r"^Sly\|", re.I), "Sly"),
    (re.compile(r"^Ath\|", re.I), "Ath"),
    (re.compile(r"^Nbe\|", re.I), "Nbe"),
]


@dataclass(frozen=True)
class GoiQuery:
    gene_id: str
    species_id: str | None  # Sly / Ath / Nbe / Cpe or None
    raw: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--orthogroups", required=True, help="OrthoFinder Orthogroups.tsv")
    p.add_argument("--genes-of-interest", required=True)
    p.add_argument("--id-map", required=True, help="Concatenated QC id map TSV")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-og-list", required=True)
    return p.parse_args()


def gene_key(seq_id: str) -> str:
    for rx in LOCUS_RES:
        m = rx.match(seq_id)
        if m:
            return m.group(1)
    return ISOFORM_SUFFIX.sub("", seq_id)


def infer_species_id(gene_id: str) -> str | None:
    for rx, sp in INFER_SPECIES:
        if rx.search(gene_id):
            return sp
    return None


def load_goi(path: str) -> list[GoiQuery]:
    queries: list[GoiQuery] = []
    with open(path) as fh:
        for line in fh:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue

            species_id = None
            gene_id = None

            # tab / comma separated: gene_id, species_id
            if "\t" in raw or ("," in raw and not raw.upper().startswith(("XP_", "NP_", "WP_", "YP_"))):
                parts = re.split(r"[\t,]+", raw)
                parts = [p.strip() for p in parts if p.strip()]
                if len(parts) >= 2 and re.fullmatch(r"[A-Za-z]{2,6}", parts[1]):
                    gene_id, species_id = parts[0], parts[1]
                elif len(parts) >= 2 and re.fullmatch(r"[A-Za-z]{2,6}", parts[0]):
                    species_id, gene_id = parts[0], parts[1]
                else:
                    gene_id = parts[0]
            # species:gene or species|gene (species code is short)
            elif re.match(r"^[A-Za-z]{2,6}[:|]", raw):
                species_id, gene_id = re.split(r"[:|]", raw, maxsplit=1)
                species_id = species_id.strip()
                gene_id = gene_id.strip()
            else:
                gene_id = raw.split()[0]

            if not gene_id:
                continue
            if species_id is None:
                species_id = infer_species_id(gene_id)
            queries.append(GoiQuery(gene_id=gene_id, species_id=species_id, raw=raw))
    return queries


def member_parts(token: str) -> tuple[str | None, str]:
    """Return (species_id, original_id) from Sly|Genus_species|orig or bare id."""
    parts = token.split("|")
    if len(parts) >= 3:
        return parts[0], parts[-1]
    if len(parts) == 2:
        return parts[0], parts[-1]
    return None, token


def ids_equivalent(query: str, candidate: str) -> bool:
    """Flexible equality between a GOI string and a proteome original id / gene key."""
    q = query.strip()
    c = candidate.strip()
    if not q or not c:
        return False
    if q == c or q.lower() == c.lower():
        return True

    qk = gene_key(q)
    ck = gene_key(c)
    if qk.lower() == ck.lower():
        return True

    # prefix: Solyc01g000570 vs Solyc01g000570.2.1
    if c.lower().startswith(q.lower()) or q.lower().startswith(c.lower()):
        # avoid tiny accidental prefixes
        shorter = min(len(q), len(c))
        if shorter >= 6:
            return True

    return False


def match_member(
    token: str,
    query: GoiQuery,
    gene_key_to_species: dict[str, set[str]],
) -> bool:
    sp, orig = member_parts(token)
    if query.species_id and sp and sp != query.species_id:
        return False

    if ids_equivalent(query.gene_id, token) or ids_equivalent(query.gene_id, orig):
        return True
    if ids_equivalent(query.gene_id, gene_key(orig)):
        return True

    # If species was inferred/specified, also require gene_key known for that species when available
    if query.species_id:
        gks = gene_key_to_species.get(gene_key(query.gene_id).lower(), set())
        if gks and query.species_id not in gks and sp != query.species_id:
            return False

    return False


def main() -> None:
    args = parse_args()
    queries = load_goi(args.genes_of_interest)
    if not queries:
        raise SystemExit("No genes of interest found in input file")

    gene_key_to_species: dict[str, set[str]] = {}
    with open(args.id_map) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gk = (row.get("gene_key") or gene_key(row["original_id"])).lower()
            gene_key_to_species.setdefault(gk, set()).add(row["species_id"])

    rows_out = []
    selected_ogs = []
    found_raw: set[str] = set()

    with open(args.orthogroups) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        species_cols = [c for c in reader.fieldnames if c != "Orthogroup"]
        for row in reader:
            og = row["Orthogroup"]
            members = []
            for col in species_cols:
                cell = row.get(col) or ""
                if not cell.strip():
                    continue
                members.extend([m.strip() for m in cell.split(",") if m.strip()])

            hits: list[tuple[GoiQuery, str]] = []
            for q in queries:
                for mem in members:
                    if match_member(mem, q, gene_key_to_species):
                        hits.append((q, mem))
                        found_raw.add(q.raw)
                        break

            if not hits:
                continue

            selected_ogs.append(og)
            hit_genes = []
            hit_members = []
            hit_species = []
            seen = set()
            for q, mem in hits:
                key = (q.gene_id, mem)
                if key in seen:
                    continue
                seen.add(key)
                hit_genes.append(q.gene_id)
                hit_members.append(mem)
                sp = q.species_id or member_parts(mem)[0] or ""
                hit_species.append(sp)

            rows_out.append(
                {
                    "Orthogroup": og,
                    "genes_of_interest": ";".join(hit_genes),
                    "query_species_ids": ";".join(hit_species),
                    "matched_members": ";".join(hit_members),
                    "n_members": len(members),
                    "n_species_with_members": sum(
                        1 for col in species_cols if (row.get(col) or "").strip()
                    ),
                    "members": ",".join(members),
                }
            )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        fieldnames = [
            "Orthogroup",
            "genes_of_interest",
            "query_species_ids",
            "matched_members",
            "n_members",
            "n_species_with_members",
            "members",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    with open(args.out_og_list, "w") as fh:
        for og in selected_ogs:
            fh.write(og + "\n")

    missing = [q.raw for q in queries if q.raw not in found_raw]
    inferred = sum(1 for q in queries if q.species_id)
    print(
        f"[select_orthogroups] selected {len(selected_ogs)} orthogroups; "
        f"{len(queries) - len(missing)}/{len(queries)} GOIs found "
        f"({inferred} with species set/inferred)"
    )
    if missing:
        print(
            f"[select_orthogroups] WARNING: {len(missing)} GOIs not found: "
            f"{', '.join(missing[:20])}"
        )
        print(
            "[select_orthogroups] Tip: for ambiguous IDs use Species:gene "
            "(e.g. Cpe:XP_023545678.1 or Ath:AT1G01010)"
        )


if __name__ == "__main__":
    main()
