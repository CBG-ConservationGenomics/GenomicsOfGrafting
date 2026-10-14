#!/usr/bin/env python3
"""Select OrthoFinder orthogroups that contain literature genes of interest.

Input file, one gene per line (tab or whitespace before the species):

  WOX4<TAB>Solanum lycopersicum
  AT1G46480<TAB>Arabidopsis thaliana
  Solyc04g078650 Nicotiana benthamiana

Column 1 is the gene symbol or locus ID of the functionally characterized gene.
Column 2 is genus and species (a short species id such as Sly is also accepted).

That protein is recorded as reference_member and is the sequence later
analyses compare against. Symbols are resolved with GFF Name / Alias /
gene_synonym attributes when --gff files are provided.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_names import (  # noqa: E402
    attach_proteins,
    index_gffs,
    infer_species_id,
    load_queries,
    match_in_species,
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

NAME_ATTRS = (
    "Name",
    "Alias",
    "alias",
    "gene",
    "gene_synonym",
    "synonym",
    "symbol",
    "locus_tag",
    "protein_id",
    "ID",
    "Note",
    "note",
    "description",
    "product",
    "gene_name",
)

# Whole-word symbols inside free-text notes (ITAG often stores WOX4 in Note).
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._-]{2,30}$")

# Literature prefixes sometimes glued onto a symbol (SlWOX4) but absent from the GFF.
SPECIES_PREFIXES = {
    "sly": ("sl",),
    "ath": ("at",),
    "nbe": ("nb",),
    "cpe": ("cp",),
}


@dataclass(frozen=True)
class GoiQuery:
    gene_name: str
    taxon: str
    raw: str
    line_no: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--orthogroups",
        required=True,
        help="OrthoFinder Orthogroups.tsv, or a Results_* directory that also contains unassigned genes",
    )
    p.add_argument("--genes-of-interest", required=True)
    p.add_argument("--id-map", required=True, help="Concatenated QC id map TSV")
    p.add_argument("--gff", nargs="*", default=[], help="Per-species GFF3 files (optional)")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-og-list", required=True)
    p.add_argument(
        "--out-crosswalk",
        default="",
        help="TSV of literature name → proteome id (default: beside --out-tsv)",
    )
    return p.parse_args()


def gene_key(seq_id: str) -> str:
    for rx in LOCUS_RES:
        m = rx.match(seq_id)
        if m:
            return m.group(1)
    return seq_id


def norm(text: str) -> str:
    return re.sub(r"[\s_]+", "", text).lower()


def id_variants(token: str) -> list[str]:
    variants = [token]
    if ":" in token:
        variants.append(token.split(":")[-1])
    return variants


def load_goi(path: str) -> list[GoiQuery]:
    queries: list[GoiQuery] = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            if "\t" in raw:
                gene, taxon = raw.split("\t", 1)
            else:
                parts = raw.split(None, 1)
                if len(parts) != 2:
                    raise SystemExit(
                        f"{path}:{line_no}: expected "
                        "'gene_name<TAB>Genus species' "
                        "(e.g. WOX4\\tSolanum lycopersicum)"
                    )
                gene, taxon = parts
            gene = gene.strip()
            taxon = taxon.strip()
            if not gene or not taxon:
                raise SystemExit(
                    f"{path}:{line_no}: gene name and genus/species are both required"
                )
            queries.append(GoiQuery(gene, taxon, raw, line_no))
    return queries


def binomial_key(taxon: str) -> str:
    cleaned = taxon.replace("_", " ").split()
    if len(cleaned) >= 2 and not re.fullmatch(r"[A-Za-z]{2,6}", taxon):
        return f"{cleaned[0]}_{cleaned[1]}"
    return taxon.replace(" ", "_")


def display_binomial(key: str) -> str:
    parts = key.split("_")
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return key.replace("_", " ")


def load_id_map(path: str) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or "new_id" not in rows[0] or "genus_species" not in rows[0]:
        raise SystemExit(f"{path} is missing new_id / genus_species columns")
    return rows


def resolve_taxon(query: GoiQuery, rows: list[dict]) -> tuple[str, str]:
    """Return (species_id, binomial display) for a query taxon."""
    genus_to_sp: dict[str, str] = {}
    sp_to_genus: dict[str, str] = {}
    for row in rows:
        genus = (row.get("genus_species") or "").strip()
        sid = (row.get("species_id") or "").strip()
        if genus and sid:
            genus_to_sp[genus.lower()] = sid
            sp_to_genus.setdefault(sid.lower(), genus)

    key = binomial_key(query.taxon).lower()
    if key in genus_to_sp:
        sid = genus_to_sp[key]
        original = next(g for g in sp_to_genus.values() if g.lower() == key)
        return sid, display_binomial(original)

    if query.taxon.lower() in sp_to_genus:
        genus = sp_to_genus[query.taxon.lower()]
        return query.taxon, display_binomial(genus)

    known = ", ".join(sorted(display_binomial(g) for g in {r["genus_species"] for r in rows if r.get("genus_species")}))
    raise SystemExit(
        f"genes of interest line {query.line_no}: unknown species '{query.taxon}'. "
        f"Use genus and species from the proteomes ({known})."
    )


def parse_attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in field.split(";"):
        if "=" not in part:
            continue
        key, val = part.split("=", 1)
        out[key.strip()] = unquote(val.strip())
    return out


def extract_names(attrs: dict[str, str]) -> set[str]:
    names: set[str] = set()
    for key in NAME_ATTRS:
        raw = attrs.get(key)
        if not raw:
            continue
        for piece in re.split(r"[,|]", raw):
            piece = piece.strip()
            if not piece:
                continue
            names.add(piece)
            if ":" in piece:
                names.add(piece.split(":")[-1])
            for token in re.split(r"[\s,;|/()\[\]:]+", piece):
                if SYMBOL_RE.match(token) and (any(c.isdigit() for c in token) or token.isupper()):
                    names.add(token)
    return {n for n in names if len(n) >= 2}


def load_alias_index(paths: list[str]) -> dict[str, set[str]]:
    """Map a normalized locus token to names that sit on that feature or its parents."""
    by_token: dict[str, set[str]] = defaultdict(set)
    for path in paths:
        file = Path(path)
        if not file.is_file() or file.stat().st_size == 0:
            continue
        features: dict[str, dict] = {}
        with file.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attrs = parse_attrs(parts[8])
                fid = attrs.get("ID") or attrs.get("protein_id") or attrs.get("Name")
                if not fid:
                    continue
                parents = [p.strip() for p in attrs.get("Parent", "").split(",") if p.strip()]
                features[fid] = {"parents": parents, "names": extract_names(attrs)}

        memo: dict[str, set[str]] = {}

        def collect(fid: str, stack: set[str] | None = None) -> set[str]:
            if fid in memo:
                return memo[fid]
            if fid not in features:
                return set()
            stack = stack or set()
            if fid in stack:
                return set()
            stack.add(fid)
            names = set(features[fid]["names"])
            for parent in features[fid]["parents"]:
                names |= collect(parent, stack)
            stack.remove(fid)
            memo[fid] = names
            return names

        for fid in features:
            names = collect(fid)
            tokens: set[str] = set()
            for variant in id_variants(fid):
                tokens.add(norm(variant))
            for name in names:
                tokens.add(norm(name))
            normalized_names = {norm(n) for n in names if len(norm(n)) >= 2}
            for tok in tokens:
                if len(tok) < 2:
                    continue
                by_token[tok].update(normalized_names)
                by_token[tok].add(tok)
    return by_token


def id_aliases(token: str) -> set[str]:
    """Forms OrthoFinder may write for one proteome accession."""
    token = token.strip()
    if not token:
        return set()
    aliases = {token, token.replace("|", "_")}
    if "|" in token:
        aliases.add(token.split("|")[-1])
    return aliases


def load_sequence_ids(root: Path) -> dict[str, str]:
    """Map OrthoFinder internal IDs (0_42) back to the original FASTA accession."""
    mapping: dict[str, str] = {}
    if not root.is_dir():
        return mapping
    for path in root.glob("**/SequenceIDs.txt"):
        for line in path.read_text(errors="replace").splitlines():
            if ":" not in line:
                continue
            internal, original = line.split(":", 1)
            original = original.strip().split()[0].rstrip(",")
            if internal.strip() and original:
                mapping[internal.strip()] = original
    return mapping


def discover_orthogroup_tables(path: Path) -> list[Path]:
    """Orthogroups.tsv plus unassigned genes, from a file or a Results directory."""
    if path.is_file():
        tables = [path]
        extra = path.with_name("Orthogroups_UnassignedGenes.tsv")
        if extra.is_file():
            tables.append(extra)
        return tables
    tables: list[Path] = []
    for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
        tables.extend(p for p in sorted(path.glob(f"**/{name}")) if p.is_file())
    if not any(p.name == "Orthogroups.tsv" for p in tables):
        raise SystemExit(f"[select_orthogroups] ERROR: Orthogroups.tsv not found under {path}")
    return tables


def translate_member(token: str, seq_ids: dict[str, str]) -> str:
    token = token.strip()
    return seq_ids.get(token, token)


def member_lookup(members: list[str]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for member in members:
        for alias in id_aliases(member):
            lookup.setdefault(alias, member)
    return lookup


def member_hit(id_row: dict, lookup: dict[str, str]) -> bool:
    for token in (id_row["new_id"], id_row["original_id"], id_row.get("gene_key") or ""):
        if any(alias in lookup for alias in id_aliases(token)):
            return True
    return False


def names_for_row(row: dict, alias_index: dict[str, set[str]]) -> set[str]:
    found = {norm(row["original_id"]), norm(row.get("gene_key") or gene_key(row["original_id"]))}
    for key in (row["original_id"], row.get("gene_key") or ""):
        if key:
            found |= alias_index.get(norm(key), set())
    return {n for n in found if n}


def strip_species_prefix(gene_name: str, species_id: str) -> str | None:
    prefixes = SPECIES_PREFIXES.get(species_id.lower(), ())
    lowered = gene_name.lower()
    for prefix in prefixes:
        if lowered.startswith(prefix) and len(gene_name) - len(prefix) >= 3:
            return gene_name[len(prefix) :]
    return None


def resolve_query(
    query: GoiQuery,
    species_id: str,
    rows: list[dict],
    alias_index: dict[str, set[str]],
) -> dict | None:
    qn = norm(query.gene_name)
    stripped = strip_species_prefix(query.gene_name, species_id)
    stripped_n = norm(stripped) if stripped else ""
    hits: list[tuple[int, int, dict]] = []
    for row in rows:
        if row.get("species_id") != species_id:
            continue
        names = names_for_row(row, alias_index)
        gene_key_n = norm(row.get("gene_key") or gene_key(row["original_id"]))
        if qn == norm(row["original_id"]) or qn == gene_key_n:
            rank = 0
        elif qn in names:
            rank = 1
        elif stripped_n and stripped_n in names and stripped_n != qn:
            rank = 2
        else:
            continue
        length = int(row["length"]) if str(row.get("length") or "").isdigit() else 0
        hits.append((rank, -length, row))
    if not hits:
        return None
    hits.sort(key=lambda item: (item[0], item[1], item[2]["new_id"]))
    best_rank = hits[0][0]
    top = [item for item in hits if item[0] == best_rank]
    if len(top) > 1:
        ids = ", ".join(item[2]["new_id"] for item in top[:6])
        print(
            f"[select_orthogroups] WARNING: {query.gene_name} matches {len(top)} "
            f"proteins in {species_id}; using {top[0][2]['new_id']}. Also: {ids}"
        )
    elif best_rank == 2:
        print(
            f"[select_orthogroups] {query.gene_name} matched as {stripped} "
            f"in {species_id} ({top[0][2]['new_id']})"
        )
    return top[0][2]


def binomial_for(species_id: str, rows: list[dict]) -> str:
    for row in rows:
        if (row.get("species_id") or "").lower() == species_id.lower() and row.get("genus_species"):
            return display_binomial(row["genus_species"])
    return species_id


def resolve_literature(queries, id_rows, gff_paths):
    """Match each literature name to one protein. Returns resolutions and the crosswalk rows."""
    index = index_gffs(gff_paths)
    views = attach_proteins(id_rows, index)
    species_ids = sorted({row["species_id"] for row in id_rows if row.get("species_id")})
    resolved = []
    missing = []
    crosswalk = []

    for query in queries:
        detail = ""
        method = ""
        view = None
        species_id = ""
        binomial = ""
        if query.taxon:
            species_id, binomial = resolve_taxon(query, id_rows)
            method, view, detail = match_in_species(query, species_id, views)
        else:
            inferred = infer_species_id(query)
            if inferred:
                species_id = inferred
                binomial = binomial_for(species_id, id_rows)
                method, view, detail = match_in_species(query, species_id, views)
            else:
                hits = []
                for sid in species_ids:
                    hit_method, hit_view, hit_detail = match_in_species(query, sid, views)
                    if hit_view is not None:
                        hits.append((sid, hit_method, hit_view, hit_detail))
                if len(hits) == 1:
                    species_id, method, view, detail = hits[0]
                    binomial = binomial_for(species_id, id_rows)
                elif len(hits) > 1:
                    method = "ambiguous"
                    detail = "; ".join(
                        f"{sid} {hit_view.row.get('original_id')} via {hit_method}"
                        for sid, hit_method, hit_view, _hit_detail in hits[:8]
                    )
        row = view.row if view is not None else None
        if row is None:
            missing.append(query)
            print(
                f"[select_orthogroups] {query.gene_name}: unmatched"
                + (f" ({method}: {detail})" if detail else " (no locus, symbol, or description)")
            )
        else:
            resolved.append((query, species_id, binomial, row))
            print(
                f"[select_orthogroups] {query.gene_name} -> {row.get('original_id')} "
                f"via {method}" + (f" [{detail}]" if method == "description" else "")
            )
        crosswalk.append(
            {
                "query_gene": query.gene_name,
                "query_taxon": query.taxon or binomial,
                "species_id": species_id,
                "match_method": method or "unmatched",
                "original_id": (row or {}).get("original_id", ""),
                "gene_key": (row or {}).get("gene_key", ""),
                "new_id": (row or {}).get("new_id", ""),
                "annotation_text": detail if method in {"description", "ambiguous"} else "",
            }
        )
    return resolved, missing, crosswalk


def main() -> None:
    args = parse_args()
    queries = load_queries(args.genes_of_interest)
    if not queries:
        raise SystemExit("No genes of interest found in input file")

    id_rows = load_id_map(args.id_map)
    resolved, missing, crosswalk = resolve_literature(queries, id_rows, args.gff)

    rows_out = []
    selected_ogs = []
    found_ids: set[str] = set()
    og_root = Path(args.orthogroups)
    seq_ids = load_sequence_ids(og_root)
    tables = discover_orthogroup_tables(og_root)
    if any(p.name.startswith("Orthogroups_Unassigned") for p in tables):
        print("[select_orthogroups] also searching Orthogroups_UnassignedGenes.tsv")

    seen_ogs: set[str] = set()
    for table in tables:
        unassigned = table.name.startswith("Orthogroups_Unassigned")
        with table.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            species_cols = [c for c in (reader.fieldnames or []) if c != "Orthogroup"]
            for row in reader:
                og = (row.get("Orthogroup") or "").strip()
                if not og:
                    continue
                if unassigned and not og.startswith("UNASSIGNED_"):
                    og = f"UNASSIGNED_{og}"
                if og in seen_ogs:
                    continue
                members: list[str] = []
                for col in species_cols:
                    cell = row.get(col) or ""
                    if cell.strip():
                        members.extend(translate_member(m, seq_ids) for m in cell.split(",") if m.strip())
                lookup = member_lookup(members)

                hits = []
                for query, sid, binomial, id_row in resolved:
                    if member_hit(id_row, lookup):
                        hits.append((query, sid, binomial, id_row["new_id"]))
                if not hits:
                    continue

                seen_ogs.add(og)
                selected_ogs.append(og)
                # File order is preserved in `resolved`; the first hit is the reference.
                ref_query, _ref_sid, ref_binomial, ref_member = hits[0]
                found_ids.add(ref_member)
                for _q, _sid, _binomial, new_id in hits:
                    found_ids.add(new_id)

                matched_ids: list[str] = []
                for _q, _s, _b, new_id in hits:
                    if new_id not in matched_ids:
                        matched_ids.append(new_id)
                alias_to_new: dict[str, str] = {}
                for new_id in matched_ids:
                    for alias in id_aliases(new_id):
                        alias_to_new[alias] = new_id
                canon_members: list[str] = []
                seen_members: set[str] = set()
                for member in members:
                    canon = member
                    for alias in id_aliases(member):
                        if alias in alias_to_new:
                            canon = alias_to_new[alias]
                            break
                    if canon not in seen_members:
                        canon_members.append(canon)
                        seen_members.add(canon)
                for new_id in matched_ids:
                    if new_id not in seen_members:
                        canon_members.append(new_id)
                        seen_members.add(new_id)
                members = canon_members
                rows_out.append(
                    {
                        "Orthogroup": og,
                        "genes_of_interest": ";".join(q.gene_name for q, _s, _b, _n in hits),
                        "query_species_ids": ";".join(sid for _q, sid, _b, _n in hits),
                        "reference_gene": ref_query.gene_name,
                        "reference_species": ref_binomial,
                        "reference_member": ref_member,
                        "matched_members": ";".join(matched_ids),
                        "n_members": len(members),
                        "n_species_with_members": sum(1 for col in species_cols if (row.get(col) or "").strip()),
                        "members": ",".join(members),
                    }
                )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Orthogroup",
        "genes_of_interest",
        "query_species_ids",
        "reference_gene",
        "reference_species",
        "reference_member",
        "matched_members",
        "n_members",
        "n_species_with_members",
        "members",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    with open(args.out_og_list, "w") as fh:
        for og in selected_ogs:
            fh.write(og + "\n")

    crosswalk_path = Path(args.out_crosswalk) if args.out_crosswalk else Path(args.out_tsv).with_name("gene_id_crosswalk.tsv")
    crosswalk_path.parent.mkdir(parents=True, exist_ok=True)
    with crosswalk_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "query_gene",
                "query_taxon",
                "species_id",
                "match_method",
                "original_id",
                "gene_key",
                "new_id",
                "annotation_text",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(crosswalk)
    print(f"[select_orthogroups] name crosswalk -> {crosswalk_path}")

    not_in_og = [q for q, _s, _b, id_row in resolved if id_row["new_id"] not in found_ids]
    print(
        f"[select_orthogroups] selected {len(selected_ogs)} orthogroups; "
        f"{len(queries) - len(missing) - len(not_in_og)}/{len(queries)} literature genes found"
    )
    if missing:
        print(
            "[select_orthogroups] WARNING: not in the proteome annotation "
            "(symbol or locus ID): "
            + ", ".join(f"{q.gene_name} ({q.taxon})" for q in missing[:20])
        )
        print(
            "[select_orthogroups] Symbols match a GFF gene symbol, or a description "
            "that names the same family and number (WOX4 matches "
            "'wuschel related homeobox 4'). Locus IDs such as Solyc04g078650 "
            "or AT1G46480 match with or without an isoform suffix. "
            "See gene_id_crosswalk.tsv for the name used."
        )
    if not_in_og:
        print(
            "[select_orthogroups] WARNING: matched a protein that is not in any orthogroup: "
            + ", ".join(q.gene_name for q in not_in_og[:20])
        )
    if not rows_out:
        unmatched = ", ".join(f"{q.gene_name} ({q.taxon})" for q in (missing + not_in_og)[:20])
        raise SystemExit(
            "[select_orthogroups] ERROR: no orthogroups contain the literature genes, "
            "so annotation was not started. Unmatched: "
            + (unmatched or "(none resolved)")
            + ". Use a locus ID or a GFF gene symbol, then genus and species "
            "(for example: WOX4<TAB>Solanum lycopersicum)."
        )


if __name__ == "__main__":
    main()
