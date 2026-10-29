#!/usr/bin/env python3
"""Normalize literature gene names onto proteome and GFF identifiers.

Annotations and papers use different names for the same gene:

  SlWOX4                         literature prefix plus symbol
  Solyc04g078650                 locus, no isoform
  Solyc04g078650.4.1             protein / mRNA id
  wuschel related homobox 4      ITAG description, symbol not stored

Locus IDs match after isoform suffixes are removed. A symbol matches when
the annotation stores that symbol, or when the description names the same
family and the same number (WOX4 matches "wuschel related homobox 4").
A trailing "-like" is not expanded to a different gene number.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote


LOCUS_RES = [
    re.compile(r"^(Solyc\d+g\d+)", re.IGNORECASE),
    re.compile(r"^(AT[1-5MC]G\d+)", re.IGNORECASE),
    re.compile(r"^(Niben\d+Scf\w+g\d+)", re.IGNORECASE),
    re.compile(r"^(Niben\w+g\d+)", re.IGNORECASE),
    re.compile(r"^(Nb\w+g\d+)", re.IGNORECASE),
    re.compile(r"^(Cp[\w.]+g\d+)", re.IGNORECASE),
    re.compile(r"^(LOC\d+)", re.IGNORECASE),
    re.compile(r"^((?:XP|NP|WP|YP)_\d+)", re.IGNORECASE),
    re.compile(r"^(VIT_\w+)", re.IGNORECASE),
    re.compile(r"^(CA\d+g\d+)", re.IGNORECASE),
]

# Species letters glued onto a symbol in papers (SlWOX4), not part of the annotation.
SPECIES_PREFIXES = {
    "sly": ("sl",),
    "ath": ("at",),
    "nbe": ("nb",),
    "cpe": ("cp",),
}
PREFIX_TO_SPECIES = {prefix: sid for sid, prefixes in SPECIES_PREFIXES.items() for prefix in prefixes}

# Expanded names used in place of the symbol. The number still has to match,
# so WOX4 does not land on "WUSCHEL-related homeobox 13".
FAMILY_SYNONYMS = {
    "wox": ("wuschel",),
    "wrky": ("wrky",),
    "erf": ("ethylene-responsive transcription factor", "ethylene responsive transcription factor"),
    "nac": ("nac domain",),
    "anac": ("nac domain", "anac"),
    "dof": ("dof",),
    "xth": ("xyloglucan",),
    "pal": ("phenylalanine ammonia",),
    "wind": ("wound induced dedifferentiation",),
    "apx": ("ascorbate peroxidase",),
    "cat": ("catalase",),
    "iaa": ("aux/iaa",),
    "gh9": ("endoglucanase", "glucanase"),
}

CLEAN_SYMBOL = re.compile(r"^([A-Za-z]{3,})(\d+(?:\.\d+)?)$")
LIKE_SUFFIX = re.compile(r"[-_ ]?like$", re.IGNORECASE)
ID_ATTRS = ("ID", "Name", "Alias", "alias", "gene", "gene_synonym", "synonym", "symbol", "locus_tag", "protein_id")
DESC_ATTRS = ("Note", "note", "description", "product")
FEATURE_TYPES = {"gene", "mrna", "cds", "protein"}


@dataclass
class NameQuery:
    gene_name: str
    taxon: str
    raw: str
    line_no: int
    pieces: tuple[str, ...]
    like: bool


@dataclass
class Bundle:
    ids: set[str] = field(default_factory=set)
    symbols: set[str] = field(default_factory=set)
    descriptions: list[str] = field(default_factory=list)
    # locus_tag / protein accession. BAC synonyms such as F8D20.60 are not primary.
    primary: set[str] = field(default_factory=set)


@dataclass
class ProteinView:
    row: dict
    ids: set[str]
    symbols: set[str]
    descriptions: list[str]
    primary: set[str] = field(default_factory=set)


@dataclass
class Resolution:
    query: NameQuery
    species_id: str
    binomial: str
    row: dict | None
    method: str
    detail: str


def norm(text: str) -> str:
    return re.sub(r"[\s_]+", "", text).lower()


def gene_key(seq_id: str) -> str:
    for rx in LOCUS_RES:
        match = rx.match(seq_id)
        if match:
            return match.group(1)
    return seq_id


def id_forms(token: str) -> set[str]:
    """Locus and protein ids, with type prefixes and isoform suffixes removed.

    Trailing ``.1`` / ``.4.1`` is an isoform or accession version only on a
    recognized locus (``AT4G35550.1``, ``Solyc04g078650.4.1``, ``NP_195280.1``).
    An Arabidopsis BAC gene name such as ``F8D20.60`` stays whole. Stripping
    ``.60`` would leave the clone name ``F8D20``, which every gene on that
    clone shares, and the longest of those proteins would be selected.
    """
    forms: set[str] = set()
    raw = token.strip()
    if not raw:
        return forms
    variants = [raw]
    if ":" in raw:
        variants.append(raw.split(":")[-1])
    if raw.startswith(("cds-", "rna-", "gene-")):
        variants.append(raw.split("-", 1)[1])
    for variant in variants:
        current = variant
        forms.add(norm(current))
        forms.add(norm(gene_key(current)))
        if gene_key(current) == current and not locus_species(current):
            continue
        while re.search(r"\.\d+$", current):
            current = re.sub(r"\.\d+$", "", current)
            forms.add(norm(current))
            forms.add(norm(gene_key(current)))
    return {form for form in forms if len(form) >= 3}


def locus_species(token: str) -> str | None:
    lowered = token.lower()
    if lowered.startswith("solyc"):
        return "sly"
    if re.match(r"at[1-5mc]g\d+", lowered):
        return "ath"
    if lowered.startswith("niben") or lowered.startswith("nb"):
        return "nbe"
    if lowered.startswith("cp") and re.match(r"cp[\w.]+g\d+", lowered):
        return "cpe"
    return None


def split_pieces(gene_name: str) -> tuple[tuple[str, ...], bool]:
    like = bool(LIKE_SUFFIX.search(gene_name))
    text = LIKE_SUFFIX.sub("", gene_name).strip()
    parts = re.split(r"[;/]|[()]", text)
    pieces = tuple(part.strip() for part in parts if part.strip())
    return pieces or (gene_name.strip(),), like


def load_queries(path: str) -> list[NameQuery]:
    queries: list[NameQuery] = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            if "\t" in raw:
                gene, taxon = raw.split("\t", 1)
                taxon = taxon.strip()
            else:
                parts = raw.split(None, 1)
                if len(parts) == 2 and " " in parts[1]:
                    gene, taxon = parts
                else:
                    gene, taxon = parts[0], ""
            gene = gene.strip()
            if not gene:
                raise SystemExit(f"{path}:{line_no}: gene name is empty")
            pieces, like = split_pieces(gene)
            queries.append(NameQuery(gene, taxon.strip(), raw, line_no, pieces, like))
    return queries


def strip_prefix(piece: str, species_id: str) -> str | None:
    if gene_key(piece) != piece and locus_species(piece):
        return None
    prefixes = SPECIES_PREFIXES.get(species_id.lower(), ())
    lowered = piece.lower()
    for prefix in prefixes:
        if lowered.startswith(prefix) and len(piece) - len(prefix) >= 3:
            rest = piece[len(prefix) :]
            if rest[0].isalpha() or rest[0].isdigit():
                return rest
    return None


def symbol_parts(piece: str) -> tuple[str, str] | None:
    match = CLEAN_SYMBOL.fullmatch(piece)
    if not match:
        return None
    return match.group(1).lower(), match.group(2)


def parse_attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in field.split(";"):
        if "=" not in part:
            continue
        key, val = part.split("=", 1)
        out[key.strip()] = unquote(val.strip())
    return out


def index_gffs(paths: list[str]) -> dict[str, Bundle]:
    """Map every locus and protein id to the symbols and descriptions on that gene."""
    features: dict[str, dict] = {}
    for path in paths:
        file = Path(path)
        if not file.is_file() or file.stat().st_size == 0:
            continue
        with file.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2].lower() not in FEATURE_TYPES:
                    continue
                attrs = parse_attrs(parts[8])
                fid = attrs.get("ID") or attrs.get("protein_id") or attrs.get("Name")
                if not fid:
                    continue
                ids: set[str] = set()
                symbols: set[str] = set()
                for key in ID_ATTRS:
                    raw = attrs.get(key)
                    if not raw:
                        continue
                    for piece in re.split(r"[,|]", raw):
                        piece = piece.strip()
                        if not piece:
                            continue
                        ids |= id_forms(piece)
                        if ":" in piece:
                            ids |= id_forms(piece.split(":")[-1])
                        if not locus_species(piece) and gene_key(piece) == piece and re.fullmatch(r"[A-Za-z][A-Za-z0-9._-]{2,30}", piece):
                            symbols.add(norm(piece))
                descriptions = []
                for key in DESC_ATTRS:
                    raw = attrs.get(key)
                    if raw:
                        descriptions.append(raw)
                # The AGI in locus_tag is the gene. A BAC name in gene_synonym is not.
                primary: set[str] = set()
                for key in ("locus_tag", "protein_id", "ID", "Name"):
                    raw = attrs.get(key) or ""
                    for piece in re.split(r"[,|]", raw):
                        piece = piece.strip()
                        if not piece:
                            continue
                        owned = (
                            key == "locus_tag"
                            or bool(locus_species(piece))
                            or gene_key(piece) != piece
                            or piece.startswith(("cds-", "rna-", "gene-"))
                            or ":" in piece
                        )
                        if owned:
                            primary |= id_forms(piece)
                dbxref = attrs.get("Dbxref") or ""
                for piece in dbxref.split(","):
                    if ":" in piece:
                        ids |= id_forms(piece.split(":")[-1])
                parents = [parent.strip() for parent in attrs.get("Parent", "").split(",") if parent.strip()]
                features[fid] = {
                    "parents": parents,
                    "ids": ids,
                    "symbols": symbols,
                    "descriptions": descriptions,
                    "primary": primary,
                }

    children: dict[str, list[str]] = defaultdict(list)
    for fid, feat in features.items():
        for parent in feat["parents"]:
            children[parent].append(fid)

    memo: dict[str, Bundle] = {}

    def bundle_for(fid: str, stack: set[str] | None = None) -> Bundle:
        if fid in memo:
            return memo[fid]
        feat = features.get(fid)
        if feat is None:
            return Bundle()
        stack = stack or set()
        if fid in stack:
            return Bundle()
        stack.add(fid)
        bundle = Bundle(
            set(feat["ids"]),
            set(feat["symbols"]),
            list(feat["descriptions"]),
            set(feat["primary"]),
        )
        for parent in feat["parents"]:
            parent_bundle = bundle_for(parent, stack)
            bundle.ids |= parent_bundle.ids
            bundle.symbols |= parent_bundle.symbols
            bundle.descriptions.extend(parent_bundle.descriptions)
            bundle.primary |= parent_bundle.primary
        for child in children.get(fid, []):
            child_bundle = bundle_for(child, stack)
            bundle.ids |= child_bundle.ids
            bundle.symbols |= child_bundle.symbols
            bundle.descriptions.extend(child_bundle.descriptions)
            bundle.primary |= child_bundle.primary
        stack.remove(fid)
        # Descriptions are repeated up and down the model; keep one copy.
        seen: set[str] = set()
        unique = []
        for text in bundle.descriptions:
            if text not in seen:
                seen.add(text)
                unique.append(text)
        bundle.descriptions = unique
        memo[fid] = bundle
        return bundle

    # The first feature to claim a token keeps it. Copying every identifier
    # across a shared token merges neighboring genes: a BAC clone name, or a
    # Note that cites another TAIR locus, would then make the longest protein
    # on that clone win.
    index: dict[str, Bundle] = {}
    for fid in features:
        bundle = bundle_for(fid)
        for key in bundle.ids:
            index.setdefault(key, bundle)
    return index


def attach_proteins(rows: list[dict], index: dict[str, Bundle]) -> list[ProteinView]:
    views = []
    for row in rows:
        ids = id_forms(row.get("original_id") or "")
        ids |= id_forms(row.get("gene_key") or "")
        primary = set(ids)
        symbols: set[str] = set()
        descriptions: list[str] = []
        seen_desc: set[str] = set()
        for key in list(ids):
            bundle = index.get(key)
            if bundle is None:
                continue
            ids |= bundle.ids
            primary |= bundle.primary
            symbols |= bundle.symbols
            for text in bundle.descriptions:
                if text not in seen_desc:
                    seen_desc.add(text)
                    descriptions.append(text)
        views.append(ProteinView(row, ids, symbols, descriptions, primary))
    return views


def _number_in(number: str, text: str) -> bool:
    bare = number.lstrip("0") or "0"
    return re.search(rf"(?<![A-Za-z0-9.])0*{re.escape(bare)}(?![A-Za-z0-9.])", text, re.IGNORECASE) is not None


def _family_in(stem: str, text: str) -> bool:
    low = text.lower()
    for synonym in FAMILY_SYNONYMS.get(stem, (stem,)):
        if " " in synonym or "-" in synonym or "/" in synonym:
            if synonym in low:
                return True
        elif re.search(rf"\b{re.escape(synonym)}\b", low):
            return True
    return False


def _query_ids(query: NameQuery, species_id: str) -> set[str]:
    found: set[str] = set()
    for piece in query.pieces:
        candidates = [piece]
        stripped = strip_prefix(piece, species_id)
        if stripped:
            candidates.append(stripped)
        for candidate in candidates:
            if locus_species(candidate) or gene_key(candidate) != candidate:
                found |= id_forms(candidate)
    return found


def _query_symbols(query: NameQuery, species_id: str) -> list[str]:
    symbols = []
    for piece in query.pieces:
        candidates = [piece]
        stripped = strip_prefix(piece, species_id)
        if stripped:
            candidates.append(stripped)
        for candidate in candidates:
            if locus_species(candidate) or gene_key(candidate) != candidate:
                continue
            symbols.append(norm(candidate))
    # The annotation may keep the prefix the literature omits: GH9B3 is ATGH9B3.
    prefixes = SPECIES_PREFIXES.get(species_id.lower(), ())
    expanded = list(symbols)
    for symbol in symbols:
        for prefix in prefixes:
            expanded.append(prefix + symbol)
    return expanded


def _description_hit(query: NameQuery, species_id: str, view: ProteinView) -> str:
    if query.like:
        return ""
    for piece in query.pieces:
        candidates = [piece]
        stripped = strip_prefix(piece, species_id)
        if stripped:
            candidates.append(stripped)
        for candidate in candidates:
            parts = symbol_parts(candidate)
            if parts is None:
                continue
            stem, number = parts
            for text in view.descriptions:
                if _family_in(stem, text) and _number_in(number, text):
                    return text
    return ""


def match_in_species(query: NameQuery, species_id: str, views: list[ProteinView]) -> tuple[str, ProteinView | None, str]:
    """Return (method, view, detail). method is ambiguous when several description hits exist."""
    pool = [view for view in views if (view.row.get("species_id") or "").lower() == species_id.lower()]
    wanted_ids = _query_ids(query, species_id)
    locus_hits = [view for view in pool if wanted_ids & view.ids]
    if locus_hits:
        # AT1G71930 is the locus_tag of VND7. A longer neighbor that only
        # shares the BAC clone name must not take that id.
        owned = [view for view in locus_hits if wanted_ids & view.primary]
        chosen = owned or locus_hits
        chosen.sort(key=lambda view: (-_length(view), view.row.get("new_id") or ""))
        return "locus", chosen[0], chosen[0].row.get("original_id") or ""

    wanted_symbols = _query_symbols(query, species_id)
    symbol_hits = [view for view in pool if any(symbol in view.symbols for symbol in wanted_symbols)]
    if symbol_hits:
        symbol_hits.sort(key=lambda view: (-_length(view), view.row.get("new_id") or ""))
        return "symbol", symbol_hits[0], wanted_symbols[0]

    described = []
    for view in pool:
        text = _description_hit(query, species_id, view)
        if text:
            described.append((view, text))
    if len(described) == 1:
        view, text = described[0]
        return "description", view, text
    if len(described) > 1:
        preview = "; ".join(
            f"{view.row.get('original_id')} ({text[:80]})" for view, text in described[:6]
        )
        return "ambiguous", None, preview
    return "", None, ""


def _length(view: ProteinView) -> int:
    raw = str(view.row.get("length") or "")
    return int(raw) if raw.isdigit() else 0


def infer_species_id(query: NameQuery) -> str | None:
    for piece in query.pieces:
        sid = locus_species(piece)
        if sid:
            return sid
        for prefix, species_id in PREFIX_TO_SPECIES.items():
            if piece.lower().startswith(prefix) and len(piece) - len(prefix) >= 3:
                rest = piece[len(prefix) :]
                if rest[0].isalpha():
                    return species_id
    return None
