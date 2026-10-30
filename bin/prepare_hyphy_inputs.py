#!/usr/bin/env python3
"""Match a codon alignment to the ASR species tree for HyPhy.

IQ-TREE names tips with species ids (Ath, Sly). Codon alignments keep the
full sequence ids (Ath|Arabidopsis_thaliana|NP_195280.1). HyPhy rejects that
pair. This expands each species tip into the codon sequences from that
species and tags foreground species with {Foreground}.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


class Node:
    def __init__(self) -> None:
        self.name: str | None = None
        self.length: str | None = None
        self.children: list[Node] = []
        self.foreground: bool = False

    def tips(self) -> list[Node]:
        if not self.children:
            return [self]
        found: list[Node] = []
        for child in self.children:
            found.extend(child.tips())
        return found


_NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
_PVALUE = re.compile(
    r'"p-value"\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)',
    re.I,
)


def parse_newick(text: str) -> Node:
    raw = text.strip()
    if not raw:
        raise ValueError("empty tree")
    if raw.endswith(";"):
        raw = raw[:-1]
    pos = 0

    def skip_ws() -> None:
        nonlocal pos
        while pos < len(raw) and raw[pos].isspace():
            pos += 1

    def parse_label() -> str:
        nonlocal pos
        skip_ws()
        if pos >= len(raw):
            return ""
        if raw[pos] in "'\"":
            quote = raw[pos]
            pos += 1
            chars: list[str] = []
            while pos < len(raw):
                if raw[pos] == quote:
                    if pos + 1 < len(raw) and raw[pos + 1] == quote:
                        chars.append(quote)
                        pos += 2
                        continue
                    pos += 1
                    break
                chars.append(raw[pos])
                pos += 1
            return "".join(chars)
        start = pos
        while pos < len(raw) and raw[pos] not in ",):;":
            pos += 1
        return raw[start:pos].strip()

    def parse_length() -> str | None:
        nonlocal pos
        skip_ws()
        if pos >= len(raw) or raw[pos] != ":":
            return None
        pos += 1
        skip_ws()
        match = _NUMBER.match(raw, pos)
        if not match:
            raise ValueError(f"missing branch length at {raw[pos:pos + 20]!r}")
        pos = match.end()
        return match.group(0)

    def parse_node() -> Node:
        nonlocal pos
        skip_ws()
        node = Node()
        if pos < len(raw) and raw[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                skip_ws()
                if pos < len(raw) and raw[pos] == ",":
                    pos += 1
                    continue
                if pos < len(raw) and raw[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"malformed newick near {raw[pos:pos + 20]!r}")
        label = parse_label()
        node.foreground = "Foreground" in re.findall(r"\{(.*?)\}", label)
        label = re.sub(r"\{.*?\}", "", label).strip()
        if node.children and re.fullmatch(r"[\d.]+", label or ""):
            label = ""
        node.name = label or None
        node.length = parse_length()
        return node

    root = parse_node()
    skip_ws()
    if pos < len(raw):
        raise ValueError(f"trailing newick text {raw[pos:pos + 20]!r}")
    if not root.tips():
        raise ValueError("tree has no tips")
    return root


def add_lengths(left: str | None, right: str | None) -> str | None:
    if left is None:
        return right
    if right is None:
        return left
    return format(float(left) + float(right), ".10g")


def collapse(node: Node) -> None:
    for child in node.children:
        collapse(child)
    changed = True
    while changed:
        changed = False
        lifted: list[Node] = []
        for child in node.children:
            if len(child.children) == 1:
                grand = child.children[0]
                grand.length = add_lengths(child.length, grand.length)
                grand.foreground = grand.foreground or child.foreground
                lifted.append(grand)
                changed = True
            else:
                lifted.append(child)
        node.children = lifted


def lift_root(root: Node) -> Node:
    while len(root.children) == 1:
        child = root.children[0]
        child.length = add_lengths(root.length, child.length)
        child.foreground = child.foreground or root.foreground
        root = child
    return root


def species_of(seq_id: str, species_ids: set[str]) -> str | None:
    if seq_id in species_ids:
        return seq_id
    head = seq_id.split("|", 1)[0]
    if head in species_ids:
        return head
    prefixed = [
        sp
        for sp in species_ids
        if seq_id.startswith(sp + "_") or seq_id.startswith(sp + "|")
    ]
    if not prefixed:
        return None
    return max(prefixed, key=len)


def expand(node: Node, groups: dict[str, list[str]], foreground: set[str]) -> Node | None:
    if not node.children:
        species = node.name or ""
        seqs = groups.get(species, [])
        if not seqs:
            return None
        if len(seqs) == 1:
            node.name = seqs[0]
            node.foreground = species in foreground
            return node
        node.name = None
        node.foreground = False
        node.children = []
        for seq_id in seqs:
            child = Node()
            child.name = seq_id
            child.length = "0.000001"
            child.foreground = species in foreground
            node.children.append(child)
        return node
    kept = []
    for child in node.children:
        expanded = expand(child, groups, foreground)
        if expanded is not None:
            kept.append(expanded)
    node.children = kept
    if not kept:
        return None
    return node


def mark_foreground_clades(node: Node) -> bool:
    if not node.children:
        return node.foreground
    flags = [mark_foreground_clades(child) for child in node.children]
    if flags and all(flags):
        node.foreground = True
    return node.foreground


def quote_name(name: str) -> str:
    if re.search(r"[(),:;\s'\[\]]", name):
        return "'" + name.replace("'", "''") + "'"
    return name


def write_newick(node: Node, is_root: bool = False) -> str:
    if node.children:
        body = ",".join(write_newick(child) for child in node.children)
        text = f"({body})"
    else:
        text = quote_name(node.name or "")
    if node.foreground and not is_root:
        text += "{Foreground}"
    if node.length is not None and not is_root:
        text += f":{node.length}"
    return text


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal name
        if name is not None:
            records.append((name, "".join(chunks)))

    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    flush()
    return records


def read_foreground(path: Path) -> set[str]:
    return {ln.strip() for ln in path.read_text(errors="replace").splitlines() if ln.strip()}


def write_status(path: Path, status: str, note: str, n_fg: int, n_seq: int) -> None:
    note = note.replace("\t", " ").replace("\n", " ").strip()
    path.write_text(f"{status}\t{note}\t{n_fg}\t{n_seq}\n")


def prepare(
    codon_path: Path,
    tree_path: Path,
    foreground_path: Path,
    out_fasta: Path,
    out_tree: Path,
    out_status: Path,
) -> None:
    try:
        records = read_fasta(codon_path)
        foreground = read_foreground(foreground_path)
        root = parse_newick(tree_path.read_text(errors="replace"))
    except Exception as exc:
        write_status(out_status, "skip", f"could not read HyPhy inputs: {exc}", 0, 0)
        return

    species_ids = {tip.name for tip in root.tips() if tip.name}
    groups: dict[str, list[str]] = {}
    order: list[str] = []
    seen: set[str] = set()
    for seq_id, _seq in records:
        if seq_id in seen:
            continue
        seen.add(seq_id)
        species = species_of(seq_id, species_ids)
        if species is None:
            continue
        groups.setdefault(species, []).append(seq_id)
        order.append(seq_id)

    expanded = expand(root, groups, foreground)
    if expanded is None:
        write_status(out_status, "skip", "no codon sequences match the species tree", 0, 0)
        return
    collapse(expanded)
    expanded = lift_root(expanded)
    tips = expanded.tips()
    n_seq = len(tips)
    n_fg = sum(1 for tip in tips if tip.foreground)
    if n_seq < 3:
        write_status(
            out_status,
            "skip",
            f"fewer than 3 codon sequences match the species tree (n={n_seq})",
            n_fg,
            n_seq,
        )
        return

    mark_foreground_clades(expanded)
    by_id = {seq_id: seq for seq_id, seq in records}
    fasta_lines = []
    for tip in tips:
        seq = by_id.get(tip.name or "", "")
        fasta_lines.append(f">{tip.name}\n{seq}\n")
    out_fasta.write_text("".join(fasta_lines))
    out_tree.write_text(write_newick(expanded, is_root=True) + ";\n")
    write_status(
        out_status,
        "ok",
        "codon sequence ids placed on the species tree",
        n_fg,
        n_seq,
    )


def pvalue_from_json(text: str) -> str:
    idx = text.lower().find('"test results"')
    if idx >= 0:
        match = _PVALUE.search(text[idx:])
        if match:
            return match.group(1)
    found = [item.group(1) for item in _PVALUE.finditer(text)]
    if not found:
        return ""
    return min(found, key=float)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    prep = sub.add_parser("prepare", help="Write a codon FASTA and labeled tree HyPhy can read")
    prep.add_argument("--codon", required=True, type=Path)
    prep.add_argument("--tree", required=True, type=Path)
    prep.add_argument("--foreground", required=True, type=Path)
    prep.add_argument("--out-fasta", required=True, type=Path)
    prep.add_argument("--out-tree", required=True, type=Path)
    prep.add_argument("--out-status", required=True, type=Path)

    pv = sub.add_parser("pvalue", help="Print the HyPhy test p-value, or nothing")
    pv.add_argument("json_path", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.cmd == "prepare":
        prepare(
            args.codon,
            args.tree,
            args.foreground,
            args.out_fasta,
            args.out_tree,
            args.out_status,
        )
        return 0
    try:
        text = args.json_path.read_text(errors="replace")
    except OSError:
        text = ""
    print(pvalue_from_json(text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
