#!/usr/bin/env python3
"""Write genome and codon-alignment HTML pages next to the family report."""

from __future__ import annotations

import argparse
import csv
import html
import re
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--provenance", nargs="*", default=[])
    p.add_argument("--busco", default="")
    p.add_argument("--codon-report", default="")
    p.add_argument("--codon-dir", default="")
    p.add_argument("--families", default="", help="gene_families_of_interest.tsv, for the protein page")
    p.add_argument("--protein-alignments", default="", help="OrthoFinder MultipleSequenceAlignments directory")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path: str | Path | None) -> list[dict]:
    if not path:
        return []
    file = Path(path)
    if not file.is_file() or file.stat().st_size == 0:
        return []
    lines = [ln for ln in file.read_text(errors="replace").splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines or "\t" not in lines[0]:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def num(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def comma(value: float | None, digits: int = 0) -> str:
    if value is None:
        return "—"
    if digits == 0:
        return f"{int(round(value)):,}"
    return f"{value:,.{digits}f}"


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def mix(a: tuple[int, int, int], b: tuple[int, int, int], t: float) -> tuple[int, int, int]:
    t = max(0.0, min(1.0, t))
    return tuple(round(a[i] + (b[i] - a[i]) * t) for i in range(3))


def style_for(rgb: tuple[int, int, int], t: float) -> str:
    color = "#f8f1e8" if t > 0.55 else "#1c1915"
    return f"background:rgb({rgb[0]},{rgb[1]},{rgb[2]});color:{color}"


def good_style(pct: float) -> str:
    t = max(0.0, min(1.0, (pct - 55) / 45))
    return style_for(mix((246, 236, 220), (27, 67, 50), t), t)


def bad_style(pct: float) -> str:
    t = max(0.0, min(1.0, pct / 18))
    return style_for(mix((246, 236, 220), (138, 62, 43), t), t)


def load_provenance(paths: list[str]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        rows.extend(read_tsv(path))
    return rows


def genome_rows(provenance: list[dict], busco: list[dict]) -> list[dict]:
    by: dict[str, dict] = {}
    for row in provenance:
        sid = (row.get("species_id") or "").strip()
        if not sid:
            continue
        item = by.setdefault(sid, {"id": sid})
        item["taxon"] = row.get("configured_taxon") or item.get("taxon") or ""
        item["proteins"] = num(row.get("n_proteins"))
        item["cds"] = num(row.get("n_cds"))
        item["expected"] = num(row.get("expected_protein_coding"))
        item["source"] = row.get("source") or ""
        item["organism"] = row.get("remote_organism") or ""
        item["warning"] = row.get("warning") or ""
    for row in busco:
        sid = (row.get("species_id") or "").strip()
        if not sid or "complete_pct" not in row:
            continue
        item = by.setdefault(sid, {"id": sid})
        item["busco"] = {
            "lineage": row.get("lineage") or "",
            "complete": num(row.get("complete_pct")),
            "single": num(row.get("single_pct")),
            "duplicated": num(row.get("duplicated_pct")),
            "fragmented": num(row.get("fragmented_pct")),
            "missing": num(row.get("missing_pct")),
            "n_markers": num(row.get("n_markers")),
            "min_complete": num(row.get("min_complete_pct")),
            "pass": (row.get("pass") or "").strip().lower(),
        }
    return [by[sid] for sid in sorted(by)]


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


SHELL = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{
    --ink: #1c1915; --muted: #5e584e; --paper: #f3efe6; --card: #fffdf8;
    --line: #e3dacb; --green: #1b4332; --rust: #8a3e2b; --sand: #f4ecdc;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; color: var(--ink); background: var(--paper);
    font-family: "Avenir Next", "Segoe UI", sans-serif; font-size: 15px; line-height: 1.45;
  }}
  .hero {{ background: var(--green); color: #f6f1e6; padding: 2.2rem 1.5rem 4rem; }}
  .wrap {{ max-width: 1180px; margin: 0 auto; }}
  .kicker {{
    margin: 0 0 0.4rem; letter-spacing: 0.14em; text-transform: uppercase;
    font-size: 0.72rem; opacity: 0.75;
  }}
  h1 {{
    margin: 0; font-family: "Iowan Old Style", Palatino, Georgia, serif;
    font-weight: 500; font-size: 2.4rem; letter-spacing: -0.02em;
  }}
  .lede {{ max-width: 42rem; margin: 0.7rem 0 0; color: #e4efe8; }}
  .nav {{ display: flex; flex-wrap: wrap; gap: 1.1rem; margin: 1.15rem 0 0; }}
  .nav a {{ color: #f6f1e6; text-decoration: none; border-bottom: 1px solid transparent; }}
  .nav a[aria-current="page"] {{ border-bottom-color: #f6f1e6; }}
  .stats {{
    display: grid; grid-template-columns: repeat(4, 1fr); gap: 0.8rem; margin-top: -2.4rem;
  }}
  .stat {{
    background: var(--card); border: 1px solid var(--line); border-radius: 10px;
    padding: 0.9rem 1rem 0.8rem; box-shadow: 0 8px 24px rgba(28, 25, 21, 0.05);
  }}
  .stat b {{
    display: block; font-family: "Iowan Old Style", Palatino, Georgia, serif;
    font-size: 1.7rem; font-weight: 500;
  }}
  .stat span {{ color: var(--muted); font-size: 0.82rem; line-height: 1.3; }}
  .panel {{
    background: var(--card); border: 1px solid var(--line); border-radius: 12px;
    padding: 1rem 1rem 0.6rem; margin-top: 1.1rem;
  }}
  h2 {{ margin: 0.2rem 0 0.35rem; font-size: 1.05rem; }}
  .note {{ color: var(--muted); font-size: 0.86rem; margin: 0 0 0.8rem; }}
  .bar-row {{ display: grid; grid-template-columns: 11rem 1fr auto; gap: 0.7rem; align-items: center; margin: 0.45rem 0; }}
  .bar-name {{ font-weight: 600; }}
  .bar-name small {{ display: block; color: var(--muted); font-weight: 400; }}
  .track {{ background: var(--sand); border-radius: 999px; height: 0.72rem; overflow: hidden; }}
  .fill {{ background: var(--green); height: 100%; border-radius: 999px; }}
  .bar-num {{ font-variant-numeric: tabular-nums; min-width: 5.5rem; text-align: right; }}
  .stack-row {{ display: grid; grid-template-columns: 11rem 1fr auto; gap: 0.7rem; align-items: center; margin: 0.55rem 0; }}
  .stack {{ display: flex; height: 0.85rem; border-radius: 999px; overflow: hidden; background: var(--sand); }}
  .stack i {{ display: block; height: 100%; }}
  .legend-pills {{ display: flex; flex-wrap: wrap; gap: 0.7rem; color: var(--muted); font-size: 0.8rem; margin: 0.2rem 0 0.6rem; }}
  .legend-pills i {{ display: inline-block; width: 0.7rem; height: 0.7rem; border-radius: 2px; margin-right: 0.25rem; vertical-align: -1px; }}
  table {{ width: 100%; border-collapse: collapse; }}
  th, td {{ text-align: left; padding: 0.55rem 0.65rem; vertical-align: top; }}
  th {{
    font-size: 0.75rem; letter-spacing: 0.04em; text-transform: uppercase; color: var(--muted);
    font-weight: 600; border-bottom: 1px solid var(--line); white-space: nowrap;
  }}
  td {{ border-top: 1px solid var(--line); }}
  .heat {{ text-align: center; font-variant-numeric: tabular-nums; border-radius: 6px; min-width: 4.2rem; }}
  .scroller {{ overflow-x: auto; }}
  .badge {{
    display: inline-block; border-radius: 999px; padding: 0.05rem 0.45rem; font-size: 0.75rem;
    background: #e7f3ec; color: var(--green);
  }}
  .badge.fail {{ background: #f3ddd4; color: var(--rust); }}
  .warn {{ color: var(--rust); font-size: 0.8rem; }}
  .aln {{ margin: 0.4rem 0 1.1rem; }}
  .block {{ margin: 0.15rem 0 0.55rem; }}
  .seq, .pos {{
    font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
    font-size: 12.5px; line-height: 1.35; white-space: pre;
  }}
  .pos {{ color: var(--muted); }}
  .lab {{ display: inline-block; width: 16rem; overflow: hidden; text-overflow: ellipsis; vertical-align: bottom; }}
  .codon {{ margin-right: 0.28em; }}
  .bA, .aaH {{ background: #d7ebdf; color: #1b4332; }}
  .bC, .aaPos {{ background: #d5e4f2; color: #1a3f66; }}
  .bG, .aaP {{ background: #f4e4bc; color: #6d4e09; }}
  .bT, .aaNeg {{ background: #f3ddd4; color: #8a3e2b; }}
  .bN, .aaG {{ background: #ece7df; color: #5e584e; }}
  .bg {{ color: #b7ad9f; }}
  .tick {{ margin-right: 0.45em; }}
  .footer-note {{ color: var(--muted); font-size: 0.82rem; padding: 0.8rem 0 1.6rem; }}
  @media (max-width: 800px) {{
    .stats {{ grid-template-columns: 1fr 1fr; }}
    .bar-row, .stack-row {{ grid-template-columns: 1fr; }}
    h1 {{ font-size: 1.8rem; }}
  }}
</style>
</head>
<body>
<header class="hero">
  <div class="wrap">
    <p class="kicker">GraftingGeneFamilies</p>
    <h1>{heading}</h1>
    <p class="lede">{lede}</p>
    <nav class="nav" aria-label="Summary pages">{nav}</nav>
  </div>
</header>
<main class="wrap">
{body}
<p class="footer-note">{footer}</p>
</main>
</body>
</html>
"""


def nav(active: str) -> str:
    items = [
        ("families", "index.html", "Families"),
        ("genomes", "genomes.html", "Genomes"),
        ("proteins", "protein_alignments.html", "Protein alignments"),
        ("codons", "codon_alignments.html", "Codon alignments"),
    ]
    parts = []
    for key, href, label in items:
        current = ' aria-current="page"' if key == active else ""
        parts.append(f'<a href="{href}"{current}>{label}</a>')
    return "".join(parts)


def page(title: str, heading: str, lede: str, active: str, body: str, footer: str) -> str:
    return SHELL.format(
        title=esc(title),
        heading=esc(heading),
        lede=esc(lede),
        nav=nav(active),
        body=body,
        footer=esc(footer),
    )


def cards(items: list[tuple[str, str]]) -> str:
    bits = []
    for value, label in items:
        bits.append(f'<article class="stat"><b>{esc(value)}</b><span>{esc(label)}</span></article>')
    return '<section class="stats">' + "".join(bits) + "</section>"


def heat_cell(pct: float | None, kind: str) -> str:
    if pct is None:
        return '<td class="heat">·</td>'
    style = good_style(pct) if kind == "good" else bad_style(pct)
    return f'<td class="heat" style="{style}">{pct:.1f}</td>'


def stacked_bar(busco: dict) -> str:
    parts = [
        ("single", busco.get("single"), "#1b4332"),
        ("duplicated", busco.get("duplicated"), "#40916c"),
        ("fragmented", busco.get("fragmented"), "#c4a15a"),
        ("missing", busco.get("missing"), "#8a3e2b"),
    ]
    total = sum(p for _n, p, _c in parts if p)
    if not total:
        return '<div class="stack"></div>'
    spans = []
    for name, pct, color in parts:
        if not pct:
            continue
        width = 100.0 * pct / total
        spans.append(
            f'<i title="{esc(name)} {pct:.1f}%" style="width:{width:.2f}%;background:{color}"></i>'
        )
    return '<div class="stack">' + "".join(spans) + "</div>"


def render_genomes(rows: list[dict]) -> str:
    proteins = [row["proteins"] for row in rows if row.get("proteins") is not None]
    cds_vals = [row["cds"] for row in rows if row.get("cds") is not None]
    busco_rows = [row for row in rows if row.get("busco")]
    complete = [row["busco"]["complete"] for row in busco_rows if row["busco"].get("complete") is not None]
    passed = sum(1 for row in busco_rows if row["busco"].get("pass") == "pass")
    threshold = next(
        (row["busco"]["min_complete"] for row in busco_rows if row["busco"].get("min_complete") is not None),
        None,
    )
    card_html = cards(
        [
            (comma(len(rows)), "Genomes"),
            (comma(sum(proteins) if proteins else None), "Proteins"),
            (comma(sum(cds_vals) if cds_vals else None), "CDS sequences"),
            (
                f"{passed} of {len(busco_rows)}" if busco_rows else "—",
                "BUSCO passed",
            ),
        ]
    )
    max_proteins = max(proteins) if proteins else 0
    bars = []
    for row in rows:
        count = row.get("proteins")
        width = 0 if not max_proteins or count is None else 100.0 * count / max_proteins
        extra = []
        if row.get("cds") is not None:
            extra.append(f"{comma(row['cds'])} CDS")
        if row.get("expected") is not None:
            extra.append(f"{comma(row['expected'])} protein-coding in the annotation")
        detail = " · ".join(extra)
        warn = f'<div class="warn">{esc(row["warning"])}</div>' if row.get("warning") else ""
        bars.append(
            f'''<div class="bar-row">
              <div class="bar-name">{esc(row["id"])}<small>{esc(row.get("taxon") or row.get("organism") or "")}</small></div>
              <div><div class="track"><div class="fill" style="width:{width:.1f}%"></div></div>{warn}</div>
              <div class="bar-num">{comma(count)}<small style="display:block;color:#5e584e;font-weight:400">{esc(detail)}</small></div>
            </div>'''
        )
    protein_panel = f'''<section class="panel">
      <h2>Proteins per genome</h2>
      <p class="note">Each bar is the number of sequences in the protein FASTA retrieved for that genome. The CDS count is the nucleotide file staged with it. The annotation count, when the source reported one, is the number of protein-coding genes it claimed.</p>
      {"".join(bars) if bars else '<p class="note">No proteome counts were recorded for this run.</p>'}
    </section>'''

    if busco_rows:
        gate = (
            f" A species passes when complete BUSCOs are at least {comma(threshold, 0)}%."
            if threshold is not None
            else ""
        )
        stacks = []
        heat_rows = []
        table_rows = []
        for row in rows:
            busco = row.get("busco")
            if not busco:
                continue
            stacks.append(
                f'''<div class="stack-row">
                  <div class="bar-name">{esc(row["id"])}<small>{esc(row.get("taxon") or "")}</small></div>
                  {stacked_bar(busco)}
                  <div class="bar-num">{comma(busco.get("complete"), 1)}%</div>
                </div>'''
            )
            heat_rows.append(
                "<tr><td>"
                + esc(row["id"])
                + "</td>"
                + heat_cell(busco.get("complete"), "good")
                + heat_cell(busco.get("single"), "good")
                + heat_cell(busco.get("duplicated"), "bad")
                + heat_cell(busco.get("fragmented"), "bad")
                + heat_cell(busco.get("missing"), "bad")
                + "</tr>"
            )
            status = busco.get("pass") or ""
            badge = ""
            if status == "pass":
                badge = '<span class="badge">Pass</span>'
            elif status == "fail":
                badge = '<span class="badge fail">Fail</span>'
            table_rows.append(
                "<tr><td>"
                + esc(row["id"])
                + "</td><td>"
                + esc(busco.get("lineage") or "")
                + "</td><td>"
                + comma(row.get("proteins"))
                + "</td><td>"
                + comma(busco.get("complete"), 2)
                + "</td><td>"
                + comma(busco.get("single"), 2)
                + "</td><td>"
                + comma(busco.get("duplicated"), 2)
                + "</td><td>"
                + comma(busco.get("fragmented"), 2)
                + "</td><td>"
                + comma(busco.get("missing"), 2)
                + "</td><td>"
                + comma(busco.get("n_markers"))
                + "</td><td>"
                + badge
                + "</td></tr>"
            )
        median_txt = comma(median(complete), 1)
        busco_panel = f'''<section class="panel">
      <h2>BUSCO</h2>
      <p class="note">Protein-mode BUSCO. Complete is single-copy plus duplicated. Median complete is {median_txt}%.{esc(gate)} The stacked bar is the share of markers in each class. In the heatmap, complete and single-copy darken toward green as they rise. Duplicated, fragmented, and missing darken toward brown, scaled so a small percentage is already visible.</p>
      <div class="legend-pills">
        <span><i style="background:#1b4332"></i>Single-copy</span>
        <span><i style="background:#40916c"></i>Duplicated</span>
        <span><i style="background:#c4a15a"></i>Fragmented</span>
        <span><i style="background:#8a3e2b"></i>Missing</span>
      </div>
      {"".join(stacks)}
      <div class="scroller">
        <table>
          <thead><tr>
            <th>Species</th><th>Complete</th><th>Single</th><th>Duplicated</th><th>Fragmented</th><th>Missing</th>
          </tr></thead>
          <tbody>{"".join(heat_rows)}</tbody>
        </table>
      </div>
    </section>
    <section class="panel">
      <h2>Counts</h2>
      <div class="scroller">
        <table>
          <thead><tr>
            <th>Species</th><th>Lineage</th><th>Proteins</th><th>Complete</th><th>Single</th>
            <th>Duplicated</th><th>Fragmented</th><th>Missing</th><th>Markers</th><th>Gate</th>
          </tr></thead>
          <tbody>{"".join(table_rows)}</tbody>
        </table>
      </div>
    </section>'''
    else:
        busco_panel = '''<section class="panel">
      <h2>BUSCO</h2>
      <p class="note">BUSCO was not recorded for this run.</p>
    </section>'''
    return page(
        "Genomes",
        "Genomes",
        "Protein and CDS counts for each genome, and BUSCO completeness where that step ran.",
        "genomes",
        card_html + protein_panel + busco_panel,
        "Darker green in the BUSCO heatmap is a higher complete or single-copy percentage. Darker brown is a higher duplicated, fragmented, or missing percentage.",
    )


def read_fasta(path: Path, nucleotide: bool = True) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                records.append((name, clean_seq("".join(chunks), nucleotide)))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        records.append((name, clean_seq("".join(chunks), nucleotide)))
    return records


def clean_seq(seq: str, nucleotide: bool) -> str:
    seq = seq.upper()
    if nucleotide:
        seq = seq.replace("U", "T")
    return seq


def short_label(seq_id: str) -> str:
    if "|" in seq_id:
        parts = seq_id.split("|")
        return f"{parts[0]} {parts[-1]}"
    match = re.match(r"^([A-Za-z0-9]+)_([A-Z][a-z]+)_([a-z]+)_(.+)$", seq_id)
    if match:
        return f"{match.group(1)} {match.group(4)}"
    return seq_id


BASE_CLASS = {"A": "bA", "C": "bC", "G": "bG", "T": "bT", "U": "bT"}


def codon_html(codon: str) -> str:
    letters = []
    for base in codon:
        if base == "-":
            letters.append('<span class="bg">-</span>')
        else:
            cls = BASE_CLASS.get(base, "bN")
            letters.append(f'<span class="{cls}">{esc(base)}</span>')
    return '<span class="codon">' + "".join(letters) + "</span>"


def alignment_blocks(records: list[tuple[str, str]], width: int = 60) -> str:
    if not records:
        return ""
    length = max(len(seq) for _name, seq in records)
    padded = [(name, seq.ljust(length, "-")) for name, seq in records]
    blocks = []
    for start in range(0, length, width):
        rows = [
            f'<div class="pos"><span class="lab"></span>{start + 1}</div>'
        ]
        for name, seq in padded:
            chunk = seq[start : start + width]
            codons = [codon_html(chunk[i : i + 3]) for i in range(0, len(chunk), 3)]
            rows.append(
                f'<div class="seq" title="{esc(name)}"><span class="lab">{esc(short_label(name))}</span>{"".join(codons)}</div>'
            )
        blocks.append('<div class="block">' + "".join(rows) + "</div>")
    return "".join(blocks)


STATUS_LABEL = {
    "ok": "Written",
    "too_few_sequences": "Too few sequences",
    "alignment_missing": "No protein alignment",
    "unassigned_no_alignment": "Unassigned, no alignment",
}


def find_alignment(codon_dir: Path | None, og: str) -> Path | None:
    if codon_dir is None or not codon_dir.is_dir():
        return None
    for name in (f"{og}.codon.fna", f"{og}.codon.fa", f"{og}.codon.fasta", f"{og}.fna", f"{og}.fa"):
        candidate = codon_dir / name
        if candidate.is_file() and candidate.stat().st_size:
            return candidate
    return None


def render_codons(report: list[dict], codon_dir: Path | None) -> str:
    written = [row for row in report if (row.get("status") or "") == "ok"]
    missing = sum(int(num(row.get("n_missing_cds")) or 0) for row in report)
    unassigned = sum(1 for row in report if (row.get("status") or "").startswith("unassigned"))
    card_html = cards(
        [
            (comma(len(written)), "Alignments written"),
            (comma(len(report)), "Orthogroups"),
            (comma(missing), "Sequences missing CDS"),
            (comma(unassigned), "Unassigned, no alignment"),
        ]
    )
    table_rows = []
    blocks = []
    for row in report:
        og = row.get("Orthogroup") or ""
        status = row.get("status") or ""
        label = STATUS_LABEL.get(status, status or "—")
        table_rows.append(
            "<tr><td><a href=\"#og-"
            + esc(og)
            + "\">"
            + esc(og)
            + "</a></td><td>"
            + comma(num(row.get("n_in_msa")))
            + "</td><td>"
            + comma(num(row.get("n_codon")))
            + "</td><td>"
            + comma(num(row.get("n_missing_cds")))
            + "</td><td>"
            + comma(num(row.get("n_project_fail")))
            + "</td><td>"
            + esc(label)
            + "</td></tr>"
        )
        path = find_alignment(codon_dir, og) if og else None
        if path:
            records = read_fasta(path)
            length = max((len(seq) for _n, seq in records), default=0)
            blocks.append(
                f'''<section class="panel" id="og-{esc(og)}">
                  <h2>{esc(og)}</h2>
                  <p class="note">{len(records)} sequences, {comma(length)} nucleotides. Codons are grouped in threes. A, C, G, and T keep their own color, and a gap is a hyphen.</p>
                  <div class="aln">{alignment_blocks(records)}</div>
                </section>'''
            )
        elif og:
            blocks.append(
                f'''<section class="panel" id="og-{esc(og)}">
                  <h2>{esc(og)}</h2>
                  <p class="note">{esc(label)}. No codon alignment was written for this orthogroup.</p>
                </section>'''
            )
    table = f'''<section class="panel">
      <h2>Codon alignment report</h2>
      <p class="note">Each alignment is the CDS projected onto the OrthoFinder protein alignment. A gap in the protein becomes three gaps in the codon alignment.</p>
      <div class="legend-pills">
        <span><i class="bA"></i>A</span>
        <span><i class="bC"></i>C</span>
        <span><i class="bG"></i>G</span>
        <span><i class="bT"></i>T</span>
        <span><i class="bg" style="background:#f4ecdc"></i>gap</span>
      </div>
      <div class="scroller">
        <table>
          <thead><tr>
            <th>Orthogroup</th><th>In protein MSA</th><th>Codon sequences</th>
            <th>Missing CDS</th><th>Projection failed</th><th>Status</th>
          </tr></thead>
          <tbody>{"".join(table_rows) if table_rows else '<tr><td colspan="6">No codon alignment report was recorded for this run.</td></tr>'}</tbody>
        </table>
      </div>
    </section>'''
    return page(
        "Codon alignments",
        "Codon alignments",
        "CDS sequences projected onto the protein alignments for the selected orthogroups.",
        "codons",
        card_html + table + "".join(blocks),
        "An unassigned orthogroup has no protein alignment, so there is no codon alignment to show.",
    )


AA_CLASS = {
    "A": "aaH",
    "V": "aaH",
    "I": "aaH",
    "L": "aaH",
    "M": "aaH",
    "F": "aaH",
    "W": "aaH",
    "P": "aaH",
    "S": "aaP",
    "T": "aaP",
    "N": "aaP",
    "Q": "aaP",
    "Y": "aaP",
    "C": "aaP",
    "K": "aaPos",
    "R": "aaPos",
    "H": "aaPos",
    "D": "aaNeg",
    "E": "aaNeg",
    "G": "aaG",
}


def residue_group(chunk: str) -> str:
    letters = []
    for aa in chunk:
        if aa == "-":
            letters.append('<span class="bg">-</span>')
        else:
            letters.append(f'<span class="{AA_CLASS.get(aa, "bN")}">{esc(aa)}</span>')
    return '<span class="tick">' + "".join(letters) + "</span>"


def protein_blocks(records: list[tuple[str, str]], width: int = 60) -> str:
    if not records:
        return ""
    length = max(len(seq) for _name, seq in records)
    padded = [(name, seq.ljust(length, "-")) for name, seq in records]
    blocks = []
    for start in range(0, length, width):
        rows = [f'<div class="pos"><span class="lab"></span>{start + 1}</div>']
        for name, seq in padded:
            chunk = seq[start : start + width]
            groups = [residue_group(chunk[i : i + 10]) for i in range(0, len(chunk), 10)]
            rows.append(
                f'<div class="seq" title="{esc(name)}"><span class="lab">{esc(short_label(name))}</span>{"".join(groups)}</div>'
            )
        blocks.append('<div class="block">' + "".join(rows) + "</div>")
    return "".join(blocks)


def find_protein_alignment(align_dir: Path | None, og: str) -> Path | None:
    if align_dir is None or not align_dir.is_dir():
        return None
    preferred = [
        align_dir / f"{og}.fa",
        align_dir / f"{og}.faa",
        align_dir / f"{og}.fasta",
        align_dir / f"{og}.aln",
        align_dir / "MultipleSequenceAlignments" / f"{og}.fa",
    ]
    for candidate in preferred:
        if candidate.is_file() and candidate.stat().st_size and ".codon" not in candidate.name:
            return candidate
    hits = [hit for hit in align_dir.rglob(f"{og}.*") if ".codon" not in hit.name]
    msa_hits = [hit for hit in hits if "MultipleSequenceAlignments" in hit.parts]
    for hit in msa_hits or hits:
        if hit.suffix.lower() in {".fa", ".faa", ".fasta", ".aln"} and hit.stat().st_size:
            return hit
    return None


def selected_orthogroups(families_path: str, codon_report: list[dict]) -> list[str]:
    ogs: list[str] = []
    for row in read_tsv(families_path):
        og = (row.get("Orthogroup") or "").strip()
        if og and og not in ogs:
            ogs.append(og)
    for row in codon_report:
        og = (row.get("Orthogroup") or "").strip()
        if og and og not in ogs:
            ogs.append(og)
    return ogs


def render_proteins(ogs: list[str], align_dir: Path | None) -> str:
    shown = 0
    n_seq = 0
    missing = 0
    table_rows = []
    blocks = []
    for og in ogs:
        path = find_protein_alignment(align_dir, og)
        if not path:
            missing += 1
            note = "Unassigned genes have no OrthoFinder alignment." if og.startswith("UNASSIGNED_") else "No protein alignment was found."
            table_rows.append(
                f"<tr><td><a href=\"#aa-{esc(og)}\">{esc(og)}</a></td><td>—</td><td>—</td><td>{esc(note)}</td></tr>"
            )
            blocks.append(
                f'''<section class="panel" id="aa-{esc(og)}">
                  <h2>{esc(og)}</h2>
                  <p class="note">{esc(note)}</p>
                </section>'''
            )
            continue
        records = read_fasta(path, nucleotide=False)
        length = max((len(seq) for _name, seq in records), default=0)
        shown += 1
        n_seq += len(records)
        table_rows.append(
            "<tr><td><a href=\"#aa-"
            + esc(og)
            + "\">"
            + esc(og)
            + "</a></td><td>"
            + comma(len(records))
            + "</td><td>"
            + comma(length)
            + "</td><td>Shown</td></tr>"
        )
        blocks.append(
            f'''<section class="panel" id="aa-{esc(og)}">
              <h2>{esc(og)}</h2>
              <p class="note">{len(records)} sequences, {comma(length)} columns. Residues are grouped in tens. Hydrophobic residues are green, polar residues gold, positive blue, negative rust, and glycine gray.</p>
              <div class="aln">{protein_blocks(records)}</div>
            </section>'''
        )
    card_html = cards(
        [
            (comma(shown), "Alignments shown"),
            (comma(len(ogs)), "Orthogroups"),
            (comma(n_seq), "Sequences"),
            (comma(missing), "Missing alignment"),
        ]
    )
    table = f'''<section class="panel">
      <h2>Protein alignments</h2>
      <p class="note">These are the OrthoFinder amino-acid alignments for the selected orthogroups. The codon page is the same alignment with each residue expanded to nucleotides.</p>
      <div class="legend-pills">
        <span><i class="aaH"></i>Hydrophobic</span>
        <span><i class="aaP"></i>Polar</span>
        <span><i class="aaPos"></i>Positive</span>
        <span><i class="aaNeg"></i>Negative</span>
        <span><i class="aaG"></i>Glycine</span>
        <span><i class="bg" style="background:#f4ecdc"></i>gap</span>
      </div>
      <div class="scroller">
        <table>
          <thead><tr><th>Orthogroup</th><th>Sequences</th><th>Columns</th><th>Status</th></tr></thead>
          <tbody>{"".join(table_rows) if table_rows else '<tr><td colspan="4">No selected orthogroups were available for this page.</td></tr>'}</tbody>
        </table>
      </div>
    </section>'''
    return page(
        "Protein alignments",
        "Protein alignments",
        "Amino-acid alignments of the selected orthogroups.",
        "proteins",
        card_html + table + "".join(blocks),
        "An unassigned orthogroup has no protein alignment.",
    )


def copy_table(src: str, dest: Path) -> None:
    rows = read_tsv(src)
    if not rows:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = load_provenance(args.provenance)
    busco = read_tsv(args.busco)
    report = read_tsv(args.codon_report)
    codon_dir = Path(args.codon_dir) if args.codon_dir else None
    if codon_dir is not None and not codon_dir.is_dir():
        codon_dir = None
    protein_dir = Path(args.protein_alignments) if args.protein_alignments else None
    if protein_dir is not None and not protein_dir.is_dir():
        protein_dir = None
    ogs = selected_orthogroups(args.families, report)
    rows = genome_rows(provenance, busco)
    (outdir / "genomes.html").write_text(render_genomes(rows))
    (outdir / "protein_alignments.html").write_text(render_proteins(ogs, protein_dir))
    (outdir / "codon_alignments.html").write_text(render_codons(report, codon_dir))
    copy_table(args.busco, outdir / "tables" / "busco_summary.tsv")
    copy_table(args.codon_report, outdir / "tables" / "codon_alignment_report.tsv")
    print(f"[write_summary_pages] {len(rows)} genomes, {len(ogs)} protein alignments, {len(report)} codon rows → {outdir}")


if __name__ == "__main__":
    main()
