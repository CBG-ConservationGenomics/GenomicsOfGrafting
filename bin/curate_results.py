#!/usr/bin/env python3
"""Write a browsable HTML report and a small folder of curated tables."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import date
from pathlib import Path


SPECIES_LABELS = {
    "Sly": "S. lycopersicum",
    "Nbe": "N. benthamiana",
    "Ath": "A. thaliana",
    "Cpe": "C. pepo",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--families", required=True, help="gene_families_of_interest.tsv")
    p.add_argument("--wide-identity", required=True, help="gene_families_pctid_by_species.tsv")
    p.add_argument("--ranking", default=None, help="graft_candidate_ranking.tsv")
    p.add_argument("--alphafold-dir", default=None, help="Directory containing og_priority.tsv")
    p.add_argument("--provenance", nargs="*", default=[], help="Per-species download provenance TSVs")
    p.add_argument("--identity-matrix-pairs", default=None)
    p.add_argument("--identity-matrix-labels", default=None)
    p.add_argument(
        "--runtime-seconds",
        type=float,
        default=None,
        help="Seconds from job start until this report is written",
    )
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def load_matrices(pairs_path: str | None, labels_path: str | None) -> dict[str, dict]:
    """Square identity matrix for each orthogroup, keyed by sequence id."""
    if not pairs_path or not labels_path:
        return {}
    pairs_file = Path(pairs_path)
    labels_file = Path(labels_path)
    if not pairs_file.is_file() or not labels_file.is_file():
        return {}
    labels_by_og: dict[str, list[dict]] = {}
    with labels_file.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            labels_by_og.setdefault(row["Orthogroup"], []).append(row)
    values: dict[tuple[str, str, str], float] = {}
    with pairs_file.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = num(row.get("pct_identity"))
            if pid is None:
                continue
            og = row["Orthogroup"]
            values[(og, row["seq_a"], row["seq_b"])] = pid
            values[(og, row["seq_b"], row["seq_a"])] = pid
    matrices: dict[str, dict] = {}
    for og, rows in labels_by_og.items():
        ids = [row["sequence_id"] for row in rows]
        grid = []
        for left in ids:
            grid.append([None if left == right else values.get((og, left, right)) for right in ids])
        matrices[og] = {
            "labels": [row["label"] for row in rows],
            "full": ids,
            "reference": [row.get("is_reference") == "yes" for row in rows],
            "values": grid,
        }
    return matrices


def read_plain_tsv(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_provenance(paths: list[str]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        file = Path(path)
        if file.is_file() and file.stat().st_size:
            rows.extend(read_plain_tsv(str(file)))
    rows.sort(key=lambda row: (row.get("species_id") or "", row.get("proteins_source") or ""))
    return rows


def read_tsv(path: str | None) -> list[dict]:
    if not path:
        return []
    file = Path(path)
    if not file.is_file() or file.stat().st_size == 0:
        return []
    with file.open() as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    if not lines or "Orthogroup" not in lines[0]:
        return []
    import io

    return list(csv.DictReader(io.StringIO("".join(lines)), delimiter="\t"))


def num(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def split_list(value: str | None, sep: str) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in str(value).split(sep) if part.strip()]


def mined_gene_count(records: list[dict]) -> int:
    """Literature genes that were placed in an orthogroup."""
    seen: set[tuple[str, str]] = set()
    for rec in records:
        genes = split_list(rec.get("genes"), ";")
        species = split_list(rec.get("queries"), ";")
        for i, gene in enumerate(genes):
            sid = species[i] if i < len(species) else ""
            seen.add((gene, sid))
    return len(seen)


def species_run_counts(provenance: list[dict], identity_species: list[dict]) -> tuple[int, int | None]:
    """Species included in the job, and how many of those were remote downloads."""
    if not provenance:
        return len(identity_species), None
    seen: set[str] = set()
    downloaded = 0
    for row in provenance:
        sid = (row.get("species_id") or "").strip()
        if not sid or sid in seen:
            continue
        seen.add(sid)
        if (row.get("source") or "").strip().lower() in {"ncbi", "url"}:
            downloaded += 1
    return len(seen), downloaded


def format_runtime(seconds: float | None) -> str:
    if seconds is None:
        return "—"
    total = max(0, int(round(float(seconds))))
    days, rem = divmod(total, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, _secs = divmod(rem, 60)
    if days:
        return f"{days}d {hours}h" if hours else f"{days}d"
    if hours:
        return f"{hours}h {minutes}m" if minutes else f"{hours}h"
    if minutes:
        return f"{minutes}m"
    return f"{total}s"


def species_from_wide(rows: list[dict]) -> list[dict]:
    if not rows:
        return []
    prefix = "best_pct_id_"
    if not any(key.startswith(prefix) for key in rows[0]):
        prefix = "mean_pct_id_"
    species = []
    for key in rows[0]:
        if key.startswith(prefix):
            sid = key[len(prefix) :]
            species.append({"id": sid, "label": SPECIES_LABELS.get(sid, sid)})
    return species


def build_records(
    families: list[dict],
    wide: list[dict],
    ranking: list[dict],
    folding: list[dict],
    species: list[dict],
) -> list[dict]:
    wide_by = {r.get("Orthogroup", ""): r for r in wide}
    rank_by = {r.get("Orthogroup", ""): r for r in ranking}
    fold_by = {r.get("Orthogroup", ""): r for r in folding}

    records = []
    for row in families:
        og = row.get("Orthogroup", "")
        wide_row = wide_by.get(og, {})
        rank = rank_by.get(og, {})
        fold = fold_by.get(og, {})
        identity = {}
        for sp in species:
            best = num(wide_row.get(f"best_pct_id_{sp['id']}"))
            if best is None:
                best = num(wide_row.get(f"mean_pct_id_{sp['id']}"))
            identity[sp["id"]] = best
        priority = num(rank.get("priority_score"))
        if priority is None:
            priority = num(fold.get("priority_score"))
        records.append(
            {
                "orthogroup": og,
                "genes": row.get("genes_of_interest", ""),
                "queries": row.get("query_species_ids", ""),
                "matched": row.get("matched_members", ""),
                "reference": row.get("reference_id", ""),
                "reference_gene": row.get("reference_gene", ""),
                "reference_species": row.get("reference_species", ""),
                "n_members": num(row.get("n_members")),
                "n_species": num(row.get("n_species_with_members")),
                "mean_id": num(row.get("mean_pct_id_to_ref")),
                "min_id": num(row.get("min_pct_id_to_ref")),
                "max_id": num(row.get("max_pct_id_to_ref")),
                "name": row.get("eggnog_preferred_name", ""),
                "description": row.get("eggnog_description", ""),
                "gos": split_list(row.get("eggnog_gos"), ","),
                "interpro": row.get("interpro_domains", ""),
                "signalp": row.get("signalp", ""),
                "deeptmhmm": row.get("deeptmhmm", ""),
                "annotation_from": row.get("annotation_source_id", ""),
                "members": split_list(row.get("member_ids"), ","),
                "identity": identity,
                "priority_score": priority,
                "n_asr": num(rank.get("n_asr_subs_on_shift_branches")),
                "shift_subs": rank.get("example_shift_subs") or fold.get("example_shift_subs") or "",
                "divergence": rank.get("divergence_candidate") or fold.get("divergence_candidate") or "",
                "mantel_p": num(rank.get("partial_mantel_p")),
                "mean_div_failure": num(rank.get("mean_div_failure")),
                "mean_div_success": num(rank.get("mean_div_success")),
                "n_selection": num(rank.get("n_selection_hits")),
                "n_convergence": num(rank.get("n_convergence_ok")),
                "fold_rank": num(fold.get("fold_rank")),
                "matrix": None,
            }
        )

    def sort_key(rec: dict) -> tuple:
        priority = rec["priority_score"]
        best_id = rec["max_id"]
        return (
            -(priority if priority is not None else -1),
            -(best_id if best_id is not None else -1),
            rec["orthogroup"],
        )

    records.sort(key=sort_key)
    return records


def write_readme(path: Path, copied: list[str]) -> None:
    lines = [
        "Curated results",
        "===============",
        "",
        "Open index.html in a browser. The page is self-contained:",
        "search, sort, and expand a row for annotation and phenotype detail.",
        "genomes.html has protein counts and BUSCO. protein_alignments.html",
        "shows the OrthoFinder amino-acid alignments. codon_alignments.html",
        "shows the CDS alignments projected onto those protein MSAs.",
        "",
        "tables/",
    ]
    for name in copied:
        lines.append(f"  {name}")
    lines.append("")
    lines.append("These tables are copies of the pipeline outputs used to build the report.")
    lines.append("Priority score, when present, is the graft-candidate rank")
    lines.append("(ASR substitutions on phenotype-shift branches, divergence beyond")
    lines.append("relatedness, HyPhy hits, and convergence).")
    lines.append("")
    path.write_text("\n".join(lines))


HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Curated gene families</title>
<style>
  :root {
    --ink: #1c1915;
    --muted: #5e584e;
    --paper: #f3efe6;
    --card: #fffdf8;
    --line: #e3dacb;
    --green: #1b4332;
    --green-2: #2d6a4f;
    --leaf: #e7f3ec;
    --rust: #8a3e2b;
    --sand: #f4ecdc;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    color: var(--ink);
    background: var(--paper);
    font-family: "Avenir Next", "Segoe UI", sans-serif;
    font-size: 15px;
    line-height: 1.45;
  }
  .hero {
    background: var(--green);
    color: #f6f1e6;
    padding: 2.4rem 1.5rem 4.2rem;
  }
  .wrap { max-width: 1180px; margin: 0 auto; }
  .kicker {
    margin: 0 0 0.4rem;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    font-size: 0.72rem;
    opacity: 0.75;
  }
  h1 {
    margin: 0;
    font-family: "Iowan Old Style", Palatino, Georgia, serif;
    font-weight: 500;
    font-size: 2.4rem;
    letter-spacing: -0.02em;
  }
  .lede { max-width: 40rem; margin: 0.7rem 0 0; color: #e4efe8; }
  .stats {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 0.8rem;
    margin-top: -2.4rem;
  }
  .stat {
    background: var(--card);
    border: 1px solid var(--line);
    border-radius: 10px;
    padding: 0.9rem 1rem 0.8rem;
    box-shadow: 0 8px 24px rgba(28, 25, 21, 0.05);
  }
  .stat b {
    display: block;
    font-family: "Iowan Old Style", Palatino, Georgia, serif;
    font-size: 1.7rem;
    font-weight: 500;
  }
  .stat span { color: var(--muted); font-size: 0.82rem; line-height: 1.3; }
  .panel {
    margin-top: 1.2rem;
    background: var(--card);
    border: 1px solid var(--line);
    border-radius: 12px;
    overflow: hidden;
  }
  h2 {
    margin: 0;
    padding: 0.9rem 1rem 0.15rem;
    font-family: "Iowan Old Style", Palatino, Georgia, serif;
    font-weight: 500;
    font-size: 1.25rem;
  }
  .source-note { margin: 0; padding: 0 1rem 0.7rem; color: var(--muted); font-size: 0.85rem; }
  tr.warn td { color: var(--rust); }
  td.url { word-break: break-all; }
  .heat-wrap { margin-top: 0.9rem; overflow-x: auto; }
  .heat-note { margin: 0 0 0.55rem; color: var(--muted); font-size: 0.85rem; max-width: 46rem; }
  table.heat { width: auto; border-collapse: separate; border-spacing: 3px; }
  table.heat th, table.heat td {
    position: static;
    padding: 0;
    border: 0;
    background: transparent;
    text-transform: none;
    letter-spacing: 0;
    font-weight: 500;
    vertical-align: middle;
  }
  table.heat th.col { height: 7.2rem; vertical-align: bottom; min-width: 2.7rem; }
  table.heat th.col span {
    display: inline-block;
    writing-mode: vertical-rl;
    transform: rotate(180deg);
    white-space: nowrap;
    font-size: 0.72rem;
    letter-spacing: 0;
    text-transform: none;
  }
  table.heat th.rowlab {
    text-align: right;
    padding-right: 0.45rem;
    white-space: nowrap;
    font-size: 0.75rem;
    color: var(--ink);
  }
  table.heat th.ref, table.heat th.ref span { font-weight: 700; }
  table.heat td.cell {
    width: 2.85rem;
    min-width: 2.85rem;
    height: 1.7rem;
    text-align: center;
    font-size: 0.68rem;
    font-variant-numeric: tabular-nums;
    border-radius: 2px;
  }
  table.heat td.diag { background: #f4efe6; }
  .heat-scale {
    display: flex;
    flex-wrap: wrap;
    gap: 0.35rem 0.45rem;
    align-items: center;
    margin-top: 0.45rem;
    color: var(--muted);
    font-size: 0.75rem;
  }
  .heat-scale i {
    width: 1.7rem;
    height: 0.7rem;
    display: inline-block;
    border-radius: 2px;
  }
  .toolbar {
    display: flex;
    flex-wrap: wrap;
    gap: 0.8rem 1.2rem;
    align-items: center;
    padding: 0.85rem 1rem;
    border-bottom: 1px solid var(--line);
    background: #fbf8f2;
  }
  input[type="search"] {
    flex: 1 1 16rem;
    border: 1px solid var(--line);
    background: white;
    border-radius: 8px;
    padding: 0.55rem 0.7rem;
    font: inherit;
    color: var(--ink);
  }
  label.check { color: var(--muted); font-size: 0.88rem; }
  .count { margin-left: auto; color: var(--muted); font-size: 0.82rem; }
  .scroller { overflow-x: auto; }
  table { width: 100%; border-collapse: collapse; }
  th, td { text-align: left; padding: 0.65rem 0.7rem; vertical-align: top; }
  th {
    position: sticky;
    top: 0;
    background: #f7f3eb;
    font-size: 0.75rem;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    color: var(--muted);
    font-weight: 600;
    border-bottom: 1px solid var(--line);
    white-space: nowrap;
  }
  th button {
    all: unset;
    cursor: pointer;
  }
  th button:hover { color: var(--ink); }
  tbody tr.family { border-top: 1px solid var(--line); cursor: pointer; }
  tbody tr.family:hover { background: #fbf8f2; }
  tbody tr.family.open { background: var(--leaf); }
  .og { font-variant-numeric: tabular-nums; font-weight: 600; }
  .name { font-family: "Iowan Old Style", Palatino, Georgia, serif; }
  .muted { color: var(--muted); }
  .idcell {
    text-align: center;
    font-variant-numeric: tabular-nums;
    border-radius: 6px;
    min-width: 3.2rem;
  }
  .badges { display: flex; gap: 0.3rem; flex-wrap: wrap; }
  .badge {
    display: inline-block;
    border-radius: 999px;
    padding: 0.1rem 0.45rem;
    font-size: 0.72rem;
    letter-spacing: 0.03em;
    background: var(--sand);
    color: var(--ink);
  }
  .badge.sp { background: #f3ddd4; color: var(--rust); }
  .badge.tm { background: #e3eaf3; color: #23415f; }
  .badge.score { background: var(--green); color: #f6f1e6; }
  .detail td { background: #fbf8f2; border-top: none; cursor: default; }
  .detail-grid {
    display: grid;
    grid-template-columns: 9rem 1fr;
    gap: 0.35rem 0.8rem;
    padding: 0.3rem 0.2rem 0.6rem;
  }
  .detail-grid dt { color: var(--muted); font-size: 0.78rem; padding-top: 0.15rem; }
  .detail-grid dd { margin: 0; }
  .chips { display: flex; flex-wrap: wrap; gap: 0.3rem; }
  .chip {
    background: white;
    border: 1px solid var(--line);
    border-radius: 999px;
    padding: 0.05rem 0.45rem;
    font-size: 0.75rem;
  }
  .detail-grid dd { min-width: 0; }
  .members { font-variant-numeric: tabular-nums; font-size: 0.86rem; overflow-wrap: anywhere; }
  .legend {
    color: var(--muted);
    font-size: 0.82rem;
    padding: 0.2rem 0 1.5rem;
  }
  .empty { padding: 1.5rem; color: var(--muted); }
  .nav { display: flex; flex-wrap: wrap; gap: 1.1rem; margin: 1.15rem 0 0; }
  .nav a { color: #f6f1e6; text-decoration: none; border-bottom: 1px solid transparent; }
  .nav a[aria-current="page"] { border-bottom-color: #f6f1e6; }
  @media (max-width: 800px) {
    .stats { grid-template-columns: 1fr 1fr; }
    h1 { font-size: 1.8rem; }
  }
</style>
</head>
<body>
<header class="hero">
  <div class="wrap">
    <p class="kicker">GraftingGeneFamilies</p>
    <h1>Curated gene families</h1>
    <p class="lede">Orthogroups that contain the genes of interest, with percent identity to the reference, functional annotation, and graft-phenotype priority when that step was run.</p>
    <nav class="nav" aria-label="Summary pages">
      <a href="index.html" aria-current="page">Families</a>
      <a href="genomes.html">Genomes</a>
      <a href="protein_alignments.html">Protein alignments</a>
      <a href="codon_alignments.html">Codon alignments</a>
    </nav>
  </div>
</header>
<main class="wrap">
  <section class="stats" id="stats"></section>
  <section class="panel" id="sources-panel">
    <h2>Proteome sources</h2>
    <p class="source-note">Each remote file is recorded with the organism and annotation the source reported, the sequence count, and a checksum. A warning means that report does not match the species configured for the run.</p>
    <div class="scroller">
      <table>
        <thead>
          <tr>
            <th>ID</th><th>Configured</th><th>Source</th><th>Remote organism</th>
            <th>Accession or URL</th><th>Annotation</th><th>Proteins</th><th>Expected</th><th>Warning</th>
          </tr>
        </thead>
        <tbody id="sources"></tbody>
      </table>
    </div>
  </section>
  <section class="panel">
    <div class="toolbar">
      <input id="q" type="search" placeholder="Search gene, orthogroup, description, domain" aria-label="Search families">
      <label class="check"><input id="signals" type="checkbox"> Signal peptide</label>
      <label class="check" id="ranked-wrap"><input id="ranked" type="checkbox"> Priority score &gt; 0</label>
      <span class="count" id="count"></span>
    </div>
    <div class="scroller">
      <table>
        <thead id="head"></thead>
        <tbody id="body"></tbody>
      </table>
    </div>
  </section>
  <p class="legend" id="legend"></p>
</main>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
const DATA = JSON.parse(document.getElementById("payload").textContent);
const state = {
  q: "",
  signals: false,
  ranked: false,
  sort: DATA.has_ranking ? "priority_score" : "max_id",
  dir: "desc",
  open: null
};

const q = document.getElementById("q");
const signals = document.getElementById("signals");
const ranked = document.getElementById("ranked");
if (!DATA.has_ranking) document.getElementById("ranked-wrap").hidden = true;

q.addEventListener("input", () => { state.q = q.value.trim().toLowerCase(); draw(); });
signals.addEventListener("change", () => { state.signals = signals.checked; draw(); });
ranked.addEventListener("change", () => { state.ranked = ranked.checked; draw(); });

function isSignal(s) {
  if (!s) return false;
  const t = s.trim().toUpperCase();
  if (!t || t === "OTHER" || t.startsWith("OTHER") || t === "NO_SP" || t === "NONE") return false;
  return t.startsWith("SP") || t.includes("SIGNAL");
}
function isTm(s) {
  if (!s) return false;
  const t = s.trim().toLowerCase();
  if (!t || t === "glob" || t === "globular" || t === "none" || t === "no") return false;
  return t.includes("tm") || t.includes("transmem");
}
function fmt(n, digits) {
  if (n === null || n === undefined || Number.isNaN(n)) return "";
  return Number(n).toLocaleString(undefined, { maximumFractionDigits: digits, minimumFractionDigits: digits });
}
function idStyle(n) {
  if (n === null || n === undefined) return "background:#efeae2;color:#8a8378";
  const t = Math.max(0, Math.min(1, (n - 35) / 65));
  const r = Math.round(244 + (27 - 244) * t);
  const g = Math.round(236 + (67 - 236) * t);
  const b = Math.round(220 + (50 - 220) * t);
  const color = t > 0.55 ? "#f6f1e6" : "#1c1915";
  return `background:rgb(${r},${g},${b});color:${color}`;
}
function haystack(row) {
  return [
    row.orthogroup, row.genes, row.name, row.description, row.interpro,
    row.signalp, row.deeptmhmm, row.members.join(" "), row.gos.join(" ")
  ].join(" ").toLowerCase();
}
function valueFor(row, key) {
  if (key.startsWith("id:")) return row.identity[key.slice(3)];
  return row[key];
}
function compare(a, b) {
  const va = valueFor(a, state.sort);
  const vb = valueFor(b, state.sort);
  const dir = state.dir === "asc" ? 1 : -1;
  const aMiss = va === null || va === undefined || va === "";
  const bMiss = vb === null || vb === undefined || vb === "";
  if (aMiss && bMiss) return a.orthogroup.localeCompare(b.orthogroup);
  if (aMiss) return 1;
  if (bMiss) return -1;
  if (typeof va === "number" && typeof vb === "number") return (va - vb) * dir;
  return String(va).localeCompare(String(vb)) * dir;
}
function filtered() {
  return DATA.families.filter(row => {
    if (state.q && !haystack(row).includes(state.q)) return false;
    if (state.signals && !isSignal(row.signalp)) return false;
    if (state.ranked && !(row.priority_score > 0)) return false;
    return true;
  }).sort(compare);
}
function columns() {
  const cols = [
    { key: "orthogroup", label: "Orthogroup" },
    { key: "genes", label: "Genes of interest" },
    { key: "name", label: "Name" }
  ];
  DATA.species.forEach(sp => cols.push({
    key: "id:" + sp.id,
    label: sp.id,
    title: "Closest protein to the reference in " + sp.label
  }));
  cols.push({ key: "max_id", label: "Best %ID", title: "Closest protein to the reference in this family" });
  cols.push({ key: "features", label: "Features", nosort: true });
  if (DATA.has_ranking) cols.push({ key: "priority_score", label: "Priority" });
  if (DATA.has_folding) cols.push({ key: "fold_rank", label: "Fold rank" });
  return cols;
}
function drawHead() {
  const head = document.getElementById("head");
  const cells = columns().map(col => {
    if (col.nosort) return `<th>${col.label}</th>`;
    const mark = state.sort === col.key ? (state.dir === "asc" ? " ↑" : " ↓") : "";
    const title = col.title ? ` title="${col.title}"` : "";
    return `<th><button type="button" data-sort="${col.key}"${title}>${col.label}${mark}</button></th>`;
  }).join("");
  head.innerHTML = `<tr><th></th>${cells}</tr>`;
  head.querySelectorAll("button").forEach(btn => {
    btn.addEventListener("click", () => {
      const key = btn.dataset.sort;
      if (state.sort === key) state.dir = state.dir === "asc" ? "desc" : "asc";
      else { state.sort = key; state.dir = key === "orthogroup" || key === "genes" || key === "name" ? "asc" : "desc"; }
      draw();
    });
  });
}
function heatColor(n) {
  const t = Math.max(0, Math.min(1, (Number(n) - 55) / 40));
  const r = Math.round(246 + (110 - 246) * t);
  const g = Math.round(228 + (42 - 228) * t);
  const b = Math.round(210 + (18 - 210) * t);
  return { background: `rgb(${r},${g},${b})`, color: t > 0.58 ? "#f8f1e8" : "#1c1915" };
}
function matrixHtml(row) {
  const m = row.matrix;
  if (!m || !m.labels || !m.labels.length) return "";
  const head = m.labels.map((label, i) => {
    const mark = m.reference[i] ? " *" : "";
    return `<th class="col${m.reference[i] ? " ref" : ""}" title="${esc(m.full[i])}"><span>${esc(label)}${mark}</span></th>`;
  }).join("");
  const body = m.labels.map((label, i) => {
    const mark = m.reference[i] ? " *" : "";
    const cells = m.values[i].map((val, j) => {
      if (i === j) return `<td class="cell diag"></td>`;
      if (val === null || val === undefined) return `<td class="cell"></td>`;
      const style = heatColor(val);
      const title = `${m.labels[i]} vs ${m.labels[j]}: ${fmt(val, 1)}%`;
      return `<td class="cell" style="background:${style.background};color:${style.color}" title="${esc(title)}">${fmt(val, 1)}</td>`;
    }).join("");
    return `<tr><th class="rowlab${m.reference[i] ? " ref" : ""}" title="${esc(m.full[i])}">${esc(label)}${mark}</th>${cells}</tr>`;
  }).join("");
  const scale = [60, 70, 80, 90, 100].map(n => {
    const style = heatColor(n);
    return `<i style="background:${style.background}"></i><span>${n}%</span>`;
  }).join("");
  return `<div class="heat-wrap">
    <p class="heat-note">Pairwise identity for every protein in this alignment. The characterized protein is marked * and listed first. Homologs from the same species each have their own row.</p>
    <table class="heat"><thead><tr><th></th>${head}</tr></thead><tbody>${body}</tbody></table>
    <div class="heat-scale">${scale}</div>
  </div>`;
}
function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g, ch => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
  }[ch]));
}
function badges(row) {
  const bits = [];
  if (isSignal(row.signalp)) bits.push('<span class="badge sp">Signal peptide</span>');
  if (isTm(row.deeptmhmm)) bits.push('<span class="badge tm">TM helix</span>');
  return bits.join("") || '<span class="muted">—</span>';
}
function detail(row) {
  const items = [];
  const add = (label, html) => { if (html) items.push(`<dt>${label}</dt><dd>${html}</dd>`); };
  const characterized = [row.reference_gene, row.reference_species].filter(Boolean).join(" · ");
  add("Characterized gene", esc(characterized));
  add("Description", esc(row.description));
  if (row.gos.length) add("GO", `<div class="chips">${row.gos.map(g => `<span class="chip">${esc(g)}</span>`).join("")}</div>`);
  add("InterPro", esc(row.interpro));
  add("SignalP", esc(row.signalp));
  add("DeepTMHMM", esc(row.deeptmhmm));
  add("Reference", esc(row.reference));
  add("Matched queries", esc(row.matched));
  if (row.n_members !== null) add("Members", `${fmt(row.n_members, 0)} proteins in ${fmt(row.n_species, 0)} species`);
  if (row.members.length) add("Protein IDs", `<div class="members">${row.members.map(esc).join(", ")}</div>`);
  add("Annotation from", esc(row.annotation_from));
  if (DATA.has_ranking || row.priority_score !== null) {
    add("Priority score", row.priority_score === null ? "" : fmt(row.priority_score, 0));
    add("ASR on shift branches", row.n_asr === null ? "" : fmt(row.n_asr, 0));
    add("Example substitutions", esc(row.shift_subs));
    add("Divergence beyond relatedness", esc(row.divergence));
    add("Partial Mantel p", row.mantel_p === null ? "" : fmt(row.mantel_p, 3));
    add("HyPhy hits", row.n_selection === null ? "" : fmt(row.n_selection, 0));
    add("Convergence tests", row.n_convergence === null ? "" : fmt(row.n_convergence, 0));
  }
  if (row.fold_rank !== null) add("AlphaFold queue", `fold rank ${fmt(row.fold_rank, 0)}`);
  const min = row.min_id === null ? "—" : fmt(row.min_id, 1);
  const max = row.max_id === null ? "—" : fmt(row.max_id, 1);
  add("% identity to reference", `Closest protein is ${max}%. The most distant protein is ${min}%. Each species column is that closest protein, and the matrix below keeps every homolog.`);
  return `<tr class="detail"><td></td><td colspan="${columns().length}"><dl class="detail-grid">${items.join("")}</dl>${matrixHtml(row)}</td></tr>`;
}
function draw() {
  drawHead();
  const rows = filtered();
  const body = document.getElementById("body");
  if (!rows.length) {
    body.innerHTML = `<tr><td class="empty" colspan="${columns().length + 1}">No families match.</td></tr>`;
  } else {
    body.innerHTML = rows.map(row => {
      const open = state.open === row.orthogroup;
      const ids = DATA.species.map(sp => {
        const n = row.identity[sp.id];
        const text = n === null || n === undefined ? "·" : fmt(n, 0);
        return `<td class="idcell" style="${idStyle(n)}" title="${esc(sp.label)}">${text}</td>`;
      }).join("");
      const priority = row.priority_score === null ? ""
        : (row.priority_score > 0
          ? `<span class="badge score">${fmt(row.priority_score, 0)}</span>`
          : `<span class="muted">${fmt(row.priority_score, 0)}</span>`);
      const fold = row.fold_rank === null ? "" : fmt(row.fold_rank, 0);
      const extra = (DATA.has_ranking ? `<td>${priority}</td>` : "") + (DATA.has_folding ? `<td>${fold}</td>` : "");
      return `<tr class="family${open ? " open" : ""}" data-og="${esc(row.orthogroup)}">
        <td>${open ? "▾" : "▸"}</td>
        <td class="og">${esc(row.orthogroup)}</td>
        <td>${esc(row.genes)}</td>
        <td class="name">${esc(row.name) || '<span class="muted">—</span>'}</td>
        ${ids}
        <td class="idcell" style="${idStyle(row.max_id)}" title="Closest protein to the reference">${row.max_id === null ? "·" : fmt(row.max_id, 0)}</td>
        <td><div class="badges">${badges(row)}</div></td>
        ${extra}
      </tr>${open ? detail(row) : ""}`;
    }).join("");
    body.querySelectorAll("tr.family").forEach(tr => {
      tr.addEventListener("click", () => {
        state.open = state.open === tr.dataset.og ? null : tr.dataset.og;
        draw();
      });
    });
  }
  document.getElementById("count").textContent = `${rows.length} of ${DATA.families.length}`;
}
const summary = DATA.summary || {};
const cards = [
  ["Genes mined", summary.genes_mined, "Literature genes placed in an orthogroup"],
  ["Species in run", summary.species_in_run, "Species whose proteomes were included in this job"],
  ["Proteomes downloaded", summary.proteomes_downloaded == null ? "—" : summary.proteomes_downloaded, "Remote NCBI or URL protein downloads. A local file is not counted."],
  ["Run time", summary.runtime || "—", "Time from job start until this report was written"]
];
document.getElementById("stats").innerHTML = cards.map(([label, value, tip]) =>
  `<article class="stat" title="${esc(tip)}"><b>${esc(value)}</b><span>${esc(label)}</span></article>`
).join("");
const bits = [`Report date ${esc(DATA.generated)}.`];
bits.push("Genes mined are the literature genes placed in an orthogroup. Species in the run are the proteomes included in this job. Proteomes downloaded counts remote NCBI or URL fetches, not a local file. Run time is from job start until this report.");
bits.push("Each species column is the closest protein from that species to the reference, and Best %ID is the single closest protein in the family. Darker green is higher identity.");
if (DATA.has_ranking) bits.push("Priority adds ASR substitutions on phenotype-shift branches, divergence beyond species relatedness, HyPhy hits, and convergence.");
if (DATA.has_folding) bits.push("Fold rank is the AlphaFold queue order (1 is first).");
bits.push("TSV copies are in the tables folder beside this file.");
document.getElementById("legend").textContent = bits.join(" ");
function renderSources() {
  const rows = DATA.provenance || [];
  const panel = document.getElementById("sources-panel");
  if (!rows.length) {
    panel.hidden = true;
    return;
  }
  document.getElementById("sources").innerHTML = rows.map(row => {
    const ann = [row.annotation_name, row.annotation_date].filter(Boolean).join(" · ");
    return `<tr class="${row.warning ? "warn" : ""}">
      <td>${esc(row.species_id)}</td>
      <td>${esc(row.configured_taxon)}</td>
      <td>${esc(row.source)}</td>
      <td>${esc(row.remote_organism) || "—"}</td>
      <td class="url">${esc(row.proteins_source)}</td>
      <td>${esc(ann) || "—"}</td>
      <td>${esc(row.n_proteins)}</td>
      <td>${esc(row.expected_protein_coding) || "—"}</td>
      <td>${esc(row.warning) || ""}</td>
    </tr>`;
  }).join("");
}
renderSources();
draw();
</script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    table_dir = outdir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    families = read_tsv(args.families)
    wide = read_tsv(args.wide_identity)
    ranking = read_tsv(args.ranking)
    folding: list[dict] = []
    af_dir = Path(args.alphafold_dir) if args.alphafold_dir else None
    if af_dir and af_dir.is_dir():
        folding = read_tsv(str(af_dir / "og_priority.tsv"))

    species = species_from_wide(wide)
    records = build_records(families, wide, ranking, folding, species)
    provenance = load_provenance(args.provenance)
    n_run, n_downloaded = species_run_counts(provenance, species)
    matrices = load_matrices(args.identity_matrix_pairs, args.identity_matrix_labels)
    for rec in records:
        rec["matrix"] = matrices.get(rec["orthogroup"])
    payload = {
        "generated": date.today().isoformat(),
        "species": species,
        "has_ranking": bool(ranking),
        "has_folding": bool(folding),
        "provenance": provenance,
        "summary": {
            "genes_mined": mined_gene_count(records),
            "species_in_run": n_run,
            "proteomes_downloaded": n_downloaded,
            "runtime": format_runtime(args.runtime_seconds),
        },
        "families": records,
    }
    html = HTML.replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False).replace("<", "\\u003c"))
    (outdir / "index.html").write_text(html)

    copied = []

    def copy_table(src: str | None, name: str, rows: list[dict]) -> None:
        if not rows or not src:
            return
        shutil.copyfile(src, table_dir / name)
        copied.append(name)

    copy_table(args.families, "gene_families_of_interest.tsv", families)
    copy_table(args.wide_identity, "gene_families_pctid_by_species.tsv", wide)
    copy_table(args.ranking, "graft_candidate_ranking.tsv", ranking)
    if folding and af_dir:
        copy_table(str(af_dir / "og_priority.tsv"), "alphafold_og_priority.tsv", folding)
    if provenance:
        dest = table_dir / "download_provenance.tsv"
        with dest.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(provenance[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(provenance)
        copied.append("download_provenance.tsv")
    copy_table(args.identity_matrix_pairs, "identity_matrix_pairs.tsv", read_plain_tsv(args.identity_matrix_pairs) if args.identity_matrix_pairs else [])
    copy_table(args.identity_matrix_labels, "identity_matrix_labels.tsv", read_plain_tsv(args.identity_matrix_labels) if args.identity_matrix_labels else [])

    write_readme(outdir / "README.txt", copied)
    print(f"[curate_results] {len(records)} families → {outdir / 'index.html'}")


if __name__ == "__main__":
    main()
