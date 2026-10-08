#!/usr/bin/env python3
"""Extract unified sequence features for structure mapping.

Inputs: cleaned protein FASTA(s), SignalP, InterProScan TSV, DeepTMHMM TSV,
optional eggNOG annotations, optional functional_feature_rules.yaml.

Output: sequence_features.tsv with 1-based sequence coordinates.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proteins", nargs="+", required=True)
    p.add_argument("--rules", default=None, help="functional_feature_rules.yaml")
    p.add_argument("--signalp", default=None)
    p.add_argument("--interproscan", default=None)
    p.add_argument("--deeptmhmm", default=None)
    p.add_argument("--eggnog", default=None)
    p.add_argument(
        "--transferred-sites",
        default=None,
        help="transferred_sites.tsv from transfer_curated_sites.py (MSA-mapped)",
    )
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-classifications", required=True)
    return p.parse_args()


def load_rules(path: str | None) -> dict:
    if not path or not Path(path).exists() or yaml is None:
        return {"classes": {}, "glycosylation_motif": "N[^P][ST]", "curated_sites": {}}
    return yaml.safe_load(Path(path).read_text()) or {}


def read_commented_tsv(path: str | None) -> list[dict]:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return []
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []
    import io

    return list(csv.DictReader(io.StringIO("\n".join(lines) + "\n"), delimiter="\t"))


def classify_protein(text: str, classes: dict) -> list[str]:
    t = (text or "").lower()
    hits = []
    for cls, spec in (classes or {}).items():
        for kw in spec.get("keywords") or []:
            if kw.lower() in t:
                hits.append(cls)
                break
    return hits


def parse_interpro_rows(path: str | None) -> list[dict]:
    """Return feature-like dicts from InterProScan TSV (headerless or headed)."""
    if not path or not Path(path).exists():
        return []
    feats = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            prot, md5, length, analysis, sig_acc, sig_desc, start, end = parts[:8]
            ipr = parts[11] if len(parts) > 11 else ""
            ipr_desc = parts[12] if len(parts) > 12 else ""
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            feats.append(
                {
                    "protein_id": prot,
                    "feature_type": "interpro_domain",
                    "feature_id": ipr or sig_acc,
                    "feature_name": ipr_desc or sig_desc,
                    "start": start_i,
                    "end": end_i,
                    "source": f"InterProScan:{analysis}",
                    "note": sig_desc,
                }
            )
            # Active site / binding site analyses often named explicitly
            blob = f"{analysis} {sig_desc} {ipr_desc}".lower()
            if "active" in blob or "binding" in blob or "catalytic" in blob:
                feats.append(
                    {
                        "protein_id": prot,
                        "feature_type": "active_site",
                        "feature_id": sig_acc,
                        "feature_name": sig_desc,
                        "start": start_i,
                        "end": end_i,
                        "source": f"InterProScan:{analysis}",
                        "note": ipr_desc,
                    }
                )
    return feats


def parse_signalp(path: str | None) -> list[dict]:
    rows = read_commented_tsv(path)
    feats = []
    for r in rows:
        pid = r.get("ID") or r.get("# ID") or r.get("name") or ""
        if not pid:
            # try first column
            pid = list(r.values())[0] if r else ""
        pred = (r.get("Prediction") or r.get("SP") or r.get("prediction") or "").lower()
        # SignalP 6 often has CS Position like "25-26" or "cleavage site: XX"
        cs = r.get("CS Position") or r.get("Cleaveage site") or r.get("cleavage_site") or ""
        end = None
        m = re.search(r"(\d+)\s*[-–]\s*(\d+)", str(cs))
        if m:
            end = int(m.group(1))
        elif pred.startswith("sp") or "signal" in pred:
            # fallback unknown length — mark 1–N unknown; skip if no CS
            continue
        if pid and end:
            feats.append(
                {
                    "protein_id": pid.split()[0],
                    "feature_type": "signal_peptide",
                    "feature_id": "SignalP",
                    "feature_name": pred or "signal_peptide",
                    "start": 1,
                    "end": end,
                    "source": "SignalP",
                    "note": str(cs),
                }
            )
    return feats


def parse_deeptmhmm(path: str | None) -> list[dict]:
    rows = read_commented_tsv(path)
    feats = []
    for r in rows:
        pid = r.get("id") or r.get("ID") or (list(r.values())[0] if r else "")
        topo = r.get("topology") or r.get("predicted_topology") or r.get("prediction") or ""
        if not pid or not topo:
            continue
        # Topology strings vary; record whole chain as annotation span unknown
        # Prefer explicit TM helix runs of 'M' if present
        for m in re.finditer(r"M+", str(topo)):
            feats.append(
                {
                    "protein_id": str(pid).split()[0],
                    "feature_type": "transmembrane",
                    "feature_id": "TM",
                    "feature_name": "TM_helix",
                    "start": m.start() + 1,
                    "end": m.end(),
                    "source": "DeepTMHMM",
                    "note": "",
                }
            )
        if not re.search(r"M+", str(topo)):
            feats.append(
                {
                    "protein_id": str(pid).split()[0],
                    "feature_type": "topology",
                    "feature_id": "DeepTMHMM",
                    "feature_name": "topology",
                    "start": 1,
                    "end": len(str(topo)) or 1,
                    "source": "DeepTMHMM",
                    "note": str(topo)[:200],
                }
            )
    return feats


def glycosylation_sequons(protein_id: str, seq: str, motif: str) -> list[dict]:
    feats = []
    for m in re.finditer(motif, seq.upper()):
        feats.append(
            {
                "protein_id": protein_id,
                "feature_type": "glycosylation_sequon",
                "feature_id": "NXS/T",
                "feature_name": m.group(0),
                "start": m.start() + 1,
                "end": m.end(),
                "source": "motif",
                "note": "N-X-S/T sequon (predicted)",
            }
        )
    return feats


def main() -> None:
    args = parse_args()
    rules = load_rules(args.rules)
    classes = rules.get("classes") or {}
    motif = rules.get("glycosylation_motif") or "N[^P][ST]"

    seqs = {}
    for fa in args.proteins:
        for rec in SeqIO.parse(fa, "fasta"):
            seqs[rec.id] = str(rec.seq).upper().rstrip("*")

    # eggNOG text per protein for classification
    eggnog_text = defaultdict(str)
    for r in read_commented_tsv(args.eggnog):
        q = r.get("query_name") or r.get("#query") or r.get("query") or ""
        if not q:
            continue
        eggnog_text[q] = " ".join(
            [
                r.get("Preferred_name") or "",
                r.get("Description") or "",
                r.get("COG_category") or "",
            ]
        )

    feats = []
    feats.extend(parse_signalp(args.signalp))
    feats.extend(parse_interpro_rows(args.interproscan))
    feats.extend(parse_deeptmhmm(args.deeptmhmm))

    # InterPro text blob per protein for classification
    ipr_text = defaultdict(str)
    for f in feats:
        if f["feature_type"] in {"interpro_domain", "active_site"}:
            ipr_text[f["protein_id"]] += " " + (f.get("feature_name") or "") + " " + (f.get("note") or "")

    for pid, seq in seqs.items():
        feats.extend(glycosylation_sequons(pid, seq, motif))

    # Curated sites on the REFERENCE only (raw numbers). Prefer MSA transfer via
    # transfer_curated_sites.py → --transferred-sites for orthologs.
    curated = rules.get("curated_sites") or []
    curated_entries: list[dict] = []
    if isinstance(curated, list):
        curated_entries = curated
    elif isinstance(curated, dict):
        for pid, sites in curated.items():
            for s in sites or []:
                e = dict(s)
                e["reference_id"] = pid
                curated_entries.append(e)
    for e in curated_entries:
        ref = e.get("reference_id")
        if not ref or ref not in seqs:
            continue
        for res in e.get("residues") or []:
            feats.append(
                {
                    "protein_id": ref,
                    "feature_type": e.get("feature", "curated_site"),
                    "feature_id": "curated",
                    "feature_name": e.get("note") or e.get("feature", "curated"),
                    "start": int(res),
                    "end": int(res),
                    "source": "curated_reference",
                    "note": (e.get("note") or "")
                    + " [reference only — use transfer_curated_sites for orthologs]",
                }
            )

    # MSA-transferred curated sites (from transfer_curated_sites.py)
    if args.transferred_sites and Path(args.transferred_sites).exists():
        with open(args.transferred_sites) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("status") != "ok":
                    continue
                tid = row.get("target_id") or ""
                try:
                    pos = int(row["target_resnum"])
                except (KeyError, ValueError, TypeError):
                    continue
                feats.append(
                    {
                        "protein_id": tid,
                        "feature_type": row.get("feature_type") or "curated_site",
                        "feature_id": "curated_msa",
                        "feature_name": row.get("note")
                        or row.get("feature_type")
                        or "curated",
                        "start": pos,
                        "end": pos,
                        "source": "curated_msa_transfer",
                        "note": (
                            f"ref={row.get('reference_id')}:"
                            f"{row.get('ref_resnum')}→col{row.get('msa_column')} "
                            f"{row.get('note') or ''}"
                        ).strip(),
                    }
                )

    # classifications
    class_rows = []
    for pid in seqs:
        blob = f"{eggnog_text.get(pid, '')} {ipr_text.get(pid, '')}"
        hits = classify_protein(blob, classes)
        class_rows.append(
            {
                "protein_id": pid,
                "classes": ";".join(hits) if hits else "",
                "n_features": sum(1 for f in feats if f["protein_id"] == pid),
            }
        )
        for cls in hits:
            feats.append(
                {
                    "protein_id": pid,
                    "feature_type": "protein_class",
                    "feature_id": cls,
                    "feature_name": (classes.get(cls) or {}).get("description") or cls,
                    "start": 1,
                    "end": len(seqs[pid]),
                    "source": "rules",
                    "note": "",
                }
            )

    fields = [
        "protein_id",
        "feature_type",
        "feature_id",
        "feature_name",
        "start",
        "end",
        "source",
        "note",
    ]
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(feats)

    with open(args.out_classifications, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["protein_id", "classes", "n_features"], delimiter="\t")
        w.writeheader()
        w.writerows(class_rows)

    print(
        f"[extract_sequence_features] proteins={len(seqs)} features={len(feats)} "
        f"classified={sum(1 for r in class_rows if r['classes'])} → {args.out_tsv}"
    )


if __name__ == "__main__":
    main()
