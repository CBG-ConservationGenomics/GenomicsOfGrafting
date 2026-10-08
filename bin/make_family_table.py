#!/usr/bin/env python3
"""Merge orthogroup selection, % identity, and annotations into a polished table."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selected-ogs", required=True)
    p.add_argument("--identity-summary", required=True)
    p.add_argument("--identity-pairwise", required=True)
    p.add_argument("--eggnog", default=None)
    p.add_argument("--interproscan", default=None)
    p.add_argument("--signalp", default=None)
    p.add_argument("--deeptmhmm", default=None)
    p.add_argument("--out-table", required=True)
    p.add_argument("--out-wide-identity", required=True)
    return p.parse_args()


def read_tsv(path: str | None) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        # skip comment lines for eggNOG etc.
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []
    # recreate reader from filtered lines
    import io

    return list(csv.DictReader(io.StringIO("".join(lines)), delimiter="\t"))


def index_eggnog(rows: list[dict]) -> dict[str, dict]:
    out = {}
    for r in rows:
        # eggNOG-mapper v2 columns vary; common: query_name / Description / COG_category / GOs
        q = r.get("query_name") or r.get("#query") or r.get("query") or ""
        if not q:
            continue
        out[q] = {
            "eggnog_description": r.get("Description") or r.get("eggNOG_OGs") or "",
            "eggnog_preferred_name": r.get("Preferred_name") or "",
            "eggnog_gos": r.get("GOs") or "",
            "eggnog_ec": r.get("EC") or "",
            "eggnog_kog": r.get("COG_category") or "",
        }
    return out


def index_interpro(rows: list[dict]) -> dict[str, str]:
    """Collapse InterProScan TSV accessions/descriptions per protein."""
    # InterProScan TSV has no header by default
    by_prot: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        # if headerless, DictReader may use first row as keys — handle both
        vals = list(r.values())
        keys = list(r.keys())
        if "Protein Accession" in r or "protein_id" in r:
            prot = r.get("Protein Accession") or r.get("protein_id")
            desc = r.get("Signature Description") or r.get("InterPro Description") or ""
            ipr = r.get("InterPro Accession") or ""
        else:
            # standard IPS TSV: 0=protein, 5=signature desc, 11=IPR, 12=IPR desc
            if len(vals) < 13:
                continue
            prot = vals[0]
            desc = vals[12] or vals[5]
            ipr = vals[11]
        if prot and (ipr or desc):
            by_prot[prot].add(f"{ipr}:{desc}" if ipr else desc)
    return {k: "; ".join(sorted(v))[:2000] for k, v in by_prot.items()}


def index_signalp(rows: list[dict]) -> dict[str, str]:
    out = {}
    for r in rows:
        q = r.get("ID") or r.get("# ID") or r.get("name") or list(r.values())[0]
        pred = r.get("Prediction") or r.get("SP") or r.get("prediction") or ""
        out[q] = pred
    return out


def index_deeptmhmm(rows: list[dict]) -> dict[str, str]:
    out = {}
    for r in rows:
        q = r.get("id") or r.get("ID") or r.get("query") or list(r.values())[0]
        # summarize topology if present
        topo = r.get("predicted_topology") or r.get("topology") or r.get("prediction") or ""
        if not topo and len(r) > 1:
            topo = list(r.values())[1]
        out[q] = topo
    return out


def annotate_for_members(members: list[str], eggnog, ipr, signalp, tmhmm) -> dict:
    """Pick annotation from first member that has any annotation, prefer reference-like order."""
    result = {
        "eggnog_description": "",
        "eggnog_preferred_name": "",
        "eggnog_gos": "",
        "interpro_domains": "",
        "signalp": "",
        "deeptmhmm": "",
        "annotation_source_id": "",
    }
    for m in members:
        hit = False
        if m in eggnog:
            result.update({k: v for k, v in eggnog[m].items() if k.startswith("eggnog")})
            hit = True
        if m in ipr and not result["interpro_domains"]:
            result["interpro_domains"] = ipr[m]
            hit = True
        if m in signalp and not result["signalp"]:
            result["signalp"] = signalp[m]
            hit = True
        if m in tmhmm and not result["deeptmhmm"]:
            result["deeptmhmm"] = tmhmm[m]
            hit = True
        if hit and not result["annotation_source_id"]:
            result["annotation_source_id"] = m
        if result["eggnog_description"] and result["interpro_domains"]:
            break
    return result


def main() -> None:
    args = parse_args()

    selected = read_tsv(args.selected_ogs)
    id_summary = {r["Orthogroup"]: r for r in read_tsv(args.identity_summary)}
    pairwise = read_tsv(args.identity_pairwise)

    eggnog = index_eggnog(read_tsv(args.eggnog))
    # InterProScan may be headerless — special-case
    ipr = {}
    if args.interproscan and Path(args.interproscan).exists():
        raw_rows = []
        with open(args.interproscan) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                raw_rows.append(
                    {
                        "Protein Accession": parts[0] if parts else "",
                        "Signature Description": parts[5] if len(parts) > 5 else "",
                        "InterPro Accession": parts[11] if len(parts) > 11 else "",
                        "InterPro Description": parts[12] if len(parts) > 12 else "",
                    }
                )
        ipr = index_interpro(raw_rows)

    signalp = index_signalp(read_tsv(args.signalp))
    tmhmm = index_deeptmhmm(read_tsv(args.deeptmhmm))

    # wide identity: OG x species mean %ID to reference
    by_og_sp: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in pairwise:
        try:
            by_og_sp[r["Orthogroup"]][r["target_species_id"]].append(float(r["pct_identity"]))
        except (KeyError, ValueError):
            continue
    species_ids = sorted({sp for og in by_og_sp.values() for sp in og})

    wide_rows = []
    polished = []

    for row in selected:
        og = row["Orthogroup"]
        members = [m.strip() for m in (row.get("members") or "").split(",") if m.strip()]
        ann = annotate_for_members(members, eggnog, ipr, signalp, tmhmm)
        idn = id_summary.get(og, {})

        polished.append(
            {
                "Orthogroup": og,
                "genes_of_interest": row.get("genes_of_interest", ""),
                "query_species_ids": row.get("query_species_ids", ""),
                "matched_members": row.get("matched_members", ""),
                "reference_id": idn.get("reference_id", ""),
                "n_members": row.get("n_members", ""),
                "n_species_with_members": row.get("n_species_with_members", ""),
                "mean_pct_id_to_ref": idn.get("mean_pct_id_to_ref", ""),
                "min_pct_id_to_ref": idn.get("min_pct_id_to_ref", ""),
                "max_pct_id_to_ref": idn.get("max_pct_id_to_ref", ""),
                **ann,
                "member_ids": row.get("members", ""),
            }
        )

        wide = {
            "Orthogroup": og,
            "genes_of_interest": row.get("genes_of_interest", ""),
            "query_species_ids": row.get("query_species_ids", ""),
            "reference_id": idn.get("reference_id", ""),
        }
        for sp in species_ids:
            vals = by_og_sp[og].get(sp, [])
            wide[f"mean_pct_id_{sp}"] = f"{(sum(vals) / len(vals)):.2f}" if vals else ""
        wide_rows.append(wide)

    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "Orthogroup",
        "genes_of_interest",
        "query_species_ids",
        "matched_members",
        "reference_id",
        "n_members",
        "n_species_with_members",
        "mean_pct_id_to_ref",
        "min_pct_id_to_ref",
        "max_pct_id_to_ref",
        "eggnog_preferred_name",
        "eggnog_description",
        "eggnog_gos",
        "interpro_domains",
        "signalp",
        "deeptmhmm",
        "annotation_source_id",
        "member_ids",
    ]
    with open(args.out_table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(polished)

    wide_fields = ["Orthogroup", "genes_of_interest", "query_species_ids", "reference_id"] + [
        f"mean_pct_id_{sp}" for sp in species_ids
    ]
    with open(args.out_wide_identity, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=wide_fields, delimiter="\t")
        w.writeheader()
        w.writerows(wide_rows)

    print(f"[make_family_table] wrote {len(polished)} orthogroup rows")


if __name__ == "__main__":
    main()
