#!/usr/bin/env python3
"""Test whether protein divergence predicts graft failure beyond species relatedness.

For each orthogroup (and pooled):
  - Point-biserial / Mann-Whitney: divergence vs success
  - Partial Mantel (permutation): phenotype distance ~ protein distance | species distance

Phenotype distance: 0 if same success class, 1 otherwise (per directed pair, symmetrized mean).
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair-divergence", required=True)
    p.add_argument("--species-distance", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--n-perm", type=int, default=999)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_species_distance(path: str) -> dict[tuple[str, str], float]:
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0][1:]
    dist = {}
    for r in rows[1:]:
        a = r[0]
        for j, b in enumerate(header):
            if j + 1 >= len(r) or r[j + 1] == "":
                continue
            dist[(a, b)] = float(r[j + 1])
    return dist


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denx = math.sqrt(sum((a - mx) ** 2 for a in x))
    deny = math.sqrt(sum((b - my) ** 2 for b in y))
    if denx == 0 or deny == 0:
        return float("nan")
    return num / (denx * deny)


def partial_mantel(
    d_y: list[float], d_x: list[float], d_z: list[float], n_perm: int, rng: random.Random
) -> tuple[float, float]:
    """Partial correlation of vectorized distances Y ~ X | Z with permutations of Y."""
    # residualize Y and X on Z
    def resid(a, b):
        r = pearson(a, b)
        if math.isnan(r):
            return a[:]
        # a_pred = mean_a + r * (sd_a/sd_b) * (b-mean_b) — use simple projection
        mb = sum(b) / len(b)
        ma = sum(a) / len(a)
        sb2 = sum((v - mb) ** 2 for v in b)
        if sb2 == 0:
            return [v - ma for v in a]
        cov = sum((u - ma) * (v - mb) for u, v in zip(a, b))
        beta = cov / sb2
        return [u - (ma + beta * (v - mb)) for u, v in zip(a, b)]

    ry = resid(d_y, d_z)
    rx = resid(d_x, d_z)
    obs = pearson(ry, rx)
    if math.isnan(obs):
        return obs, float("nan")

    count = 0
    for _ in range(n_perm):
        shuf = ry[:]
        rng.shuffle(shuf)
        r = pearson(shuf, rx)
        if not math.isnan(r) and abs(r) >= abs(obs):
            count += 1
    pval = (count + 1) / (n_perm + 1)
    return obs, pval


def mann_whitney_u(a: list[float], b: list[float]) -> float:
    """Two-sided asymptotic p-value (normal approx); returns p or nan."""
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    combined = [(v, 0) for v in a] + [(v, 1) for v in b]
    combined.sort(key=lambda t: t[0])
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg = (i + j + 1) / 2.0  # 1-based avg rank
        for k in range(i, j):
            ranks[k] = avg
        i = j
    r_a = sum(ranks[k] for k, (_v, g) in enumerate(combined) if g == 0)
    n1, n2 = len(a), len(b)
    u1 = r_a - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    mu = n1 * n2 / 2.0
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sigma == 0:
        return float("nan")
    z = (u - mu) / sigma
    # two-sided from normal
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return p


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)
    sp_dist = load_species_distance(args.species_distance)

    by_og: dict[str, list[dict]] = defaultdict(list)
    with open(args.pair_divergence) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") != "ok" or not row.get("pct_identity"):
                continue
            by_og[row["Orthogroup"]].append(row)

    out_rows = []

    def analyze(og: str, rows: list[dict]) -> dict:
        # divergence as (100 - pct_id) / 100
        succ = []
        fail = []
        for r in rows:
            div = (100.0 - float(r["pct_identity"])) / 100.0
            if r["success"] == "yes":
                succ.append(div)
            else:
                fail.append(div)

        # vectorize unique undirected pairs for Mantel-like test
        # use directed pairs as observations for simple correlation with species distance
        y = []  # phenotype distance: 0 success, 1 failure
        x = []  # protein divergence
        z = []  # species distance
        for r in rows:
            sc, rs = r["scion"], r["rootstock"]
            if (sc, rs) not in sp_dist and (rs, sc) not in sp_dist:
                continue
            sd = sp_dist.get((sc, rs), sp_dist.get((rs, sc)))
            y.append(0.0 if r["success"] == "yes" else 1.0)
            x.append((100.0 - float(r["pct_identity"])) / 100.0)
            z.append(sd)

        r_xy = pearson(x, y)
        r_zy = pearson(z, y)
        partial_r, partial_p = (
            partial_mantel(y, x, z, args.n_perm, rng) if len(x) >= 4 else (float("nan"), float("nan"))
        )
        mw_p = mann_whitney_u(fail, succ)

        return {
            "Orthogroup": og,
            "n_pairs": len(rows),
            "n_success": len(succ),
            "n_failure": len(fail),
            "mean_div_success": f"{(sum(succ)/len(succ)):.4f}" if succ else "",
            "mean_div_failure": f"{(sum(fail)/len(fail)):.4f}" if fail else "",
            "mannwhitney_p": f"{mw_p:.4g}" if not math.isnan(mw_p) else "",
            "pearson_div_vs_fail": f"{r_xy:.4f}" if not math.isnan(r_xy) else "",
            "pearson_spdist_vs_fail": f"{r_zy:.4f}" if not math.isnan(r_zy) else "",
            "partial_mantel_r": f"{partial_r:.4f}" if not math.isnan(partial_r) else "",
            "partial_mantel_p": f"{partial_p:.4g}" if not math.isnan(partial_p) else "",
            "candidate_beyond_relatedness": (
                "yes"
                if (not math.isnan(partial_p) and partial_p < 0.05 and not math.isnan(partial_r) and partial_r > 0)
                else "no"
            ),
        }

    for og, rows in sorted(by_og.items()):
        out_rows.append(analyze(og, rows))

    # pooled
    all_rows = [r for rows in by_og.values() for r in rows]
    if all_rows:
        pooled = analyze("POOLED", all_rows)
        out_rows.append(pooled)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "Orthogroup",
        "n_pairs",
        "n_success",
        "n_failure",
        "mean_div_success",
        "mean_div_failure",
        "mannwhitney_p",
        "pearson_div_vs_fail",
        "pearson_spdist_vs_fail",
        "partial_mantel_r",
        "partial_mantel_p",
        "candidate_beyond_relatedness",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    n_cand = sum(1 for r in out_rows if r["candidate_beyond_relatedness"] == "yes")
    print(f"[test_divergence_vs_phenotype] ogs={len(out_rows)} candidates={n_cand} → {args.out_tsv}")


if __name__ == "__main__":
    main()
