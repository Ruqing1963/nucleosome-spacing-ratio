"""
analyze_robustness.py
---------------------
Robustness of the shuffle-invariant (scalar, non-repulsive) verdict for the
S. cerevisiae unique map (reproduces Table 2):
  * excluding large nucleosome-free-region gaps to isolate the bulk array
    (gaps > 300 / 250 / 220 bp);
  * randomly thinning the array by keeping only the highest-confidence
    nucleosomes (top 50 / 25 / 10% by NCP-score-to-noise ratio).
In every case <r>_obs ~= <r>_shuf: the VALUE moves with the marginal, but the
scalar (non-repulsive) character does not.

Usage:
    python analyze_robustness.py --map data/raw/cerevisiae_unique_map.txt
"""
import argparse
import csv
import os
import numpy as np

from stats_common import (load_map, spacings_by_chrom, r_pooled, r_shuffled,
                          SC_CHROMS)


def _sub(rows, gapcap=None, conf_pct=None):
    if conf_pct is not None:
        sn = np.array([r[3] for r in rows])
        thr = np.percentile(sn, conf_pct)
        rows = [r for r in rows if r[3] >= thr]
    spac = spacings_by_chrom(rows, SC_CHROMS, gapcap=gapcap)
    r_obs, n = r_pooled(spac)
    r_shuf, _ = r_shuffled(spac, n_perm=120)
    return n, r_obs, r_shuf


def main(map_path, outdir):
    rows = load_map(map_path)
    table = [("all unique (baseline)", _sub(rows))]
    for cap in (300, 250, 220):
        table.append((f"bulk array, exclude gaps >{cap} bp", _sub(rows, gapcap=cap)))
    for pct, lab in ((50, "top 50% confidence (thinned)"),
                     (75, "top 25% confidence (thinned)"),
                     (90, "top 10% confidence (thinned)")):
        table.append((lab, _sub(rows, conf_pct=pct)))

    os.makedirs(os.path.join(outdir, "results"), exist_ok=True)
    out = os.path.join(outdir, "results", "robustness.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["treatment", "n_pairs", "r_obs", "r_shuf"])
        for lab, (n, ro, rs) in table:
            w.writerow([lab, n, round(ro, 4), round(rs, 4)])
            print(f"{lab:40s} N={n:6d}  <r>_obs={ro:.4f}  <r>_shuf={rs:.4f}")
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", default="data/raw/cerevisiae_unique_map.txt")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    main(a.map, a.outdir)
