"""
analyze_crossspecies.py
-----------------------
Cross-species comparison of S. cerevisiae and S. pombe unique maps:
  * writes results/r_summary.csv and results/r_per_chromosome.csv (Table 1);
  * writes data/spacing_histogram_pombe.csv;
  * reproduces Figure 2 (both species are shuffle-invariant; the tighter
    S. pombe marginal yields the higher <r>).

Usage:
    python analyze_crossspecies.py \
        --cerevisiae data/raw/cerevisiae_unique_map.txt \
        --pombe      data/raw/pombe_unique_map.txt
"""
import argparse
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from stats_common import (load_map, spacings_by_chrom, r_pooled, r_shuffled,
                          lag1_autocorr, all_spacings,
                          R_POISSON, R_GOE, R_GUE, SC_CHROMS, SP_CHROMS)


def summarize(rows, order):
    spac = spacings_by_chrom(rows, order)
    al = all_spacings(spac)
    r_obs, n = r_pooled(spac)
    r_shuf, r_sd = r_shuffled(spac)
    ac = lag1_autocorr(spac)
    return dict(spac=spac, all=al, n=n, r_obs=r_obs, r_shuf=r_shuf,
                r_sd=r_sd, ac=ac,
                median=int(np.median(al)), cv=float(al.std() / al.mean()))


def main(sc_path, sp_path, outdir):
    SC = summarize(load_map(sc_path), SC_CHROMS)
    SP = summarize(load_map(sp_path), SP_CHROMS)

    os.makedirs(os.path.join(outdir, "results"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "data"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "figures"), exist_ok=True)

    with open(os.path.join(outdir, "results", "r_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["species", "n_pairs", "r_obs", "r_shuf", "r_shuf_sd",
                    "lag1_autocorr", "median_spacing_bp", "CV"])
        for name, D in [("S_cerevisiae", SC), ("S_pombe", SP)]:
            w.writerow([name, D["n"], round(D["r_obs"], 4), round(D["r_shuf"], 4),
                        round(D["r_sd"], 4), round(D["ac"], 4),
                        D["median"], round(D["cv"], 3)])

    with open(os.path.join(outdir, "results", "r_per_chromosome.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["species", "chromosome", "n_pairs", "r_obs"])
        for name, D in [("S_cerevisiae", SC), ("S_pombe", SP)]:
            for c, d in D["spac"].items():
                ro, n = r_pooled({c: d})
                w.writerow([name, c, n, round(ro, 4)])

    h, e = np.histogram(SP["all"], bins=np.arange(80, 601, 2))
    ctr = (0.5 * (e[1:] + e[:-1])).astype(int)
    with open(os.path.join(outdir, "data", "spacing_histogram_pombe.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["spacing_bp_bin_center", "count"])
        for b, c in zip(ctr, h):
            w.writerow([int(b), int(c)])

    # ---- Figure 2 ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for v, lab, col in [(R_POISSON, f'Poisson {R_POISSON}', '#C0392B'),
                        (R_GOE, f'GOE {R_GOE}', '#2471A3'),
                        (R_GUE, f'GUE {R_GUE}', '#27AE60')]:
        ax[0].axvline(v, color=col, ls='--', lw=1.2, label=lab)
    ax[0].scatter([SC["r_obs"]], [1], s=120, color='k', zorder=5,
                  label=f'S. cerevisiae obs {SC["r_obs"]:.3f}')
    ax[0].scatter([SC["r_shuf"]], [1], s=60, facecolor='none', edgecolor='orange',
                  lw=2, zorder=6, label='cerevisiae shuffled')
    ax[0].scatter([SP["r_obs"]], [0.6], s=120, color='#7D3C98', zorder=5,
                  label=f'S. pombe obs {SP["r_obs"]:.3f}')
    ax[0].scatter([SP["r_shuf"]], [0.6], s=60, facecolor='none', edgecolor='orange',
                  lw=2, zorder=6, label='pombe shuffled')
    ax[0].set_xlim(0.35, 1.0); ax[0].set_ylim(0, 1.6); ax[0].set_yticks([])
    ax[0].set_xlabel(r'$\langle r\rangle$'); ax[0].legend(fontsize=7.5, loc='upper left')
    ax[0].set_title('(a) both species: obs = shuffled = scalar, NOT GOE')
    b = np.arange(80, 400, 3)
    ax[1].hist(SC["all"], bins=b, density=True, histtype='step', color='k', lw=1.4,
               label=f'S. cerevisiae (med {SC["median"]} bp)')
    ax[1].hist(SP["all"], bins=b, density=True, histtype='step', color='#7D3C98', lw=1.4,
               label=f'S. pombe (med {SP["median"]} bp)')
    ax[1].axvline(147, color='#2471A3', ls=':', lw=1, label='147 bp core (hard core)')
    ax[1].set_xlabel('nucleosome center-to-center spacing (bp)'); ax[1].set_ylabel('density')
    ax[1].legend(fontsize=8)
    ax[1].set_title('(b) pombe spacing shorter/tighter -> higher <r>\n'
                    '(scalar marginal tracks species NRL)')
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, "figures", f"nuc_crossspecies.{ext}"),
                    dpi=120, bbox_inches="tight")

    for name, D in [("S. cerevisiae", SC), ("S. pombe", SP)]:
        print(f"{name:14s} N={D['n']:6d}  <r>_obs={D['r_obs']:.4f}  "
              f"<r>_shuf={D['r_shuf']:.4f}  lag1={D['ac']:+.4f}  "
              f"median={D['median']}bp  CV={D['cv']:.3f}")
    print("wrote results/r_summary.csv, results/r_per_chromosome.csv, "
          "data/spacing_histogram_pombe.csv, figures/nuc_crossspecies.{pdf,png}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cerevisiae", default="data/raw/cerevisiae_unique_map.txt")
    ap.add_argument("--pombe", default="data/raw/pombe_unique_map.txt")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    main(a.cerevisiae, a.pombe, a.outdir)
