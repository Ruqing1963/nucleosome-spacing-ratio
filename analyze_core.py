"""
analyze_core.py
---------------
Core consecutive-spacing-ratio analysis for the S. cerevisiae unique map:
  * pooled and per-chromosome <r>, the permutation (shuffle) null, lag-1
    autocorrelation, median spacing and CV (Table 1, cerevisiae row);
  * the spacing distribution (164 bp mode, 147 bp hard core, ~10 bp comb);
  * the superposition demonstration: overlaying k independent arrays drives
    <r> toward the Poisson value.
Reproduces Figure 1 and writes results/superposition.csv and
data/spacing_histogram_cerevisiae.csv.

Usage:
    python analyze_core.py --map data/raw/cerevisiae_unique_map.txt
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
                          R_POISSON, R_GOE, R_GUE, SC_CHROMS)


def _rmean(s):
    s = np.asarray(s, float); s = s[s > 0]
    return float((np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])).mean())


def main(map_path, outdir):
    rows = load_map(map_path)
    spac = spacings_by_chrom(rows, SC_CHROMS)
    alls = all_spacings(spac)

    r_obs, n = r_pooled(spac)
    r_shuf, r_sd = r_shuffled(spac)
    ac = lag1_autocorr(spac)
    print(f"S. cerevisiae  N={n}  <r>_obs={r_obs:.4f}  <r>_shuf={r_shuf:.4f}"
          f"+/-{r_sd:.4f}  lag1={ac:+.4f}  median={np.median(alls):.0f}bp"
          f"  CV={alls.std()/alls.mean():.3f}")

    # ---- superposition: overlay k independent arrays -> Poisson ----
    posd = {c: np.concatenate([[0], np.cumsum(d)]) for c, d in spac.items()}
    chrl = list(posd.keys())
    rng = np.random.default_rng(1)
    ks = [1, 2, 3, 4, 6, 8, 12, 16]
    rk, rk_iid = [], []
    for k in ks:
        vr = []
        for _ in range(40):
            pick = rng.choice(chrl, size=k, replace=False)
            m = np.sort(np.concatenate([posd[c] for c in pick]))
            d = np.diff(m); d = d[d > 0]
            if len(d) > 3:
                vr.append(_rmean(d))
        vs = []
        for _ in range(40):
            seqs = [np.cumsum(rng.choice(alls, size=4000)) for _ in range(k)]
            m = np.sort(np.concatenate(seqs)); d = np.diff(m); d = d[d > 0]
            vs.append(_rmean(d))
        rk.append(np.mean(vr)); rk_iid.append(np.mean(vs))

    os.makedirs(os.path.join(outdir, "results"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "data"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "figures"), exist_ok=True)
    with open(os.path.join(outdir, "results", "superposition.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k_superposed", "r_overlay_real_chromosomes",
                    "r_overlay_iid_marginal"])
        for k, a, b in zip(ks, rk, rk_iid):
            w.writerow([k, round(float(a), 4), round(float(b), 4)])

    # ---- derived spacing histogram (copyright-safe) ----
    h, e = np.histogram(alls, bins=np.arange(80, 601, 2))
    ctr = (0.5 * (e[1:] + e[:-1])).astype(int)
    with open(os.path.join(outdir, "data", "spacing_histogram_cerevisiae.csv"),
              "w", newline="") as f:
        w = csv.writer(f); w.writerow(["spacing_bp_bin_center", "count"])
        for b, c in zip(ctr, h):
            w.writerow([int(b), int(c)])

    # ---- Figure 1 ----
    sb = alls[(alls >= 120) & (alls <= 420)]
    hist, edges = np.histogram(sb, bins=np.arange(120, 421, 2))
    cc = 0.5 * (edges[1:] + edges[:-1])
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    ax[0].axvline(R_POISSON, color='#C0392B', ls='--', lw=1.3, label=f'Poisson {R_POISSON}')
    ax[0].axvline(R_GOE, color='#2471A3', ls='--', lw=1.3, label=f'GOE {R_GOE}')
    ax[0].axvline(R_GUE, color='#27AE60', ls='--', lw=1.3, label=f'GUE {R_GUE}')
    ax[0].axvline(r_obs, color='k', lw=2.4, label=f'observed {r_obs:.3f}')
    ax[0].axvline(r_shuf, color='orange', ls=':', lw=2, label=f'shuffled {r_shuf:.3f}')
    ax[0].set_xlim(0.35, 1.0); ax[0].set_yticks([]); ax[0].set_xlabel(r'$\langle r\rangle$')
    ax[0].legend(fontsize=8, loc='upper left')
    ax[0].set_title('(a) <r>: observed = shuffled\n= scalar/steric, NOT GOE')
    ax[1].plot(cc, hist, color='#555', lw=1)
    for m in range(130, 420, 10):
        ax[1].axvline(m, color='#ddd', lw=.4, zorder=0)
    ax[1].axvline(164, color='#C0392B', ls='--', lw=1, label='164 bp (repeat)')
    ax[1].axvline(147, color='#2471A3', ls=':', lw=1, label='147 bp (core/hard-core)')
    ax[1].set_xlabel('nucleosome center-to-center spacing (bp)'); ax[1].set_ylabel('count')
    ax[1].legend(fontsize=8)
    ax[1].set_title('(b) spacing distribution\n(~164 bp mode, hard core, 10-bp comb)')
    ax[2].plot(ks, rk, 'o-', color='#7D3C98', label='overlay real chromosomes')
    ax[2].plot(ks, rk_iid, 's--', color='#999', label='overlay i.i.d. (empirical marginal)')
    ax[2].axhline(R_POISSON, color='#C0392B', ls='--', lw=1, label=f'Poisson {R_POISSON}')
    ax[2].axhline(r_obs, color='k', ls=':', lw=1, label=f'single chain {r_obs:.3f}')
    ax[2].set_xlabel('k = number of independent sequences superposed')
    ax[2].set_ylabel(r'$\langle r\rangle$')
    ax[2].legend(fontsize=8)
    ax[2].set_title('(c) superposition of independent arrays -> Poisson')
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, "figures", f"nuc_fig.{ext}"),
                    dpi=120, bbox_inches="tight")
    print("wrote results/superposition.csv, data/spacing_histogram_cerevisiae.csv, figures/nuc_fig.{pdf,png}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", default="data/raw/cerevisiae_unique_map.txt",
                    help="path to the S. cerevisiae unique chemical map "
                         "(Brogaard 2012; download per README)")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    main(a.map, a.outdir)
