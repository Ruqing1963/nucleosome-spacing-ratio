"""
stats_common.py
---------------
Shared helpers for the consecutive level-spacing ratio <r> analysis of
base-pair-resolution nucleosome dyad maps.

The maps themselves are third-party data and are NOT redistributed here; download
them from the original accessions (see README) into data/raw/ and pass their
paths to the analysis scripts.

Expected map format (whitespace-delimited, one nucleosome per line):
    chromosome   dyad_position(bp)   NCP_score   NCP_score_to_noise
"""
import collections
import numpy as np

# --- RMT reference values for the consecutive spacing ratio <r> ---
R_POISSON = 0.386   # uncorrelated / superposition limit
R_GOE     = 0.531   # Gaussian orthogonal ensemble (level repulsion)
R_GUE     = 0.603   # Gaussian unitary ensemble (stronger repulsion)
#  (a single rigidly regular "picket-fence" sequence tends to <r> -> 1)

# --- standard chromosome orders ---
SC_CHROMS = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII',
             'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII',
             'chrXIV', 'chrXV', 'chrXVI']          # S. cerevisiae (16)
SP_CHROMS = ['chrI', 'chrII', 'chrIII']            # S. pombe (3)


def load_map(path):
    """Read a unique chemical map. Returns list of
    (chrom:str, pos:int, ncp_score:float, ncp_to_noise:float)."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) >= 4:
                try:
                    rows.append((p[0], int(p[1]), float(p[2]), float(p[3])))
                except ValueError:
                    continue
    return rows


def spacings_by_chrom(rows, chrom_order=None, gapcap=None):
    """Consecutive center-to-center spacings within each chromosome.

    gapcap : if given, discard spacings larger than this many bp (used to
             excise large nucleosome-free-region gaps and isolate the bulk
             array)."""
    by = collections.defaultdict(list)
    for c, pos, _s, _sn in rows:
        by[c].append(pos)
    order = chrom_order or sorted(by)
    out = {}
    for c in order:
        if c in by:
            d = np.diff(np.sort(np.array(by[c])))
            if gapcap is not None:
                d = d[d <= gapcap]
            out[c] = d[d > 0]
    return out


def r_pooled(spac):
    """Mean consecutive spacing ratio r_i = min(s_i,s_{i+1})/max(s_i,s_{i+1}),
    pooled within chromosomes (never across chromosome boundaries).
    Returns (mean_r, n_pairs)."""
    rs = []
    for d in spac.values():
        if len(d) >= 3:
            rs.append(np.minimum(d[:-1], d[1:]) / np.maximum(d[:-1], d[1:]))
    rr = np.concatenate(rs)
    return float(rr.mean()), int(len(rr))


def r_shuffled(spac, n_perm=200, seed=0):
    """Permutation (shuffle) null: randomly permute the ORDER of the spacings
    within each chromosome. This preserves the marginal spacing distribution
    exactly while destroying any sequential correlation; it is equivalent to an
    i.i.d. renewal drawn from the empirical linker distribution.
    Returns (mean_r, std_r) over n_perm permutations."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_perm):
        sd = {}
        for c, d in spac.items():
            dd = d.copy()
            rng.shuffle(dd)
            sd[c] = dd
        vals.append(r_pooled(sd)[0])
    return float(np.mean(vals)), float(np.std(vals))


def lag1_autocorr(spac):
    """Mean lag-1 autocorrelation of the spacings (clearly negative for a rigid
    / level-repulsive spectrum, ~0 for a renewal process)."""
    a = []
    for d in spac.values():
        x = d.astype(float) - d.mean()
        a.append(np.sum(x[:-1] * x[1:]) / np.sum(x * x))
    return float(np.mean(a))


def all_spacings(spac):
    """Flatten the per-chromosome spacings into one array."""
    return np.concatenate(list(spac.values()))
