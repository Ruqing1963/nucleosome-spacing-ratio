# nucleosome-spacing-ratio

**The consecutive spacing-ratio of nucleosome arrays: high regularity from steric renewal, not level repulsion**

Ruqing Chen — GUT Geoservice Inc., Montréal, Québec, Canada

This repository contains the analysis code, derived data and figures for a
random-matrix-theory (RMT) **consecutive level-spacing ratio ⟨r⟩** analysis of
base-pair-resolution nucleosome dyad maps in two yeasts.

---

## Summary

The spacings between consecutive nucleosome dyads are summarised with the
unfold-free ratio

```
r_i = min(s_i, s_{i+1}) / max(s_i, s_{i+1})
```

whose mean ⟨r⟩ distinguishes uncorrelated/superposed sequences (Poisson,
⟨r⟩ = 0.386), level-repulsive Wigner–Dyson spectra (GOE 0.531 / GUE 0.603) and
rigidly regular sequences (⟨r⟩ → 1), with **no spectral unfolding**.

Applied to the chemical maps of *Saccharomyces cerevisiae* and
*Schizosaccharomyces pombe*, the dyad arrays give **⟨r⟩ ≈ 0.78 and 0.81** — far
above every RMT ensemble, superficially suggesting strong order. A
**permutation (shuffle) test** overturns this reading: randomising the order of
the consecutive linker lengths leaves ⟨r⟩ unchanged (within 0.001), and the
lag-1 spacing autocorrelation is ≈ 0. The high ⟨r⟩ is therefore **entirely a
property of the marginal linker-length distribution** — a steric hard-core
(Tonks-gas-like) renewal — with **no sequential level repulsion**.

The verdict is robust to confidence thresholds and to nucleosome-free-region
trimming, and holds in both species, with ⟨r⟩ tracking each species' repeat
length. Superposing independent arrays, or randomly thinning a single array,
both drive ⟨r⟩ back toward the Poisson value, confirming that the order is local
and steric. **The permutation test is a necessary control** when ⟨r⟩ is applied
to genomic or spatial point patterns.

### Key numbers

| Species | N pairs | ⟨r⟩_obs | ⟨r⟩_shuf | lag-1 autocorr. | median spacing (CV) |
|---|---:|---:|---:|---:|---|
| *S. cerevisiae* | 67,516 | 0.783 | 0.782 | +0.05 | 164 bp (0.61) |
| *S. pombe*      | 75,822 | 0.811 | 0.811 | +0.07 | 151 bp (0.47) |

(Reference: Poisson 0.386, GOE 0.531, GUE 0.603, regular → 1.)

---

## Repository layout

```
nucleosome-spacing-ratio/
├── README.md
├── LICENSE                     # CC BY 4.0 (raw third-party maps NOT included)
├── requirements.txt            # numpy, matplotlib
├── code/
│   ├── stats_common.py         # <r>, shuffle null, lag-1 autocorr, map parser
│   ├── analyze_core.py         # core <r> + spacing dist + superposition (Fig 1)
│   ├── analyze_robustness.py   # gap-trim + confidence thinning (Table 2)
│   └── analyze_crossspecies.py # both species (Table 1, Fig 2)
├── data/                       # DERIVED data only (copyright-safe)
│   ├── spacing_histogram_cerevisiae.csv
│   └── spacing_histogram_pombe.csv
├── results/
│   ├── r_summary.csv           # Table 1
│   ├── r_per_chromosome.csv    # per-chromosome <r>, both species
│   ├── robustness.csv          # Table 2
│   └── superposition.csv       # <r> vs k independent arrays superposed
├── figures/
│   ├── nuc_fig.pdf / .png             # Figure 1
│   └── nuc_crossspecies.pdf / .png    # Figure 2
└── paper/
    ├── nuc_paper.pdf
    ├── nuc_paper.tex
    └── (figure pdfs for compilation)
```

---

## Data sources (not redistributed here)

The raw base-pair-resolution **unique** chemical maps are third-party data and
are **not** included in this repository. Download them and place them under
`data/raw/`:

- **S. cerevisiae** — Brogaard, Xi, Wang & Widom (2012), *Nature* **486**, 496–501,
  doi:10.1038/nature11142. GEO accession **GSE36063**; the unique map is
  Supplementary Table 2 (Nature SI). 67,548 dyads, 16 chromosomes.
- **S. pombe** — Moyle-Heyrman *et al.* (2013), *PNAS* **110**, 20158–20163,
  doi:10.1073/pnas.1315809110. GEO accession **GSE46975**; the unique map is
  **Supporting Information Dataset S1** (`sd01.txt`). 75,828 dyads, 3 chromosomes.

Use the **unique** (distinct-nucleosome) maps, **not** the redundant maps: the
latter enumerate overlapping candidate positions per genomic nucleosome (median
centre-to-centre distance ~12 bp, below the 147 bp footprint) and are unsuitable
for inter-nucleosome spacing.

Expected file format (whitespace-delimited, one nucleosome per line):
`chromosome  dyad_position(bp)  NCP_score  NCP_score_to_noise`.

---

## Reproducing the analysis

```bash
pip install -r requirements.txt

# place the downloaded unique maps here:
#   data/raw/cerevisiae_unique_map.txt
#   data/raw/pombe_unique_map.txt

cd code
python analyze_core.py        --map data/raw/cerevisiae_unique_map.txt --outdir ..
python analyze_robustness.py  --map data/raw/cerevisiae_unique_map.txt --outdir ..
python analyze_crossspecies.py --cerevisiae data/raw/cerevisiae_unique_map.txt \
                               --pombe      data/raw/pombe_unique_map.txt --outdir ..
```

This regenerates everything in `results/`, the derived histograms in `data/`,
and the figures in `figures/`. (Paths are relative to `code/`; adjust `--outdir`
as needed.)

---

## Method notes

- **⟨r⟩ is two-ended.** The Poisson value 0.386 is the *uncorrelated/superposition*
  limit, **not** the value of a single regular sequence (which tends to 1). A
  large ⟨r⟩ alone does not imply level repulsion.
- **The shuffle test is the decisive control.** Permuting the order of the
  spacings preserves the marginal but destroys sequential correlation; if
  ⟨r⟩_obs ≈ ⟨r⟩_shuf the value is scalar (set by the marginal), not repulsive.
- **Discretisation robustness.** The ~10 bp rotational comb and integer-bp
  resolution make the spacings partly discretised, but ⟨r⟩_obs and ⟨r⟩_shuf are
  computed on the identically discretised spacings, so the shuffle test is
  unaffected.

---

## Citation

If you use this code or data, please cite the paper (`paper/nuc_paper.pdf`) and
this repository. A Zenodo DOI will be assigned on release and added here.

## License

Code, derived data, figures and paper: **CC BY 4.0** (see `LICENSE`).
The raw chemical maps are **not** covered by this license and must be obtained
from the original sources above.
