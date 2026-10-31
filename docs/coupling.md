# How coupling is measured

[`coupling.measure`](api.md#citeqc.coupling.measure) asks, for each antibody:
*within a homogeneous group of cells, does the antibody signal track its own
transcript?*

## Inputs

| Table | Rows | Columns |
|---|---|---|
| `adt` | cells | antibodies, normalised (for example with [`clr`](api.md#citeqc.coupling.clr)) |
| `rna` | cells | genes, normalised (for example with [`log_normalise`](api.md#citeqc.coupling.log_normalise)) |
| `meta` | cells | the stratum columns (default `donor`, `celltype`) and the raw library sizes `nCount_ADT`, `nCount_RNA` |

Cells are matched by barcode (the row index); cells missing from any table
are dropped. `pairs` maps each antibody to its cognate gene, usually from
[`mapping.map_antibodies`](api.md#citeqc.mapping.map_antibodies).

## Three constraints the function enforces

1. **Stratify by donor x cell type.** Differences in cell-type composition move
   protein and mRNA together, so an unstratified correlation mostly measures
   composition and is strongly inflated.
2. **Use equal-sized pseudo-groups.** A mean over more cells is less noisy, so
   groups of unequal size introduce a size-driven gradient across every
   feature (see [`size_confounding`](diagnostics.md#size_confounding)).
3. **Regress out library depth** in both modalities before correlating.

## The algorithm

For every stratum (each combination of the `strata` columns) with at least
`group_size x min_groups` cells (10 x 25 = 250 by default):

1. Shuffle the cells and split them into `n // group_size` pseudo-groups of
   exactly `group_size` cells; leftover cells are not used in that draw.
2. Average the antibody and RNA values within each pseudo-group, and average
   `log10(depth + 1)` for each modality.
3. For each antibody, regress its pseudo-group means on the mean ADT depth
   (ordinary least squares with an intercept), and its cognate gene's means on
   the mean RNA depth. Keep the residuals.
4. Take the Spearman correlation of the two residual vectors across
   pseudo-groups. Skip it if fewer than `min_groups` pseudo-groups have finite
   values (for example because the antibody was not measured in that donor),
   or if either residual has zero variance.
5. Repeat for `n_draw` random partitions (10 by default).

Each stratum x draw gives one observation. The summary has one row per
antibody:

| Column | Meaning |
|---|---|
| `median_rho` | median of all observations: the headline number |
| `q25`, `q75` | quartiles across strata and draws |
| `n_obs` | number of observations; antibodies with fewer than `min_obs` (10) are left out |
| `n_strata` | number of strata that contributed |

`measure(..., return_raw=True)` also returns every observation with its
stratum labels and the residual standard deviations. Observations are **not
independent**: draws within a stratum share cells, and strata within a donor
share the donor. For confidence intervals, resample donors, recompute the
per-marker medians, and ask how often a marker keeps its rank.

## How to read the result

Absolute correlations depend on how finely cell types are split, so they are
**not comparable between datasets**. Compare
[`rank_percentile`](api.md#citeqc.coupling.rank_percentile) instead. Within a
dataset, calibrate against a marker that should couple well in the cells you
are looking at (CD8/*CD8A* in CD8 T cells, CD19/*CD19* in B cells). That value is
your ceiling; markers far below it are not readable from RNA.

The median runs over all strata, including strata where the marker is not
expressed. A lineage marker such as CD8 is therefore diluted by the CD4 T-cell
strata. Use `return_raw=True` for a per-cell-type breakdown.

## Is a low value specific to that marker?

A correlation of 0.02 is hard to read on its own: it may be a decoupled marker,
or a marker the assay cannot measure.
[`coupling.specificity`](api.md#citeqc.coupling.specificity) correlates every
antibody with every gene in a background set, using the same stratification,
grouping and depth adjustment, and reports where the cognate gene falls in that
antibody's row (`spec_percentile`, `rank_of_cognate`, `z_vs_background`).

`spec_percentile` of 1.0 means the protein tracks its own transcript better than
any other gene; about 0.5 means the cognate gene is indistinguishable from the
background. The background is descriptive, not a null: genes that truly co-vary
with the antibody sit in it. A marker can be low on `measure` yet specific
here, so weak coupling is not the same as no coupling. With the same `seed`,
the cognate entries equal `measure`'s medians.

## Normalisation helpers

- [`clr`](api.md#citeqc.coupling.clr) computes `log1p(x)` minus the per-cell
  mean of `log1p(x)` across antibodies (missing values are ignored in the mean).
  Other tools implement CLR slightly differently; if you need to match another
  pipeline exactly, normalise there and pass the normalised tables.
- [`log_normalise`](api.md#citeqc.coupling.log_normalise) computes
  `log1p(counts / total x 10,000)`. Pass `total=meta["nCount_RNA"]` when your
  table holds a subset of genes; by default the row sums are used.

Both are compositional: each cell is scaled by a quantity computed from all its
features. On a small antibody panel, a strongly coupled antibody can therefore
leak a little of its transcript's signal into every other antibody of the
same cell. The [quick start](quickstart.md#reading-the-result) shows this
for CD45RA and CD45RO, which share the cognate gene *PTPRC*.

## The antibody-to-gene table

145 entries covering common TotalSeq naming, with an annotation where a
documented post-transcriptional mechanism exists: `splice_isoform` (CD45RA,
CD45RO), `vesicular` (for example CD152/CTLA-4 and CD154/CD40L) and `IRE_3UTR`
(CD71/*TFRC*). Inspect it with `citeqc map` or
[`mapping.load_map`](api.md#citeqc.mapping.load_map).

[`map_antibodies`](api.md#citeqc.mapping.map_antibodies) tries the raw name,
then a normalised form (`CD3--UCHT1-TSA` becomes `CD3`), then that form with a
trailing replicate index removed (`CD38-1` becomes `CD38`). Antibodies it
cannot place are **returned as unmapped, never guessed**. Extend it with
`map_antibodies(..., extra={"MyAb": "GENE"})`.

Three caveats: the mapping is one-to-one, so it cannot express multi-gene or
post-translationally modified epitopes; `HLA-DR` maps to *HLA-DRA* alone; and
CD45RA and CD45RO both map to *PTPRC*, which is the point rather than a
limitation.
