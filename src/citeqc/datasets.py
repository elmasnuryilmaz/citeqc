"""A small synthetic CITE-seq dataset with a known answer, for trying citeqc offline.

**Nothing here is real data.** The counts are simulated, so that every
antibody's behaviour is known by construction and you can check that citeqc
recovers it. Use it to learn the workflow, not to draw biological conclusions.

Design (see `DESIGN`):

* 4 donors x 2 T-cell types (``CD4 T``, ``CD8 T``), 320-479 cells per stratum.
* 15 antibodies and 25 genes, as raw UMI counts. The RNA table is a subset of
  the transcriptome: ``meta.nCount_RNA`` is the whole-library size, larger than
  the row sums of `rna`.
* Within each stratum, every gene has a per-cell latent state. *Coupled*
  antibodies are driven by the state of their cognate gene plus
  protein-specific noise; *decoupled* antibodies (CD45RA, CD152, CD154) vary
  independently of their transcript; the B-cell markers CD19 and CD20 are absent
  from T cells and carry only ambient antibody background, so they serve as
  negative controls; IgG1 is an isotype control with no cognate gene.
* Donor ``D4`` was stained with a panel that lacks CD154 (NaN in `adt`), so
  `diagnostics.panel_consistency` has something to find.
* ``BG01``-``BG12`` are background genes with no antibody, for
  `coupling.specificity`.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

__all__ = ["synthetic", "pseudoclonotypes", "write_example", "DESIGN"]

CELLTYPES = ("CD4 T", "CD8 T")
DONORS = ("D1", "D2", "D3", "D4")

#: antibody -> (cognate gene or None, design). The design is the ground truth.
DESIGN = {
    "CD3": ("CD3E", "coupled"),
    "CD2": ("CD2", "coupled"),
    "CD5": ("CD5", "coupled"),
    "CD7": ("CD7", "coupled"),
    "CD27": ("CD27", "coupled"),
    "CD28": ("CD28", "coupled"),
    "CD62L": ("SELL", "coupled"),
    "CD127": ("IL7R", "coupled"),
    "CD45RO": ("PTPRC", "coupled"),
    "CD45RA": ("PTPRC", "decoupled"),
    "CD152": ("CTLA4", "decoupled"),
    "CD154": ("CD40LG", "decoupled"),
    "CD19": ("CD19", "absent"),
    "CD20": ("MS4A1", "absent"),
    "IgG1": (None, "isotype"),
}

# expected UMIs per 10,000 per cell, (CD4 T, CD8 T)
_RNA_LEVEL = {
    "CD3E": (12, 12), "CD2": (6, 8), "CD5": (4, 3), "CD7": (3, 6), "CD27": (3, 3),
    "CD28": (3, 1.5), "SELL": (5, 4), "IL7R": (6, 3), "PTPRC": (25, 25),
    "CTLA4": (1.0, 0.4), "CD40LG": (1.5, 0.3), "CD19": (0.01, 0.01), "MS4A1": (0.01, 0.01),
}
# antibody counts per cell at unit depth, (CD4 T, CD8 T); absent markers are 0
_ADT_LEVEL = {
    "CD3": (90, 90), "CD2": (40, 50), "CD5": (50, 35), "CD7": (20, 40), "CD27": (30, 30),
    "CD28": (35, 20), "CD62L": (30, 25), "CD127": (25, 12), "CD45RO": (40, 30),
    "CD45RA": (35, 50), "CD152": (6, 3), "CD154": (8, 2), "CD19": (0, 0), "CD20": (0, 0),
    "IgG1": (0, 0),
}
_AMBIENT = 3.0          # antibody background counts per cell at unit depth
_LATENT_SD = 0.6        # within-stratum biological variation of each gene (log scale)


def synthetic(seed: int = 0):
    """Simulate the example dataset.

    The same `seed` gives the same tables on the same NumPy version.

    Returns
    -------
    adt : DataFrame
        Cells x 15 antibodies, raw counts (float, NaN where not measured).
    rna : DataFrame
        Cells x 25 genes, raw counts.
    meta : DataFrame
        Cells x [donor, celltype, nCount_ADT, nCount_RNA].
    """
    rng = np.random.default_rng(seed)
    genes = list(_RNA_LEVEL) + [f"BG{i:02d}" for i in range(1, 13)]
    bg_level = rng.lognormal(np.log(3), 1.0, 12)
    level = {g: _RNA_LEVEL.get(g) or (bg_level[i - 13],) * 2 for i, g in enumerate(genes)}
    donor_shift = {d: rng.normal(0, 0.2, len(genes)) for d in DONORS}
    antibodies = list(DESIGN)
    beta = {a: rng.uniform(0.7, 1.0) for a in antibodies}        # protein ~ transcript slope
    noise = {a: rng.uniform(0.3, 0.5) for a in antibodies}       # protein-specific noise

    adt_parts, rna_parts, meta_parts = [], [], []
    for d in DONORS:
        for t, ct in enumerate(CELLTYPES):
            n = int(rng.integers(320, 480))
            z = rng.normal(0, _LATENT_SD, (n, len(genes)))          # per-cell gene state
            mu = np.array([level[g][t] for g in genes]) * np.exp(donor_shift[d])
            depth_rna = np.round(rng.lognormal(np.log(6000), 0.45, n))
            lam = depth_rna[:, None] / 1e4 * mu * np.exp(z - _LATENT_SD ** 2 / 2)
            rna_parts.append(rng.poisson(lam))
            depth_adt = rng.lognormal(0, 0.5, n)
            counts = np.empty((n, len(antibodies)))
            for j, a in enumerate(antibodies):
                gene, design = DESIGN[a]
                eps = rng.normal(0, noise[a], n)
                if design == "coupled":
                    state = beta[a] * z[:, genes.index(gene)] + eps
                else:                     # independent of the cognate transcript
                    state = rng.normal(0, _LATENT_SD, n) + eps
                sd = np.sqrt(np.var(state)) if design == "coupled" else np.hypot(_LATENT_SD, noise[a])
                expr = _ADT_LEVEL[a][t] * np.exp(state - sd ** 2 / 2)
                counts[:, j] = rng.poisson(depth_adt * (expr + _AMBIENT))
            if d == "D4":                                         # panel without CD154
                counts[:, antibodies.index("CD154")] = np.nan
            cells = [f"{d}_{ct.replace(' ', '')}_{i:04d}" for i in range(n)]
            adt_parts.append(pd.DataFrame(counts, index=cells, columns=antibodies))
            meta_parts.append(pd.DataFrame({
                "donor": d, "celltype": ct,
                "nCount_ADT": np.nansum(counts, axis=1).astype(int),
                "nCount_RNA": depth_rna.astype(int)}, index=cells))
    adt = pd.concat(adt_parts)
    rna = pd.DataFrame(np.vstack(rna_parts), index=adt.index, columns=genes)
    meta = pd.concat(meta_parts)
    for df in (adt, rna, meta):
        df.index.name = "cell"
    return adt, rna, meta


def pseudoclonotypes(rna: pd.DataFrame, meta: pd.DataFrame, n_groups: int = 120,
                     min_size: int = 5, max_size: int = 300, seed: int = 0):
    """Groups of unequal size drawn from the same cells: a null for group-size effects.

    Each group is a random sample of cells (sizes log-uniform between `min_size`
    and `max_size`), so no gene truly depends on group size. The value per group
    is a pseudobulk, ``log1p`` of the mean counts per 10,000 (`meta.nCount_RNA`
    as library size). Any association with size that
    `diagnostics.size_confounding` finds is therefore the artefact.

    Returns
    -------
    values : DataFrame
        Groups x genes.
    sizes : Series
        Cells per group.
    members : dict
        ``{group: array of cell ids}``, the input `equal_n_resample` takes.
    """
    rng = np.random.default_rng(seed)
    cp10k = rna.div(meta.loc[rna.index, "nCount_RNA"], axis=0) * 1e4
    sizes = np.round(10 ** rng.uniform(np.log10(min_size), np.log10(max_size), n_groups)).astype(int)
    names = [f"clone{i:03d}" for i in range(n_groups)]
    members = {g: rng.choice(rna.index.to_numpy(), s, replace=False) for g, s in zip(names, sizes)}
    values = pd.DataFrame({g: np.log1p(cp10k.loc[m].mean()) for g, m in members.items()}).T
    values.index.name = "group"
    return values, pd.Series(sizes, index=names, name="n_cells"), members


def write_example(directory, seed: int = 0) -> dict:
    """Write the example dataset as CSV files for the command line.

    Files: ``adt.csv``, ``rna.csv``, ``meta.csv`` (cells in rows), ``groups.csv``
    (pseudo-clonotype pseudobulk plus an ``n_cells`` column, for
    ``citeqc diagnose --groups``) and ``design.csv`` (the ground truth).
    Returns {name: path}.
    """
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    adt, rna, meta = synthetic(seed)
    values, sizes, _ = pseudoclonotypes(rna, meta, seed=seed)
    design = pd.DataFrame([(a, g or "", d) for a, (g, d) in DESIGN.items()],
                          columns=["antibody", "gene", "design"])
    paths = {k: out / f"{k}.csv" for k in ("adt", "rna", "meta", "groups", "design")}
    adt.to_csv(paths["adt"])
    rna.to_csv(paths["rna"])
    meta.to_csv(paths["meta"])
    values.join(sizes).to_csv(paths["groups"])
    design.to_csv(paths["design"], index=False)
    return paths
