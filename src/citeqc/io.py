"""Readers. CSV/TSV always, parquet with pyarrow, AnnData/MuData with the `anndata` extra."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

__all__ = ["read_matrix", "from_mudata"]

_COMPRESSION = {".gz", ".bz2", ".xz", ".zip", ".zst"}


def _separator(p: Path) -> str:
    """Tab for .tsv/.txt (optionally compressed, e.g. .tsv.gz), comma otherwise."""
    suffixes = [s.lower() for s in p.suffixes]
    if suffixes and suffixes[-1] in _COMPRESSION:
        suffixes = suffixes[:-1]
    return "\t" if suffixes and suffixes[-1] in (".tsv", ".txt") else ","


def read_matrix(path, transpose: bool = False) -> pd.DataFrame:
    """Read a table whose first column is the row index (cell barcodes or group names).

    Tab-separated for `.tsv`/`.txt`, comma-separated otherwise; gzip, bz2, xz,
    zip and zstd compression are detected from the extension. `.parquet` files
    need `pyarrow` (or `fastparquet`). Pass `transpose=True` if the file is
    features x cells.
    """
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        df = pd.read_parquet(p)
    else:
        df = pd.read_csv(p, sep=_separator(p), index_col=0)
    if transpose:
        df = df.T
    df.index = df.index.astype(str)
    return df


def _row_sums(adata) -> pd.Series:
    """Per-cell total of `adata.X`, dense or sparse, indexed by cell barcode."""
    tot = np.asarray(adata.X.sum(axis=1), dtype=float).ravel()
    return pd.Series(tot, index=adata.obs_names)


def from_mudata(path, adt_key: str = "prot", rna_key: str = "rna",
                strata=("donor", "celltype")):
    """Read `(adt, rna, meta)` from a `.h5mu` file. Requires the `anndata` extra.

    `adt` and `rna` are the `.X` matrices of the two modalities as cells x
    features DataFrames, and `meta` is the global `.obs`, which must contain
    the `strata` columns. Cells missing from either modality are dropped from
    all three tables.

    If `.obs` has no `nCount_ADT` / `nCount_RNA` column, they are computed as
    the per-cell sum of each modality's `.X`. That is only a library size if
    `.X` holds raw counts: if your `.X` is already normalised, add the raw
    totals to `.obs` before calling this.
    """
    try:
        import mudata
    except ImportError as e:                                  # pragma: no cover
        raise ImportError("reading .h5mu needs the optional dependency: "
                          "pip install 'citeqc[anndata]'") from e
    md = mudata.read_h5mu(str(path))
    for k in (adt_key, rna_key):
        if k not in md.mod:
            raise KeyError(f"modality {k!r} not in file; found {list(md.mod)}")
    a_mod, r_mod = md.mod[adt_key], md.mod[rna_key]
    for c in strata:
        if c not in md.obs.columns:
            raise KeyError(f"obs is missing stratum column {c!r}")
    cells = a_mod.obs_names.intersection(r_mod.obs_names).intersection(md.obs_names)
    if len(cells) == 0:
        raise ValueError(f"modalities {adt_key!r} and {rna_key!r} share no cells")
    meta = md.obs.loc[cells].copy()
    if "nCount_ADT" not in meta:
        meta["nCount_ADT"] = _row_sums(a_mod).reindex(cells).to_numpy()
    if "nCount_RNA" not in meta:
        meta["nCount_RNA"] = _row_sums(r_mod).reindex(cells).to_numpy()
    adt = a_mod.to_df().loc[cells]
    rna = r_mod.to_df().loc[cells]
    return adt, rna, meta
