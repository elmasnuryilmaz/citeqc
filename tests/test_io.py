import gzip
import numpy as np
import pandas as pd
import pytest
from citeqc.io import read_matrix, from_mudata


@pytest.fixture
def table():
    return pd.DataFrame({"CD3": [1.0, 2.0], "CD19": [0.0, 5.0]},
                        index=pd.Index(["c1", "c2"], name="cell"))


@pytest.mark.parametrize("name,sep", [
    ("adt.csv", ","), ("adt.tsv", "\t"), ("adt.txt", "\t"),
    ("adt.csv.gz", ","), ("adt.tsv.gz", "\t"),
    ("sample.v2.tsv.gz", "\t"),        # extra dots used to fall back to commas
    ("ADT.TSV", "\t"),
])
def test_read_matrix_detects_separator(tmp_path, table, name, sep):
    path = tmp_path / name
    text = table.to_csv(sep=sep)
    if name.endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    out = read_matrix(path)
    pd.testing.assert_frame_equal(out, table, check_names=False)


def test_read_matrix_transposes_features_by_cells(tmp_path, table):
    path = tmp_path / "adt.csv"
    table.T.to_csv(path)
    out = read_matrix(path, transpose=True)
    assert list(out.index) == ["c1", "c2"] and list(out.columns) == ["CD3", "CD19"]


# ------------------------------------------------------- the anndata extra
def _write_h5mu(path, sparse=True, drop_rna_cell=False):
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    sp = pytest.importorskip("scipy.sparse")
    rng = np.random.default_rng(0)
    cells = [f"c{i}" for i in range(6)]
    adt_x = rng.poisson(5, (6, 3)).astype(np.float32)
    rna_x = rng.poisson(2, (6, 4)).astype(np.float32)
    prot = ad.AnnData(sp.csr_matrix(adt_x) if sparse else adt_x,
                      obs=pd.DataFrame(index=cells),
                      var=pd.DataFrame(index=["CD3", "CD4", "CD8"]))
    rna_cells = cells[:-1] if drop_rna_cell else cells
    rna = ad.AnnData(sp.csr_matrix(rna_x[:len(rna_cells)]) if sparse else rna_x[:len(rna_cells)],
                     obs=pd.DataFrame(index=rna_cells),
                     var=pd.DataFrame(index=["CD3E", "CD4", "CD8A", "ACTB"]))
    mdata = md.MuData({"prot": prot, "rna": rna})
    mdata.obs["donor"] = "d1"
    mdata.obs["celltype"] = ["T"] * 3 + ["B"] * 3
    mdata.write(str(path))
    return adt_x, rna_x, cells


@pytest.mark.parametrize("sparse", [True, False])
def test_from_mudata_reads_both_modalities(tmp_path, sparse):
    adt_x, rna_x, cells = _write_h5mu(tmp_path / "x.h5mu", sparse=sparse)
    adt, rna, meta = from_mudata(tmp_path / "x.h5mu")
    assert list(adt.index) == cells and list(adt.columns) == ["CD3", "CD4", "CD8"]
    assert rna.shape == (6, 4)
    assert {"donor", "celltype", "nCount_ADT", "nCount_RNA"} <= set(meta.columns)
    np.testing.assert_allclose(meta["nCount_ADT"].to_numpy(), adt_x.sum(axis=1))
    np.testing.assert_allclose(meta["nCount_RNA"].to_numpy(), rna_x.sum(axis=1))


def test_from_mudata_keeps_cells_present_in_both_modalities(tmp_path):
    _write_h5mu(tmp_path / "x.h5mu", drop_rna_cell=True)
    adt, rna, meta = from_mudata(tmp_path / "x.h5mu")
    assert len(adt) == len(rna) == len(meta) == 5
    assert list(adt.index) == list(rna.index) == list(meta.index)
    assert meta["nCount_RNA"].notna().all()


def test_from_mudata_errors(tmp_path):
    _write_h5mu(tmp_path / "x.h5mu")
    with pytest.raises(KeyError, match="modality 'adt'"):
        from_mudata(tmp_path / "x.h5mu", adt_key="adt")
    with pytest.raises(KeyError, match="stratum column 'batch'"):
        from_mudata(tmp_path / "x.h5mu", strata=("donor", "batch"))
