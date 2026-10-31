import numpy as np
import pandas as pd
import pytest
from citeqc import coupling, datasets, mapping, diagnostics as dx


@pytest.fixture(scope="module")
def example():
    return datasets.synthetic(seed=0)


@pytest.fixture(scope="module")
def normalised(example):
    adt, rna, meta = example
    A = pd.DataFrame(coupling.clr(adt.values), index=adt.index, columns=adt.columns)
    R = coupling.log_normalise(rna, total=meta.nCount_RNA)
    return A, R, meta


def test_shapes_and_labels(example):
    adt, rna, meta = example
    assert list(adt.columns) == list(datasets.DESIGN)
    assert rna.shape[1] == 25 and len(adt) == len(rna) == len(meta)
    assert meta.groupby(["donor", "celltype"]).size().min() >= 250   # enough for defaults
    assert (meta.nCount_RNA > rna.sum(axis=1)).all()                 # rna is a gene subset
    assert adt.loc[meta.donor == "D4", "CD154"].isna().all()
    assert adt.drop(columns="CD154").notna().all().all()


def test_same_seed_same_data(example):
    adt, rna, meta = datasets.synthetic(seed=0)
    pd.testing.assert_frame_equal(adt, example[0])
    pd.testing.assert_frame_equal(rna, example[1])


def test_every_design_antibody_maps_as_designed(example):
    adt, rna, _ = example
    pairs, unmapped = mapping.map_antibodies(adt.columns, rna.columns)
    assert unmapped == ["IgG1"]
    for ab, (gene, _) in datasets.DESIGN.items():
        assert pairs.get(ab) == gene


def test_measure_recovers_the_design(normalised):
    A, R, meta = normalised
    pairs, _ = mapping.map_antibodies(A.columns, R.columns)
    rho = coupling.measure(A, R, meta, pairs).set_index("marker").median_rho
    design = {a: d for a, (_, d) in datasets.DESIGN.items() if a in rho.index}
    coupled = [a for a, d in design.items() if d == "coupled"]
    other = [a for a, d in design.items() if d != "coupled"]
    assert rho[coupled].min() > 0.25
    assert rho[other].abs().max() < 0.15
    assert rho[coupled].min() > rho[other].max()


def test_specificity_ranks_coupled_markers_first(normalised):
    A, R, meta = normalised
    pairs, _ = mapping.map_antibodies(A.columns, R.columns)
    _, S = coupling.specificity(A, R, meta, pairs, genes=list(R.columns), n_draw=3)
    s = S.set_index("marker")
    coupled = [a for a, (_, d) in datasets.DESIGN.items() if d == "coupled"]
    assert (s.loc[coupled, "rank_of_cognate"] == 1).all()


def test_diagnostics_find_what_was_planted(example):
    adt, rna, meta = example
    panel = dx.panel_consistency(adt, meta, by="donor")
    assert not panel.passed and "CD154" not in panel.detail["core_panel"]
    values, sizes, members = datasets.pseudoclonotypes(rna, meta)
    assert not dx.size_confounding(values, sizes).passed
    cp10k = rna.div(meta.nCount_RNA, axis=0) * 1e4
    draw = next(dx.equal_n_resample(members, n=5, n_draw=1))
    equal = pd.DataFrame({g: np.log1p(cp10k.loc[c].mean()) for g, c in draw.items()}).T
    assert dx.size_confounding(equal, sizes).passed


def test_write_example(tmp_path):
    paths = datasets.write_example(tmp_path / "ex")
    assert set(paths) == {"adt", "rna", "meta", "groups", "design"}
    assert all(p.exists() for p in paths.values())
    groups = pd.read_csv(paths["groups"], index_col=0)
    assert "n_cells" in groups.columns and groups.index[0].startswith("clone")
