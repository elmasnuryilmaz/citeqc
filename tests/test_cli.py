import pandas as pd
import pytest
from citeqc.cli import main


def _write(tmp_path, synthetic):
    adt, rna, meta = synthetic
    adt.to_csv(tmp_path / "adt.csv"); rna.to_csv(tmp_path / "rna.csv")
    meta.to_csv(tmp_path / "meta.csv")
    return [str(tmp_path / f) for f in ("adt.csv", "rna.csv", "meta.csv")]


def test_couple_end_to_end(tmp_path, synthetic, capsys):
    a, r, m = _write(tmp_path, synthetic)
    out = tmp_path / "res.csv"
    assert main(["couple", "--adt", a, "--rna", r, "--meta", m,
                 "--out", str(out), "--draws", "3"]) == 0
    res = pd.read_csv(out).set_index("marker")
    assert res.loc["CD8", "median_rho"] > res.loc["CD152", "median_rho"]
    assert {"pct_rank", "mechanism"} <= set(res.columns)
    assert res.loc["CD152", "mechanism"] == "vesicular"


def test_diagnose_exit_code_signals_flags(tmp_path, synthetic):
    adt, _, meta = synthetic
    adt = adt.copy(); adt.loc[meta.donor == "d2", "CD19"] = None
    adt.to_csv(tmp_path / "adt.csv"); meta.to_csv(tmp_path / "meta.csv")
    rc = main(["diagnose", "--adt", str(tmp_path / "adt.csv"),
               "--meta", str(tmp_path / "meta.csv")])
    assert rc == 1                                   # a flag must be non-zero


def test_diagnose_without_inputs_errors():
    with pytest.raises(SystemExit):
        main(["diagnose"])


def test_map_query(capsys):
    assert main(["map", "--query", "ctla"]) == 0
    assert "CTLA4" in capsys.readouterr().out


# ----------------------------------------------- example data + regressions
@pytest.fixture
def example_dir(tmp_path):
    assert main(["example", str(tmp_path / "ex")]) == 0
    return tmp_path / "ex"


def test_example_then_couple_with_normalise(example_dir, tmp_path):
    out = tmp_path / "coupling.csv"
    ex = example_dir
    assert main(["couple", "--adt", str(ex / "adt.csv"), "--rna", str(ex / "rna.csv"),
                 "--meta", str(ex / "meta.csv"), "--normalise", "--draws", "3",
                 "--out", str(out)]) == 0
    res = pd.read_csv(out).set_index("marker")
    assert res.loc["CD3", "median_rho"] > 0.3 > abs(res.loc["CD45RA", "median_rho"])
    assert res.loc["CD45RA", "mechanism"] == "splice_isoform"


def test_normalise_scales_rna_by_library_size(monkeypatch, example_dir, tmp_path):
    """log-CP10K must divide by meta's nCount_RNA, not by the row sum of a gene subset."""
    from citeqc import coupling
    seen = {}
    real = coupling.log_normalise

    def spy(counts, total=None, scale=1e4):
        seen["total"] = total
        return real(counts, total=total, scale=scale)
    monkeypatch.setattr(coupling, "log_normalise", spy)
    ex = example_dir
    assert main(["couple", "--adt", str(ex / "adt.csv"), "--rna", str(ex / "rna.csv"),
                 "--meta", str(ex / "meta.csv"), "--normalise", "--draws", "2",
                 "--out", str(tmp_path / "o.csv")]) == 0
    meta = pd.read_csv(ex / "meta.csv", index_col=0)
    rna = pd.read_csv(ex / "rna.csv", index_col=0)
    assert (seen["total"].reindex(meta.index) == meta.nCount_RNA).all()
    assert (meta.nCount_RNA > rna.sum(axis=1)).all()      # the old row-sum scaling differs


def test_diagnose_groups_with_named_groups(example_dir, capsys):
    """Group names in the index used to make size_confounding raise."""
    rc = main(["diagnose", "--groups", str(example_dir / "groups.csv")])
    assert rc == 1
    assert "[FLAG] size_confounding" in capsys.readouterr().out


def test_diagnose_flags_the_example_panel(example_dir, capsys):
    rc = main(["diagnose", "--adt", str(example_dir / "adt.csv"),
               "--meta", str(example_dir / "meta.csv")])
    assert rc == 1 and "core panel = 14 markers" in capsys.readouterr().out


def test_diagnose_adt_without_meta_is_an_error(example_dir):
    with pytest.raises(SystemExit, match="both --adt and --meta"):
        main(["diagnose", "--adt", str(example_dir / "adt.csv")])


def test_version(capsys):
    from citeqc import __version__
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0 and __version__ in capsys.readouterr().out
