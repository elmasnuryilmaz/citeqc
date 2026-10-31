# Installation

citeqc needs Python 3.10 or newer and depends only on NumPy, pandas and SciPy.

```bash
pip install citeqc
```

Check that the command-line tool is on your path:

```bash
citeqc --version
```

## Optional extras

| Extra | Adds | Needed for |
|---|---|---|
| `anndata` | `anndata`, `mudata` | reading `.h5mu` files with [`io.from_mudata`](api.md#citeqc.io.from_mudata) |
| `test` | `pytest` | running the test suite |
| `docs` | MkDocs, Material for MkDocs, mkdocstrings | building this site |

```bash
pip install 'citeqc[anndata]'
```

Reading `.parquet` tables with [`io.read_matrix`](api.md#citeqc.io.read_matrix)
also needs `pyarrow` (`pip install pyarrow`).

## From source

```bash
git clone https://github.com/elmasnuryilmaz/citeqc
cd citeqc
pip install -e '.[test,anndata]'
pytest
```

To build the documentation locally:

```bash
pip install -e '.[docs]'
mkdocs serve
```

## Supported versions

Tests run on Python 3.10 to 3.14, and once on Python 3.10 with the oldest
dependency releases citeqc claims to support (NumPy 1.23, pandas 1.5,
SciPy 1.10).
