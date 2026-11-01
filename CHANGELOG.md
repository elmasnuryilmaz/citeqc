# Changelog

All notable changes to citeqc are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [0.2.0] - 2026-11-01

### Fixed

- `diagnostics.size_confounding` paired each group with a size by index label
  instead of by row position. Groups with any index other than `0..n-1` were
  affected: group names (what `citeqc diagnose --groups` always reads) raised
  "no feature had variance to correlate", and a shuffled integer index paired
  groups with the wrong sizes and could turn a flag into a pass. Sizes now
  follow the row order of `values`; a Series of sizes with the same labels is
  matched by label.
- `diagnostics.panel_consistency` compared panel *sizes*, so donors measuring
  different antibodies passed as long as they measured the same number. It now
  compares the measured sets (new `detail["n_panels"]`). It also raises a clear
  error when `adt` and `meta` share no cells, instead of an `IndexError`.
- `io.from_mudata` built `nCount_ADT` / `nCount_RNA` from a 2-D matrix when
  `.X` was sparse, and misaligned them when the two modalities held different
  cells. It now keeps the cells present in both modalities and computes the
  totals per cell.
- `io.read_matrix` read files such as `sample.v2.tsv.gz` as comma-separated.
  Tab detection now looks at the last suffix before any compression suffix,
  case-insensitively.
- `citeqc couple --normalise` scaled RNA by the row sums of the input table,
  which is not counts per 10,000 of the library when the table holds a subset
  of genes. It now uses `nCount_RNA` from `--meta`.
- `citeqc diagnose` skipped `panel_consistency` without a word when only one
  of `--adt` / `--meta` was given; this is now an error.
- The declared minimum SciPy version was too low: SciPy 1.9 has no
  `spearmanr(...).statistic`. The requirement is now `scipy>=1.10`, and CI
  tests the minimum versions.

### Added

- `citeqc.datasets`: a small synthetic CITE-seq dataset whose design is known
  by construction (`synthetic`, `DESIGN`), a null for group-size effects
  (`pseudoclonotypes`), and `write_example`; `citeqc example` writes it as CSV.
- `examples/quickstart.py`, run in CI.
- `coupling.measure(..., min_obs=10)`: the minimum number of observations per
  marker was hard-coded; the error now says how many were reached.
- `citeqc --version`, and help text for every command-line option.
- Documentation site (MkDocs Material) with a quick start, CLI and API
  reference, and an explanation of each diagnostic.
- Release workflow publishing to PyPI with Trusted Publishing on version tags
  (TestPyPI on manual runs), and a workflow deploying the docs to GitHub Pages.

### Changed

- Requires Python 3.10 or newer; Python 3.9 reached end of life in October
  2025. Tested on Python 3.10 to 3.14 and with the minimum dependency versions.
- Package metadata: SPDX licence expression with the licence file, project
  URLs, Python version classifiers. The version is set in one place,
  `citeqc.__version__`.
- `tests/validate_on_real_data.py` is now in English.

## [0.1.1] - 2026-09-30

### Added

- `coupling.specificity()`, `measure(..., return_raw=True)` and
  `diagnostics.control_check` (`ambient_check` kept as an alias).

## [0.1.0] - 2026-09-27

- Initial release (as `citecheck`, renamed to `citeqc` the same day).

[0.2.0]: https://github.com/elmasnuryilmaz/citeqc/compare/7d57fd7...HEAD
[0.1.1]: https://github.com/elmasnuryilmaz/citeqc/commit/efbaaad
[0.1.0]: https://github.com/elmasnuryilmaz/citeqc/commit/e7b4850
