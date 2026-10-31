# CLI reference

`citeqc` has four subcommands:

| Command | Does |
|---|---|
| `citeqc couple` | measure protein~mRNA coupling per antibody and write a CSV |
| `citeqc diagnose` | run the [diagnostics](diagnostics.md) whose inputs you pass; exit code 1 if any flags |
| `citeqc map` | search the bundled antibody-to-gene table |
| `citeqc example` | write the [synthetic example dataset](quickstart.md) as CSV files |

**Input tables** are CSV (or TSV for `.tsv`/`.txt`; `.gz` and other
compression is detected) with the row labels in the first column. For
`couple`, `adt`, `rna` and `meta` have **cells in rows**; pass `--transpose` if
`adt` and `rna` are features x cells. `meta` needs the stratum columns and the
raw library sizes `nCount_ADT` and `nCount_RNA`.

For `diagnose`:

| Option | Table |
|---|---|
| `--adt` + `--meta` | as for `couple`; missing antibodies must be empty (NaN), not 0 |
| `--stat` | markers in rows, one column: your per-marker statistic |
| `--groups` | groups in rows, features in columns, plus a size column (`--size-col`, default `n_cells`) |
| `--pvalues` | one column of p-values |
| `--bounded` | features x samples, values in [0, 1] |

The help text below is generated from the code when the site is built.

<!-- citeqc-cli-help -->
