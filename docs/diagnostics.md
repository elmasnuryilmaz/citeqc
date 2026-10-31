# Diagnostics

Five checks for artefacts that produce confident, biologically plausible false
positives. Each is built around the same question: *what would this analysis
report if there were no signal at all?* None of them needs a ground truth;
they use negative controls that the data already contain.

A check **flags** a problem; it does not say where the problem comes from.

| Check | What goes wrong | How you can tell |
|---|---|---|
| [`panel_consistency`](#panel_consistency) | Antibody panels differ between donors; "absent" and "not measured" get confused | Donors measure different sets of antibodies; some donors are unstained |
| [`control_check`](#control_check) | The statistic carries a marker-independent gradient | Lineage-impossible markers (CD19 on T cells) score like genuine ones |
| [`size_confounding`](#size_confounding) | Group means drift with group size, inflating low-expressed features | Almost every feature moves the same direction with size |
| [`pvalue_sanity`](#pvalue_sanity) | Constant or saturated features flatten the null | p-values *depleted* near zero, or a spike in the interior |
| [`saturation_check`](#saturation_check) | Bounded features pinned at 0/1 carry no information but enter FDR | A large fraction is pinned across all samples |

Every check returns a [`Finding`](api.md#citeqc.diagnostics.Finding) with
`.passed`, a printable `.summary`, and `.detail`, a dict with the numbers
behind the decision.

```python
from citeqc import diagnostics as dx

dx.panel_consistency(adt, meta, by="donor")
dx.control_check(my_zscores, control_markers=["CD19", "CD20"])
dx.size_confounding(pseudobulk, sizes=n_cells_per_group)
dx.pvalue_sanity(my_pvalues)
dx.saturation_check(psi_matrix)
```

From the shell, `citeqc diagnose` runs every check whose inputs you pass and
exits with code 1 if any check flags, so it can stop a pipeline:

```bash
citeqc diagnose --adt adt.csv --meta meta.csv --pvalues p.csv
```

## `panel_consistency`

**Artefact.** Multi-site and multi-batch CITE-seq often uses different antibody
panels. When the tables are merged, an antibody that a donor was never stained
for looks like an antibody the donor does not express, and "marker absent"
silently becomes "marker not measured".

**How it decides.** Within each group of `by` (default `donor`), an antibody
counts as measured if more than `min_detect` (default 0.5) of the group's cells
have a non-missing value. Groups that measured nothing are reported separately
as unstained. The check passes only if every stained group measured the same
set of antibodies and no group is unstained.

**What to do.** Restrict cross-donor analyses to `detail["core_panel"]`, the
antibodies measured in every stained group, and exclude unstained groups
(`detail["empty_groups"]`).

!!! note
    Antibodies that were not measured must be **missing (NaN)**, not zero.
    Merging tools often fill missing columns with 0, which makes them
    indistinguishable from real zeros and invisible to this check.

## `control_check`

**Artefact.** A per-marker statistic (a Z-score, -log10 p, an effect size)
contains a gradient that has nothing to do with the markers, so markers that
cannot be expressed in the cells analysed score as if they were.

**How it decides.** You pass your statistic and a list of `control_markers`
that are truly absent from the cells analysed (CD19 and CD20 on T cells, CD3
on monocytes). The check computes the mean absolute statistic of the controls
divided by that of the other markers (or of `signal_markers`, if given), and
flags when this ratio reaches `ratio_threshold` (default 0.5).

Choose controls that are really absent: CD16, CD56 and CD161 occur on CD8
T-cell subsets and make poor controls there.

**What to do.** The check diagnoses; it does not explain. Test group size first
([`size_confounding`](#size_confounding), `equal_n_resample`); ambient signal
is the second suspect. If you try an ambient correction, judge it on control
markers *held out* of the correction. In the analysis this package was
developed on, resampling to equal N removed the gradient and regressing out an
ambient score did not (see the
[README](https://github.com/elmasnuryilmaz/citeqc#check-your-analysis)).

`ambient_check` is an alias of `control_check`, kept for existing code.

## `size_confounding`

**Artefact.** A per-group statistic (clonotype pseudobulk, cluster means)
depends on how many cells each group contains. A mean over more cells is less
noisy, and for skewed, zero-inflated features its typical value also shifts
with the number of cells averaged. The result looks like a gene programme that
scales with clone size, but it is produced by the averaging.

**How it decides.** For every feature, the Spearman correlation (by default)
between its value and `log10(size + 1)` across groups. The signature of the
artefact is a *pervasive one-directional* association: real biology does not
move almost every feature the same way. The check flags when the larger of the
fraction of positive and the fraction of negative correlations reaches
`skew_threshold` (default 0.75). `detail["direction"]` says which way.

`sizes` must be in the same row order as `values`. A pandas Series indexed by
the same group labels is matched by label.

**What to do.** Resample every group to the same number of cells with
[`equal_n_resample`](api.md#citeqc.diagnostics.equal_n_resample), recompute the
statistic, and run the check again. Groups smaller than `n` are dropped: that
exclusion is the price of removing the artefact, so report it. The
[quick start](quickstart.md) shows the flag and the fix on data where no gene
depends on group size.

## `pvalue_sanity`

**Artefact.** A test that cannot work produces a p-value distribution no null
or signal would. Commonly many features are constant or saturated, give p
close to 1 or pile up at one value, and dilute the multiple-testing correction.

**How it decides.** Needs at least 50 p-values in [0, 1]. For each `alpha`
(0.001, 0.01, 0.05) it computes observed / expected counts below `alpha`. At
the smallest `alpha`, a ratio below 0.8 is *depleted* (flag), above 1.5
*enriched* (signal), and anything in between *uniform* (no signal; passes). It
also estimates the null fraction `pi0` (Storey's estimator at lambda = 0.5)
and flags an interior spike: the tallest of 20 histogram bins more than 3x the
uniform height and not at either end.

**What to do.** A working test is never depleted near zero. Drop constant and
saturated features (see [`saturation_check`](#saturation_check)) and test
again.

## `saturation_check`

**Artefact.** Bounded values such as PSI or cell-type fractions pinned at 0 or
1 in every sample carry no information, but they still enter
multiple-testing correction and flatten the p-value distribution.

**How it decides.** `values` is features x samples; only features complete
across samples are considered. A feature is saturated if all its values are at
most `lo` (0.02) or all at least `hi` (0.98), and testable if it is not
saturated and its range exceeds `min_range` (0.05). The check flags when more
than `max_saturated` (25%) of features are saturated.

**What to do.** Test only `detail["testable"]`.
