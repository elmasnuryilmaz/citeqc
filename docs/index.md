# citeqc

**Protein-mRNA coupling and artefact diagnostics for CITE-seq data.**

citeqc answers two questions about a CITE-seq experiment:

1. **Can this surface marker be read off the transcriptome?** For each
   antibody, [`coupling.measure`](api.md#citeqc.coupling.measure) measures how
   well its signal tracks its own (cognate) transcript within homogeneous
   groups of cells, and [`coupling.specificity`](api.md#citeqc.coupling.specificity)
   asks whether the cognate gene stands out among other genes.
2. **Is my analysis producing confident but false positives?** Five
   [diagnostics](diagnostics.md) each use a negative control the data already
   contains to ask *what would this analysis report if there were no signal?*

Why the first question matters: some surface proteins are set by mechanisms
that the transcript cannot see. CD45RA is an alternative splice product of
*PTPRC*, so an antibody separates CD45RA from CD45RO while total *PTPRC* mRNA
cannot. CTLA-4 (CD152) and CD40L (CD154) are held in intracellular vesicles
until the cell is triggered, so their surface level reflects trafficking
rather than transcription. The bundled antibody-to-gene table annotates
markers with such documented mechanisms.

## Where to start

- [Installation](installation.md): `pip install citeqc`.
- [Quick start](quickstart.md): run the whole workflow in seconds on a bundled
  synthetic dataset whose answer is known.
- [How coupling is measured](coupling.md): the three design constraints the
  measurement enforces, and how to read the numbers.
- [Diagnostics](diagnostics.md): what each check guards against and how it
  decides.
- [CLI reference](cli.md) and [Python API](api.md).

## What citeqc does not do

It does not denoise antibody counts (see totalVI or dsb), correct ambient RNA
(see SoupX or CellBender, neither of which addresses the antibody signal), or
normalise beyond two simple helpers. It measures and it checks. Run the
measurement on the normalisation you trust, and compare ranks between
normalisations.

## Citing

If you use citeqc, please cite it using the metadata in
[`CITATION.cff`](https://github.com/elmasnuryilmaz/citeqc/blob/main/CITATION.cff)
(GitHub's "Cite this repository" button reads it), the accompanying paper (in
preparation), and the datasets you analyse.
