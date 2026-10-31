"""citeqc quick start on the bundled synthetic dataset (simulated, not real data).

Runs in a few seconds:

    python examples/quickstart.py

Every antibody's behaviour is known by construction (see citeqc.datasets.DESIGN),
so you can check each result against the truth.
"""
import numpy as np
import pandas as pd
from citeqc import coupling, datasets, mapping, diagnostics as dx

pd.set_option("display.width", 120)

# 1. Data: raw counts, cells in rows. meta holds donor, celltype and library sizes.
adt, rna, meta = datasets.synthetic(seed=0)
print(f"{len(adt)} cells, {adt.shape[1]} antibodies, {rna.shape[1]} genes, "
      f"{meta.donor.nunique()} donors x {meta.celltype.nunique()} cell types")

# 2. Is the same antibody panel measured in every donor? Check before anything else.
panel = dx.panel_consistency(adt, meta, by="donor")
print("\n" + str(panel))
print("  panel sizes per donor:", panel.detail["panel_sizes"])

# 3. Normalise: CLR for antibodies, log counts per 10k of the library for RNA.
adt_n = pd.DataFrame(coupling.clr(adt.values), index=adt.index, columns=adt.columns)
rna_n = coupling.log_normalise(rna, total=meta["nCount_RNA"])

# 4. Pair antibodies with their cognate genes. Unknown names are reported, never guessed.
pairs, unmapped = mapping.map_antibodies(adt_n.columns, rna_n.columns)
print(f"\nmapped {len(pairs)}/{adt.shape[1]} antibodies; unmapped: {unmapped}")

# 5. Coupling: within donor x cell type, equal-sized pseudo-groups, depth regressed out.
res = coupling.measure(adt_n, rna_n, meta, pairs)
res["pct_rank"] = coupling.rank_percentile(res)
res["mechanism"] = [mapping.mechanism_of(m) for m in res.marker]
res["design"] = [datasets.DESIGN[m][1] for m in res.marker]
print("\nprotein~mRNA coupling (design = the truth the data were simulated from):")
print(res[["marker", "gene", "median_rho", "q25", "q75", "n_strata", "pct_rank",
           "mechanism", "design"]].round(3).to_string(index=False))

# 6. Specificity: does each protein track its own transcript better than other genes?
_, spec = coupling.specificity(adt_n, rna_n, meta, pairs, genes=list(rna_n.columns))
print("\nwhere the cognate gene ranks among all 25 genes:")
print(spec[["marker", "cognate_rho", "spec_percentile", "rank_of_cognate"]]
      .round(3).to_string(index=False))

# 7. Negative controls: CD19 and CD20 are absent from T cells and must not "couple".
print("\n" + str(dx.control_check(res.set_index("marker").median_rho,
                                   control_markers=["CD19", "CD20"])))

# 8. Group size: pseudo-clonotypes of 5-300 cells drawn from the same cells, so no
#    gene truly depends on group size. The pseudobulk still drifts with size ...
values, sizes, members = datasets.pseudoclonotypes(rna, meta)
print("\n" + str(dx.size_confounding(values, sizes)))

# ... and averaging exactly 5 cells per group removes the drift.
cp10k = rna.div(meta["nCount_RNA"], axis=0) * 1e4
draw = next(dx.equal_n_resample(members, n=5, n_draw=1, seed=0))
equal_n = pd.DataFrame({g: np.log1p(cp10k.loc[cells].mean()) for g, cells in draw.items()}).T
print(dx.size_confounding(equal_n, sizes.loc[equal_n.index]))
