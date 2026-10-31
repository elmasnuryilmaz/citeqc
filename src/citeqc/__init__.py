"""citeqc: protein/mRNA coupling and artefact diagnostics for CITE-seq data.

Two entry points:

    from citeqc import coupling, diagnostics

`coupling.measure()` quantifies, for each antibody, how well its signal is
predicted by its cognate transcript, answering the question "can I read this
marker off the transcriptome?".

`diagnostics.*` implement five checks for artefacts that generate large,
plausible-looking false positives in multimodal single-cell analysis.
Each check is built around a negative control rather than a threshold.

`datasets.synthetic()` simulates a small CITE-seq dataset with a known answer,
for trying the package without downloading anything.
"""
__version__ = "0.2.0"

from .mapping import load_map, map_antibodies, MECHANISMS
from .coupling import measure, specificity, clr
from . import diagnostics, datasets

__all__ = ["measure", "specificity", "clr", "load_map", "map_antibodies", "MECHANISMS",
           "diagnostics", "datasets", "__version__"]
