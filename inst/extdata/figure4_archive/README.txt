ProtVis archived B73-Y12 GSEA reproduction
=================================================

Purpose
-------
This archive supports the dedicated "B73-Y12 Figure reproduction" mode in
ProtVis GSEA. It is not a pre-rendered figure and it is not a precomputed
GSEA result table.

Bundled content
---------------
Root_VE_phenylpropanoid.pvg stores a compact representation of:
1. the frozen Root_VE B73-minus-Y12 ranked protein metric used to reproduce
   the supplied phenylpropanoid GSEA panel;
2. the positions and identifiers of Phenylpropanoid biosynthesis (KEGG 00940)
   proteins in that ranked list.

Runtime calculation
-------------------
ProtVis reconstructs the ranked vector, recalculates the weighted running
enrichment score (weight exponent p = 1), and runs fgseaMultilevel at runtime
with seed 20260920, minSize = 5, maxSize = 500, and eps = 0.

Reproduction audit
------------------
The supplied panel geometry is reproduced by the frozen full retained-protein
B73-minus-Y12 log2FC ranking:
- ranked proteins: 11,049
- pathway hits: 196
- positive ES maximum: 0.2769240517
- negative ES minimum / absolute extreme: -0.3064468069
- absolute extreme position: 9,893
- rank range: 11.324951 to -7.920156

This is intentionally kept separate from the regular ProtVis GSEA workflow,
which remains available as "Standard GSEA".
