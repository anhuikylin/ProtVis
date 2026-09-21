ProtVis GSEA reproduction benchmark
=====================================

Two reproduction paths are intentionally available.

1. Exact supplied Root_VE figure (default)
------------------------------------------
The supplied Phenylpropanoid biosynthesis panel is matched by the bundled
historical Root_VE full retained-protein ranking:
- 11,049 ranked proteins;
- 196 map00940 hits;
- historical B73-minus-Y12 log2FC rank metric;
- no truncation of the ranked list to KEGG-annotated proteins;
- weighted running enrichment score, p = 1;
- seed 20260920, minSize = 5, maxSize = 500, eps = 0;
- absolute ES extreme at position 9,893, ES = -0.3064468069.

This archive is not a static plot or a precomputed GSEA table. ProtVis
reconstructs the ranked metric and hit positions and recalculates the running
ES and fgseaMultilevel result at runtime.

2. Recalculate from loaded DEP result
-------------------------------------
This path follows the current GSEA_phenylpropanoid_B73_Y12.R script logic:
- load the complete tested-protein DEP table;
- use moderated limma t by default;
- no significance pre-filter by default;
- intersect the ranked list with the KEGG-annotated tested background;
- run Phenylpropanoid biosynthesis (map00940) weighted GSEA.

The two paths are kept separate because the supplied historical panel geometry
does not match the currently loaded recommended DEP result or the current
script's annotated-background/t-ranking combination.

Multiple comparisons
--------------------
The interface can load and analyse several previous DEP comparisons at once.
Root_VE is the default focus. In exact-figure mode, Root_VE uses the historical
figure archive while other selected comparisons are recalculated from the
loaded DEP tables.
