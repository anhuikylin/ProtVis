ProtVis GSEA reproduction benchmark
=====================================

Purpose
-------
The active "B73-Y12 Figure reproduction" workflow now loads the complete
previous DEP result tables produced by ProtVis Differential Analysis. Multiple
B73-vs-Y12 comparisons can be loaded and analysed together; Root_VE is selected
as the default focus comparison when it is present.

The source GSEA_phenylpropanoid_B73_Y12.R workflow is reproduced by using:
- the complete tested-protein DEP result for each comparison;
- moderated limma t statistic as the default ranking metric;
- no significance pre-filter by default;
- comparison-specific intersection with the bundled Enrichmentdb2 annotation;
- Phenylpropanoid biosynthesis (KEGG map00940);
- weighted running enrichment score, p = 1;
- fgseaMultilevel with seed 20260920, minSize = 5, maxSize = 500, eps = 0.

Users can change the ranking metric or apply custom P/FDR/log2FC inclusion
thresholds, but these settings intentionally move away from the exact source
script. Exact figure reproduction additionally requires the loaded DEP result
to come from the archived historical limma workflow.

Legacy visual benchmark
-----------------------
Root_VE_phenylpropanoid.pvg is retained only as a compact regression benchmark
for the previously supplied figure geometry. It is not the active GSEA data
source and it is not a precomputed GSEA result table or static image.
