Built-in maize-teosinte KEGG Figure 3C-D reproduction data
============================================================

File
----
maize_teosinte_kegg_figure3_cd.csv

Purpose
-------
This table is the built-in numeric source for the Enrichment analysis >
KEGG Enrichment Analysis > Maize-teosinte reproduction tab.

Archived reference
------------------
03.Maize_Teosinte_Jul02_2024/04.result/01.Publish_figures/KEGG enrichment.png

Numeric reconstruction / verification source
--------------------------------------------
The displayed pathway-stage points were reconstructed from the archived
directional DEP + KEGG workflow and checked against the archived Figure 3C-D
image. The archived workflow uses:
- five B73 versus Y12 developmental comparisons;
- significant DEPs defined by adj.P.Val < 0.05 and abs(logFC) > 1;
- positive logFC as B73-higher and negative logFC as Y12-higher;
- clusterProfiler::compareCluster(fun = "enricher");
- KEGG TERM2GENE/TERM2NAME from Enrichmentdb2.xlsx;
- pvalueCutoff = 0.05 and qvalueCutoff = 1;
- raw enrichment p value as the dot colour and gene Count as dot size.

Integrity checks
----------------
- Total displayed enrichment points: 29
- Panel C (Zea mays ssp. mays enriched): 18
- Panel D (Zea mays ssp. mexicana enriched): 11
- Panel D intentionally has no Root_V1.V2 point because no pathway from that
  directional comparison is displayed after the archived enrichment filtering.

The application validates these counts, pathway labels, keys, and numeric
values at load time; regression tests are in:
tests/testthat/test-enrichment-maize-teosinte-reproduction.R
