Built-in maize-teosinte enrichment background

Source: Enrichmentdb2.xlsx supplied with the ProtVis maize-teosinte data.
Its bytes were verified against the copy in 02.MaizeTeosintePro.zip.

The bundled tables are deliberately in the same table format produced by
Toolkits > Background Make:

  * GO_background: TERM, GENE, NAME
  * KEGG_background: TERM, GENE, NAME

They are stored as compact tab-separated .xz package data so the installed
package remains small. ProtVis reads them as the two standard tables above;
users can still upload or download Background Make workbooks in .xlsx form.
They were made by joining the source TERM-to-GENE sheets (t2g.go and
t2g.kegg) to the first non-empty pathway/term name in the corresponding t2n
sheets. Rows without a term, gene, or name are not included.

No precomputed enrichment result or static plot is included. Directional KEGG
is calculated from the current DEP result.

"Figure 3 reproduction (compareCluster)" reconstructs the historical workflow:
archived Step6 limma, the Step4 any-detected protein filter, BH < 0.05 and
|log2FC| > 1 directional lists, compareCluster/enricher with pvalueCutoff =
0.05 and qvalueCutoff = 1, no user-supplied universe, raw pvalue colour, and
clusterProfiler-compatible dotplot selection. If the current DEP was produced
with the recommended workflow, Directional KEGG automatically rebuilds the
archived DEP from Step4_data_transformed.rda + Step6_data_normalization.rda
without overwriting the user's current DEP results.

"Standard ORA" remains available for general datasets and uses each
comparison's retained tested proteins as the enrichment universe with BH FDR
filtering.
