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

No precomputed KEGG enrichment result or static plot is included. Directional
KEGG is recalculated at runtime.

"Figure 3 reproduction (archived DEP)" uses the five significant DEP protein
lists written by the original 02.MaizeTeosintePro/01.src/02.depforseedling.R
workflow to 03.progress/03.dep (Root_VE, Root_V1.V2, Root_V4,
Leaf_VE.V1.V2, and Leaf_V4.V6.V8). The lists are bundled in a compact
membership representation and internally audited against the historical
direction counts. ProtVis then reruns compareCluster/enricher with the bundled
Enrichmentdb2 background, pvalueCutoff = 0.05, qvalueCutoff = 1, no
user-supplied universe, raw pvalue colour, and clusterProfiler-compatible
dotplot selection. This avoids treating a modern Step4/Step6 reconstruction as
identical to the frozen historical DEP analysis while still calculating KEGG
from protein lists rather than bundling a precomputed enrichment table or
static plot.

"Standard ORA" remains available"Standard ORA" remains available for general datasets and uses each
comparison's retained tested proteins as the enrichment universe with BH FDR
filtering.
