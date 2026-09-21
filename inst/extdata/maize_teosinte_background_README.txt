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
uses the current DEP result and the proteins actually tested in each comparison
as that comparison's enrichment universe.
