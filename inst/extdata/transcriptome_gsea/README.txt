ProtVis transcriptome GSEA assets
================================

Runtime KEGG background
-----------------------
ProtVis reuses inst/extdata/maize_teosinte_KEGG_background.tsv.xz. That file
was generated from 02.MaizeTeosintePro/02.data/Enrichmentdb2.xlsx and its
source bytes were previously verified against 02.MaizeTeosintePro.zip.

Historical method
-----------------
The transcriptome GSEA implementation follows:
  02.MaizeTeosintePro/01.src/RNAseq.R

The archived code ranks the complete DESeq2 B73-vs-Y12 result by
log2FoldChange and passes that complete ranked universe directly to
clusterProfiler::GSEA. KEGG annotation defines pathway membership; genes
without KEGG annotation remain in the ranked list. The historical settings are
minGSSize = 5, maxGSSize = 500, pvalueCutoff = 1,
TERM2GENE = t2g.kegg and TERM2NAME = t2n.kegg.

Regression reference
--------------------
RNAseq_GSEA_reference_00940.tsv contains the historical Phenylpropanoid
biosynthesis (00940) result extracted from:
  02.MaizeTeosintePro/03.progress/04.RNAseq/KEGGenrich.xlsx

This reference is for validation only; it is not used as the active result.

Archive audit
-------------
03.Maize_Teosinte_Jul02_2024.zip was also inspected. It contains the later
five-comparison proteomics workflow and updated B73/Y12 naming context, but no
RNA-seq GSEA script or transcriptome expression/group table. Therefore no
transcriptome GSEA method is falsely attributed to that archive.

Expression/group data
---------------------
No Gene Expression Data, group-information table, or sample-information table
is bundled for this toolkit. Users upload those files at runtime.
