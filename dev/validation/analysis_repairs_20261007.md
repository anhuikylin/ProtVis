# Analysis repairs and executed validation

Base: anhuikylin/ProtVis dev, 59fe3d073dc98eb53d9719168df87f6eb01f748a.

## Changes

- Apply the requested post-search cutoff to rank-one target PSMs at spectrum,
  peptide and protein q-value levels. LFQ rows must pass their own q-value
  cutoff and match an accepted peptide/protein assignment. Retain unfiltered
  tables and source files. Recovered runs are filtered using the same rules.
- Sum positive LFQ intensities from unique target peptides only. Keep shared
  peptides in retained LFQ tables without duplicating them across proteins.
  All-zero or missing protein/sample observations remain missing. Require
  unambiguous file-to-sample mapping. This is unique-peptide sum quantification,
  not MaxLFQ or parsimony-based protein-group inference.
- Use Group1 minus Group2 for programmatic log2FC, consistent with the existing
  Shiny contrasts. Check log2 input scale before programmatic and matrix-based
  Shiny differential analyses. Record scale through preprocessing and historical
  matrix activation; reject raw, ln/log10, standardized and VSN values.
- Real Sage execution exposed single-file JSON array serialization and optional
  modification-map issues, also repaired. Distinguish Trypsin from Trypsin/P.
  The runner explicitly uses Sage's documented telemetry-disable CLI flag.
- Real downstream execution exposed a mismatch between the native annotation
  field names and the pipeline enrichment reader. Support GO_annotation and
  KEGG_annotation, combine mappings, and restrict selected IDs to the annotated
  universe.
- Correct two pre-existing test expectations: limma's archived logFC sorting is
  by absolute effect size, and anyDuplicated returns integer zero. Both failures
  were reproduced before those test changes; archived analysis code is unchanged.

## Executed checks

Environment: Ubuntu 24.04, R 4.3.3, limma 3.58.1, statmod 1.5.0.
Core R files were sourced directly; this was not an installed-package R CMD
check or an interactive Shiny browser test.

- All package R files and the reusable validation script parse successfully.
- 30 new regression assertions pass, with no warnings or skips.
- 87 existing relevant DEP/schema-v4 assertions pass, with no warnings or skips.
- Real B73/Y12 quantitative matrix: 12,300 rows after contamination/site/reverse
  flag filtering, 30 samples; 12,115 rows remain after the validation missingness
  filter. Executed transformation, imputation, normalization, PCA, five
  differential comparisons, enrichment, correlation network, exports and
  checkpoint restoration. Restoration passed through the public checkpoint API.
- Recommended DEP was tested separately on observed log2 values before
  imputation, sample-median centered, minimum two observations in each group,
  robust/trended limma, BH <= 0.05 and absolute log2FC > 1.
- On identical complete observed input, programmatic mean differences and limma
  effects agree to a maximum absolute discrepancy of 1.11e-15.

| Comparison (B73 vs Y12) | Quantitatively tested | Significant quantitative proteins | Single-group detections |
|---|---:|---:|---:|
| Root_VE | 7,575 | 87 | 1,948 |
| Root_V1.V2 | 7,578 | 150 | 1,948 |
| Root_V4 | 7,554 | 263 | 1,956 |
| Leaf_VE.V1.V2 | 7,567 | 221 | 1,957 |
| Leaf_V4.V6.V8 | 7,546 | 107 | 1,964 |

Single-group detections are descriptive presence/absence evidence, not proteins
passing the quantitative BH test. The generic pipeline retains its Welch-test
engine and imputed-matrix input; its significant counts can differ from
Recommended DEP. This repair aligns direction and scale validation, not the
statistical engines. Zero selected IDs legitimately produce a skipped ORA result.

Recommended DEP lists also reached native ORA with comparison-specific tested
backgrounds: all five comparisons returned 430 mapped terms with finite P
values. Aggregated results are in recommended_enrichment_real_data_summary.csv.
These validation settings/results do not replace manuscript-specific analyses.

## Sage execution and LFQ validation

- Official Sage v0.14.7 Linux GNU release asset executed successfully against
  Sage's real single-scan PXD016766 fixture (Q99536 FASTA). The asset reports
  version 0.14.6 in results.json. One raw PSM was written, zero passed the 1%
  confidence gates, and zero LFQ rows were retained. No pseudo-count protein
  matrix was invented. The optional-modifications-disabled configuration also
  executed successfully.
- The initial run's default telemetry attempt was blocked by automatic review.
  Subsequent runs explicitly disabled telemetry using Sage's native CLI option.
- A real mspms Sage LFQ output contained 1,555 peptide rows and 12 samples.
  969 rows passed LFQ q_value <= 0.01 and produced 222 protein rows under the
  unique-peptide sum policy. This tests roll-up of genuine LFQ output; it does
  not validate its upstream spectrum/peptide/protein q-values without PSM data.
- Sage 0.14.7 internally seeds LFQ tracing from peptide_q <= 0.01. Raising the
  post-search cutoff cannot recover features the engine never traced.

Sources:
https://github.com/lazear/sage/tree/v0.14.7/tests
https://github.com/lazear/sage/blob/v0.14.7/crates/sage/src/lfq.rs
https://github.com/lazear/sage/blob/v0.14.7/crates/sage-cli/src/main.rs
https://github.com/baynec2/mspms/blob/master/inst/extdata/sage_lfq.tsv

Input MD5:
- verified B73/Y12 matrix: c526b70b12c5be107ffcd34d60b7b52a
- Enrichmentdb2.xlsx: 990912ef1f09d90279e29b43dee0c288
- real Sage LFQ table: 8019d6e66c5e85a2cb272fc20d4ce5c1

## Outstanding full raw-data validation

A continuous replicated mzML -> Sage -> filtered LFQ -> protein matrix ->
differential/enrichment run has NOT been completed. The available mzML fixture
contains only one scan and produces no high-confidence LFQ values. Completing
that validation requires matched mzML files from a replicated two-group
experiment, the matching FASTA, and file-to-sample/group metadata. Synthetic
replicate files or fabricated identification confidence were not substituted.

For the completed B73/Y12 matrix workflow, rerun:

```sh
Rscript dev/validate_analysis_repairs.R matrix.csv terms.csv output_directory
```

terms.csv must contain protein_id and term, one mapping per row. The script
accepts the verified matrix's Protein_ID, filter flags and B73/Y12 sample names.
