test_that("Sage filters target rank-one identifications and LFQ confidence", {
  psms <- data.frame(peptide = c("AA", "BB", "CC", "DD", "EE"),
    proteins = c("P1", "P2", "P3", "rev_P4", "P5"),
    label = c(1, 1, 1, -1, 1), rank = c(1, 1, 2, 1, 1),
    spectrum_q = c(.001, .02, .001, .001, .001),
    peptide_q = .001, protein_q = c(.001, .001, .001, .001, .03))
  lfq <- data.frame(peptide = psms$peptide, proteins = psms$proteins,
    q_value = .001, s.mzML = 10)
  bundle <- .protvis_sage_filter_bundle(list(psms = psms, lfq_table = lfq), .01)
  expect_identical(bundle$psms$peptide, "AA")
  expect_identical(bundle$lfq_table$peptide, "AA")
  expect_equal(nrow(bundle$psms_raw), 5)
  relaxed <- .protvis_sage_filter_bundle(bundle, .05)
  expect_equal(nrow(relaxed$psms), 3)
  expect_error(.protvis_sage_filter_bundle(list(psms = psms[, -6], lfq_table = lfq)), "missing")
  lfq$q_value[1] <- NA
  expect_equal(nrow(.protvis_sage_filter_bundle(list(psms = psms, lfq_table = lfq))$lfq_table), 0)
})

test_that("Sage excludes shared evidence, preserves missingness, and maps samples", {
  lfq <- data.frame(proteins = c("P1", "P1", "P1;P2", "P2", "rev_P3"),
    a.mzML = c(10, 20, 1000, 0, 90), b.mzML = c(0, NA, 2000, 40, 50))
  info <- data.frame(sample_id = c("A", "B"), mzml_file = c("a.mzML", "b.mzML"))
  x <- .protvis_sage_lfq_protein_matrix(lfq, info, c("a.mzML", "b.mzML"))
  expect_equal(x["P1", "A"], 30)
  expect_true(is.na(x["P1", "B"]))
  expect_true(is.na(x["P2", "A"]))
  expect_equal(x["P2", "B"], 40)
  expect_equal(rownames(x), c("P1", "P2"))
  expect_error(.protvis_sage_lfq_protein_matrix(lfq, info, "missing.mzML"), "match")
})

test_that("pipeline log2FC has the same Group1-minus-Group2 sign as Shiny", {
  x <- rbind(P1 = c(12, 13, 11, 8, 9, 7), P2 = c(8, 7, 9, 12, 11, 13))
  colnames(x) <- paste0("s", 1:6)
  info <- data.frame(sample_id = colnames(x), group = rep(c("A", "B"), each = 3))
  dataset <- create_protvis_dataset(x, sample_info = info, metadata = list(expression_scale = "log2"))
  result <- .protvis_differential_analysis(dataset, list(group1 = "A", group2 = "B"))
  table <- result$analysis_results$differential_analysis$table
  expect_equal(table$log2FC[match("P1", table$protein_id)], 4)
  expect_equal(table$log2FC[match("P2", table$protein_id)], -4)
  for (scale in c("raw", "ln", "log10", "standardized", "vsn")) {
    bad <- .protvis_set_expression_scale(dataset, scale)
    expect_error(.protvis_differential_analysis(bad, list(input_scale = "log2")), "requires log2")
  }
  unknown <- dataset
  unknown$metadata$expression_scale <- NULL
  expect_error(.protvis_require_log2_scale(unknown), "requires log2")
  expect_invisible(.protvis_require_log2_scale(unknown, "log2"))
})

test_that("native GO/KEGG annotations reach pipeline enrichment", {
  x <- matrix(1:12, 3, dimnames=list(c('P1','P2','P3'), paste0('s',1:4)))
  dataset <- create_protvis_dataset(x, annotation=list(
    GO_annotation=data.frame(protein_id=c('P1','P2'),term='GO:1'),
    KEGG_annotation=data.frame(protein_id=c('P1','P3'),term='map1')))
  result <- .protvis_enrichment(dataset, list(selected_ids=c('P1','unmapped')))
  expect_equal(nrow(result$analysis_results$enrichment$table), 2)
  expect_identical(result$analysis_results$enrichment$selected, 'P1')
  expect_true(all(is.finite(result$analysis_results$enrichment$table$p_value)))
})

test_that("Sage distinguishes Trypsin/P and validates the requested FDR", {
  folder <- tempfile(); dir.create(folder)
  fasta <- file.path(folder,'p.fasta'); mzml <- file.path(folder,'s.mzML')
  file.create(fasta,mzml)
  trypsin <- protvis_sage_build_config(fasta,mzml,folder,list(enzyme='Trypsin'))
  unrestricted <- protvis_sage_build_config(fasta,mzml,folder,list(enzyme='Trypsin/P'))
  expect_identical(trypsin$database$enzyme$restrict, 'P')
  expect_null(unrestricted$database$enzyme$restrict)
  expect_error(protvis_sage_build_config(fasta,mzml,folder,list(fdr=NA)), 'FDR')
  expect_error(protvis_sage_build_config(fasta,mzml,folder,list(fdr=0)), 'FDR')
})

test_that("reactivating a raw matrix cannot inherit a later log2 scale", {
  x <- matrix(1:12,3,dimnames=list(paste0('P',1:3),paste0('s',1:4)))
  dataset <- create_protvis_dataset(x, metadata=list(expression_scale='raw'))
  initial <- dataset$workflow$active_matrix_run_id
  dataset <- .protvis_set_expression_scale(dataset, 'log2')
  activated <- protvis_activate_result(dataset, initial)
  expect_identical(activated$metadata$expression_scale, 'raw')
  expect_error(.protvis_require_log2_scale(activated), 'requires log2')
})
