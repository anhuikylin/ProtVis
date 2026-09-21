testthat::test_that("transcriptome GSEA background contains phenylpropanoid pathway", {
  bg <- ProtVis:::.protvis_deg_load_gsea_background()
  testthat::expect_true(nrow(bg$TERM2GENE) > 50000)
  testthat::expect_true(nrow(bg$TERM2NAME) > 150)
  row <- bg$TERM2NAME[bg$TERM2NAME$TERM == "00940", , drop = FALSE]
  testthat::expect_equal(nrow(row), 1L)
  testthat::expect_identical(
    row$NAME[[1L]],
    "Phenylpropanoid biosynthesis"
  )
})


testthat::test_that("transcriptome GSEA rank follows RNAseq.R log2FC default", {
  x <- data.frame(
    GeneID = c("G1", "G2", "G3", "G4"),
    log2FoldChange = c(-2, 3, 0.5, 1),
    stat = c(-4, 2, 1, 3),
    pvalue = c(0.01, 0.02, 0.4, 0.05),
    stringsAsFactors = FALSE
  )
  ranks <- ProtVis:::.protvis_deg_prepare_gsea_rank(
    x,
    metric = "log2FoldChange"
  )
  testthat::expect_identical(
    names(ranks),
    c("G2", "G4", "G3", "G1")
  )
  testthat::expect_equal(
    unname(ranks),
    c(3, 1, 0.5, -2)
  )
})


testthat::test_that("historical transcriptome GSEA 00940 reference is bundled", {
  ref <- ProtVis:::.protvis_deg_load_gsea_reference()
  testthat::expect_identical(ref$ID[[1L]], "00940")
  testthat::expect_identical(
    ref$Description[[1L]],
    "Phenylpropanoid biosynthesis"
  )
  testthat::expect_equal(ref$setSize[[1L]], 328)
  testthat::expect_equal(
    ref$NES[[1L]],
    -1.360826149505172,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    ref$p.adjust[[1L]],
    0.1009765165855614,
    tolerance = 1e-12
  )
})


testthat::test_that("transcriptome GSEA assets do not bundle expression or group tables", {
  path <- system.file(
    "extdata",
    "transcriptome_gsea",
    package = "ProtVis"
  )
  if (!nzchar(path)) {
    path <- file.path(
      "inst",
      "extdata",
      "transcriptome_gsea"
    )
  }
  files <- list.files(path)
  testthat::expect_false(any(grepl(
    "expression|group|sample",
    files,
    ignore.case = TRUE
  )))
})


testthat::test_that("transcriptome GSEA keeps the complete DESeq2 ranked universe", {
  source <- paste(
    deparse(ProtVis:::.protvis_deg_run_gsea),
    collapse = "\n"
  )
  testthat::expect_match(
    source,
    "complete",
    ignore.case = TRUE
  )
  testthat::expect_false(grepl(
    "ranks <- ranks\\[annotated\\]",
    source
  ))
})
