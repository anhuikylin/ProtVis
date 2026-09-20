testthat::test_that("archived DEP defaults match 03.Maize_Teosinte_Jul02_2024", {
  defaults <- ProtVis:::.protvis_dep_archived_defaults()

  testthat::expect_equal(defaults$logfc, 1)
  testthat::expect_equal(defaults$fdr, 0.05)
  testthat::expect_identical(defaults$p_metric, "adj.P.Val")
  testthat::expect_identical(defaults$adjust_method, "BH")
  testthat::expect_identical(defaults$sort_by, "logFC")
  testthat::expect_true(defaults$matrix_shift)
  testthat::expect_identical(defaults$protein_universe, "shared_preknn")
  testthat::expect_equal(defaults$min_detected, 2L)
})

testthat::test_that("DEP classification uses adjusted P values and strict archived cutoffs", {
  result <- data.frame(
    logFC = c(1.1, -1.2, 1.0, 2, -2),
    P.Value = c(0.001, 0.001, 0.001, 0.001, 0.001),
    adj.P.Val = c(0.04, 0.03, 0.01, 0.05, 0.051),
    stringsAsFactors = FALSE
  )

  classified <- ProtVis:::.protvis_dep_classify(
    result, logfc = 1, cutoff = 0.05, p_metric = "adj.P.Val"
  )

  testthat::expect_identical(
    classified$regulation,
    c(
      "Upregulated",
      "Downregulated",
      "Not significant",
      "Not significant",
      "Not significant"
    )
  )
})

testthat::test_that("shared pre-KNN protein universe requires detection in both genotypes", {
  mat <- matrix(
    c(
      1, 1, NA, 2, 2, NA,
      1, NA, NA, 2, 2, 2,
      1, 1, 1, 2, NA, NA,
      1, 1, 1, 2, 2, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("P", 1:4),
      c("B73_1", "B73_2", "B73_3", "Y12_1", "Y12_2", "Y12_3")
    )
  )

  ids <- ProtVis:::.protvis_dep_shared_ids(
    mat,
    c("B73_1", "B73_2", "B73_3"),
    c("Y12_1", "Y12_2", "Y12_3"),
    min_detected = 2L
  )

  testthat::expect_identical(ids, c("P1", "P4"))
})

testthat::test_that("DEP summary counts up, down and not significant proteins", {
  results <- list(
    B73_Root_VE_vs_Y12_Root_VE = data.frame(
      regulation = c(
        "Upregulated", "Upregulated",
        "Downregulated", "Not significant"
      ),
      stringsAsFactors = FALSE
    )
  )
  out <- ProtVis:::.protvis_dep_summary_table(results)

  testthat::expect_equal(sum(out$Protein_number), 4)
  testthat::expect_equal(
    out$Protein_number[out$Direction == "Upregulated"], 2
  )
  testthat::expect_equal(
    out$Protein_number[out$Direction == "Downregulated"], 1
  )
})

testthat::test_that("archived limma helper returns BH adjusted table", {
  testthat::skip_if_not_installed("limma")

  mat <- matrix(
    c(
      5, 5.2, 5.1, 3, 3.1, 3.2,
      4, 4.1, 4.2, 4, 4.1, 4.2,
      2, 2.2, 2.1, 5, 5.1, 5.2,
      6, 6.1, 6.2, 6, 6.2, 6.1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("P", 1:4),
      c("B1", "B2", "B3", "Y1", "Y2", "Y3")
    )
  )

  out <- ProtVis:::.protvis_dep_run_limma_archived(
    mat,
    c("B1", "B2", "B3"),
    c("Y1", "Y2", "Y3"),
    "B73_Root_VE",
    "Y12_Root_VE",
    adjust_method = "BH",
    sort_by = "logFC",
    matrix_shift = TRUE
  )

  testthat::expect_true(all(c(
    "ID", "logFC", "P.Value", "adj.P.Val", "FC", "Group1", "Group2"
  ) %in% names(out)))
  testthat::expect_true(all(diff(out$logFC) <= 0))
})
