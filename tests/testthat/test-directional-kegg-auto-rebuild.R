testthat::test_that("directional KEGG can auto-rebuild archived DEP", {
  testthat::skip_if_not_installed("limma")

  samples <- c(
    "B73_1", "B73_2", "B73_3",
    "Y12_1", "Y12_2", "Y12_3"
  )
  sample_info <- data.frame(
    sample_id = samples,
    group = c(
      rep("B73_Root_VE", 3),
      rep("Y12_Root_VE", 3)
    ),
    stringsAsFactors = FALSE
  )

  set.seed(1)
  step6 <- matrix(
    rnorm(60, 5, 0.2),
    nrow = 10,
    dimnames = list(
      paste0("P", 1:10),
      samples
    )
  )
  step6[1:3, 1:3] <- step6[1:3, 1:3] + 3
  step6[4:6, 4:6] <- step6[4:6, 4:6] + 3

  step4 <- step6
  step4[10, ] <- NA_real_

  current_dep <- list(
    B73_Root_VE_vs_Y12_Root_VE = data.frame(
      ID = paste0("P", 1:10),
      Group1 = "B73_Root_VE",
      Group2 = "Y12_Root_VE",
      analysis_mode = "recommended",
      protein_universe = "both_genotypes",
      matrix_source = "Step4_data_transformed",
      test_method = "ebayes_robust",
      stringsAsFactors = FALSE
    )
  )

  rebuilt <- ProtVis:::.protvis_directional_rebuild_archived_dep(
    dep_results = current_dep,
    comparisons = names(current_dep),
    normalized_matrix = step6,
    detection_matrix = step4,
    sample_info = sample_info
  )

  testthat::expect_named(rebuilt, names(current_dep))
  result <- rebuilt[[1]]
  testthat::expect_true(all(c(
    "ID", "logFC", "P.Value", "adj.P.Val",
    "regulation", "analysis_mode",
    "protein_universe", "matrix_source",
    "test_method", "reproduction_source"
  ) %in% names(result)))
  testthat::expect_true(all(result$analysis_mode == "archived"))
  testthat::expect_true(all(
    result$protein_universe == "archived_any_detected"
  ))
  testthat::expect_true(all(
    result$matrix_source == "Step6_data_normalization"
  ))
  testthat::expect_true(all(
    result$test_method == "archived_eBayes"
  ))
  testthat::expect_false("P10" %in% result$ID)

  status <- ProtVis:::.protvis_directional_archived_compatibility(
    rebuilt,
    names(rebuilt)
  )
  testthat::expect_true(status$ok)
})


testthat::test_that("auto-rebuild does not overwrite current DEP", {
  testthat::skip_if_not_installed("limma")

  samples <- c("A1", "A2", "A3", "B1", "B2", "B3")
  info <- data.frame(
    sample_id = samples,
    group = c(rep("A", 3), rep("B", 3)),
    stringsAsFactors = FALSE
  )
  step6 <- matrix(
    seq_len(48),
    nrow = 8,
    dimnames = list(paste0("P", 1:8), samples)
  )
  step4 <- step6
  current <- list(
    A_vs_B = data.frame(
      ID = paste0("P", 1:8),
      Group1 = "A",
      Group2 = "B",
      analysis_mode = "recommended",
      protein_universe = "both_genotypes",
      matrix_source = "Step4_data_transformed",
      test_method = "ebayes_robust",
      stringsAsFactors = FALSE
    )
  )
  before <- current

  invisible(
    ProtVis:::.protvis_directional_rebuild_archived_dep(
      current,
      "A_vs_B",
      step6,
      step4,
      info
    )
  )

  testthat::expect_identical(current, before)
})
