testthat::test_that("Recommended DEP defaults avoid historical Step6 artifacts", {
  defaults <- ProtVis:::.protvis_dep_recommended_defaults()

  testthat::expect_identical(defaults$mode, "recommended")
  testthat::expect_identical(defaults$matrix_source, "Step4_data_transformed")
  testthat::expect_false(defaults$matrix_shift)
  testthat::expect_true(defaults$center_samples)
  testthat::expect_identical(defaults$protein_universe, "both_genotypes")
  testthat::expect_equal(defaults$min_detected, 2L)
  testthat::expect_identical(defaults$test_method, "ebayes_robust")
  testthat::expect_identical(defaults$p_metric, "adj.P.Val")
  testthat::expect_identical(defaults$volcano_p_metric, "P.Value")
})

testthat::test_that("Recommended DEP median-centres observed values and preserves missingness", {
  mat <- matrix(
    c(
      10, 20, 30,
      12, NA, 32,
      14, 24, 34
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      paste0("P", 1:3),
      c("S1", "S2", "S3")
    )
  )

  centered <- ProtVis:::.protvis_dep_prepare_recommended_matrix(
    mat,
    center_samples = TRUE
  )

  testthat::expect_equal(
    apply(centered, 2, stats::median, na.rm = TRUE),
    c(S1 = 0, S2 = 0, S3 = 0)
  )
  testthat::expect_true(is.na(centered["P2", "S2"]))
  testthat::expect_true(any(centered < 0, na.rm = TRUE))
  testthat::expect_false(all(centered[is.finite(centered)] > 0))
})

testthat::test_that("Recommended detection filter requires real observations in both groups", {
  mat <- matrix(
    c(
      1, 1, NA, 2, 2, NA,
      1, 1, NA, NA, NA, NA,
      NA, NA, NA, 2, 2, 2,
      1, NA, NA, 2, NA, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("P", 1:4),
      c("A1", "A2", "A3", "B1", "B2", "B3")
    )
  )

  both <- ProtVis:::.protvis_dep_shared_ids(
    mat,
    c("A1", "A2", "A3"),
    c("B1", "B2", "B3"),
    mode = "both_genotypes",
    min_detected = 2L
  )
  either <- ProtVis:::.protvis_dep_shared_ids(
    mat,
    c("A1", "A2", "A3"),
    c("B1", "B2", "B3"),
    mode = "either_genotype",
    min_detected = 2L
  )

  testthat::expect_identical(both, "P1")
  testthat::expect_identical(either, c("P1", "P2", "P3"))
})

testthat::test_that("Presence-absence proteins are reported separately instead of imputed into limma", {
  mat <- matrix(
    c(
      1, 1, 1, NA, NA, NA,
      NA, NA, NA, 2, 2, 2,
      1, 1, NA, 2, 2, NA,
      1, NA, NA, 2, NA, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("P", 1:4),
      c("A1", "A2", "A3", "B1", "B2", "B3")
    )
  )

  out <- ProtVis:::.protvis_dep_presence_absence(
    mat,
    c("A1", "A2", "A3"),
    c("B1", "B2", "B3"),
    min_detected = 2L
  )

  testthat::expect_identical(out$ID, c("P1", "P2"))
  testthat::expect_identical(
    out$Pattern,
    c("Detected in Group1 only", "Detected in Group2 only")
  )
})

testthat::test_that("Recommended limma uses robust trend eBayes without positive shifting", {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("statmod")

  set.seed(42)
  mat <- matrix(
    stats::rnorm(100 * 6, mean = 6, sd = 0.35),
    nrow = 100,
    dimnames = list(
      paste0("P", seq_len(100)),
      c("A1", "A2", "A3", "B1", "B2", "B3")
    )
  )
  mat[1:10, 1:3] <- mat[1:10, 1:3] + 2

  out <- ProtVis:::.protvis_dep_run_limma_recommended(
    mat,
    c("A1", "A2", "A3"),
    c("B1", "B2", "B3"),
    "A",
    "B",
    adjust_method = "BH",
    sort_by = "P",
    test_method = "ebayes_robust",
    lfc = 1
  )

  testthat::expect_true(all(c(
    "ID", "logFC", "P.Value", "adj.P.Val", "FC",
    "Group1", "Group2", "fc_threshold_tested"
  ) %in% names(out)))
  testthat::expect_false(any(out$fc_threshold_tested))
  testthat::expect_true(all(diff(out$P.Value) >= 0))
  testthat::expect_gt(stats::median(out$logFC[out$ID %in% paste0("P", 1:10)]), 1)
})


testthat::test_that("parallel evidence summary separates quantitative and presence-absence evidence", {
  results <- list(
    A_vs_B = data.frame(
      ID = c("P1", "P2", "P3", "P4"),
      regulation = c(
        "Upregulated",
        "Downregulated",
        "Not significant",
        "Upregulated"
      ),
      stringsAsFactors = FALSE
    )
  )
  presence <- list(
    A_vs_B = data.frame(
      ID = c("P5", "P6", "P7"),
      Pattern = c(
        "Detected in Group1 only",
        "Detected in Group2 only",
        "Detected in Group1 only"
      ),
      stringsAsFactors = FALSE
    )
  )

  out <- ProtVis:::.protvis_dep_evidence_summary(
    results,
    presence,
    "A_vs_B"
  )

  testthat::expect_equal(
    out$Protein_number[
      out$Evidence == "Quantitative DEP" &
        out$Direction == "Upregulated"
    ],
    2L
  )
  testthat::expect_equal(
    out$Protein_number[
      out$Evidence == "Quantitative DEP" &
        out$Direction == "Downregulated"
    ],
    1L
  )
  testthat::expect_equal(
    out$Protein_number[
      out$Evidence == "Presence/absence" &
        out$Direction == "Detected in Group1 only"
    ],
    2L
  )
  testthat::expect_equal(
    out$Protein_number[
      out$Evidence == "Presence/absence" &
        out$Direction == "Detected in Group2 only"
    ],
    1L
  )
  total <- out[
    out$Evidence == "All retained differential evidence",
    , drop = FALSE
  ]
  testthat::expect_equal(
    total$Protein_number[
      total$Direction == "Quantitative DEP tested"
    ],
    4L
  )
  testthat::expect_equal(total$Total_proteins, rep(7L, 3L))
})

testthat::test_that("presence-absence evidence retains observed-intensity context", {
  mat <- matrix(
    c(
      5, 6, 7, NA, NA, NA,
      NA, NA, NA, 8, 9, 10
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("P1", "P2"),
      c("A1", "A2", "A3", "B1", "B2", "B3")
    )
  )

  out <- ProtVis:::.protvis_dep_presence_absence(
    mat,
    c("A1", "A2", "A3"),
    c("B1", "B2", "B3"),
    min_detected = 2L
  )

  testthat::expect_true(all(c(
    "Group1_median_observed",
    "Group2_median_observed",
    "Detection_difference"
  ) %in% names(out)))
  testthat::expect_equal(
    out$Group1_median_observed[out$ID == "P1"],
    6
  )
  testthat::expect_true(
    is.na(out$Group2_median_observed[out$ID == "P1"])
  )
  testthat::expect_equal(
    out$Detection_difference[out$ID == "P1"],
    3
  )
})

testthat::test_that("Evidence 2 count overview uses the candidate table", {
  presence <- data.frame(
    ID = paste0("P", 1:5),
    Pattern = c(
      "Detected in Group1 only",
      "Detected in Group2 only",
      "Detected in Group1 only",
      "Detected in Group1 only",
      "Detected in Group2 only"
    ),
    stringsAsFactors = FALSE
  )

  counts <- ProtVis:::.protvis_dep_presence_count_data(
    presence, "B73_Root_VE", "Y12_Root_VE"
  )

  testthat::expect_identical(
    as.character(counts$Direction),
    c(
      "B73_Root_VE detected / Y12_Root_VE not detected",
      "Y12_Root_VE detected / B73_Root_VE not detected"
    )
  )
  testthat::expect_equal(counts$Protein_number, c(3L, 2L))

  empty <- ProtVis:::.protvis_dep_presence_count_data(
    data.frame(), "B73", "Y12"
  )
  testthat::expect_equal(empty$Protein_number, c(0L, 0L))
})
