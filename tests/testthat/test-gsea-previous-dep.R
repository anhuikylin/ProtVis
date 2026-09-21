testthat::test_that("previous DEP GSEA defaults to moderated t ranking", {
  dep <- data.frame(
    ID = paste0("P", 1:8),
    logFC = c(1, -1, 0.5, -0.5, 2, -2, 0.2, -0.2),
    t = c(8, -7, 6, -5, 4, -3, 2, -1),
    P.Value = c(1e-5, 1e-4, 0.001, 0.002, 0.01, 0.02, 0.2, 0.3),
    adj.P.Val = c(1e-4, 1e-3, 0.01, 0.02, 0.05, 0.08, 0.4, 0.5),
    Group1 = "B73_Root_VE",
    Group2 = "Y12_Root_VE",
    stringsAsFactors = FALSE
  )

  ranked <- ProtVis:::.protvis_gsea_prepare_dep_rank(
    dep,
    rank_metric = "t",
    filter_mode = "all"
  )

  testthat::expect_identical(
    names(ranked$rank_vector),
    paste0("P", c(1, 3, 5, 7, 8, 6, 4, 2))
  )
  testthat::expect_equal(
    unname(ranked$rank_vector),
    c(8, 6, 4, 2, -1, -3, -5, -7)
  )
  testthat::expect_equal(ranked$n_after_filter, 8L)
})


testthat::test_that("custom previous DEP filters are optional and explicit", {
  dep <- data.frame(
    ID = paste0("P", 1:5),
    logFC = c(2, 1.5, 0.2, -1.8, -0.1),
    t = c(5, 4, 3, -6, -2),
    P.Value = c(0.001, 0.02, 0.5, 0.003, 0.7),
    adj.P.Val = c(0.01, 0.04, 0.8, 0.02, 0.9),
    stringsAsFactors = FALSE
  )

  all_ranked <- ProtVis:::.protvis_gsea_prepare_dep_rank(
    dep,
    rank_metric = "t",
    filter_mode = "all"
  )
  filtered <- ProtVis:::.protvis_gsea_prepare_dep_rank(
    dep,
    rank_metric = "t",
    filter_mode = "custom",
    p_metric = "adj.P.Val",
    p_cutoff = 0.05,
    abs_logfc = 1
  )

  testthat::expect_equal(length(all_ranked$rank_vector), 5L)
  testthat::expect_identical(
    sort(names(filtered$rank_vector)),
    c("P1", "P2", "P4")
  )
})


testthat::test_that("Root_VE is the default focus among multiple comparisons", {
  comparisons <- c(
    "B73_Root_V4_vs_Y12_Root_V4",
    "B73_Leaf_VE.V1.V2_vs_Y12_Leaf_VE.V1.V2",
    "B73_Root_VE_vs_Y12_Root_VE",
    "B73_Root_V1.V2_vs_Y12_Root_V1.V2"
  )
  testthat::expect_identical(
    ProtVis:::.protvis_gsea_rootve_comparison(comparisons),
    "B73_Root_VE_vs_Y12_Root_VE"
  )
})


testthat::test_that("archived DEP provenance is detected for exact reproduction", {
  dep <- data.frame(
    ID = "P1",
    logFC = 1,
    t = 2,
    P.Value = 0.01,
    adj.P.Val = 0.02,
    analysis_mode = "archived",
    protein_universe = "archived_any_detected",
    matrix_source = "Step6_data_normalization",
    test_method = "archived_eBayes",
    stringsAsFactors = FALSE
  )
  status <- ProtVis:::.protvis_gsea_dep_compatibility(dep)
  testthat::expect_true(status$exact)

  dep$analysis_mode <- "recommended"
  status2 <- ProtVis:::.protvis_gsea_dep_compatibility(dep)
  testthat::expect_false(status2$exact)
})
