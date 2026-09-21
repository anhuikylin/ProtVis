test_that("archived directional KEGG reproduces compareCluster semantics", {
  skip_if_not_installed("clusterProfiler")

  ids <- paste0("P", seq_len(100))
  make_result <- function(stage) {
    p <- rep(0.8, 100)
    p[1:10] <- 1e-08
    p[21:30] <- 1e-08
    logfc <- rep(0, 100)
    logfc[1:10] <- 2
    logfc[21:30] <- -2

    data.frame(
      ID = ids,
      logFC = logfc,
      P.Value = p,
      adj.P.Val = p.adjust(p, method = "BH"),
      regulation = "Not significant",
      Group1 = paste0("B73_", stage),
      Group2 = paste0("Y12_", stage),
      analysis_mode = "archived",
      protein_universe = "archived_any_detected",
      matrix_source = "Step6_data_normalization",
      test_method = "archived_eBayes",
      stringsAsFactors = FALSE
    )
  }

  dep <- list(
    B73_Root_VE_vs_Y12_Root_VE =
      make_result("Root_VE"),
    B73_Root_V4_vs_Y12_Root_V4 =
      make_result("Root_V4")
  )

  background <- rbind(
    data.frame(
      TERM = "T_UP",
      GENE = paste0("P", 1:20),
      NAME = "Up pathway"
    ),
    data.frame(
      TERM = "T_DOWN",
      GENE = paste0("P", 21:40),
      NAME = "Down pathway"
    ),
    data.frame(
      TERM = "T_OTHER",
      GENE = paste0("P", 41:100),
      NAME = "Other pathway"
    )
  )

  result <- ProtVis:::.protvis_directional_kegg_data(
    dep,
    background,
    comparisons = names(dep),
    top_n = 10,
    method = "archived_comparecluster",
    evidence_mode = "quantitative",
    pvalue_cutoff = 0.05
  )

  expect_gt(nrow(result), 0)
  expect_identical(
    attr(result, "analysis_method"),
    "archived_comparecluster"
  )
  expect_true(
    all(c("Upregulated", "Downregulated") %in%
          unique(result$Direction))
  )
  expect_true(any(
    result$Direction == "Upregulated" &
      result$Description == "Up pathway"
  ))
  expect_true(any(
    result$Direction == "Downregulated" &
      result$Description == "Down pathway"
  ))
  expect_true(all(
    result$Analysis_mode ==
      "Figure 3 compareCluster reproduction"
  ))

  objects <- attr(result, "comparecluster_objects")
  expect_true(is.list(objects))
  expect_s4_class(objects$Upregulated, "compareClusterResult")
  expect_s4_class(objects$Downregulated, "compareClusterResult")

  # The archived engine must rebuild the directional list from P.Value/logFC;
  # it must not trust the current regulation labels.
  expect_true(all(dep[[1]]$regulation == "Not significant"))
  expect_true(nrow(result) > 0)
})


test_that("archived reproduction rejects non-archived DEP provenance", {
  dep <- list(
    demo = data.frame(
      ID = paste0("P", 1:20),
      logFC = rep(2, 20),
      P.Value = rep(1e-06, 20),
      Group1 = "B73_Root_VE",
      Group2 = "Y12_Root_VE",
      analysis_mode = "recommended",
      protein_universe = "both_genotypes",
      matrix_source = "Step4_data_transformed",
      test_method = "ebayes_robust",
      stringsAsFactors = FALSE
    )
  )

  status <- ProtVis:::.protvis_directional_archived_compatibility(
    dep,
    "demo"
  )
  expect_false(status$ok)
  expect_match(
    status$message,
    "Figure 3 reproduction requires DEP > Archived reproduction",
    fixed = TRUE
  )
})


test_that("archived background keeps first TERM per pathway NAME", {
  background <- data.frame(
    TERM = c("T1", "T2", "T1", "T3"),
    GENE = c("P1", "P2", "P3", "P4"),
    NAME = c("Same name", "Same name", "Same name", "Other"),
    stringsAsFactors = FALSE
  )

  archived <- ProtVis:::.protvis_directional_background(
    background,
    archived = TRUE
  )

  expect_equal(
    archived$TERM2NAME$TERM,
    c("T1", "T3")
  )
  expect_true(all(
    archived$TERM2GENE$TERM %in% c("T1", "T3")
  ))
  expect_false("T2" %in% archived$TERM2GENE$TERM)
})
