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


test_that("archived compareCluster resolves enricher without attaching clusterProfiler", {
  skip_if_not_installed("clusterProfiler")

  ids <- paste0("P", seq_len(80))
  p <- rep(0.8, 80)
  p[1:12] <- 1e-08
  p[21:32] <- 1e-08
  logfc <- rep(0, 80)
  logfc[1:12] <- 2
  logfc[21:32] <- -2

  dep <- list(
    B73_Root_VE_vs_Y12_Root_VE = data.frame(
      ID = ids,
      logFC = logfc,
      P.Value = p,
      adj.P.Val = p.adjust(p, method = "BH"),
      regulation = "Not significant",
      Group1 = "B73_Root_VE",
      Group2 = "Y12_Root_VE",
      analysis_mode = "archived",
      protein_universe = "archived_any_detected",
      matrix_source = "Step6_data_normalization",
      test_method = "archived_eBayes",
      stringsAsFactors = FALSE
    )
  )

  background <- rbind(
    data.frame(
      TERM = "T_UP",
      GENE = paste0("P", 1:18),
      NAME = "Up pathway"
    ),
    data.frame(
      TERM = "T_DOWN",
      GENE = paste0("P", 21:38),
      NAME = "Down pathway"
    ),
    data.frame(
      TERM = "T_OTHER",
      GENE = paste0("P", 39:80),
      NAME = "Other pathway"
    )
  )

  # Reproduction must work through namespace imports alone; users should not
  # need library(clusterProfiler) in the Shiny session.
  result <- ProtVis:::.protvis_directional_kegg_data(
    dep_results = dep,
    kegg_background = background,
    comparisons = names(dep),
    method = "archived_comparecluster",
    evidence_mode = "quantitative",
    top_n = 10,
    pvalue_cutoff = 0.05
  )

  expect_gt(nrow(result), 0)
  expect_identical(
    attr(result, "analysis_method"),
    "archived_comparecluster"
  )
})


test_that("bundled Figure 3 DEP membership reproduces frozen direction counts", {
  lists <- ProtVis:::.protvis_directional_load_figure3_gene_lists()
  expect_equal(
    lists$counts$B73_higher,
    c(2287L, 2167L, 2442L, 2291L, 2283L)
  )
  expect_equal(
    lists$counts$Y12_higher,
    c(1315L, 1477L, 1350L, 1425L, 1419L)
  )
  expect_identical(
    names(lists$group1),
    c(
      "Root_VE", "Root_V1.V2", "Root_V4",
      "Leaf_VE.V1.V2", "Leaf_V4.V6.V8"
    )
  )
  expect_true(all(grepl(
    "^(Zm00001d|PZ00001a)",
    unlist(c(lists$group1, lists$group2))
  )))
})


test_that("exact Figure 3 archive reproduces the historical KEGG pathway sets", {
  skip_if_not_installed("clusterProfiler")

  background <- ProtVis:::.protvis_load_builtin_enrichment_background()
  result <- ProtVis:::.protvis_directional_kegg_data(
    dep_results = NULL,
    kegg_background = background$KEGG_background,
    comparisons =
      ProtVis:::.protvis_directional_figure3_spec()$Comparison,
    method = "figure3_archive",
    evidence_mode = "quantitative",
    top_n = 10,
    pvalue_cutoff = 0.05
  )

  expect_gt(nrow(result), 0)
  expect_identical(
    attr(result, "analysis_method"),
    "figure3_archive"
  )

  up <- sort(unique(trimws(
    result$Description[result$Direction == "Upregulated"]
  )))
  down <- sort(unique(trimws(
    result$Description[result$Direction == "Downregulated"]
  )))

  expect_equal(
    up,
    sort(c(
      "Pentose and glucuronate interconversions",
      "Phenylpropanoid biosynthesis",
      "Galactose metabolism",
      "Lipid biosynthesis proteins",
      "Fatty acid biosynthesis"
    ))
  )
  expect_equal(
    down,
    sort(c(
      "Phenylpropanoid biosynthesis",
      "Metabolism of terpenoids and polyketides",
      "Monoterpenoid biosynthesis",
      "Glutathione metabolism",
      "Photosynthesis"
    ))
  )
})


test_that("directional color controls default to blue without changing the statistic", {
  expect_identical(
    ProtVis:::.protvis_directional_plot_colour(
      "#2C7FB8", "#000000"
    ),
    "#2C7FB8"
  )
  expect_identical(
    ProtVis:::.protvis_directional_plot_colour(
      "not-a-colour", "#F7FBFF"
    ),
    "#F7FBFF"
  )

  source <- paste(
    deparse(ProtVis:::.protvis_directional_kegg_panel),
    collapse = "\n"
  )
  expect_match(source, 'color = "pvalue"', fixed = TRUE)
  expect_match(source, 'name = "pvalue"', fixed = TRUE)
  expect_match(source, 'fill = p.adjust', fixed = TRUE)
  expect_match(source, 'name = "BH-adjusted P"', fixed = TRUE)
})
