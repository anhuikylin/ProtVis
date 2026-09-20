testthat::test_that("built-in maize-teosinte KEGG Figure 3C-D data pass integrity checks", {
  df <- ProtVis:::.protvis_maize_teosinte_kegg_data()

  testthat::expect_invisible(
    ProtVis:::.protvis_validate_maize_teosinte_kegg_data(df)
  )
  testthat::expect_equal(nrow(df), 29L)
  testthat::expect_equal(sum(df$Panel == "C"), 18L)
  testthat::expect_equal(sum(df$Panel == "D"), 11L)
  testthat::expect_false(any(
    df$Panel == "D" & df$Cluster == "Root_V1.V2"
  ))
})

testthat::test_that("maize-teosinte KEGG reproduction keeps key archived counts", {
  df <- ProtVis:::.protvis_maize_teosinte_kegg_data()

  pick <- function(panel, cluster, description) {
    df[
      df$Panel == panel &
        df$Cluster == cluster &
        df$Description == description,
      ,
      drop = FALSE
    ]
  }

  root_v12 <- pick(
    "C",
    "Root_V1.V2",
    "Phenylpropanoid biosynthesis"
  )
  testthat::expect_equal(nrow(root_v12), 1L)
  testthat::expect_equal(root_v12$Count, 45)
  testthat::expect_equal(
    root_v12$pvalue,
    8.37198056160672e-08,
    tolerance = 1e-15
  )

  leaf_v4v8 <- pick(
    "D",
    "Leaf_V4.V6.V8",
    "Photosynthesis"
  )
  testthat::expect_equal(nrow(leaf_v4v8), 1L)
  testthat::expect_equal(leaf_v4v8$Count, 12)
  testthat::expect_equal(
    leaf_v4v8$pvalue,
    0.001225228463218804,
    tolerance = 1e-15
  )
})

testthat::test_that("maize-teosinte KEGG panels preserve publication categories", {
  df <- ProtVis:::.protvis_maize_teosinte_kegg_data()
  p_c <- ProtVis:::.protvis_maize_teosinte_kegg_panel(df, "C")
  p_d <- ProtVis:::.protvis_maize_teosinte_kegg_panel(df, "D")
  combined <- ProtVis:::plot_maize_teosinte_kegg_reproduction(df)

  testthat::expect_s3_class(p_c, "ggplot")
  testthat::expect_s3_class(p_d, "ggplot")
  testthat::expect_true(inherits(combined, "patchwork"))

  c_x <- levels(p_c$data$Cluster_label)
  d_x <- levels(p_d$data$Cluster_label)
  testthat::expect_identical(
    c_x,
    c("Root_VE", "Root_V1.V2", "Root_V4", "Leaf_VE-V2", "Leaf_V4-V8")
  )
  testthat::expect_identical(
    d_x,
    c("Root_VE", "Root_V4", "Leaf_VE-V2", "Leaf_V4-V8")
  )
})
