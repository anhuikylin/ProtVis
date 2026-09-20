testthat::test_that("triplicate grouping keeps three replicates together", {
  samples <- c(
    "B73_Root_VE_1", "B73_Root_VE_2", "B73_Root_VE_3",
    "Y12_Root_VE_1", "Y12_Root_VE_2", "Y12_Root_VE_3"
  )
  groups <- ProtVis:::.protvis_sample_group_values(
    sample_info = NULL,
    sample_ids = samples,
    mode = "triplicate"
  )
  testthat::expect_identical(
    groups,
    c(
      rep("B73_Root_VE", 3),
      rep("Y12_Root_VE", 3)
    )
  )
})

testthat::test_that("species grouping follows encoded B73 and Y12 identities", {
  samples <- c("B73_Root_VE_1", "Y12_Root_VE_1")
  deliberately_wrong <- data.frame(
    sample_id = samples,
    species = c("Root", "Shoot"),
    stringsAsFactors = FALSE
  )
  species <- ProtVis:::.protvis_sample_group_values(
    deliberately_wrong, samples, mode = "species"
  )
  testthat::expect_identical(
    species,
    c("Zea mays ssp. mays", "Zea mays ssp. mexicana")
  )
})

testthat::test_that("tissue grouping follows encoded sample names", {
  samples <- c("B73_Root_VE_1", "B73_Leaf_VE.V1.V2_1")
  deliberately_wrong <- data.frame(
    sample_id = samples,
    tissue = c("Shoot", "Root"),
    stringsAsFactors = FALSE
  )
  tissue <- ProtVis:::.protvis_sample_group_values(
    deliberately_wrong, samples, mode = "tissue"
  )
  testthat::expect_identical(tissue, c("Root", "Shoot"))
})

testthat::test_that("manual shape palette supports more than six groups", {
  groups <- paste0("Group", 1:10)
  shapes <- ProtVis:::.protvis_shape_palette(groups)
  testthat::expect_equal(length(shapes), 10L)
  testthat::expect_equal(length(unique(shapes)), 10L)
  testthat::expect_identical(names(shapes), groups)
})

testthat::test_that("group palette is stable and named", {
  groups <- c("A", "B", "C")
  first <- ProtVis:::.protvis_group_palette(groups)
  second <- ProtVis:::.protvis_group_palette(groups)
  testthat::expect_identical(first, second)
  testthat::expect_identical(names(first), groups)
})
