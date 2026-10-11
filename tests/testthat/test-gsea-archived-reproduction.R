test_that("bundled Root_VE GSEA matches the moderated-t KEGG background", {
  archive_path <- ProtVis:::.protvis_gsea_archive_path()
  expect_true(file.exists(archive_path))

  prepared <- ProtVis:::.protvis_gsea_prepare_archived()

  expect_equal(nrow(prepared$rank_tbl), 6905L)
  expect_equal(sum(prepared$rank_tbl$in_pathway), 196L)
  expect_equal(prepared$extreme_position, 6202L)
  expect_equal(
    prepared$extreme_es,
    -0.3984104,
    tolerance = 1e-6
  )
  expect_equal(
    max(prepared$rank_tbl$running_ES),
    0.2562033,
    tolerance = 1e-6
  )
  expect_equal(
    max(prepared$rank_tbl$rank_metric),
    154.317739,
    tolerance = 1e-6
  )
  expect_equal(
    min(prepared$rank_tbl$rank_metric),
    -175.843961,
    tolerance = 1e-6
  )
})


test_that("bundled Root_VE GSEA uses the moderated-t KEGG ranking", {
  archive <- ProtVis:::.protvis_gsea_load_archived_rootve()

  expect_identical(archive$stage, "Root_VE")
  expect_identical(
    archive$pathway,
    "Phenylpropanoid biosynthesis"
  )
  expect_identical(archive$rank_metric_name, "limma moderated t statistic")
  expect_equal(length(archive$rank_metric), 6905L)
  expect_true(all(diff(archive$rank_metric) <= 1e-12))
})
