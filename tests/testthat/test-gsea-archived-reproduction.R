test_that("archived Root_VE GSEA reproduction matches the supplied panel geometry", {
  archive_path <- ProtVis:::.protvis_gsea_archive_path()
  expect_true(file.exists(archive_path))

  prepared <- ProtVis:::.protvis_gsea_prepare_archived()

  expect_equal(nrow(prepared$rank_tbl), 11049L)
  expect_equal(sum(prepared$rank_tbl$in_pathway), 196L)
  expect_equal(prepared$extreme_position, 9893L)
  expect_equal(
    prepared$extreme_es,
    -0.3064468,
    tolerance = 1e-6
  )
  expect_equal(
    max(prepared$rank_tbl$running_ES),
    0.2769241,
    tolerance = 1e-6
  )
  expect_equal(
    max(prepared$rank_tbl$rank_metric),
    11.324951,
    tolerance = 1e-6
  )
  expect_equal(
    min(prepared$rank_tbl$rank_metric),
    -7.920156,
    tolerance = 1e-6
  )
})


test_that("archived Root_VE GSEA uses logFC full retained ranking", {
  archive <- ProtVis:::.protvis_gsea_load_archived_rootve()

  expect_identical(archive$stage, "Root_VE")
  expect_identical(
    archive$pathway,
    "Phenylpropanoid biosynthesis"
  )
  expect_identical(archive$rank_metric_name, "logFC")
  expect_equal(length(archive$rank_metric), 11049L)
  expect_true(all(diff(archive$rank_metric) <= 1e-12))
})
