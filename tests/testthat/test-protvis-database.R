test_that("ProtVisDatabase supplies every declared local data resource", {
  manifest <- ProtVisDatabase::protvis_database_manifest()
  expect_gte(nrow(manifest), 20L)
  paths <- vapply(manifest$path, function(path) {
    ProtVisDatabase::protvis_database_path(
      "extdata", path,
      must_work = TRUE
    )
  }, character(1))
  expect_true(all(file.exists(paths)))
  expect_true(all(nzchar(manifest$sha256)))
})
