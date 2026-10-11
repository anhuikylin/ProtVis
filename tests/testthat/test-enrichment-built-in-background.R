test_that("built-in maize-teosinte background uses Background Make format", {
  background <- ProtVis:::.protvis_load_builtin_enrichment_background()

  expect_named(background, c("GO_background", "KEGG_background"))
  for (table in background) {
    expect_true(is.data.frame(table))
    expect_true(all(c("TERM", "GENE", "NAME") %in% names(table)))
    expect_gt(nrow(table), 0)
    expect_false(anyNA(table[, c("TERM", "GENE", "NAME")]))
  }
})
