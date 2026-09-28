test_that("subcellular localization normalizes common provider labels", {
  expect_equal(.subcellular_normalize_locations("plasma membrane"), "Cell membrane")
  expect_equal(.subcellular_normalize_locations("chloroplastic protein"), "Chloroplast")
  expect_equal(.subcellular_normalize_locations(c("nuclear", "cytoplasmic")),
               c("Nucleus", "Cytoplasm"))
})

test_that("WoLF PSORT parser selects the highest scoring location", {
  x <- .subcellular_parse_wolfpsort(
    "ZmProtein_demo chlo: 14, cyto: 8, nucl: 3, mito: 2, plas: 1"
  )
  expect_equal(x$prediction, "Chloroplast")
  expect_true(is.data.frame(x$scores))
  expect_equal(x$scores$location[[1]], "Chloroplast")
})

test_that("consensus tolerates unavailable providers", {
  providers <- list(
    .subcellular_provider_result("Plant-mPLoc", "unavailable",
                                 error = "server unavailable"),
    .subcellular_provider_result("WoLF PSORT", "success", "Chloroplast"),
    .subcellular_provider_result("BUSCA", "success", "Chloroplast")
  )
  x <- .subcellular_consensus(providers)
  expect_equal(x$prediction, "Chloroplast")
  expect_equal(x$confidence, "Moderate")
  expect_equal(x$successful, 2L)
  expect_equal(x$requested, 3L)
})

test_that("consensus reports disagreement without hiding it", {
  providers <- list(
    .subcellular_provider_result("Plant-mPLoc", "success", "Nucleus"),
    .subcellular_provider_result("WoLF PSORT", "success", "Chloroplast"),
    .subcellular_provider_result("BUSCA", "unavailable", error = "timeout")
  )
  x <- .subcellular_consensus(providers)
  expect_setequal(x$prediction, c("Nucleus", "Chloroplast"))
  expect_equal(x$confidence, "Low")
})

test_that("application uses the new localization module", {
  ui <- paste(readLines(testthat::test_path("..", "..", "R", "app_ui.R")), collapse = "\n")
  server <- paste(readLines(testthat::test_path("..", "..", "R", "app_server.R")), collapse = "\n")
  module <- paste(readLines(testthat::test_path("..", "..", "R", "subcellular_localization.R")), collapse = "\n")

  expect_match(ui, 'nav_panel\\("Subcellular localization"')
  expect_match(ui, 'subcellular_localization_ui\\("subcellular_localization"\\)')
  expect_match(server, 'subcellular_localization_server\\("subcellular_localization"\\)')
  expect_match(module, 'checkboxGroupInput')
  expect_match(module, 'Plant-mPLoc')
  expect_match(module, 'WoLF PSORT')
  expect_match(module, 'BUSCA')
  expect_match(module, 'Completed with warnings')
})
