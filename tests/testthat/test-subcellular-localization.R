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

test_that("effective Toolkits UI uses the same localization module as the server", {
  active_ui <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzz_protein_workbench_app_ui.R")),
    collapse = "\n"
  )
  server <- paste(
    readLines(testthat::test_path("..", "..", "R", "app_server.R")),
    collapse = "\n"
  )
  expect_true(grepl(
    'subcellular_localization_ui("subcellular_localization")',
    active_ui, fixed = TRUE
  ))
  expect_true(grepl(
    'subcellular_localization_server("subcellular_localization"',
    server, fixed = TRUE
  ))
  expect_false(grepl('plant_mploc_ui("plant_mploc")', active_ui, fixed = TRUE))
})

test_that("DeepLoc web result can rescue unavailable automatic sources", {
  provider <- .subcellular_import_deeploc("chloroplast; Lysosome/Vacuole")
  expect_equal(provider$status, "success")
  expect_equal(provider$prediction, c("Chloroplast", "Lysosome/Vacuole"))
  expect_match(provider$details, "user", ignore.case = TRUE)

  providers <- list(
    .subcellular_provider_result("Plant-mPLoc", "unavailable", error = "timeout"),
    "DeepLoc 2.1" = provider
  )
  summary <- .subcellular_consensus(providers)
  expect_equal(summary$successful, 1L)
  expect_setequal(summary$prediction, c("Chloroplast", "Lysosome/Vacuole"))
  expect_error(.subcellular_import_deeploc("unknown compartment"), "DeepLoc")
})

test_that("DeepLoc FASTA export validates sequence length and protein ID", {
  fasta <- .subcellular_deeploc_fasta("ACDEFGHIKLMN", "protein 1")
  expect_true(startsWith(fasta, ">protein_1\n"))
  expect_error(.subcellular_deeploc_fasta("ACD", "short"), "at least 10")
})

test_that("web refresh target and BUSCA result table match live page formats", {
  expect_equal(
    .subcellular_parse_refresh("5; URL=/results/job.html"),
    "/results/job.html"
  )
  expect_equal(.subcellular_parse_refresh("5"), "")
  expect_equal(.subcellular_parse_refresh("0; url=javascript:alert(1)"), "")

  page <- rvest::read_html(paste0(
    '<table id="resultdata"><thead><tr><th></th>',
    '<th>Protein Accession/ID</th><th>GO-id</th><th>GO-term</th>',
    '</tr></thead><tbody><tr><td></td><td>test_protein</td>',
    '<td>GO:0005615</td><td>C:extracellular space</td>',
    '</tr></tbody></table>'
  ))
  expect_equal(.subcellular_parse_busca_page(page), "Extracellular")
  expect_equal(
    .subcellular_parse_busca_json(list(data = list(
      list("test_protein", "GO:0005615", "C:extracellular space", 1)
    ))),
    "Extracellular"
  )
  expect_length(.subcellular_parse_busca_page(rvest::read_html("<p>Queued</p>")), 0)
})

