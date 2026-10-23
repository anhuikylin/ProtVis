test_that("subcellular localization normalizes common provider labels", {
  expect_equal(.subcellular_normalize_locations("plasma membrane"), "Cell membrane")
  expect_equal(.subcellular_normalize_locations("Endoplasmic reticulum"), "Endoplasmic reticulum")
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


test_that("DeepLoc full screenshot example retains all scores and metadata", {
  x <- .subcellular_import_deeploc(
    "Endoplasmic reticulum", "Soluble", "Signal peptide",
    .subcellular_deeploc_example()
  )
  expect_equal(x$prediction, "Endoplasmic reticulum")
  expect_equal(x$membrane_types, "Soluble")
  expect_equal(x$signals, "Signal peptide")
  expect_equal(nrow(x$score_table), 14L)
  expect_equal(x$score_table$probability[x$score_table$label == "Plastid"], 0.0019)
  plastid <- .subcellular_import_deeploc("Plastid")
  expect_equal(plastid$deeploc_locations, "Plastid")
  expect_equal(plastid$prediction, "Chloroplast")
  expect_equal(.subcellular_deeploc_probabilities(
    sub("Plastid,0.0019", "Chloroplast,0.0019", .subcellular_deeploc_example(),
        fixed = TRUE)
  )$label[[6L]], "Plastid")
  expect_equal(x$score_table$probability[x$score_table$label == "Endoplasmic reticulum"], 0.7058)
  expect_error(.subcellular_deeploc_probabilities("Soluble,1.2"), "Probability")
})

test_that("DeepLoc summary matches protein ID and importance validates residue positions", {
  summary_file <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame("Protein ID" = c("A", "B"),
    "Predicted localizations" = c("Nucleus", "Cytoplasm"),
    "Predicted membrane association" = c("Soluble", "Peripheral"),
    check.names = FALSE), summary_file, row.names = FALSE)
  summary <- .subcellular_deeploc_summary(list(datapath = summary_file), "B")
  expect_equal(summary$location, "Cytoplasm")
  expect_error(.subcellular_deeploc_summary(list(datapath = summary_file), "C"), "matching")

  importance_file <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(position = c(1, 3), residue = c("A", "D"),
                              importance = c(0.7, 0.1)), importance_file, row.names = FALSE)
  importance <- .subcellular_deeploc_importance(list(datapath = importance_file))
  expect_equal(nrow(importance), 2L)
  expect_equal(.subcellular_import_deeploc("Nucleus", importance = importance,
                                          sequence = "ACDEFGHIKLMN")$sorting_importance,
               importance)
  expect_error(.subcellular_import_deeploc("Nucleus", importance = importance,
                                          sequence = "ACCEFGHIKLMN"), "do not match")
})

test_that("FASTA upload selects one protein and rejects multiple records", {
  f <- tempfile(fileext = ".fasta")
  writeLines(c(">sp|P07221|CASQ1_RABIT example", "MNAADRMGARVALLLLLVLGSPQSGVHGEEGLD"), f)
  x <- .subcellular_deeploc_read_fasta(list(datapath = f, name = "example.fasta"))
  expect_equal(x$id, "sp_P07221_CASQ1_RABIT")
  expect_true(startsWith(x$sequence, "MNAADR"))
  writeLines(c(">one", "ACDEFGHIKLMN", ">two", "ACDEFGHIKLMN"), f)
  expect_error(.subcellular_deeploc_read_fasta(list(datapath = f, name = "two.fasta")),
               "one protein")
})

test_that("official DeepLoc example summary has full predictions for each sequence", {
  path <- testthat::test_path("..", "..", "inst", "extdata",
                              "deeploc_2_1_example_summary.csv")
  x <- .subcellular_deeploc_summary(
    list(datapath = path), "sp_P07221_CASQ1_RABIT"
  )
  expect_equal(x$location, "Extracellular|Endoplasmic reticulum")
  imported <- .subcellular_import_deeploc(x$location, x$membrane, x$signals,
                                           x$probabilities)
  expect_equal(imported$prediction, c("Extracellular", "Endoplasmic reticulum"))
  expect_equal(imported$membrane_types, "Soluble")
  expect_equal(imported$signals, "Signal peptide")
  expect_equal(nrow(imported$score_table), 14L)
  expect_equal(imported$score_table$probability[
    imported$score_table$label == "Extracellular"
  ], 0.707099974155426)
  expect_error(.subcellular_deeploc_summary(list(datapath = path), "wrong_id"),
               "matching")
  expect_equal(.subcellular_deeploc_job_id(
    "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6ABB678C002764C16CB5A81F&wait=20"
  ), "6ABB678C002764C16CB5A81F")
})

test_that("official DeepLoc attention CSV infers residue positions", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(AA = c("M", "N", "A"),
                              Alpha = c(0.00448785, 0.004924497, 0.0057411827)),
                   path, row.names = FALSE)
  d <- .subcellular_deeploc_importance(list(datapath = path))
  expect_equal(d$position, 1:3)
  expect_equal(d$residue, c("M", "N", "A"))
  expect_equal(d$importance[[1L]], 0.00448785)
  expect_equal(.subcellular_import_deeploc("Endoplasmic reticulum",
                                            importance = d,
                                            sequence = "MNAADRMGARVALLLLLVLG")$sorting_importance, d)
})

test_that("DeepLoc failed job is not mistaken for a completed prediction", {
  failed_page <- paste0(
    "<html><h1>Failed run</h1><p>Unfortunately, your job failed.</p>",
    "<div>Download prediction results: CSV Summary</div></html>"
  )
  expect_equal(.subcellular_deeploc_page_state(failed_page), "failed")
  expect_equal(.subcellular_deeploc_page_state(
    "<html><h1>Queued run</h1><p>Please wait</p></html>"
  ), "queued")
  expect_equal(.subcellular_deeploc_page_state(
    "<html><h1>Finished run</h1><a>Download prediction results</a></html>"
  ), "complete")
})
