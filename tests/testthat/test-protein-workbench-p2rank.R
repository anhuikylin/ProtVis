test_that("Protein Workbench binding pockets default to P2Rank", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_p2rank_test"))

  expect_match(html, "P2Rank (recommended)", fixed = TRUE)
  expect_match(html, "fpocket (advanced / optional)", fixed = TRUE)
  expect_match(html, "P2Rank folder or executable", fixed = TRUE)
  expect_match(html, "Download P2Rank", fixed = TRUE)
  expect_match(html, "Java 17+ required", fixed = TRUE)
  expect_match(html, "RUN POCKET PREDICTION", fixed = TRUE)

  for (id in c(
    "pocket_backend",
    "p2rank_path",
    "p2rank_profile",
    "p2rank_threads",
    "p2rank_status",
    "download_pocket_results"
  )) {
    expect_match(html, id, fixed = TRUE)
  }
})

test_that("P2Rank executable detection accepts an extracted folder", {
  root <- tempfile("p2rank_")
  dir.create(root)
  exe <- file.path(root, if (.Platform$OS.type == "windows") "prank.bat" else "prank")
  writeLines("echo test", exe)

  detected <- ProtVis:::.protvis_pw_p2rank_executable(root)
  expect_true(nzchar(detected))
  expect_equal(
    normalizePath(detected, winslash = "/", mustWork = TRUE),
    normalizePath(exe, winslash = "/", mustWork = TRUE)
  )
})

test_that("P2Rank Windows Java classpath includes the bundled runtime libraries", {
  classpath <- ProtVis:::.protvis_pw_p2rank_classpath("P2Rank install")

  expect_identical(
    classpath,
    paste(
      file.path("P2Rank install", "bin", "p2rank.jar"),
      file.path("P2Rank install", "bin", "lib", "*"),
      sep = .Platform$path.sep
    )
  )
})

test_that("P2Rank prediction CSV is normalized for ProtVis", {
  pdb <- tempfile(fileext = ".pdb")
  writeLines(c(
    "ATOM      1  N   GLY A  10      11.104  13.207  10.111  1.00 20.00           N",
    "ATOM      2  CA  GLY A  10      12.000  13.500  10.500  1.00 20.00           C",
    "ATOM      3  N   SER A  11      13.104  14.207  11.111  1.00 20.00           N",
    "ATOM      4  CA  SER A  11      14.000  14.500  11.500  1.00 20.00           C",
    "END"
  ), pdb)

  out <- tempfile("p2out_")
  dir.create(out)

  pred <- data.frame(
    name = "pocket1",
    rank = 1L,
    score = 8.5,
    probability = 0.79,
    sas_points = 17L,
    surf_atoms = 10L,
    center_x = 1,
    center_y = 2,
    center_z = 3,
    residue_ids = "A_10 A_11",
    stringsAsFactors = FALSE
  )

  norm <- ProtVis:::.protvis_pw_p2rank_normalize(pred, pdb, out)
  expect_equal(nrow(norm), 1L)
  expect_equal(norm$pocket[[1L]], 1L)
  expect_equal(norm$probability[[1L]], 0.79)
  expect_equal(norm$residue_count[[1L]], 2L)
  expect_true(file.exists(norm$pocket_file[[1L]]))
})

test_that("P2Rank residue IDs preserve chain and residue number", {
  ids <- ProtVis:::.protvis_pw_p2rank_residue_ids("A_103 A_180 B_42")

  expect_equal(nrow(ids), 3L)
  expect_equal(ids$chain, c("A", "A", "B"))
  expect_equal(ids$resi, c(103L, 180L, 42L))
})

test_that("P2Rank backend supports AlphaFold profile and standard outputs", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzzz_protein_workbench_p2rank.R")),
    collapse = "\n"
  )

  expect_match(source, '"predict"', fixed = TRUE)
  expect_match(source, '"-c", "alphafold"', fixed = TRUE)
  expect_match(source, "_predictions\\\\.csv$", fixed = TRUE)
  expect_match(source, "_residues\\\\.csv$", fixed = TRUE)
  expect_false(grepl("-export_pocket_descriptors", source, fixed = TRUE))
  expect_false(grepl("-pocket_grid_format", source, fixed = TRUE))
  expect_match(source, "pocket_backend = rv_pocket$result$backend", fixed = TRUE)
  expect_match(source, ".protvis_pw_pocket_view", fixed = TRUE)
  expect_match(source, "cz.siret.prank.program.Main", fixed = TRUE)
  expect_match(source, ".protvis_pw_p2rank_classpath", fixed = TRUE)
  expect_false(grepl(".protvis_pw_p2rank_batch_command", source, fixed = TRUE))
})


test_that("P2Rank path cleaning handles quoted Windows paths", {
  cleaned <- ProtVis:::.protvis_pw_p2rank_clean_path(
    '"E:\\ount\\p2rank\\p2rank_2.5.1"'
  )
  expect_equal(cleaned, "E:/ount/p2rank/p2rank_2.5.1")
})

test_that("P2Rank diagnostics report folder and prank.bat presence", {
  root <- tempfile("p2rank_diag_")
  dir.create(root)
  bat <- file.path(root, "prank.bat")
  writeLines("@echo off", bat)

  diag <- ProtVis:::.protvis_pw_p2rank_path_diagnostics(root)
  expect_true(diag$directory_exists)
  expect_true(diag$candidate_exists[["prank_bat"]])
  expect_true(diag$found)
})


test_that("Binding pocket buttons keep icons fully visible", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_icon_test"))
  expect_match(html, "binding_pocket_root", fixed = TRUE)
  expect_match(html, "overflow:visible", fixed = TRUE)
  expect_match(html, "display:inline-flex", fixed = TRUE)
})

test_that("Built-in 1HEL structure is previewed before P2Rank starts", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzzz_protein_workbench_p2rank.R")),
    collapse = "\\n"
  )

  demo_start <- regexpr(
    "shiny::observeEvent(input$use_pocket_demo",
    source,
    fixed = TRUE
  )[[1L]]
  click_start <- regexpr(
    'shinyjs::click(session$ns("run_pocket"))',
    source,
    fixed = TRUE
  )[[1L]]

  expect_gt(demo_start, 0L)
  expect_gt(click_start, demo_start)
  demo_setup <- substr(source, demo_start, click_start)
  expect_match(demo_setup, "rv_pocket$pdb_path <- pocket_demo_path", fixed = TRUE)
  expect_match(demo_setup, "rv_pocket$pdb_text <- base::paste(", fixed = TRUE)
  expect_match(demo_setup, "rv_pocket$result <- NULL", fixed = TRUE)
  expect_match(demo_setup, "Built-in 1HEL structure loaded.", fixed = TRUE)
})

test_that("1HEL preview flushes before the P2Rank run is triggered", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzzz_protein_workbench_p2rank.R")),
    collapse = "\\n"
  )

  preview_message <- regexpr(
    "Built-in 1HEL structure loaded.",
    source,
    fixed = TRUE
  )[[1L]]
  flush_callback <- regexpr(
    "session$onFlushed(",
    source,
    fixed = TRUE
  )[[1L]]
  run_click <- regexpr(
    'shinyjs::click(session$ns("run_pocket"))',
    source,
    fixed = TRUE
  )[[1L]]

  expect_gt(preview_message, 0L)
  expect_gt(flush_callback, preview_message)
  expect_gt(run_click, flush_callback)
})
