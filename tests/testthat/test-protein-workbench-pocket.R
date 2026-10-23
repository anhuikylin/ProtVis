test_that("Protein Workbench fpocket parser reads ranked pocket descriptors", {
  info <- tempfile(fileext = "_info.txt")
  writeLines(c(
    "Pocket 1 :",
    "Score : 42.500",
    "Druggability Score : 0.812",
    "Number of Alpha Spheres : 35",
    "Total SASA : 511.2",
    "Polar SASA : 122.4",
    "Apolar SASA : 388.8",
    "Volume : 603.7",
    "Hydrophobicity score : 18.2",
    "",
    "Pocket 2 :",
    "Score : 21.100",
    "Druggability Score : 0.401",
    "Number of Alpha Spheres : 18",
    "Volume : 301.5"
  ), info)

  parsed <- ProtVis:::.protvis_pw_fpocket_parse_info(info)
  expect_equal(nrow(parsed), 2L)
  expect_equal(parsed$pocket[[1L]], 1L)
  expect_equal(parsed$score[[1L]], 42.5)
  expect_equal(parsed$druggability_score[[1L]], 0.812)
  expect_equal(parsed$volume[[1L]], 603.7)
})

test_that("Protein Workbench pocket residue parser extracts PDB residue identifiers", {
  pdb <- tempfile(fileext = ".pdb")
  writeLines(c(
    "ATOM      1  N   GLY A  10      11.104  13.207  10.111  1.00 20.00           N",
    "ATOM      2  CA  GLY A  10      12.000  13.500  10.500  1.00 20.00           C",
    "ATOM      3  N   SER A  11      13.104  14.207  11.111  1.00 20.00           N"
  ), pdb)

  residues <- ProtVis:::.protvis_pw_pdb_residues(pdb)
  expect_equal(nrow(residues), 2L)
  expect_equal(residues$resi, c(10L, 11L))
  expect_equal(residues$chain, c("A", "A"))
})

test_that("Protein Workbench exposes binding pocket prediction controls", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_pocket_test"))

  expect_match(html, "Binding pocket prediction", fixed = TRUE)
  expect_match(html, "RUN POCKET PREDICTION", fixed = TRUE)
  expect_match(html, "Current AlphaFold structure", fixed = TRUE)
  expect_match(html, "Upload PDB", fixed = TRUE)
  expect_match(html, "fpocket executable", fixed = TRUE)

  for (id in c(
    "pocket_structure_source",
    "pocket_pdb",
    "fpocket_path",
    "run_pocket",
    "pocket_selected",
    "download_pocket_table",
    "download_pocket_pdb",
    "download_fpocket_results",
    "pocket_view"
  )) {
    expect_match(html, id, fixed = TRUE)
  }
})

test_that("Protein Workbench fpocket integration keeps execution and visualization local", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzz_protein_workbench_pocket.R")),
    collapse = "\n"
  )

  expect_match(source, 'Sys.which("fpocket")', fixed = TRUE)
  expect_match(source, "system2(", fixed = TRUE)
  expect_match(source, 'args = c("-f"', fixed = TRUE)
  expect_match(source, "m_style_sphere", fixed = TRUE)
  expect_match(source, "druggability_score", fixed = TRUE)
  expect_match(source, "Full fpocket results", fixed = TRUE)
})


test_that("Protein Workbench binding pocket includes built-in 1HEL demo", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_pocket_demo_test"))

  expect_match(html, "Built-in demo · 1HEL lysozyme", fixed = TRUE)
  expect_match(html, "USE BUILT-IN DEMO", fixed = TRUE)
  expect_match(html, "Download demo PDB", fixed = TRUE)
  expect_match(html, "hen egg white lysozyme", fixed = TRUE)
  expect_match(html, "1.70 Å", fixed = TRUE)

  for (id in c(
    "use_pocket_demo",
    "download_pocket_demo"
  )) {
    expect_match(html, id, fixed = TRUE)
  }
})

test_that("Protein Workbench pocket demo reuses bundled ProtVisDatabase structure", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzz_protein_workbench_pocket.R")),
    collapse = "\n"
  )

  expect_match(source, '.protvis_pw_pocket_demo_path', fixed = TRUE)
  expect_match(source, '"extdata", "structure", "1hel.pdb"', fixed = TRUE)
  expect_match(source, 'package = "ProtVisDatabase"', fixed = TRUE)
  expect_match(source, 'selected = "demo"', fixed = TRUE)
  expect_match(source, 'shinyjs::click(session$ns("run_pocket"))', fixed = TRUE)
  expect_match(source, 'ProtVis_binding_pocket_demo_1HEL.pdb', fixed = TRUE)
})


test_that("Protein Workbench fpocket backend supports Docker fallback", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzz_protein_workbench_pocket.R")),
    collapse = "\n"
  )

  expect_match(source, ".protvis_pw_docker_executable", fixed = TRUE)
  expect_match(source, "official fpocket/fpocket image", fixed = TRUE)
  expect_match(source, '"fpocket/fpocket"', fixed = TRUE)
  expect_match(source, '"run", "--rm"', fixed = TRUE)
  expect_match(source, "Docker Desktop is the recommended fallback on Windows.", fixed = TRUE)
  expect_match(source, "fpocket_backend = rv_pocket$result$backend", fixed = TRUE)
})


test_that("Protein Workbench distinguishes Docker CLI from running daemon", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzz_protein_workbench_pocket.R")),
    collapse = "\n"
  )

  expect_match(source, ".protvis_pw_docker_status", fixed = TRUE)
  expect_match(source, '"info", "--format", "{{.ServerVersion}}"', fixed = TRUE)
  expect_match(source, 'type = "docker_stopped"', fixed = TRUE)
  expect_match(source, "Docker Desktop is installed but its engine is not running.", fixed = TRUE)
  expect_match(source, "Start Docker Desktop and wait until it is fully ready", fixed = TRUE)
})


test_that("Protein Workbench diagnoses missing Docker image and registry failures", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzzz_protein_workbench_pocket.R")),
    collapse = "\n"
  )

  expect_match(source, ".protvis_pw_docker_image_status", fixed = TRUE)
  expect_match(source, ".protvis_pw_docker_pull_image", fixed = TRUE)
  expect_match(source, '"image", "inspect"', fixed = TRUE)
  expect_match(source, '"pull", image', fixed = TRUE)
  expect_match(source, "registry-1\\\\.docker\\\\.io", fixed = TRUE)
  expect_match(source, "Docker is running, but the image could not be downloaded from the registry.", fixed = TRUE)
  expect_match(source, "fpocket_docker_image", fixed = TRUE)
  expect_match(source, "PULL / CHECK FPOCKET IMAGE", fixed = TRUE)
})
