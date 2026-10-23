test_that("Protein Workbench evolutionary analysis exposes common plant species", {
  species <- ProtVis:::.protvis_pw_plant_species()

  expect_equal(species[["Zea mays (maize)"]], "4577")
  expect_equal(species[["Oryza sativa (rice)"]], "4530")
  expect_equal(species[["Arabidopsis thaliana"]], "3702")
  expect_equal(species[["Glycine max (soybean)"]], "3847")
  expect_true(length(species) >= 10L)
})

test_that("Protein Workbench evolutionary candidate selection favors comparable reviewed proteins", {
  hits <- data.frame(
    accession = c("A", "B", "C"),
    length = c(1000L, 510L, 495L),
    entry_type = c("TrEMBL", "UniProtKB reviewed (Swiss-Prot)", "TrEMBL"),
    stringsAsFactors = FALSE
  )
  picked <- ProtVis:::.protvis_pw_evolution_pick_candidate(hits, query_length = 500L)
  expect_equal(picked$accession[[1L]], "B")
})

test_that("Protein Workbench includes plant Evolutionary analysis controls and downloads", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_evolution_test"))

  expect_match(html, "Evolutionary analysis", fixed = TRUE)
  expect_match(html, "Common plant species", fixed = TRUE)
  expect_match(html, "RUN EVOLUTIONARY ANALYSIS", fixed = TRUE)
  expect_match(html, "Neighbor joining", fixed = TRUE)
  expect_match(html, "Maximum likelihood", fixed = TRUE)

  expected_ids <- c(
    "evo_species",
    "evo_method",
    "run_evolution",
    "download_evo_table",
    "download_evo_alignment",
    "download_evo_tree",
    "download_evo_png",
    "download_evo_svg",
    "download_evo_pdf"
  )
  for (id in expected_ids) {
    expect_match(html, id, fixed = TRUE)
  }
})

test_that("Protein Workbench evolution uses DECIPHER alignment and phylogeny", {
  source <- paste(
    readLines(testthat::test_path("..", "..", "R", "zzzzzz_protein_workbench_evolution.R")),
    collapse = "\n"
  )

  expect_match(source, "DECIPHER::AlignSeqs", fixed = TRUE)
  expect_match(source, "DECIPHER::Treeline", fixed = TRUE)
  expect_match(source, "ape::as.phylo", fixed = TRUE)
  expect_match(source, "Current Protein Workbench query", fixed = TRUE)
  expect_match(source, "UniProtKB", fixed = TRUE)
})
