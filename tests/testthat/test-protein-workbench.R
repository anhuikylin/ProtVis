test_that("Protein Workbench sequence helpers return stable summaries", {
  seq <- "VLSPADKTNVKAAWAKVGNHAADFGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"

  stats <- ProtVis:::.protvis_pw_sequence_stats(seq)
  expect_equal(stats$value[stats$metric == "Length"], "141 aa")
  expect_true(grepl("Da$", stats$value[stats$metric == "Molecular weight"]))

  composition <- ProtVis:::.protvis_pw_composition(seq)
  expect_equal(sum(composition$Count), nchar(seq))
  expect_equal(round(sum(composition$Fraction), 8), 1)

  hydropathy <- ProtVis:::.protvis_pw_hydropathy(seq, 9)
  expect_equal(nrow(hydropathy), nchar(seq))
})

test_that("Protein Workbench parses UniProt-style feature records", {
  entry <- list(
    primaryAccession = "PTEST1",
    uniProtkbId = "PTEST1_TEST",
    proteinDescription = list(
      recommendedName = list(fullName = list(value = "Test protein"))
    ),
    genes = list(list(geneName = list(value = "TEST"))),
    organism = list(scientificName = "Test species", taxonId = 1),
    sequence = list(value = "ACDEFGHIKLMNPQRSTVWY", length = 20, molWeight = 2200),
    features = list(
      list(
        type = "Modified residue",
        description = "Phosphoserine",
        location = list(start = list(value = 5), end = list(value = 5))
      ),
      list(
        type = "Domain",
        description = "Example domain",
        location = list(start = list(value = 2), end = list(value = 18))
      )
    ),
    uniProtKBCrossReferences = list(
      list(database = "PDB", id = "1ABC", properties = list()),
      list(database = "InterPro", id = "IPR000001", properties = list())
    )
  )

  features <- ProtVis:::.protvis_pw_features_table(entry)
  expect_equal(nrow(features), 2)
  expect_equal(features$start, c(5L, 2L))

  ptm <- ProtVis:::.protvis_pw_ptm_table(entry)
  expect_equal(nrow(ptm), 1)
  expect_match(ptm$description, "Phosphoserine")

  xrefs <- ProtVis:::.protvis_pw_xrefs_table(entry)
  expect_true(all(c("PDB", "InterPro") %in% xrefs$database))
})

test_that("Protein Workbench is exposed as a separate Toolkits module", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_test"))
  expect_match(html, "Protein Workbench", fixed = TRUE)
  expect_match(html, "PTM &amp; sites")
  expect_match(html, "AlphaFold structure", fixed = TRUE)
  expect_match(html, "Cross-references", fixed = TRUE)

  app_html <- as.character(ProtVis::app_ui(NULL))
  expect_match(app_html, "Protein Workbench", fixed = TRUE)
  expect_match(app_html, "Protein Extract", fixed = TRUE)
  expect_match(app_html, "Plant-mPLoc", fixed = TRUE)
  expect_match(app_html, "swissmodel", fixed = TRUE)
})


test_that("Protein Workbench static plots expose PNG, SVG and PDF downloads", {
  html <- as.character(ProtVis::protein_workbench_ui("pw_download_test"))
  expected_ids <- c(
    "download_composition_plot_png",
    "download_composition_plot_svg",
    "download_composition_plot_pdf",
    "download_hydropathy_plot_png",
    "download_hydropathy_plot_svg",
    "download_hydropathy_plot_pdf",
    "download_domain_plot_png",
    "download_domain_plot_svg",
    "download_domain_plot_pdf",
    "download_interaction_plot_png",
    "download_interaction_plot_svg",
    "download_interaction_plot_pdf"
  )
  for (id in expected_ids) {
    expect_match(html, id, fixed = TRUE)
  }

  base_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  interaction_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "zzzz_protein_workbench_interaction_diagnostics.R")),
    collapse = "\n"
  )
  expect_match(base_source, "ProtVis_residue_composition", fixed = TRUE)
  expect_match(base_source, "ProtVis_Kyte_Doolittle_hydropathy", fixed = TRUE)
  expect_match(base_source, "ProtVis_domain_architecture", fixed = TRUE)
  expect_match(interaction_source, "ProtVis_STRING_interaction_network", fixed = TRUE)
})


test_that("Protein Workbench exposes AlphaFold PDB download support", {
  module_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  expect_match(module_source, "DOWNLOAD PDB", fixed = TRUE)
  expect_match(module_source, "download_alphafold_pdb", fixed = TRUE)
  expect_match(module_source, "_AlphaFold.pdb", fixed = TRUE)
  expect_match(module_source, 'contentType = "chemical/x-pdb"', fixed = TRUE)
})
