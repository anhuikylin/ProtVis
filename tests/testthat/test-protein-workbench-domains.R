test_that("Protein Workbench domain tracks use labels and adaptive height", {
  domains <- data.frame(
    accession = c("IPR0001", "IPR0002", "IPR0001"),
    name = c(
      "ATPase catalytic domain",
      "A very long transmembrane-associated domain annotation used for display testing",
      "ATPase catalytic domain"
    ),
    start = c(12L, 155L, 420L),
    end = c(120L, 255L, 590L),
    stringsAsFactors = FALSE
  )

  plot_data <- ProtVis:::.protvis_pw_domain_plot_data(domains)
  expect_equal(nrow(plot_data), 3)
  expect_equal(plot_data$label[[1]], "ATPase catalytic domain")
  expect_true(all(nchar(plot_data$label_display) <= 52L))
  expect_equal(length(unique(plot_data$track)), 3)
  expect_gt(
    ProtVis:::.protvis_pw_domain_plot_height(nrow(plot_data)),
    ProtVis:::.protvis_pw_domain_plot_height(0L)
  )
  expect_equal(ProtVis:::.protvis_pw_domain_plot_height(100L), 1600L)
})

test_that("Protein Workbench selects species-aware annotation resources", {
  maize_entry <- list(
    organism = list(scientificName = "Zea mays"),
    genes = list(list(geneName = list(value = "Zm00001eb000210")))
  )
  maize_links <- ProtVis:::.protvis_pw_external_links("A0A1D6JJK6", maize_entry)
  expect_true(all(c("Ensembl", "NCBI_Gene", "KEGG_Genes", "Plant_Reactome", "MaizeGDB") %in% names(maize_links)))
  expect_match(maize_links$MaizeGDB, "Zm00001eb000210", fixed = TRUE)

  human_entry <- list(
    organism = list(scientificName = "Homo sapiens"),
    genes = list(list(geneName = list(value = "HBB")))
  )
  human_links <- ProtVis:::.protvis_pw_external_links("P68871", human_entry)
  expect_false("MaizeGDB" %in% names(human_links))
  expect_match(human_links$Ensembl, "HBB", fixed = TRUE)
})

test_that("Protein Workbench uses a server namespace for adaptive domain output", {
  module_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  expect_match(module_source, 'plotOutput\\(session\\$ns\\("domain_plot"\\)')
  expect_match(module_source, "protvis_pw_empty_plot", fixed = TRUE)
})
