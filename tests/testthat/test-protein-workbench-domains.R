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

test_that("Protein Workbench parses and exports batch protein sequences", {
  expect_equal(
    ProtVis:::.protvis_pw_parse_identifiers("Zm00001eb000210, P68871\nZm00001eb000210"),
    c("Zm00001eb000210", "P68871")
  )
  expect_error(ProtVis:::.protvis_pw_parse_identifiers("bad id!"), "Identifiers may contain")

  fasta <- paste(
    ">sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2",
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKK",
    sep = "\n"
  )
  records <- ProtVis:::.protvis_pw_parse_fasta_records(fasta)
  expect_equal(records$accession, "P68871")
  expect_equal(records$gene, "HBB")
  expect_equal(records$taxon_id, "9606")
  expect_gt(records$length, 50)

  records$input_id <- "P68871"
  records$status <- "Retrieved"
  exported <- ProtVis:::.protvis_pw_batch_fasta(records)
  expect_match(exported, ">P68871|HBB", fixed = TRUE)
  expect_match(exported, "MVHLTPEE", fixed = TRUE)
})

test_that("Protein Workbench offers common species and a batch retrieval tab", {
  species <- ProtVis:::.protvis_pw_common_species()
  maize_versions <- ProtVis:::.protvis_pw_maize_versions()
  expect_equal(species[["Zea mays (maize)"]], "4577")
  expect_equal(species[["Homo sapiens"]], "9606")
  expect_equal(maize_versions[["B73 RefGen_v5 / NAM v5 (Zm00001eb...)"]], "b73_v5")
  expect_equal(maize_versions[["B73 RefGen_v4 (Zm00001d...)"]], "b73_v4")
  expect_equal(maize_versions[["B73 RefGen_v3 (GRMZM2G...)"]], "b73_v3")
  expect_equal(
    ProtVis:::.protvis_pw_maize_version_label("b73_v4"),
    "B73 RefGen_v4 (Zm00001d...)"
  )
  html <- as.character(ProtVis::protein_workbench_ui("pw_batch_test"))
  expect_match(html, "Batch sequences", fixed = TRUE)
  expect_match(html, "RETRIEVE SEQUENCES", fixed = TRUE)
  expect_match(html, "Maize reference version", fixed = TRUE)
  expect_match(tolower(html), "no fasta upload is needed", fixed = TRUE)
})

test_that("Protein Workbench places batch retrieval in the widened sidebar", {
  module_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  expect_match(module_source, "width = 550", fixed = TRUE)
  expect_equal(
    length(regmatches(module_source, gregexpr('ns\\("batch_ids"\\)', module_source))[[1]]),
    1L
  )
})

test_that("Protein Workbench batch retrieval protects exact and fallback matching", {
  module_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  expect_match(module_source, "gene_exact:", fixed = TRUE)
  expect_match(module_source, "protvis_pw_search_uniprot(input_id", fixed = TRUE)
  expect_match(module_source, "uniprotkb/stream", fixed = TRUE)
})

test_that("Protein Workbench confirms completed batch retrieval", {
  module_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  expect_match(module_source, "Sequence retrieval completed:", fixed = TRUE)
  expect_match(module_source, 'type = if (retrieved > 0L) "message" else "warning"', fixed = TRUE)
  expect_match(module_source, "pw-batch-status-success", fixed = TRUE)
})

test_that("Protein Workbench shows batch completion directly below the action", {
  module_source <- base::paste(
    base::readLines(testthat::test_path("..", "..", "R", "protein_workbench.R")),
    collapse = "\n"
  )
  button_at <- regexpr('ns\\("batch_run"\\)', module_source)[[1]]
  status_at <- regexpr('shiny::uiOutput\\(ns\\("batch_status"\\)\\)', module_source)[[1]]
  download_at <- regexpr('ns\\("download_batch_fasta"\\)', module_source)[[1]]
  expect_gt(status_at, button_at)
  expect_lt(status_at, download_at)
})

test_that("Protein Workbench provides maize fallback resource links", {
  links <- ProtVis:::.protvis_pw_maize_fallback_links("Zm00001eb000210")
  expect_match(links, "MaizeGDB", fixed = TRUE)
  expect_match(links, "Phytozome", fixed = TRUE)
  expect_match(links, "Zm00001eb000210", fixed = TRUE)
})
