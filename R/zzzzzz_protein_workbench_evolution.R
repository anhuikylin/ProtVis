# Protein Workbench evolutionary analysis -------------------------------
#
# Additive, protein-centric evolutionary analysis for common plant species.
# Homolog candidates are resolved from UniProt using the current protein's
# gene/protein annotation, aligned locally with DECIPHER, and converted to a
# phylogenetic tree in the Workbench. Existing Protein Workbench functionality
# is preserved unchanged.

.protvis_pw_pre_evolution_ui <- protein_workbench_ui
.protvis_pw_pre_evolution_server <- protein_workbench_server

.protvis_pw_plant_species <- function() {
  c(
    "Zea mays (maize)" = "4577",
    "Oryza sativa (rice)" = "4530",
    "Sorghum bicolor" = "4558",
    "Setaria italica" = "4555",
    "Triticum aestivum (wheat)" = "4565",
    "Hordeum vulgare (barley)" = "4513",
    "Brachypodium distachyon" = "15368",
    "Arabidopsis thaliana" = "3702",
    "Brassica napus" = "3708",
    "Glycine max (soybean)" = "3847",
    "Medicago truncatula" = "3880",
    "Solanum lycopersicum (tomato)" = "4081",
    "Solanum tuberosum (potato)" = "4113",
    "Gossypium hirsutum (cotton)" = "3635"
  )
}

.protvis_pw_primary_gene <- function(entry) {
  if (base::is.null(entry)) return("")
  genes <- entry$genes %||% base::list()
  if (!base::length(genes)) return("")
  gene <- genes[[1L]]$geneName$value %||% ""
  base::trimws(base::as.character(gene))
}

.protvis_pw_evolution_protein_name <- function(entry) {
  if (base::is.null(entry)) return("")
  name <- base::tryCatch(.protvis_pw_protein_name(entry), error = function(e) "")
  name <- base::trimws(base::as.character(name %||% ""))
  if (!base::nzchar(name)) return("")
  # UniProt descriptions can include organism-specific qualifiers; removing
  # parenthetical suffixes improves cross-species retrieval while retaining
  # the biologically meaningful protein name.
  name <- base::sub("\\s*\\([^()]+\\)\\s*$", "", name)
  base::trimws(name)
}

.protvis_pw_escape_uniprot_phrase <- function(x) {
  x <- base::as.character(x %||% "")
  x <- base::gsub("\\\\", "\\\\\\\\", x)
  base::gsub('"', '\\"', x, fixed = TRUE)
}

.protvis_pw_evolution_candidate_search <- function(entry, taxon_id, size = 5L) {
  if (base::is.null(entry) || !base::nzchar(base::as.character(taxon_id %||% ""))) {
    return(base::data.frame())
  }
  gene <- .protvis_pw_primary_gene(entry)
  protein_name <- .protvis_pw_evolution_protein_name(entry)

  queries <- base::character()
  if (base::nzchar(gene)) {
    queries <- c(queries, base::paste0("gene_exact:", gene))
  }
  if (base::nzchar(protein_name)) {
    queries <- c(
      queries,
      base::paste0('protein_name:"', .protvis_pw_escape_uniprot_phrase(protein_name), '"')
    )
  }
  if (!base::length(queries)) return(base::data.frame())

  for (query in queries) {
    hits <- base::tryCatch(
      .protvis_pw_search_uniprot(
        query = query,
        organism_id = base::as.character(taxon_id),
        size = base::as.integer(size)
      ),
      error = function(e) base::data.frame()
    )
    if (base::nrow(hits)) {
      hits$retrieval_query <- query
      return(hits)
    }
  }
  base::data.frame()
}

.protvis_pw_evolution_pick_candidate <- function(hits, query_length) {
  if (base::is.null(hits) || !base::nrow(hits)) return(NULL)
  qlen <- base::as.numeric(query_length %||% NA_real_)
  if (!base::is.finite(qlen) || qlen <= 0) qlen <- stats::median(hits$length, na.rm = TRUE)

  len <- base::as.numeric(hits$length)
  length_penalty <- base::abs(len - qlen) / base::max(1, qlen)
  length_penalty[!base::is.finite(length_penalty)] <- 10

  reviewed <- base::grepl("reviewed|Swiss-Prot", hits$entry_type, ignore.case = TRUE)
  reviewed_bonus <- ifelse(reviewed, -0.10, 0)
  score <- length_penalty + reviewed_bonus

  hits[base::which.min(score), , drop = FALSE]
}

.protvis_pw_evolution_fetch <- function(entry, selected_taxa) {
  if (base::is.null(entry)) base::stop("Resolve a UniProt protein before evolutionary analysis.", call. = FALSE)
  selected_taxa <- base::unique(base::as.character(selected_taxa %||% base::character()))
  selected_taxa <- selected_taxa[base::nzchar(selected_taxa)]
  if (base::length(selected_taxa) < 2L) {
    base::stop("Select at least two plant species.", call. = FALSE)
  }

  query_accession <- base::as.character(entry$primaryAccession %||% "")
  query_sequence <- .protvis_pw_sequence(entry)
  query_taxon <- base::as.character(entry$organism$taxonId %||% "")
  query_species <- base::as.character(entry$organism$scientificName %||% "")
  query_gene <- .protvis_pw_primary_gene(entry)
  query_name <- .protvis_pw_evolution_protein_name(entry)
  query_length <- base::nchar(query_sequence)

  species_map <- .protvis_pw_plant_species()
  taxon_to_label <- stats::setNames(base::names(species_map), base::unname(species_map))

  rows <- base::list()
  row_id <- 0L

  for (taxon in selected_taxa) {
    species_label <- taxon_to_label[[taxon]] %||% base::paste0("Taxon ", taxon)

    if (base::nzchar(query_taxon) && identical(taxon, query_taxon) &&
        base::nzchar(query_sequence)) {
      row_id <- row_id + 1L
      rows[[row_id]] <- base::data.frame(
        accession = query_accession,
        gene = query_gene,
        protein_name = query_name,
        species = query_species,
        taxon_id = query_taxon,
        length = query_length,
        sequence = query_sequence,
        source = "Current Protein Workbench query",
        status = "Query",
        stringsAsFactors = FALSE
      )
      next
    }

    hits <- .protvis_pw_evolution_candidate_search(entry, taxon, size = 5L)
    if (!base::nrow(hits)) {
      row_id <- row_id + 1L
      rows[[row_id]] <- base::data.frame(
        accession = "",
        gene = "",
        protein_name = "",
        species = species_label,
        taxon_id = taxon,
        length = NA_integer_,
        sequence = "",
        source = "UniProtKB",
        status = "Not found",
        stringsAsFactors = FALSE
      )
      next
    }

    candidate <- .protvis_pw_evolution_pick_candidate(hits, query_length)
    accession <- base::as.character(candidate$accession[[1L]] %||% "")
    candidate_entry <- base::tryCatch(
      .protvis_pw_get_uniprot(accession),
      error = function(e) NULL
    )
    seq <- if (base::is.null(candidate_entry)) "" else .protvis_pw_sequence(candidate_entry)
    gene <- if (base::is.null(candidate_entry)) {
      base::as.character(candidate$genes[[1L]] %||% "")
    } else {
      .protvis_pw_primary_gene(candidate_entry)
    }
    pname <- if (base::is.null(candidate_entry)) {
      base::as.character(candidate$protein_name[[1L]] %||% "")
    } else {
      .protvis_pw_evolution_protein_name(candidate_entry)
    }

    row_id <- row_id + 1L
    rows[[row_id]] <- base::data.frame(
      accession = accession,
      gene = gene,
      protein_name = pname,
      species = base::as.character(candidate$organism[[1L]] %||% species_label),
      taxon_id = taxon,
      length = if (base::nzchar(seq)) base::nchar(seq) else base::as.integer(candidate$length[[1L]]),
      sequence = seq,
      source = base::paste0("UniProtKB · ", candidate$retrieval_query[[1L]] %||% "annotation match"),
      status = if (base::nzchar(seq)) "Retrieved" else "Sequence unavailable",
      stringsAsFactors = FALSE
    )
  }

  table <- base::do.call(base::rbind, rows)
  base::rownames(table) <- NULL

  # Ensure the current query remains represented even if its species was not
  # selected explicitly. This also provides a stable root/reference label.
  if (base::nzchar(query_sequence) &&
      !base::any(table$status == "Query")) {
    query_row <- base::data.frame(
      accession = query_accession,
      gene = query_gene,
      protein_name = query_name,
      species = query_species,
      taxon_id = query_taxon,
      length = query_length,
      sequence = query_sequence,
      source = "Current Protein Workbench query",
      status = "Query",
      stringsAsFactors = FALSE
    )
    table <- base::rbind(query_row, table)
  }

  table
}

.protvis_pw_evolution_tip_label <- function(row, query = FALSE) {
  gene <- base::trimws(base::as.character(row$gene %||% ""))
  accession <- base::trimws(base::as.character(row$accession %||% ""))
  species <- base::trimws(base::as.character(row$species %||% ""))
  id <- if (base::nzchar(gene)) gene else accession
  if (!base::nzchar(id)) id <- "protein"
  prefix <- if (isTRUE(query)) "QUERY_" else ""
  base::paste0(prefix, id, " | ", species)
}

.protvis_pw_evolution_alignment <- function(table, processors = 1L) {
  if (!base::requireNamespace("DECIPHER", quietly = TRUE)) {
    base::stop("DECIPHER is required for Protein Workbench evolutionary analysis.", call. = FALSE)
  }
  keep <- table$status %in% c("Query", "Retrieved") & base::nzchar(table$sequence)
  data <- table[keep, , drop = FALSE]
  if (base::nrow(data) < 3L) {
    base::stop("At least three protein sequences are required to construct a phylogenetic tree.", call. = FALSE)
  }

  labels <- base::vapply(base::seq_len(base::nrow(data)), function(i) {
    .protvis_pw_evolution_tip_label(data[i, , drop = FALSE], query = identical(data$status[[i]], "Query"))
  }, character(1))
  labels <- base::make.unique(labels)

  aa <- Biostrings::AAStringSet(data$sequence)
  base::names(aa) <- labels
  DECIPHER::AlignSeqs(
    aa,
    processors = base::max(1L, base::as.integer(processors)),
    verbose = FALSE
  )
}

.protvis_pw_evolution_tree <- function(alignment, method = "NJ", processors = 1L) {
  if (!base::requireNamespace("DECIPHER", quietly = TRUE)) {
    base::stop("DECIPHER is required for Protein Workbench evolutionary analysis.", call. = FALSE)
  }
  method <- base::toupper(base::as.character(method %||% "NJ"))
  if (!method %in% c("NJ", "ME", "ML")) method <- "NJ"

  dend <- DECIPHER::Treeline(
    myXStringSet = alignment,
    method = method,
    type = "dendrogram",
    processors = base::max(1L, base::as.integer(processors)),
    verbose = FALSE
  )

  hc <- stats::as.hclust(dend)
  phy <- ape::as.phylo(hc)
  phy
}

.protvis_pw_evolution_ui <- function(id) {
  ns <- shiny::NS(id)
  species <- .protvis_pw_plant_species()
  default_species <- species[c(
    "Zea mays (maize)",
    "Oryza sativa (rice)",
    "Sorghum bicolor",
    "Setaria italica",
    "Arabidopsis thaliana",
    "Glycine max (soybean)",
    "Solanum lycopersicum (tomato)"
  )]

  shiny::div(
    style = "padding:0 16px 24px;",
    bslib::card(
      bslib::card_header(
        shiny::div(
          shiny::tags$strong("Evolutionary analysis"),
          shiny::div(
            "Plant-focused homolog retrieval, multiple-sequence alignment and phylogenetic reconstruction for the currently resolved protein.",
            style = "font-size:12px;color:#657789;font-weight:400;margin-top:2px;"
          )
        )
      ),
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(4, 8),
          bslib::card(
            bslib::card_header("Plant homolog settings"),
            shiny::checkboxGroupInput(
              ns("evo_species"),
              "Common plant species",
              choices = species,
              selected = base::unname(default_species)
            ),
            shiny::selectInput(
              ns("evo_method"),
              "Phylogenetic method",
              choices = c(
                "Neighbor joining (fast)" = "NJ",
                "Minimum evolution" = "ME",
                "Maximum likelihood" = "ML"
              ),
              selected = "NJ"
            ),
            shiny::numericInput(
              ns("evo_processors"),
              "Processors",
              value = 1,
              min = 1,
              max = 16,
              step = 1
            ),
            shiny::actionButton(
              ns("run_evolution"),
              "RUN EVOLUTIONARY ANALYSIS",
              icon = bsicons::bs_icon("diagram-2"),
              class = "btn-primary"
            ),
            shiny::uiOutput(ns("evo_status")),
            shiny::hr(),
            shiny::downloadButton(ns("download_evo_table"), "Homolog table", class = "btn-sm btn-outline-primary"),
            shiny::downloadButton(ns("download_evo_alignment"), "Aligned FASTA", class = "btn-sm btn-outline-primary"),
            shiny::downloadButton(ns("download_evo_tree"), "Newick tree", class = "btn-sm btn-outline-primary")
          ),
          shiny::div(
            bslib::card(
              full_screen = TRUE,
              bslib::card_header(
                shiny::div(
                  style = "display:flex;justify-content:space-between;align-items:center;gap:10px;width:100%;",
                  shiny::span("Plant phylogenetic tree"),
                  shiny::div(
                    shiny::downloadButton(ns("download_evo_png"), "PNG", class = "btn-sm btn-outline-primary"),
                    shiny::downloadButton(ns("download_evo_svg"), "SVG", class = "btn-sm btn-outline-primary"),
                    shiny::downloadButton(ns("download_evo_pdf"), "PDF", class = "btn-sm btn-outline-primary")
                  )
                )
              ),
              shiny::plotOutput(ns("evo_tree_plot"), height = "560px")
            )
          )
        ),
        shiny::br(),
        bslib::card(
          bslib::card_header("Retrieved plant homolog candidates"),
          DT::DTOutput(ns("evo_table"))
        )
      )
    )
  )
}

#' Protein Workbench UI
#'
#' Extends Protein Workbench with plant-focused evolutionary analysis.
#'
#' @param id Shiny module id.
#' @return Shiny UI.
#' @export
protein_workbench_ui <- function(id) {
  shiny::tagList(
    .protvis_pw_pre_evolution_ui(id),
    .protvis_pw_evolution_ui(id)
  )
}

#' Protein Workbench server
#'
#' Extends Protein Workbench with plant homolog retrieval, DECIPHER alignment
#' and phylogenetic reconstruction.
#'
#' @param id Shiny module id.
#' @param shared_state Optional shared ProtVis reactive state.
#' @return Shiny module server.
#' @export
protein_workbench_server <- function(id, shared_state = NULL) {
  .protvis_pw_pre_evolution_server(id, shared_state = shared_state)

  shiny::moduleServer(id, function(input, output, session) {
    rv_evo <- shiny::reactiveValues(
      entry = NULL,
      homologs = base::data.frame(),
      alignment = NULL,
      tree = NULL,
      message = "Resolve a UniProt protein above, select plant species, then run evolutionary analysis."
    )

    current_accession_evo <- shiny::reactive({
      base::trimws(base::as.character(input$accession %||% ""))
    })

    load_entry <- function() {
      accession <- current_accession_evo()
      if (!base::nzchar(accession)) return(NULL)
      .protvis_pw_get_uniprot(accession)
    }

    shiny::observeEvent(input$accession, {
      accession <- current_accession_evo()
      rv_evo$homologs <- base::data.frame()
      rv_evo$alignment <- NULL
      rv_evo$tree <- NULL
      if (!base::nzchar(accession)) {
        rv_evo$entry <- NULL
        rv_evo$message <- "Resolve a UniProt protein above before evolutionary analysis."
      } else {
        rv_evo$entry <- base::tryCatch(load_entry(), error = function(e) NULL)
        rv_evo$message <- if (base::is.null(rv_evo$entry)) {
          base::paste("Could not load UniProt record for", accession)
        } else {
          base::paste("Ready:", accession)
        }
      }
    }, ignoreInit = FALSE)

    shiny::observeEvent(input$clear, {
      rv_evo$entry <- NULL
      rv_evo$homologs <- base::data.frame()
      rv_evo$alignment <- NULL
      rv_evo$tree <- NULL
      rv_evo$message <- "Evolutionary analysis cleared."
    })

    shiny::observeEvent(input$run_evolution, {
      base::tryCatch({
        entry <- rv_evo$entry
        if (base::is.null(entry)) entry <- load_entry()
        if (base::is.null(entry)) base::stop("Resolve a valid UniProt protein first.", call. = FALSE)

        selected <- input$evo_species %||% base::character()
        if (base::length(selected) < 2L) base::stop("Select at least two common plant species.", call. = FALSE)

        method <- base::as.character(input$evo_method %||% "NJ")
        processors <- base::as.integer(input$evo_processors %||% 1L)

        shiny::withProgress(message = "Protein Workbench evolutionary analysis", value = 0.05, {
          shiny::setProgress(0.12, detail = "Retrieving plant homolog candidates from UniProt")
          homologs <- .protvis_pw_evolution_fetch(entry, selected)
          rv_evo$homologs <- homologs

          available <- base::sum(homologs$status %in% c("Query", "Retrieved"))
          if (available < 3L) {
            base::stop(
              base::paste0("Only ", available, " usable protein sequences were retrieved. Select more plant species or use a better-annotated UniProt protein."),
              call. = FALSE
            )
          }

          shiny::setProgress(0.55, detail = "Multiple-sequence alignment")
          alignment <- .protvis_pw_evolution_alignment(homologs, processors = processors)
          rv_evo$alignment <- alignment

          shiny::setProgress(0.82, detail = base::paste("Building", method, "phylogenetic tree"))
          tree <- .protvis_pw_evolution_tree(alignment, method = method, processors = processors)
          rv_evo$tree <- tree

          shiny::setProgress(1, detail = "Evolutionary analysis completed")
        })

        rv_evo$message <- base::paste0(
          "Completed · ",
          base::sum(rv_evo$homologs$status %in% c("Query", "Retrieved")),
          " sequences · ",
          base::toupper(method),
          " tree"
        )

        .protvis_record_shared_run(
          shared_state,
          module = "protein_workbench_evolution",
          method = base::paste0("UniProt_candidates_DECIPHER_", base::toupper(method)),
          category = "toolkits",
          parameters = base::list(
            plant_taxa = selected,
            phylogenetic_method = base::toupper(method),
            processors = processors
          ),
          tables = base::list(
            plant_homolog_candidates = rv_evo$homologs
          ),
          statistics = base::list(
            aligned_sequences = base::length(rv_evo$alignment),
            tree_tips = base::length(rv_evo$tree$tip.label)
          )
        )

        shiny::showNotification("Plant evolutionary analysis completed.", type = "message", duration = 4)
      }, error = function(e) {
        rv_evo$message <- base::paste("Evolution:", base::conditionMessage(e))
        shiny::showNotification(rv_evo$message, type = "error", duration = 8)
      })
    })

    output$evo_status <- shiny::renderUI({
      ok <- !base::is.null(rv_evo$tree)
      shiny::div(
        rv_evo$message,
        style = base::paste0(
          "margin-top:10px;font-size:12px;padding:9px 10px;border-radius:10px;",
          if (ok) "background:#edf9f3;color:#286749;" else "background:#f7f9fb;color:#657789;"
        )
      )
    })

    output$evo_table <- DT::renderDT({
      table <- rv_evo$homologs
      if (!base::nrow(table)) {
        table <- base::data.frame(Message = "Run evolutionary analysis to retrieve common-plant homolog candidates.")
      } else {
        table$sequence <- NULL
      }
      DT::datatable(
        table,
        rownames = FALSE,
        filter = if (base::nrow(table) > 1L) "top" else "none",
        options = base::list(pageLength = 15, scrollX = TRUE, autoWidth = TRUE)
      )
    })

    draw_tree <- function() {
      tree <- rv_evo$tree
      if (base::is.null(tree)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Run evolutionary analysis to display the plant phylogenetic tree.", col = "#657789")
        return(invisible(NULL))
      }
      tip_col <- ifelse(base::grepl("^QUERY_", tree$tip.label), "#d73027", "#1f4e79")
      labels <- base::sub("^QUERY_", "", tree$tip.label)
      tree_plot <- tree
      tree_plot$tip.label <- labels
      ape::plot.phylo(
        tree_plot,
        type = "phylogram",
        direction = "rightwards",
        cex = 0.82,
        tip.color = tip_col,
        edge.width = 1.2,
        no.margin = FALSE
      )
      graphics::title(main = "Plant protein phylogeny")
      graphics::mtext("Query protein highlighted in red", side = 1, line = 2.2, cex = 0.75, col = "#657789")
      invisible(NULL)
    }

    output$evo_tree_plot <- shiny::renderPlot({
      draw_tree()
    })

    output$download_evo_table <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_plant_homologs_", base::Sys.Date(), ".csv"),
      content = function(file) utils::write.csv(rv_evo$homologs, file, row.names = FALSE)
    )

    output$download_evo_alignment <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_plant_alignment_", base::Sys.Date(), ".fasta"),
      content = function(file) {
        shiny::req(rv_evo$alignment)
        Biostrings::writeXStringSet(rv_evo$alignment, filepath = file, format = "fasta")
      }
    )

    output$download_evo_tree <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_plant_phylogeny_", base::Sys.Date(), ".nwk"),
      content = function(file) {
        shiny::req(rv_evo$tree)
        ape::write.tree(rv_evo$tree, file = file)
      }
    )

    output$download_evo_png <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_plant_phylogeny_", base::Sys.Date(), ".png"),
      content = function(file) {
        grDevices::png(file, width = 2400, height = 1600, res = 300)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_tree()
      }
    )

    output$download_evo_pdf <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_plant_phylogeny_", base::Sys.Date(), ".pdf"),
      content = function(file) {
        grDevices::pdf(file, width = 10, height = 7)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_tree()
      }
    )

    output$download_evo_svg <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_plant_phylogeny_", base::Sys.Date(), ".svg"),
      content = function(file) {
        grDevices::svg(file, width = 10, height = 7)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_tree()
      }
    )
  })
}
