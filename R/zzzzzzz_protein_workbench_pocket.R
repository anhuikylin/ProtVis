# Protein Workbench binding-pocket prediction ----------------------------
#
# fpocket-backed structure analysis for the currently resolved protein.
# Uses the Protein Workbench AlphaFold structure by default and also accepts a
# user-supplied PDB. The external fpocket binary is discovered on PATH or can
# be supplied explicitly by the user.

.protvis_pw_pre_pocket_ui <- protein_workbench_ui
.protvis_pw_pre_pocket_server <- protein_workbench_server

.protvis_pw_fpocket_executable <- function(path = "") {
  path <- base::trimws(base::as.character(path %||% ""))
  if (base::nzchar(path) && base::file.exists(path)) {
    return(base::normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  exe <- base::Sys.which("fpocket")
  if (base::nzchar(exe)) {
    return(base::normalizePath(exe, winslash = "/", mustWork = TRUE))
  }
  ""
}

.protvis_pw_fpocket_info_file <- function(out_dir, stem) {
  candidates <- c(
    base::file.path(out_dir, base::paste0(stem, "_info.txt")),
    base::file.path(out_dir, base::paste0(stem, "_out_info.txt"))
  )
  hit <- candidates[base::file.exists(candidates)]
  if (base::length(hit)) return(hit[[1L]])
  found <- base::list.files(out_dir, pattern = "_info\\.txt$", full.names = TRUE)
  if (base::length(found)) found[[1L]] else ""
}

.protvis_pw_fpocket_parse_info <- function(path) {
  if (!base::nzchar(path) || !base::file.exists(path)) return(base::data.frame())
  lines <- base::readLines(path, warn = FALSE)
  rows <- base::list()
  current <- NULL

  flush_current <- function() {
    if (base::is.null(current)) return(NULL)
    keys <- base::names(current)
    values <- base::unlist(current, use.names = FALSE)
    base::data.frame(
      pocket = base::as.integer(current[["Pocket"]] %||% NA_integer_),
      score = base::as.numeric(current[["Score"]] %||% NA_real_),
      druggability_score = base::as.numeric(current[["Druggability Score"]] %||% NA_real_),
      number_of_alpha_spheres = base::as.numeric(current[["Number of Alpha Spheres"]] %||% NA_real_),
      total_sasa = base::as.numeric(current[["Total SASA"]] %||% NA_real_),
      polar_sasa = base::as.numeric(current[["Polar SASA"]] %||% NA_real_),
      apolar_sasa = base::as.numeric(current[["Apolar SASA"]] %||% NA_real_),
      volume = base::as.numeric(current[["Volume"]] %||% NA_real_),
      mean_local_hydrophobic_density = base::as.numeric(current[["Mean local hydrophobic density"]] %||% NA_real_),
      hydrophobicity_score = base::as.numeric(current[["Hydrophobicity score"]] %||% NA_real_),
      polarity_score = base::as.numeric(current[["Polarity score"]] %||% NA_real_),
      charge_score = base::as.numeric(current[["Charge score"]] %||% NA_real_),
      proportion_of_polar_atoms = base::as.numeric(current[["Proportion of polar atoms"]] %||% NA_real_),
      alpha_sphere_density = base::as.numeric(current[["Alpha sphere density"]] %||% NA_real_),
      center_of_mass_alpha_sphere_max_dist = base::as.numeric(current[["Center of mass - Alpha Sphere max dist"]] %||% NA_real_),
      flex = I(base::list(stats::setNames(values, keys))),
      stringsAsFactors = FALSE
    )
  }

  for (line in lines) {
    line <- base::trimws(line)
    if (!base::nzchar(line)) next

    m <- base::regexec("^Pocket\\s+([0-9]+)\\s*:?$", line, ignore.case = TRUE)
    g <- base::regmatches(line, m)[[1L]]
    if (base::length(g)) {
      if (!base::is.null(current)) rows[[base::length(rows) + 1L]] <- flush_current()
      current <- base::list(Pocket = base::as.integer(g[[2L]]))
      next
    }

    if (base::is.null(current) || !base::grepl(":", line, fixed = TRUE)) next
    parts <- base::strsplit(line, ":", fixed = TRUE)[[1L]]
    if (base::length(parts) < 2L) next
    key <- base::trimws(parts[[1L]])
    value <- base::trimws(base::paste(parts[-1L], collapse = ":"))
    numeric_value <- base::suppressWarnings(base::as.numeric(value))
    current[[key]] <- if (base::is.finite(numeric_value)) numeric_value else value
  }
  if (!base::is.null(current)) rows[[base::length(rows) + 1L]] <- flush_current()

  if (!base::length(rows)) return(base::data.frame())
  out <- base::do.call(base::rbind, rows)
  out <- out[base::order(out$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.protvis_pw_fpocket_pocket_file <- function(out_dir, pocket_number) {
  pockets_dir <- base::file.path(out_dir, "pockets")
  candidates <- c(
    base::file.path(pockets_dir, base::paste0("pocket", pocket_number, "_atm.pdb")),
    base::file.path(pockets_dir, base::paste0("pocket", pocket_number, "_vert.pqr"))
  )
  hit <- candidates[base::file.exists(candidates)]
  if (base::length(hit)) hit[[1L]] else ""
}

.protvis_pw_pdb_residues <- function(path) {
  if (!base::nzchar(path) || !base::file.exists(path)) return(base::data.frame())
  lines <- base::readLines(path, warn = FALSE)
  lines <- lines[base::grepl("^(ATOM  |HETATM)", lines)]
  if (!base::length(lines)) return(base::data.frame())

  pad <- function(x, n) base::sprintf(base::paste0("%-", n, "s"), x)
  lines <- base::vapply(lines, function(x) if (base::nchar(x) < 80L) pad(x, 80L) else x, character(1))

  out <- base::data.frame(
    chain = base::trimws(base::substring(lines, 22L, 22L)),
    resi = base::suppressWarnings(base::as.integer(base::trimws(base::substring(lines, 23L, 26L)))),
    residue = base::trimws(base::substring(lines, 18L, 20L)),
    stringsAsFactors = FALSE
  )
  out <- out[base::is.finite(out$resi), , drop = FALSE]
  base::unique(out)
}

.protvis_pw_run_fpocket <- function(pdb_path, fpocket_path = "") {
  if (!base::file.exists(pdb_path)) base::stop("PDB structure file was not found.", call. = FALSE)
  exe <- .protvis_pw_fpocket_executable(fpocket_path)
  if (!base::nzchar(exe)) {
    base::stop(
      "fpocket executable was not found. Install fpocket and add it to PATH, or provide the fpocket executable path.",
      call. = FALSE
    )
  }

  work_dir <- base::tempfile("protvis_fpocket_")
  base::dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  local_pdb <- base::file.path(work_dir, "protvis_structure.pdb")
  if (!base::file.copy(pdb_path, local_pdb, overwrite = TRUE)) {
    base::stop("Could not prepare the PDB file for fpocket.", call. = FALSE)
  }

  output <- base::system2(
    exe,
    args = c("-f", base::shQuote(local_pdb)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- base::attr(output, "status") %||% 0L
  out_dir <- base::file.path(work_dir, "protvis_structure_out")
  info_file <- .protvis_pw_fpocket_info_file(out_dir, "protvis_structure")

  if (!identical(base::as.integer(status), 0L) || !base::dir.exists(out_dir) || !base::nzchar(info_file)) {
    base::stop(
      base::paste0("fpocket failed. ", base::paste(utils::tail(output, 8L), collapse = " | ")),
      call. = FALSE
    )
  }

  pockets <- .protvis_pw_fpocket_parse_info(info_file)
  if (!base::nrow(pockets)) {
    base::stop("fpocket completed but no pocket descriptors could be parsed.", call. = FALSE)
  }

  pockets$pocket_file <- base::vapply(
    pockets$pocket,
    function(i) .protvis_pw_fpocket_pocket_file(out_dir, i),
    character(1)
  )
  pockets$residue_count <- base::vapply(
    pockets$pocket_file,
    function(path) base::nrow(.protvis_pw_pdb_residues(path)),
    integer(1)
  )

  base::list(
    executable = exe,
    work_dir = work_dir,
    out_dir = out_dir,
    info_file = info_file,
    stdout = output,
    pockets = pockets,
    input_pdb = local_pdb
  )
}

.protvis_pw_pocket_view <- function(pdb_text, pocket_file = "") {
  if (base::is.null(pdb_text) || !base::nzchar(pdb_text)) return(r3dmol::r3dmol())
  viewer <- r3dmol::r3dmol() |>
    r3dmol::m_add_model(data = pdb_text, format = "pdb") |>
    r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "spectrum"))

  residues <- .protvis_pw_pdb_residues(pocket_file)
  if (base::nrow(residues)) {
    chains <- base::unique(residues$chain)
    for (chain in chains) {
      chain_resi <- base::sort(base::unique(residues$resi[residues$chain == chain]))
      selection <- if (base::nzchar(chain)) {
        r3dmol::m_sel(chain = chain, resi = chain_resi)
      } else {
        r3dmol::m_sel(resi = chain_resi)
      }
      viewer <- r3dmol::m_add_style(
        viewer,
        sel = selection,
        style = r3dmol::m_style_sphere(scale = 0.45, color = "#d73027")
      )
    }
  }
  viewer |> r3dmol::m_zoom_to()
}

.protvis_pw_pocket_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding:0 16px 24px;",
    bslib::card(
      bslib::card_header(
        shiny::div(
          shiny::tags$strong("Binding pocket prediction"),
          shiny::div(
            "Detect and rank structure-based binding pockets with fpocket; predicted pocket residues are highlighted on the protein structure.",
            style = "font-size:12px;color:#657789;font-weight:400;margin-top:2px;"
          )
        )
      ),
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(4, 8),
          bslib::card(
            bslib::card_header("Pocket settings"),
            shiny::radioButtons(
              ns("pocket_structure_source"),
              "Structure source",
              choices = c(
                "Current AlphaFold structure" = "alphafold",
                "Upload PDB" = "upload"
              ),
              selected = "alphafold"
            ),
            shiny::conditionalPanel(
              condition = base::sprintf("input['%s'] === 'upload'", ns("pocket_structure_source")),
              shiny::fileInput(ns("pocket_pdb"), "PDB structure", accept = ".pdb")
            ),
            shiny::textInput(
              ns("fpocket_path"),
              "fpocket executable (optional)",
              value = "",
              placeholder = "Auto-detect from PATH"
            ),
            shiny::actionButton(
              ns("run_pocket"),
              "RUN POCKET PREDICTION",
              icon = bsicons::bs_icon("bullseye"),
              class = "btn-primary"
            ),
            shiny::uiOutput(ns("pocket_status")),
            shiny::hr(),
            shiny::selectInput(ns("pocket_selected"), "Displayed pocket", choices = character()),
            shiny::downloadButton(ns("download_pocket_table"), "Pocket table", class = "btn-sm btn-outline-primary"),
            shiny::downloadButton(ns("download_pocket_pdb"), "Pocket PDB", class = "btn-sm btn-outline-primary"),
            shiny::downloadButton(ns("download_fpocket_results"), "Full fpocket results", class = "btn-sm btn-outline-primary")
          ),
          bslib::card(
            full_screen = TRUE,
            bslib::card_header("3D pocket view"),
            r3dmol::r3dmolOutput(ns("pocket_view"), height = "560px")
          )
        ),
        shiny::br(),
        bslib::card(
          bslib::card_header("Predicted binding pockets"),
          DT::DTOutput(ns("pocket_table"))
        )
      )
    )
  )
}

#' Protein Workbench UI
#'
#' Extends Protein Workbench with fpocket binding-pocket prediction.
#'
#' @param id Shiny module id.
#' @return Shiny UI.
#' @export
protein_workbench_ui <- function(id) {
  shiny::tagList(
    .protvis_pw_pre_pocket_ui(id),
    .protvis_pw_pocket_ui(id)
  )
}

#' Protein Workbench server
#'
#' Extends Protein Workbench with fpocket execution, descriptor parsing and
#' interactive 3D pocket visualization.
#'
#' @param id Shiny module id.
#' @param shared_state Optional shared ProtVis reactive state.
#' @return Shiny module server.
#' @export
protein_workbench_server <- function(id, shared_state = NULL) {
  .protvis_pw_pre_pocket_server(id, shared_state = shared_state)

  shiny::moduleServer(id, function(input, output, session) {
    rv_pocket <- shiny::reactiveValues(
      pdb_text = "",
      pdb_path = "",
      result = NULL,
      message = "Use the current AlphaFold structure or upload a PDB, then run fpocket."
    )

    current_accession_pocket <- shiny::reactive({
      base::trimws(base::as.character(input$accession %||% ""))
    })

    load_alphafold_pdb <- function() {
      accession <- current_accession_pocket()
      if (!base::nzchar(accession)) base::stop("Resolve a UniProt accession first.", call. = FALSE)
      af <- .protvis_pw_alpha_record(.protvis_pw_alphafold(accession))
      pdb_url <- af$pdbUrl %||% af$pdb_url %||% ""
      if (!base::nzchar(base::as.character(pdb_url))) {
        base::stop("No AlphaFold PDB structure is available for the current protein.", call. = FALSE)
      }
      pdb_text <- .protvis_pw_http_text(base::as.character(pdb_url), timeout = 90, not_found = NULL)
      if (base::is.null(pdb_text) || !base::nzchar(pdb_text)) {
        base::stop("Could not download the AlphaFold PDB structure.", call. = FALSE)
      }
      path <- base::tempfile(fileext = ".pdb")
      base::writeLines(pdb_text, path, useBytes = TRUE)
      base::list(text = pdb_text, path = path)
    }

    selected_pocket_row <- shiny::reactive({
      if (base::is.null(rv_pocket$result) || !base::nrow(rv_pocket$result$pockets)) return(NULL)
      selected <- base::suppressWarnings(base::as.integer(input$pocket_selected %||% NA_integer_))
      table <- rv_pocket$result$pockets
      if (!base::is.finite(selected)) return(table[1L, , drop = FALSE])
      hit <- table[table$pocket == selected, , drop = FALSE]
      if (base::nrow(hit)) hit[1L, , drop = FALSE] else table[1L, , drop = FALSE]
    })

    shiny::observeEvent(input$run_pocket, {
      base::tryCatch({
        source <- input$pocket_structure_source %||% "alphafold"

        shiny::withProgress(message = "Binding pocket prediction", value = 0.05, {
          if (identical(source, "upload")) {
            shiny::req(input$pocket_pdb)
            if (!base::grepl("\\.pdb$", input$pocket_pdb$name, ignore.case = TRUE)) {
              base::stop("Upload a .pdb structure file.", call. = FALSE)
            }
            rv_pocket$pdb_path <- input$pocket_pdb$datapath
            rv_pocket$pdb_text <- base::paste(base::readLines(rv_pocket$pdb_path, warn = FALSE), collapse = "\n")
          } else {
            shiny::setProgress(0.18, detail = "Loading AlphaFold structure")
            structure <- load_alphafold_pdb()
            rv_pocket$pdb_path <- structure$path
            rv_pocket$pdb_text <- structure$text
          }

          shiny::setProgress(0.42, detail = "Running fpocket")
          result <- .protvis_pw_run_fpocket(
            rv_pocket$pdb_path,
            fpocket_path = input$fpocket_path %||% ""
          )
          rv_pocket$result <- result

          shiny::setProgress(0.90, detail = "Preparing pocket descriptors and 3D view")
          choices <- stats::setNames(
            base::as.character(result$pockets$pocket),
            base::paste0(
              "Pocket ", result$pockets$pocket,
              " · score ", base::format(base::round(result$pockets$score, 3), trim = TRUE),
              " · druggability ", base::format(base::round(result$pockets$druggability_score, 3), trim = TRUE)
            )
          )
          shiny::updateSelectInput(
            session,
            "pocket_selected",
            choices = choices,
            selected = base::as.character(result$pockets$pocket[[1L]])
          )

          shiny::setProgress(1, detail = "Pocket prediction completed")
        })

        rv_pocket$message <- base::paste0(
          "Completed · ", base::nrow(rv_pocket$result$pockets),
          " pockets detected · fpocket: ", base::basename(rv_pocket$result$executable)
        )

        clean_table <- rv_pocket$result$pockets
        clean_table$flex <- NULL
        clean_table$pocket_file <- base::basename(clean_table$pocket_file)

        .protvis_record_shared_run(
          shared_state,
          module = "protein_workbench_binding_pocket",
          method = "fpocket",
          category = "toolkits",
          parameters = base::list(
            structure_source = source,
            fpocket_executable = rv_pocket$result$executable
          ),
          tables = base::list(binding_pockets = clean_table),
          files = base::list(
            input_pdb = rv_pocket$result$input_pdb,
            fpocket_info = rv_pocket$result$info_file
          )
        )

        shiny::showNotification(
          base::paste("fpocket detected", base::nrow(rv_pocket$result$pockets), "binding pockets."),
          type = "message",
          duration = 4
        )
      }, error = function(e) {
        rv_pocket$result <- NULL
        rv_pocket$message <- base::paste("Pocket prediction:", base::conditionMessage(e))
        shiny::showNotification(rv_pocket$message, type = "error", duration = 8)
      })
    })

    output$pocket_status <- shiny::renderUI({
      ready <- !base::is.null(rv_pocket$result)
      executable <- .protvis_pw_fpocket_executable(input$fpocket_path %||% "")
      shiny::div(
        shiny::div(rv_pocket$message),
        shiny::div(
          if (base::nzchar(executable)) {
            base::paste0("fpocket detected: ", executable)
          } else {
            "fpocket not detected on PATH. Install fpocket or provide its executable path."
          },
          style = "font-size:11px;margin-top:4px;"
        ),
        style = base::paste0(
          "margin-top:10px;font-size:12px;padding:9px 10px;border-radius:10px;",
          if (ready) "background:#edf9f3;color:#286749;" else "background:#f7f9fb;color:#657789;"
        )
      )
    })

    output$pocket_table <- DT::renderDT({
      if (base::is.null(rv_pocket$result) || !base::nrow(rv_pocket$result$pockets)) {
        return(DT::datatable(
          base::data.frame(Message = "Run fpocket to detect structure-based binding pockets."),
          rownames = FALSE,
          options = base::list(dom = "t", paging = FALSE)
        ))
      }
      table <- rv_pocket$result$pockets
      table$flex <- NULL
      table$pocket_file <- base::basename(table$pocket_file)
      DT::datatable(
        table,
        rownames = FALSE,
        filter = "top",
        options = base::list(pageLength = 12, scrollX = TRUE, autoWidth = TRUE)
      )
    })

    output$pocket_view <- r3dmol::renderR3dmol({
      row <- selected_pocket_row()
      pocket_file <- if (base::is.null(row)) "" else base::as.character(row$pocket_file[[1L]])
      .protvis_pw_pocket_view(rv_pocket$pdb_text, pocket_file)
    })

    output$download_pocket_table <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_fpocket_", base::Sys.Date(), ".csv"),
      content = function(file) {
        shiny::req(rv_pocket$result)
        table <- rv_pocket$result$pockets
        table$flex <- NULL
        utils::write.csv(table, file, row.names = FALSE)
      }
    )

    output$download_pocket_pdb <- shiny::downloadHandler(
      filename = function() {
        row <- selected_pocket_row()
        shiny::req(row)
        base::paste0("ProtVis_pocket_", row$pocket[[1L]], ".pdb")
      },
      content = function(file) {
        row <- selected_pocket_row()
        shiny::req(row)
        src <- base::as.character(row$pocket_file[[1L]])
        if (!base::nzchar(src) || !base::file.exists(src)) base::stop("Pocket PDB is unavailable.")
        base::file.copy(src, file, overwrite = TRUE)
      }
    )

    output$download_fpocket_results <- shiny::downloadHandler(
      filename = function() base::paste0("ProtVis_fpocket_results_", base::Sys.Date(), ".zip"),
      content = function(file) {
        shiny::req(rv_pocket$result)
        old <- base::getwd()
        on.exit(base::setwd(old), add = TRUE)
        base::setwd(rv_pocket$result$out_dir)
        files <- base::list.files(".", recursive = TRUE, all.files = FALSE)
        utils::zip(file, files = files)
      }
    )
  })
}
