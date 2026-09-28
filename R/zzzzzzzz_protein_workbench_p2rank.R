# Protein Workbench P2Rank binding-pocket backend -----------------------
#
# P2Rank is the default binding-pocket predictor because its released binary
# distribution runs directly on Windows, Linux and macOS (Java 17+). The
# previous fpocket implementation remains available as an optional backend.

.protvis_pw_p2rank_clean_path <- function(path = "") {
  path <- base::trimws(base::as.character(path %||% ""))
  if (!base::nzchar(path)) return("")
  path <- base::gsub("\\\\", "/", path)
  path <- base::sub("^\"", "", path)
  path <- base::sub("\"$", "", path)
  path <- base::sub("^'", "", path)
  path <- base::sub("'$", "", path)
  path
}

.protvis_pw_p2rank_executable <- function(path = "") {
  path <- .protvis_pw_p2rank_clean_path(path)
  candidates <- base::character()

  if (base::nzchar(path)) {
    if (base::dir.exists(path)) {
      candidates <- c(
        base::file.path(path, "prank.bat"),
        base::file.path(path, "prank.cmd"),
        base::file.path(path, "prank"),
        base::file.path(path, "prank.sh")
      )
    } else if (base::file.exists(path)) {
      candidates <- c(candidates, path)
    }
  }

  candidates <- c(
    candidates,
    base::Sys.which("prank.bat"),
    base::Sys.which("prank.cmd"),
    base::Sys.which("prank"),
    base::Sys.which("prank.sh")
  )
  candidates <- base::unique(candidates[base::nzchar(candidates)])
  candidates <- candidates[base::file.exists(candidates)]

  if (!base::length(candidates)) return("")
  base::normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}

.protvis_pw_p2rank_path_diagnostics <- function(path = "") {
  cleaned <- .protvis_pw_p2rank_clean_path(path)
  directory_exists <- base::nzchar(cleaned) && base::dir.exists(cleaned)
  input_file_exists <- base::nzchar(cleaned) && base::file.exists(cleaned)

  candidates <- if (directory_exists) {
    c(
      prank_bat = base::file.path(cleaned, "prank.bat"),
      prank_cmd = base::file.path(cleaned, "prank.cmd"),
      prank = base::file.path(cleaned, "prank"),
      prank_sh = base::file.path(cleaned, "prank.sh")
    )
  } else {
    c(input = cleaned)
  }

  candidate_exists <- base::vapply(
    candidates,
    base::file.exists,
    logical(1)
  )
  executable <- .protvis_pw_p2rank_executable(cleaned)

  base::list(
    input = path,
    cleaned = cleaned,
    directory_exists = directory_exists,
    input_file_exists = input_file_exists,
    candidates = candidates,
    candidate_exists = candidate_exists,
    executable = executable,
    found = base::nzchar(executable)
  )
}

.protvis_pw_java_status <- function() {
  java <- base::Sys.which("java")
  if (!base::nzchar(java)) {
    return(base::list(
      available = FALSE,
      version = NA_integer_,
      executable = "",
      message = "Java was not found. P2Rank requires Java 17 or later."
    ))
  }

  output <- base::tryCatch(
    base::system2(java, "-version", stdout = TRUE, stderr = TRUE, timeout = 10),
    error = function(e) character()
  )
  txt <- base::paste(output, collapse = " ")
  match <- base::regexec('version "?([0-9]+)', txt)
  parts <- base::regmatches(txt, match)[[1L]]
  version <- if (base::length(parts) >= 2L) {
    base::suppressWarnings(base::as.integer(parts[[2L]]))
  } else {
    NA_integer_
  }

  base::list(
    available = TRUE,
    version = version,
    executable = base::normalizePath(java, winslash = "/", mustWork = TRUE),
    message = if (base::is.finite(version) && version >= 17L) {
      base::paste0("Java ", version, " ready")
    } else if (base::is.finite(version)) {
      base::paste0("Java ", version, " detected; P2Rank requires Java 17 or later.")
    } else {
      "Java detected; version could not be parsed."
    }
  )
}

.protvis_pw_p2rank_status <- function(path = "") {
  diagnostics <- .protvis_pw_p2rank_path_diagnostics(path)
  exe <- diagnostics$executable
  java <- .protvis_pw_java_status()

  if (!base::nzchar(exe)) {
    detail <- if (base::nzchar(diagnostics$cleaned)) {
      if (!isTRUE(diagnostics$directory_exists) &&
          !isTRUE(diagnostics$input_file_exists)) {
        base::paste0(" Path not found: ", diagnostics$cleaned)
      } else if (isTRUE(diagnostics$directory_exists)) {
        if (isTRUE(diagnostics$candidate_exists[["prank_bat"]])) {
          " prank.bat exists but could not be resolved."
        } else {
          " Folder exists, but prank.bat/prank was not found inside it."
        }
      } else {
        " Input file exists but is not a recognized P2Rank launcher."
      }
    } else {
      ""
    }

    return(base::list(
      ready = FALSE,
      executable = "",
      java = java,
      diagnostics = diagnostics,
      message = base::paste0("P2Rank was not found.", detail)
    ))
  }

  if (!isTRUE(java$available) || (base::is.finite(java$version) && java$version < 17L)) {
    return(base::list(
      ready = FALSE,
      executable = exe,
      java = java,
      message = java$message
    ))
  }

  base::list(
    ready = TRUE,
    executable = exe,
    java = java,
    diagnostics = diagnostics,
    message = base::paste0("P2Rank ready · ", exe, " · ", java$message)
  )
}

.protvis_pw_system2_p2rank <- function(executable, args, timeout = 600) {
  executable <- base::normalizePath(
    executable,
    winslash = "/",
    mustWork = TRUE
  )
  old_wd <- base::getwd()
  on.exit(base::setwd(old_wd), add = TRUE)
  base::setwd(base::dirname(executable))

  command <- base::basename(executable)
  is_batch <- base::grepl("\\.(bat|cmd)$", command, ignore.case = TRUE)
  is_windows <- identical(.Platform$OS.type, "windows")

  if (is_windows && is_batch) {
    comspec <- base::Sys.getenv("COMSPEC", unset = "cmd.exe")
    return(base::system2(
      comspec,
      args = c("/d", "/s", "/c", base::shQuote(command), args),
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout
    ))
  }

  if (!is_windows && !base::grepl("^/", command)) {
    command <- base::paste0("./", command)
  }
  base::system2(
    command,
    args = args,
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout
  )
}

.protvis_pw_p2rank_read_csv <- function(path) {
  if (!base::nzchar(path) || !base::file.exists(path)) return(base::data.frame())
  x <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )
  base::names(x) <- base::trimws(base::names(x))
  x
}

.protvis_pw_p2rank_predictions_file <- function(out_dir) {
  files <- base::list.files(
    out_dir,
    pattern = "_predictions\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (base::length(files)) files[[1L]] else ""
}

.protvis_pw_p2rank_residues_file <- function(out_dir) {
  files <- base::list.files(
    out_dir,
    pattern = "_residues\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (base::length(files)) files[[1L]] else ""
}

.protvis_pw_p2rank_residue_ids <- function(x) {
  x <- base::trimws(base::as.character(x %||% ""))
  if (!base::nzchar(x)) return(base::data.frame())

  ids <- base::strsplit(x, "\\s+")[[1L]]
  rows <- base::lapply(ids, function(id) {
    parts <- base::strsplit(id, "_", fixed = TRUE)[[1L]]
    if (base::length(parts) < 2L) return(NULL)
    chain <- base::paste(parts[-base::length(parts)], collapse = "_")
    label <- parts[[base::length(parts)]]
    resi <- base::suppressWarnings(base::as.integer(base::sub("[^0-9-].*$", "", label)))
    if (!base::is.finite(resi)) return(NULL)
    base::data.frame(chain = chain, resi = resi, stringsAsFactors = FALSE)
  })
  rows <- rows[!base::vapply(rows, base::is.null, logical(1))]
  if (!base::length(rows)) return(base::data.frame())
  base::unique(base::do.call(base::rbind, rows))
}

.protvis_pw_pocket_pdb_from_residue_ids <- function(input_pdb, residue_ids, output_pdb) {
  residues <- .protvis_pw_p2rank_residue_ids(residue_ids)
  if (!base::nrow(residues)) {
    base::writeLines(character(), output_pdb)
    return(output_pdb)
  }

  lines <- base::readLines(input_pdb, warn = FALSE)
  atom <- lines[base::grepl("^(ATOM  |HETATM)", lines)]
  if (!base::length(atom)) {
    base::writeLines(character(), output_pdb)
    return(output_pdb)
  }

  padded <- base::vapply(
    atom,
    function(x) if (base::nchar(x) < 80L) base::sprintf("%-80s", x) else x,
    character(1)
  )
  chain <- base::trimws(base::substring(padded, 22L, 22L))
  resi <- base::suppressWarnings(base::as.integer(
    base::trimws(base::substring(padded, 23L, 26L))
  ))

  keep <- base::logical(base::length(atom))
  for (i in base::seq_len(base::nrow(residues))) {
    keep <- keep | (
      chain == residues$chain[[i]] &
        resi == residues$resi[[i]]
    )
  }
  base::writeLines(c(atom[keep], "END"), output_pdb)
  output_pdb
}

.protvis_pw_p2rank_normalize <- function(predictions, input_pdb, out_dir) {
  if (!base::nrow(predictions)) return(base::data.frame())

  get_col <- function(name, default = NA) {
    if (name %in% base::names(predictions)) predictions[[name]]
    else base::rep(default, base::nrow(predictions))
  }

  pocket <- base::suppressWarnings(base::as.integer(get_col("rank")))
  if (base::all(base::is.na(pocket))) pocket <- base::seq_len(base::nrow(predictions))

  table <- base::data.frame(
    pocket = pocket,
    name = base::as.character(get_col("name", "")),
    score = base::suppressWarnings(base::as.numeric(get_col("score"))),
    probability = base::suppressWarnings(base::as.numeric(get_col("probability"))),
    center_x = base::suppressWarnings(base::as.numeric(get_col("center_x"))),
    center_y = base::suppressWarnings(base::as.numeric(get_col("center_y"))),
    center_z = base::suppressWarnings(base::as.numeric(get_col("center_z"))),
    sas_points = base::suppressWarnings(base::as.numeric(get_col("sas_points"))),
    surf_atoms = base::suppressWarnings(base::as.numeric(get_col("surf_atoms"))),
    residue_ids = base::as.character(get_col("residue_ids", "")),
    stringsAsFactors = FALSE
  )

  pocket_dir <- base::file.path(out_dir, "protvis_pockets")
  base::dir.create(pocket_dir, recursive = TRUE, showWarnings = FALSE)
  table$pocket_file <- base::vapply(
    base::seq_len(base::nrow(table)),
    function(i) {
      path <- base::file.path(
        pocket_dir,
        base::paste0("pocket", table$pocket[[i]], "_atm.pdb")
      )
      .protvis_pw_pocket_pdb_from_residue_ids(
        input_pdb,
        table$residue_ids[[i]],
        path
      )
    },
    character(1)
  )
  table$residue_count <- base::vapply(
    table$residue_ids,
    function(x) base::nrow(.protvis_pw_p2rank_residue_ids(x)),
    integer(1)
  )

  table <- table[base::order(table$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  base::rownames(table) <- NULL
  table
}

.protvis_pw_run_p2rank <- function(
  pdb_path,
  p2rank_path = "",
  profile = c("default", "alphafold"),
  threads = 1L
) {
  if (!base::file.exists(pdb_path)) {
    base::stop("PDB structure file was not found.", call. = FALSE)
  }

  status <- .protvis_pw_p2rank_status(p2rank_path)
  if (!isTRUE(status$ready)) base::stop(status$message, call. = FALSE)

  profile <- base::match.arg(profile)
  threads <- base::max(1L, base::as.integer(threads %||% 1L))

  work_dir <- base::tempfile("protvis_p2rank_")
  base::dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  local_pdb <- base::file.path(work_dir, "protvis_structure.pdb")
  if (!base::file.copy(pdb_path, local_pdb, overwrite = TRUE)) {
    base::stop("Could not prepare the PDB file for P2Rank.", call. = FALSE)
  }
  out_dir <- base::file.path(work_dir, "p2rank_output")
  base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    "predict",
    "-f", base::shQuote(local_pdb),
    "-o", base::shQuote(out_dir),
    "-threads", base::as.character(threads),
    "-visualizations", "0",
    "-export_pocket_descriptors", "1",
    "-pocket_grid_format", "csv"
  )
  if (identical(profile, "alphafold")) {
    args <- c("predict", "-c", "alphafold", args[-1L])
  }

  output <- .protvis_pw_system2_p2rank(
    status$executable,
    args = args,
    timeout = 900
  )
  exit_status <- base::attr(output, "status") %||% 0L

  predictions_file <- .protvis_pw_p2rank_predictions_file(out_dir)
  residues_file <- .protvis_pw_p2rank_residues_file(out_dir)

  if (!identical(base::as.integer(exit_status), 0L) ||
      !base::nzchar(predictions_file)) {
    base::stop(
      base::paste0(
        "P2Rank failed. ",
        base::paste(utils::tail(output, 12L), collapse = " | ")
      ),
      call. = FALSE
    )
  }

  predictions <- .protvis_pw_p2rank_read_csv(predictions_file)
  if (!base::nrow(predictions)) {
    base::stop(
      "P2Rank completed successfully but did not predict any binding pockets.",
      call. = FALSE
    )
  }

  pockets <- .protvis_pw_p2rank_normalize(predictions, local_pdb, out_dir)

  base::list(
    backend = "p2rank",
    backend_label = base::paste0("P2Rank · ", status$executable),
    executable = status$executable,
    profile = profile,
    work_dir = work_dir,
    out_dir = out_dir,
    predictions_file = predictions_file,
    residues_file = residues_file,
    stdout = output,
    pockets = pockets,
    input_pdb = local_pdb
  )
}

.protvis_pw_binding_backend_ui <- function(ns) {
  shiny::tagList(
    shiny::selectInput(
      ns("pocket_backend"),
      "Prediction backend",
      choices = c(
        "P2Rank (recommended)" = "p2rank",
        "fpocket (advanced / optional)" = "fpocket"
      ),
      selected = "p2rank"
    ),
    shiny::conditionalPanel(
      condition = base::sprintf("input['%s'] === 'p2rank'", ns("pocket_backend")),
      shiny::textInput(
        ns("p2rank_path"),
        "P2Rank folder or executable",
        value = "",
        placeholder = "Auto-detect prank.bat / prank; or select P2Rank folder"
      ),
      shiny::selectInput(
        ns("p2rank_profile"),
        "P2Rank profile",
        choices = c(
          "Default · experimental PDB structures" = "default",
          "AlphaFold / predicted structures" = "alphafold"
        ),
        selected = "default"
      ),
      shiny::numericInput(
        ns("p2rank_threads"),
        "Threads",
        value = 2,
        min = 1,
        max = 16,
        step = 1
      ),
      shiny::div(
        shiny::tags$a(
          href = "https://github.com/rdk/p2rank/releases",
          target = "_blank",
          "Download P2Rank"
        ),
        shiny::span(" · unzip only; Java 17+ required."),
        style = "font-size:11px;color:#657789;margin-bottom:8px;"
      ),
      shiny::uiOutput(ns("p2rank_status"))
    ),
    shiny::conditionalPanel(
      condition = base::sprintf("input['%s'] === 'fpocket'", ns("pocket_backend")),
      shiny::textInput(
        ns("fpocket_path"),
        "fpocket executable (optional)",
        value = "",
        placeholder = "Auto-detect native fpocket; otherwise use Docker"
      ),
      shiny::textInput(
        ns("fpocket_docker_image"),
        "Docker image",
        value = "fpocket/fpocket:latest",
        placeholder = "fpocket/fpocket:latest"
      ),
      shiny::actionButton(
        ns("pull_fpocket_image"),
        "PULL / CHECK FPOCKET IMAGE",
        icon = bsicons::bs_icon("cloud-download"),
        class = "btn-outline-primary"
      ),
      shiny::uiOutput(ns("docker_image_status"))
    )
  )
}

.protvis_pw_pocket_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::tags$style(shiny::HTML(base::paste0(
      "#", ns("binding_pocket_root"), " .btn{",
      "display:inline-flex;align-items:center;justify-content:center;",
      "gap:6px;min-height:38px;line-height:1.2;overflow:visible;",
      "white-space:normal;text-align:center;padding-top:7px;padding-bottom:7px;}",
      "#", ns("binding_pocket_root"), " .btn svg,",
      "#", ns("binding_pocket_root"), " .btn .bi{",
      "display:inline-block!important;flex:0 0 auto;overflow:visible!important;",
      "width:1em!important;height:1em!important;min-width:1em;}",
      "#", ns("binding_pocket_root"), " .download-button{overflow:visible!important;}"
    ))),
    shiny::div(
      id = ns("binding_pocket_root"),
    style = "padding:0 16px 24px;",
    bslib::card(
      bslib::card_header(
        shiny::div(
          shiny::tags$strong("Binding pocket prediction"),
          shiny::div(
            "P2Rank is the default cross-platform predictor; fpocket remains available as an optional advanced backend.",
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
                "Built-in demo · 1HEL lysozyme" = "demo",
                "Upload PDB" = "upload"
              ),
              selected = "alphafold"
            ),
            shiny::div(
              style = "margin-bottom:10px;",
              shiny::actionButton(
                ns("use_pocket_demo"),
                "USE BUILT-IN DEMO",
                icon = bsicons::bs_icon("database"),
                class = "btn-outline-primary"
              ),
              shiny::downloadButton(
                ns("download_pocket_demo"),
                "Download demo PDB",
                class = "btn-sm btn-outline-secondary"
              )
            ),
            shiny::div(
              "Demo: PDB 1HEL · hen egg white lysozyme · X-ray structure (1.70 Å).",
              style = "font-size:11px;color:#657789;margin-bottom:10px;"
            ),
            shiny::conditionalPanel(
              condition = base::sprintf("input['%s'] === 'upload'", ns("pocket_structure_source")),
              shiny::fileInput(ns("pocket_pdb"), "PDB structure", accept = ".pdb")
            ),
            .protvis_pw_binding_backend_ui(ns),
            shiny::actionButton(
              ns("run_pocket"),
              "RUN POCKET PREDICTION",
              icon = bsicons::bs_icon("bullseye"),
              class = "btn-primary"
            ),
            shiny::uiOutput(ns("pocket_status")),
            shiny::hr(),
            shiny::selectInput(
              ns("pocket_selected"),
              "Displayed pocket",
              choices = character()
            ),
            shiny::downloadButton(
              ns("download_pocket_table"),
              "Pocket table",
              class = "btn-sm btn-outline-primary"
            ),
            shiny::downloadButton(
              ns("download_pocket_pdb"),
              "Pocket PDB",
              class = "btn-sm btn-outline-primary"
            ),
            shiny::downloadButton(
              ns("download_pocket_results"),
              "Full pocket results",
              class = "btn-sm btn-outline-primary"
            )
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
  )
}

# Override the previous fpocket-first UI with P2Rank-first binding pockets.
protein_workbench_ui <- function(id) {
  shiny::tagList(
    .protvis_pw_pre_pocket_ui(id),
    .protvis_pw_pocket_ui(id)
  )
}

protein_workbench_server <- function(id, shared_state = NULL) {
  .protvis_pw_pre_pocket_server(id, shared_state = shared_state)

  shiny::moduleServer(id, function(input, output, session) {
    rv_pocket <- shiny::reactiveValues(
      pdb_text = "",
      pdb_path = "",
      result = NULL,
      message = "Choose a structure and run P2Rank binding-pocket prediction."
    )

    pocket_demo_path <- .protvis_pw_pocket_demo_path()

    current_accession_pocket <- shiny::reactive({
      base::trimws(base::as.character(input$accession %||% ""))
    })

    load_alphafold_pdb <- function() {
      accession <- current_accession_pocket()
      if (!base::nzchar(accession)) {
        base::stop("Resolve a UniProt accession first.", call. = FALSE)
      }
      af <- .protvis_pw_alpha_record(.protvis_pw_alphafold(accession))
      pdb_url <- af$pdbUrl %||% af$pdb_url %||% ""
      if (!base::nzchar(base::as.character(pdb_url))) {
        base::stop(
          "No AlphaFold PDB structure is available for the current protein.",
          call. = FALSE
        )
      }
      pdb_text <- .protvis_pw_http_text(
        base::as.character(pdb_url),
        timeout = 90,
        not_found = NULL
      )
      if (base::is.null(pdb_text) || !base::nzchar(pdb_text)) {
        base::stop("Could not download the AlphaFold PDB structure.", call. = FALSE)
      }
      path <- base::tempfile(fileext = ".pdb")
      base::writeLines(pdb_text, path, useBytes = TRUE)
      base::list(text = pdb_text, path = path)
    }

    selected_pocket_row <- shiny::reactive({
      result <- rv_pocket$result
      if (base::is.null(result) || !base::nrow(result$pockets)) return(NULL)
      selected <- base::suppressWarnings(
        base::as.integer(input$pocket_selected %||% NA_integer_)
      )
      table <- result$pockets
      if (!base::is.finite(selected)) return(table[1L, , drop = FALSE])
      hit <- table[table$pocket == selected, , drop = FALSE]
      if (base::nrow(hit)) hit[1L, , drop = FALSE] else table[1L, , drop = FALSE]
    })

    shiny::observeEvent(input$use_pocket_demo, {
      if (!base::nzchar(pocket_demo_path) || !base::file.exists(pocket_demo_path)) {
        shiny::showNotification(
          "Built-in 1HEL PDB was not found. Reinstall ProtVisDatabase.",
          type = "error",
          duration = 6
        )
        return(NULL)
      }
      shiny::updateRadioButtons(
        session,
        "pocket_structure_source",
        selected = "demo"
      )
      shiny::updateSelectInput(
        session,
        "pocket_backend",
        selected = "p2rank"
      )
      shiny::updateSelectInput(
        session,
        "p2rank_profile",
        selected = "default"
      )
      rv_pocket$message <- "Built-in 1HEL demo selected. Running P2Rank..."
      shinyjs::click(session$ns("run_pocket"))
    })

    output$download_pocket_demo <- shiny::downloadHandler(
      filename = function() "ProtVis_binding_pocket_demo_1HEL.pdb",
      content = function(file) {
        if (!base::nzchar(pocket_demo_path) ||
            !base::file.exists(pocket_demo_path)) {
          base::stop("Built-in 1HEL PDB was not found.", call. = FALSE)
        }
        base::file.copy(pocket_demo_path, file, overwrite = TRUE)
      },
      contentType = "chemical/x-pdb"
    )

    output$p2rank_status <- shiny::renderUI({
      status <- .protvis_pw_p2rank_status(input$p2rank_path %||% "")
      diagnostics <- status$diagnostics %||%
        .protvis_pw_p2rank_path_diagnostics(input$p2rank_path %||% "")

      details <- if (base::nzchar(diagnostics$cleaned)) {
        shiny::div(
          shiny::div(
            base::paste0(
              "Path: ", diagnostics$cleaned,
              " · folder: ",
              if (isTRUE(diagnostics$directory_exists)) "YES" else "NO"
            )
          ),
          if (isTRUE(diagnostics$directory_exists)) {
            shiny::div(
              base::paste0(
                "prank.bat: ",
                if (isTRUE(diagnostics$candidate_exists[["prank_bat"]])) {
                  "YES"
                } else {
                  "NO"
                }
              )
            )
          },
          style = "font-size:10px;margin-top:4px;opacity:0.9;"
        )
      } else {
        NULL
      }

      shiny::div(
        shiny::div(status$message),
        details,
        style = base::paste0(
          "font-size:11px;margin:6px 0 10px;padding:7px 9px;border-radius:8px;",
          if (isTRUE(status$ready)) {
            "background:#edf9f3;color:#286749;"
          } else {
            "background:#fff7e6;color:#8a6116;"
          }
        )
      )
    })

    output$docker_image_status <- shiny::renderUI({
      backend <- .protvis_pw_fpocket_backend(input$fpocket_path %||% "")
      if (!identical(backend$type, "docker")) return(NULL)
      status <- .protvis_pw_docker_image_status(
        input$fpocket_docker_image %||% "fpocket/fpocket:latest"
      )
      shiny::div(
        status$message,
        style = base::paste0(
          "font-size:11px;margin:6px 0 10px;padding:7px 9px;border-radius:8px;",
          if (isTRUE(status$present)) {
            "background:#edf9f3;color:#286749;"
          } else {
            "background:#fff7e6;color:#8a6116;"
          }
        )
      )
    })

    shiny::observeEvent(input$pull_fpocket_image, {
      base::tryCatch({
        image <- input$fpocket_docker_image %||% "fpocket/fpocket:latest"
        status <- .protvis_pw_docker_image_status(image)
        if (!isTRUE(status$present)) .protvis_pw_docker_pull_image(image)
        shiny::showNotification(
          base::paste("Docker image ready:", image),
          type = "message",
          duration = 4
        )
      }, error = function(e) {
        shiny::showNotification(
          base::conditionMessage(e),
          type = "error",
          duration = 12
        )
      })
    })

    shiny::observeEvent(input$run_pocket, {
      base::tryCatch({
        source <- input$pocket_structure_source %||% "alphafold"
        backend <- input$pocket_backend %||% "p2rank"

        shiny::withProgress(
          message = "Binding pocket prediction",
          value = 0.05,
          {
            if (identical(source, "upload")) {
              shiny::req(input$pocket_pdb)
              if (!base::grepl("\\.pdb$", input$pocket_pdb$name, ignore.case = TRUE)) {
                base::stop("Upload a .pdb structure file.", call. = FALSE)
              }
              rv_pocket$pdb_path <- input$pocket_pdb$datapath
              rv_pocket$pdb_text <- base::paste(
                base::readLines(rv_pocket$pdb_path, warn = FALSE),
                collapse = "\n"
              )
            } else if (identical(source, "demo")) {
              if (!base::nzchar(pocket_demo_path) ||
                  !base::file.exists(pocket_demo_path)) {
                base::stop(
                  "Built-in 1HEL PDB was not found. Reinstall ProtVisDatabase.",
                  call. = FALSE
                )
              }
              rv_pocket$pdb_path <- pocket_demo_path
              rv_pocket$pdb_text <- base::paste(
                base::readLines(pocket_demo_path, warn = FALSE),
                collapse = "\n"
              )
            } else {
              shiny::setProgress(0.15, detail = "Loading AlphaFold structure")
              structure <- load_alphafold_pdb()
              rv_pocket$pdb_path <- structure$path
              rv_pocket$pdb_text <- structure$text
            }

            shiny::setProgress(
              0.35,
              detail = if (identical(backend, "p2rank")) {
                "Running P2Rank"
              } else {
                "Running fpocket"
              }
            )

            if (identical(backend, "p2rank")) {
              profile <- input$p2rank_profile %||% "default"
              if (identical(source, "alphafold")) profile <- "alphafold"
              result <- .protvis_pw_run_p2rank(
                rv_pocket$pdb_path,
                p2rank_path = input$p2rank_path %||% "",
                profile = profile,
                threads = input$p2rank_threads %||% 2L
              )
            } else {
              result <- .protvis_pw_run_fpocket(
                rv_pocket$pdb_path,
                fpocket_path = input$fpocket_path %||% "",
                docker_image = input$fpocket_docker_image %||%
                  "fpocket/fpocket:latest"
              )
            }

            rv_pocket$result <- result
            shiny::setProgress(0.88, detail = "Preparing pocket table and 3D view")

            probability <- if ("probability" %in% base::names(result$pockets)) {
              result$pockets$probability
            } else {
              base::rep(NA_real_, base::nrow(result$pockets))
            }
            choices <- stats::setNames(
              base::as.character(result$pockets$pocket),
              base::paste0(
                "Pocket ", result$pockets$pocket,
                " · score ", base::format(
                  base::round(result$pockets$score, 3),
                  trim = TRUE
                ),
                ifelse(
                  base::is.finite(probability),
                  base::paste0(
                    " · probability ",
                    base::format(base::round(probability, 3), trim = TRUE)
                  ),
                  ""
                )
              )
            )
            shiny::updateSelectInput(
              session,
              "pocket_selected",
              choices = choices,
              selected = base::as.character(result$pockets$pocket[[1L]])
            )
            shiny::setProgress(1, detail = "Pocket prediction completed")
          }
        )

        rv_pocket$message <- base::paste0(
          "Completed · ",
          base::nrow(rv_pocket$result$pockets),
          " pockets · ",
          rv_pocket$result$backend_label
        )

        clean_table <- rv_pocket$result$pockets
        clean_table$pocket_file <- base::basename(clean_table$pocket_file)
        if ("flex" %in% base::names(clean_table)) clean_table$flex <- NULL

        files <- base::list(input_pdb = rv_pocket$result$input_pdb)
        if (identical(rv_pocket$result$backend, "p2rank")) {
          files$predictions <- rv_pocket$result$predictions_file
          if (base::nzchar(rv_pocket$result$residues_file)) {
            files$residues <- rv_pocket$result$residues_file
          }
        } else if (!base::is.null(rv_pocket$result$info_file)) {
          files$fpocket_info <- rv_pocket$result$info_file
        }

        .protvis_record_shared_run(
          shared_state,
          module = "protein_workbench_binding_pocket",
          method = rv_pocket$result$backend,
          category = "toolkits",
          parameters = base::list(
            structure_source = source,
            pocket_backend = rv_pocket$result$backend,
            backend_label = rv_pocket$result$backend_label,
            p2rank_profile = rv_pocket$result$profile %||% NULL
          ),
          tables = base::list(binding_pockets = clean_table),
          files = files
        )

        shiny::showNotification(
          base::paste(
            rv_pocket$result$backend_label,
            "detected",
            base::nrow(rv_pocket$result$pockets),
            "binding pockets."
          ),
          type = "message",
          duration = 4
        )
      }, error = function(e) {
        rv_pocket$result <- NULL
        rv_pocket$message <- base::paste(
          "Pocket prediction:",
          base::conditionMessage(e)
        )
        shiny::showNotification(
          rv_pocket$message,
          type = "error",
          duration = 10
        )
      })
    })

    output$pocket_status <- shiny::renderUI({
      ready <- !base::is.null(rv_pocket$result)
      shiny::div(
        rv_pocket$message,
        style = base::paste0(
          "margin-top:10px;font-size:12px;padding:9px 10px;border-radius:10px;",
          if (ready) {
            "background:#edf9f3;color:#286749;"
          } else {
            "background:#f7f9fb;color:#657789;"
          }
        )
      )
    })

    output$pocket_table <- DT::renderDT({
      result <- rv_pocket$result
      if (base::is.null(result) || !base::nrow(result$pockets)) {
        return(DT::datatable(
          base::data.frame(
            Message = "Run P2Rank or fpocket to detect structure-based binding pockets."
          ),
          rownames = FALSE,
          options = base::list(dom = "t", paging = FALSE)
        ))
      }

      table <- result$pockets
      if ("flex" %in% base::names(table)) table$flex <- NULL
      table$pocket_file <- base::basename(table$pocket_file)
      DT::datatable(
        table,
        rownames = FALSE,
        filter = "top",
        options = base::list(
          pageLength = 12,
          scrollX = TRUE,
          autoWidth = TRUE
        )
      )
    })

    output$pocket_view <- r3dmol::renderR3dmol({
      row <- selected_pocket_row()
      pocket_file <- if (base::is.null(row)) {
        ""
      } else {
        base::as.character(row$pocket_file[[1L]])
      }
      .protvis_pw_pocket_view(rv_pocket$pdb_text, pocket_file)
    })

    output$download_pocket_table <- shiny::downloadHandler(
      filename = function() {
        backend <- rv_pocket$result$backend %||% "pocket"
        base::paste0("ProtVis_", backend, "_pockets_", base::Sys.Date(), ".csv")
      },
      content = function(file) {
        shiny::req(rv_pocket$result)
        table <- rv_pocket$result$pockets
        if ("flex" %in% base::names(table)) table$flex <- NULL
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
        if (!base::nzchar(src) || !base::file.exists(src)) {
          base::stop("Pocket PDB is unavailable.", call. = FALSE)
        }
        base::file.copy(src, file, overwrite = TRUE)
      },
      contentType = "chemical/x-pdb"
    )

    output$download_pocket_results <- shiny::downloadHandler(
      filename = function() {
        backend <- rv_pocket$result$backend %||% "pocket"
        base::paste0(
          "ProtVis_", backend, "_binding_pocket_results_",
          base::Sys.Date(), ".zip"
        )
      },
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
