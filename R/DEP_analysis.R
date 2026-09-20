#' Differential protein analysis UI
#'
#' The default limma workflow reproduces the archived
#' 03.Maize_Teosinte_Jul02_2024 settings: Step6 normalized intensities,
#' B73-minus-Y12 contrasts, historical positive shift, BH adjustment,
#' adjusted-P < 0.05 and |log2FC| > 1.
#'
#' @param id Module namespace.
#' @return Shiny UI.
#' @export
DEP_analysis_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 350,
        shiny::actionButton(
          ns("load_data"), "LOAD DATA",
          class = "btn btn-light fw-bold pv-run-button"
        ),
        shiny::uiOutput(ns("load_status_panel")),
        shinyWidgets::progressBar(
          id = ns("load_progress"),
          value = 0,
          total = 100,
          display_pct = TRUE,
          striped = TRUE,
          status = "success",
          title = "Load progress"
        ),
        shiny::hr(),

        bslib::accordion(
          open = "Analysis parameters",
          bslib::accordion_panel(
            title = "Analysis parameters",
            icon = bsicons::bs_icon("sliders"),
            shiny::tags$small(
              "Defaults reproduce the archived 03.Maize_Teosinte_Jul02_2024 limma workflow.",
              style = "color:#6c757d;"
            ),
            shiny::numericInput(
              ns("dep_logfc"), "|log2FC| threshold",
              value = 1, min = 0, max = 10, step = 0.1
            ),
            shiny::numericInput(
              ns("dep_fdr"), "P-value/FDR threshold",
              value = 0.05, min = 0, max = 1, step = 0.01
            ),
            shiny::selectInput(
              ns("dep_p_metric"), "Significance P-value",
              choices = c(
                "Adjusted P-value / FDR (archived default)" = "adj.P.Val",
                "Raw P-value" = "P.Value"
              ),
              selected = "adj.P.Val"
            ),
            shiny::selectInput(
              ns("dep_adjust_method"), "Multiple-testing adjustment",
              choices = c(
                "BH (archived default)" = "BH",
                "BY" = "BY",
                "Bonferroni" = "bonferroni",
                "Holm" = "holm",
                "None" = "none"
              ),
              selected = "BH"
            ),
            shiny::selectInput(
              ns("dep_sort_by"), "topTable sorting",
              choices = c(
                "logFC (archived default)" = "logFC",
                "P-value" = "P",
                "B-statistic" = "B",
                "None" = "none"
              ),
              selected = "logFC"
            ),
            shiny::checkboxInput(
              ns("dep_matrix_shift"),
              "Apply historical x + abs(min(x)) shift",
              value = TRUE
            ),
            shiny::selectInput(
              ns("dep_protein_universe"), "Protein universe",
              choices = c(
                "Archived default: detected in any of the 6 comparison samples" = "archived_any_detected",
                "Strict: detected in both genotypes" = "both_genotypes",
                "All proteins in Step6 normalized matrix" = "all"
              ),
              selected = "archived_any_detected"
            ),
            shiny::conditionalPanel(
              condition = paste0(
                "input['", ns("dep_protein_universe"),
                "'] == 'both_genotypes'"
              ),
              shiny::numericInput(
                ns("dep_min_detected"),
                "Minimum detected replicates per genotype",
                value = 2, min = 1, max = 3, step = 1
              )
            )
          )
        ),

        shiny::hr(),
        shiny::h5("Demo CompareGroup"),
        shiny::tags$small(
          "The built-in demo uses the five archived B73 vs Y12 developmental comparisons.",
          style = "color:#6c757d;"
        ),
        DT::DTOutput(ns("demo_compare_table")),
        shiny::actionButton(
          ns("use_demo_data"), "USE DEMO DATA",
          class = "btn btn-light fw-bold"
        ),

        shiny::hr(),
        shinyWidgets::switchInput(
          inputId = ns("input_mode"),
          label = "CompareGroup",
          value = TRUE,
          onLabel = "Upload",
          offLabel = "Paste",
          width = "100%"
        ),
        shiny::conditionalPanel(
          condition = paste0("input['", ns("input_mode"), "'] == true"),
          shiny::fileInput(
            ns("compare_file"), NULL,
            accept = c(".csv", ".xlsx", ".xls")
          )
        ),
        shiny::conditionalPanel(
          condition = paste0("input['", ns("input_mode"), "'] == false"),
          rhandsontable::rHandsontableOutput(ns("hot_compare")),
          shiny::textAreaInput(
            ns("paste_data"), NULL,
            placeholder = "Group1\tGroup2", rows = 5
          ),
          shiny::actionButton(
            ns("apply_paste"), "APPLY PASTE",
            class = "btn btn-light fw-bold"
          )
        ),

        shiny::hr(),
        shiny::actionButton(
          ns("run_dep"), "RUN DEP",
          class = "btn btn-primary fw-bold pv-run-button"
        ),
        shiny::br(), shiny::br(),
        shinyWidgets::progressBar(
          id = ns("dep_progress"),
          value = 0,
          total = 100,
          display_pct = TRUE,
          striped = TRUE,
          status = "warning",
          title = "DEP progress"
        )
      ),

      bslib::card(
        height = "760px",
        bslib::card_header("Data Preview"),
        bslib::navset_card_tab(
          full_screen = TRUE,
          bslib::nav_panel(
            "Sample Info",
            shiny::div(
              style = "height:640px;overflow:auto;",
              DT::DTOutput(ns("sample_info"))
            )
          ),
          bslib::nav_panel(
            "Normalized Data",
            shiny::div(
              style = "height:640px;overflow:auto;",
              DT::DTOutput(ns("normalized_data"))
            )
          ),
          bslib::nav_panel(
            "Group Comparison",
            shiny::div(
              style = "height:640px;overflow:auto;",
              DT::DTOutput(ns("group_comparison"))
            )
          ),
          bslib::nav_panel(
            "DEP result",
            shiny::uiOutput(ns("dynamic_dep_tabs"))
          ),
          bslib::nav_panel(
            "DEP summary",
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                width = 245,
                colourpicker::colourInput(
                  ns("summary_down"), "Down", value = "#90EE90"
                ),
                colourpicker::colourInput(
                  ns("summary_ns"), "Not significant", value = "#B3B3B3"
                ),
                colourpicker::colourInput(
                  ns("summary_up"), "Up", value = "#FA8072"
                ),
                shiny::numericInput(
                  ns("summary_width"), "PDF width (inch)",
                  value = 7, min = 4, max = 20
                ),
                shiny::numericInput(
                  ns("summary_height"), "PDF height (inch)",
                  value = 5.5, min = 4, max = 20
                ),
                shiny::downloadButton(
                  ns("download_dep_summary"), "DOWNLOAD SUMMARY PDF"
                ),
                shiny::downloadButton(
                  ns("download_dep_counts"), "DOWNLOAD COUNTS CSV"
                )
              ),
              shiny::div(
                style = paste0(
                  "max-width:760px;height:450px;margin:18px auto 0;",
                  "padding:0 10px;"
                ),
                shiny::plotOutput(
                  ns("dep_summary_plot"),
                  height = "430px",
                  width = "100%"
                )
              )
            )
          )
        )
      )
    )
  )
}

.protvis_dep_archived_defaults <- function() {
  list(
    logfc = 1,
    fdr = 0.05,
    p_metric = "adj.P.Val",
    adjust_method = "BH",
    sort_by = "logFC",
    matrix_shift = TRUE,
    protein_universe = "archived_any_detected",
    min_detected = 2L
  )
}

.protvis_dep_stage_label <- function(group1, group2 = NULL) {
  x <- as.character(group1)
  x <- sub("^(B73|Y12)_", "", x)
  x <- sub("_[123]$", "", x)
  x <- gsub("Leaf_VE\\.V1\\.V2", "Leaf_VE-V2", x)
  x <- gsub("Leaf_V4\\.V6\\.V8", "Leaf_V4-V8", x)
  x
}

.protvis_dep_classify <- function(result, logfc = 1, cutoff = 0.05,
                                  p_metric = "adj.P.Val") {
  result <- as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
  if (!p_metric %in% names(result)) {
    stop("Significance column not found: ", p_metric, call. = FALSE)
  }
  p <- suppressWarnings(as.numeric(result[[p_metric]]))
  lfc <- suppressWarnings(as.numeric(result$logFC))
  result$regulation <- ifelse(
    is.finite(p) & is.finite(lfc) & p < cutoff & lfc > logfc,
    "Upregulated",
    ifelse(
      is.finite(p) & is.finite(lfc) & p < cutoff & lfc < -logfc,
      "Downregulated",
      "Not significant"
    )
  )
  result
}

.protvis_dep_shared_ids <- function(
    pre_knn, group1_samples, group2_samples,
    mode = "archived_any_detected", min_detected = 2L) {
  if (is.null(pre_knn)) return(NULL)
  mat <- as.matrix(pre_knn)
  storage.mode(mat) <- "numeric"
  required <- c(group1_samples, group2_samples)
  if (!all(required %in% colnames(mat))) return(NULL)

  selected <- mat[, required, drop = FALSE]
  mode <- as.character(mode %||% "archived_any_detected")

  if (identical(mode, "archived_any_detected")) {
    # Exact historical 02.depforseedling.R / reproduction logic:
    # mat_step2[, six samples] -> replace NA with 0 -> rowSums(.) > 0.
    # Step4_data_transformed is the log2(mat_step2 * 1e7) snapshot before KNN,
    # so any finite value here is equivalent to a positive mat_step2 value.
    keep <- rowSums(is.finite(selected)) > 0L
    return(rownames(mat)[keep])
  }

  min_detected <- as.integer(min_detected)
  min_detected <- max(1L, min_detected)
  g1 <- rowSums(is.finite(
    mat[, group1_samples, drop = FALSE]
  )) >= min_detected
  g2 <- rowSums(is.finite(
    mat[, group2_samples, drop = FALSE]
  )) >= min_detected
  rownames(mat)[g1 & g2]
}

.protvis_dep_run_limma_archived <- function(
    matrix, group1_samples, group2_samples, group1, group2,
    adjust_method = "BH", sort_by = "logFC", matrix_shift = TRUE) {
  x <- as.matrix(matrix[, c(group1_samples, group2_samples), drop = FALSE])
  storage.mode(x) <- "numeric"

  if (isTRUE(matrix_shift)) {
    minimum <- suppressWarnings(min(x, na.rm = TRUE))
    if (is.finite(minimum)) x <- x + abs(minimum)
  }

  groups <- factor(
    c(
      rep("G1", length(group1_samples)),
      rep("G2", length(group2_samples))
    ),
    levels = c("G1", "G2")
  )
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- c("G1", "G2")
  rownames(design) <- colnames(x)

  fit <- limma::lmFit(x, design)
  contrast <- limma::makeContrasts(G1 - G2, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit)

  out <- limma::topTable(
    fit,
    coef = 1,
    n = Inf,
    adjust.method = adjust_method,
    sort.by = sort_by
  )
  out <- stats::na.omit(out)
  out <- tibble::rownames_to_column(out, "ID")
  out$FC <- 2 ^ out$logFC
  out$Group1 <- group1
  out$Group2 <- group2
  out
}

.protvis_dep_summary_table <- function(results, comparison_order = names(results)) {
  if (!length(results)) return(data.frame())
  rows <- lapply(comparison_order, function(key) {
    df <- results[[key]]
    if (is.null(df)) return(NULL)
    parts <- strsplit(key, "_vs_", fixed = TRUE)[[1L]]
    g1 <- if (length(parts)) parts[[1L]] else key
    stage <- .protvis_dep_stage_label(g1)
    counts <- table(factor(
      df$regulation,
      levels = c("Downregulated", "Not significant", "Upregulated")
    ))
    data.frame(
      Comparison = key,
      Stage = stage,
      Direction = names(counts),
      Protein_number = as.integer(counts),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}

.protvis_dep_summary_plot <- function(counts, up = "#FA8072",
                                      ns = "#B3B3B3",
                                      down = "#90EE90") {
  if (is.null(counts) || !nrow(counts)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text", x = 0, y = 0,
          label = "Run DEP to display the summary."
        )
    )
  }

  stage_order <- unique(counts$Stage)
  counts$Stage <- factor(counts$Stage, levels = rev(stage_order))
  counts$Direction <- factor(
    counts$Direction,
    levels = c("Upregulated", "Not significant", "Downregulated"),
    labels = c("up", "not significant", "down")
  )

  title_expr <- expression(
    italic("Zea mays ssp. mays") ~ "vs" ~ italic("Zea mays ssp. mexicana")
  )

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = Stage, y = Protein_number, fill = Direction)
  ) +
    ggplot2::geom_col(
      colour = "black",
      linewidth = 0.3,
      width = 0.9,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "down" = down,
        "not significant" = ns,
        "up" = up
      ),
      breaks = c("down", "not significant", "up"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title_expr,
      x = NULL,
      y = "Protein number",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 11
      ),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.title = ggplot2::element_text(size = 10, colour = "black"),
      panel.border = ggplot2::element_rect(
        colour = "black", linewidth = 0.8
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 9)
    )
}

.protvis_dep_count_plot <- function(result, title = NULL,
                                    up = "#FA8072",
                                    ns = "#B3B3B3",
                                    down = "#90EE90") {
  counts <- data.frame(
    Direction = factor(
      c("Upregulated", "Not significant", "Downregulated"),
      levels = c("Upregulated", "Not significant", "Downregulated")
    ),
    Protein_number = as.integer(table(factor(
      result$regulation,
      levels = c("Upregulated", "Not significant", "Downregulated")
    )))
  )

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = "", y = Protein_number, fill = Direction)
  ) +
    ggplot2::geom_col(colour = "black", linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(
      "Upregulated" = up,
      "Not significant" = ns,
      "Downregulated" = down
    )) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = NULL,
      y = "Protein number",
      fill = NULL,
      title = title
    )
}

.protvis_dep_volcano_plot <- function(result, params, up, down, ns) {
  p_metric <- params$p_metric
  p <- suppressWarnings(as.numeric(result[[p_metric]]))
  p[p <= 0 | !is.finite(p)] <- .Machine$double.xmin
  plot_df <- result
  plot_df$.plot_p <- p

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = logFC, y = -log10(.plot_p), colour = regulation)
  ) +
    ggplot2::geom_point(alpha = 0.75, size = 1.4) +
    ggplot2::scale_colour_manual(
      values = c(
        "Upregulated" = up,
        "Downregulated" = down,
        "Not significant" = ns
      )
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(params$fdr),
      linetype = "dashed", colour = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-params$logfc, params$logfc),
      linetype = "dashed", colour = "black"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Log2 Fold Change",
      y = if (identical(p_metric, "adj.P.Val")) {
        "-Log10(adjusted P-value)"
      } else {
        "-Log10(P-value)"
      },
      colour = NULL
    ) +
    ggplot2::theme(legend.position = "top")
}

#' Differential protein analysis server
#'
#' @param id Module namespace.
#' @param shared_state Shared ProtVis reactive state.
#' @return Server module.
#' @export
DEP_analysis_server <- function(id, shared_state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    numeric_matrix <- function(x) {
      mat <- as.matrix(x)
      storage.mode(mat) <- "numeric"
      mat
    }

    rv <- shiny::reactiveValues(
      sample_info = NULL,
      normalized_matrix = NULL,
      pre_knn_matrix = NULL,
      compare_data = NULL,
      dep_results = list(),
      dep_summary = data.frame(),
      dep_params = .protvis_dep_archived_defaults(),
      load_success = FALSE,
      dep_ready = FALSE,
      dep_has_run = FALSE,
      volcano_baseline = list()
    )

    reset_dep <- function() {
      rv$dep_results <- list()
      rv$dep_summary <- data.frame()
      rv$dep_ready <- FALSE
      rv$dep_has_run <- FALSE
      rv$volcano_baseline <- list()
      shinyWidgets::updateProgressBar(
        session = session,
        id = "dep_progress",
        value = 0,
        total = 100
      )
    }

    demo_compare_data <- shiny::reactive({
      info <- rv$sample_info
      if (is.null(info) || !"group" %in% names(info)) {
        return(data.frame(
          Group1 = character(),
          Group2 = character(),
          stringsAsFactors = FALSE
        ))
      }
      groups <- unique(as.character(info$group))
      groups <- groups[
        !is.na(groups) & nzchar(groups) & groups != "Unassigned"
      ]

      archived_stages <- c(
        "Root_VE", "Root_V1.V2", "Root_V4",
        "Leaf_VE.V1.V2", "Leaf_V4.V6.V8"
      )
      archived <- data.frame(
        Group1 = paste0("B73_", archived_stages),
        Group2 = paste0("Y12_", archived_stages),
        stringsAsFactors = FALSE
      )
      keep <- archived$Group1 %in% groups & archived$Group2 %in% groups
      if (all(keep)) return(archived)

      b73 <- groups[grepl("^B73_", groups)]
      y12 <- groups[grepl("^Y12_", groups)]
      matched <- b73[paste0("Y12_", sub("^B73_", "", b73)) %in% y12]
      if (length(matched)) {
        return(data.frame(
          Group1 = matched,
          Group2 = paste0("Y12_", sub("^B73_", "", matched)),
          stringsAsFactors = FALSE
        ))
      }
      if (length(groups) >= 2L) {
        return(data.frame(
          Group1 = groups[[1L]],
          Group2 = groups[[2L]],
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        Group1 = character(),
        Group2 = character(),
        stringsAsFactors = FALSE
      )
    })

    output$demo_compare_table <- DT::renderDT({
      DT::datatable(
        demo_compare_data(),
        rownames = FALSE,
        options = list(
          dom = "t", paging = FALSE, searching = FALSE,
          ordering = FALSE, info = FALSE, scrollX = TRUE
        )
      )
    })

    shiny::observeEvent(input$load_data, {
      tryCatch({
        step6 <- NULL
        step4 <- NULL
        workdir <- as.character(shared_state$workdir %||% "")
        if (nzchar(workdir)) {
          p6 <- file.path(workdir, "Step6_data_normalization.rda")
          p4 <- file.path(workdir, "Step4_data_transformed.rda")
          if (file.exists(p6)) {
            step6 <- .protvis_load_stage_dataset(
              p6, expression_names = "normalized_data"
            )
          }
          if (file.exists(p4)) {
            step4 <- .protvis_load_stage_dataset(
              p4, expression_names = "transformed"
            )
          }
        }

        if (is.null(step6) && inherits(shared_state$dataset, "ProtVis_dataset")) {
          step6 <- shared_state$dataset
        }
        if (is.null(step6)) {
          stop("Step6 normalized ProtVis_dataset is unavailable.")
        }

        rv$sample_info <- step6$sample_info
        rv$normalized_matrix <- numeric_matrix(step6$expression_data)
        rv$pre_knn_matrix <- if (!is.null(step4)) {
          numeric_matrix(step4$expression_data)
        } else {
          NULL
        }
        rv$compare_data <- demo_compare_data()
        rv$load_success <- TRUE
        reset_dep()

        shared_state$dataset <- step6
        shinyWidgets::updateProgressBar(
          session = session, id = "load_progress", value = 100, total = 100
        )
        shiny::showNotification(
          if (is.null(rv$pre_knn_matrix)) {
            "✅ Step6 loaded. Step4 was not found; archived protein-universe filtering will fall back to all Step6 proteins."
          } else {
            "✅ Step6 normalized and Step4 pre-KNN matrices loaded; archived detected-protein universe can be reproduced."
          },
          type = "message",
          duration = 6
        )
      }, error = function(e) {
        rv$load_success <- FALSE
        shinyWidgets::updateProgressBar(
          session = session, id = "load_progress", value = 0, total = 100
        )
        shiny::showNotification(
          paste("Unable to load DEP data:", conditionMessage(e)),
          type = "error", duration = NULL
        )
      })
    }, ignoreInit = TRUE)

    output$load_status_panel <- shiny::renderUI({
      if (!isTRUE(rv$load_success)) {
        return(shiny::span("❌ Data not loaded", style = "color:red;"))
      }
      shiny::div(
        shiny::span("✅ Data loaded", style = "color:green;"),
        shiny::br(),
        paste(
          nrow(rv$normalized_matrix), "proteins ×",
          ncol(rv$normalized_matrix), "samples"
        ),
        if (!is.null(rv$pre_knn_matrix)) {
          shiny::tagList(
            shiny::br(),
            shiny::tags$small("Step4 available for the archived six-sample detected-protein universe.")
          )
        }
      )
    })

    shiny::observeEvent(input$use_demo_data, {
      rv$compare_data <- demo_compare_data()
      reset_dep()
    }, ignoreInit = TRUE)

    output$hot_compare <- rhandsontable::renderRHandsontable({
      rhandsontable::rhandsontable(
        rv$compare_data %||% demo_compare_data(),
        stretchH = "all"
      )
    })

    shiny::observeEvent(input$compare_file, {
      shiny::req(input$compare_file)
      ext <- tolower(tools::file_ext(input$compare_file$name))
      df <- if (ext == "csv") {
        utils::read.csv(
          input$compare_file$datapath,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      } else {
        as.data.frame(readxl::read_excel(input$compare_file$datapath))
      }
      rv$compare_data <- df
      reset_dep()
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$apply_paste, {
      shiny::req(nzchar(input$paste_data %||% ""))
      rv$compare_data <- utils::read.table(
        text = input$paste_data,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      reset_dep()
    }, ignoreInit = TRUE)

    output$sample_info <- DT::renderDT({
      shiny::req(rv$sample_info)
      DT::datatable(
        rv$sample_info,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })

    output$normalized_data <- DT::renderDT({
      shiny::req(rv$normalized_matrix)
      df <- data.frame(
        Protein_ID = rownames(rv$normalized_matrix),
        rv$normalized_matrix,
        check.names = FALSE
      )
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })

    output$group_comparison <- DT::renderDT({
      DT::datatable(
        rv$compare_data %||% data.frame(),
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE, scrollX = TRUE)
      )
    })

    shiny::observeEvent(input$run_dep, {
      shiny::req(
        isTRUE(rv$load_success),
        rv$normalized_matrix,
        rv$sample_info,
        rv$compare_data
      )

      comparisons <- as.data.frame(
        rv$compare_data,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      if (!all(c("Group1", "Group2") %in% names(comparisons))) {
        shiny::showNotification(
          "CompareGroup must contain Group1 and Group2.",
          type = "error"
        )
        return(invisible(NULL))
      }
      comparisons <- comparisons[
        stats::complete.cases(comparisons[, c("Group1", "Group2")]),
        ,
        drop = FALSE
      ]
      if (!nrow(comparisons)) {
        shiny::showNotification("No valid comparisons.", type = "error")
        return(invisible(NULL))
      }

      params <- list(
        logfc = as.numeric(input$dep_logfc %||% 1),
        fdr = as.numeric(input$dep_fdr %||% 0.05),
        p_metric = input$dep_p_metric %||% "adj.P.Val",
        adjust_method = input$dep_adjust_method %||% "BH",
        sort_by = input$dep_sort_by %||% "logFC",
        matrix_shift = isTRUE(input$dep_matrix_shift),
        protein_universe = input$dep_protein_universe %||% "archived_any_detected",
        min_detected = as.integer(input$dep_min_detected %||% 2L)
      )
      rv$dep_params <- params
      rv$dep_results <- list()
      rv$volcano_baseline <- list()
      rv$dep_ready <- FALSE
      rv$dep_has_run <- TRUE

      info <- rv$sample_info
      if (!"sample_id" %in% names(info) || !"group" %in% names(info)) {
        shiny::showNotification(
          "sample_info must contain sample_id and group columns.",
          type = "error"
        )
        return(invisible(NULL))
      }

      shiny::withProgress(message = "Running archived-compatible limma", value = 0, {
        for (i in seq_len(nrow(comparisons))) {
          g1 <- as.character(comparisons$Group1[[i]])
          g2 <- as.character(comparisons$Group2[[i]])
          shiny::setProgress(
            (i - 1) / nrow(comparisons),
            detail = paste(g1, "vs", g2)
          )

          s1 <- as.character(info$sample_id[as.character(info$group) == g1])
          s2 <- as.character(info$sample_id[as.character(info$group) == g2])
          s1 <- s1[s1 %in% colnames(rv$normalized_matrix)]
          s2 <- s2[s2 %in% colnames(rv$normalized_matrix)]
          if (!length(s1) || !length(s2)) next

          matrix_use <- rv$normalized_matrix
          effective_universe <- "all"
          if (!identical(params$protein_universe, "all")) {
            ids <- .protvis_dep_shared_ids(
              rv$pre_knn_matrix,
              s1,
              s2,
              mode = params$protein_universe,
              min_detected = params$min_detected
            )
            if (!is.null(ids) && length(ids)) {
              ids <- intersect(ids, rownames(matrix_use))
              matrix_use <- matrix_use[ids, , drop = FALSE]
              effective_universe <- params$protein_universe
            }
          }

          result <- .protvis_dep_run_limma_archived(
            matrix_use,
            s1, s2, g1, g2,
            adjust_method = params$adjust_method,
            sort_by = params$sort_by,
            matrix_shift = params$matrix_shift
          )
          result <- .protvis_dep_classify(
            result,
            logfc = params$logfc,
            cutoff = params$fdr,
            p_metric = params$p_metric
          )
          result$protein_universe <- effective_universe

          key <- paste0(g1, "_vs_", g2)
          rv$dep_results[[key]] <- result
          rv$volcano_baseline[[paste0("show_volcano_", i)]] <-
            input[[paste0("show_volcano_", i)]] %||% 0
        }
        shiny::setProgress(1)
      })

      comparison_order <- paste0(
        comparisons$Group1, "_vs_", comparisons$Group2
      )
      rv$dep_summary <- .protvis_dep_summary_table(
        rv$dep_results, comparison_order
      )
      rv$dep_ready <- length(rv$dep_results) > 0

      if (inherits(shared_state$dataset, "ProtVis_dataset") && rv$dep_ready) {
        dataset <- shared_state$dataset
        dataset$analysis_results$DEP <- list(
          results = rv$dep_results,
          summary = rv$dep_summary,
          parameters = rv$dep_params,
          comparisons = comparisons
        )
        dataset <- .protvis_append_process(
          dataset,
          "differential_analysis",
          status = "success",
          parameters = rv$dep_params
        )
        .protvis_ui_sync_state(dataset, shared_state)
        workdir <- as.character(shared_state$workdir %||% "")
        if (nzchar(workdir) && dir.exists(workdir)) {
          .protvis_save_stage_dataset(
            dataset,
            file.path(workdir, "Step7_differential_analysis.rda")
          )
        }
      }

      shiny::showNotification(
        if (rv$dep_ready) {
          "✅ DEP analysis completed. Volcano plots remain unloaded until SHOW VOLCANO is clicked."
        } else {
          "No comparison produced a valid DEP result."
        },
        type = if (rv$dep_ready) "message" else "warning",
        duration = 6
      )
    }, ignoreInit = TRUE)

    output$dep_summary_plot <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          isTRUE(rv$dep_ready) && nrow(rv$dep_summary) > 0,
          "Run DEP to display the summary."
        )
      )
      print(.protvis_dep_summary_plot(
        rv$dep_summary,
        up = input$summary_up %||% "#FA8072",
        ns = input$summary_ns %||% "#B3B3B3",
        down = input$summary_down %||% "#90EE90"
      ))
    })

    output$download_dep_summary <- shiny::downloadHandler(
      filename = function() paste0("DEP_summary_", Sys.Date(), ".pdf"),
      content = function(file) {
        shiny::req(nrow(rv$dep_summary) > 0)
        plot <- .protvis_dep_summary_plot(
          rv$dep_summary,
          up = input$summary_up %||% "#FA8072",
          ns = input$summary_ns %||% "#B3B3B3",
          down = input$summary_down %||% "#90EE90"
        )
        ggplot2::ggsave(
          file, plot = plot, device = "pdf",
          width = input$summary_width %||% 7,
          height = input$summary_height %||% 5.5,
          units = "in"
        )
      }
    )

    output$download_dep_counts <- shiny::downloadHandler(
      filename = function() paste0("DEP_counts_", Sys.Date(), ".csv"),
      content = function(file) {
        shiny::req(nrow(rv$dep_summary) > 0)
        utils::write.csv(rv$dep_summary, file, row.names = FALSE)
      }
    )

    output$dynamic_dep_tabs <- shiny::renderUI({
      if (!isTRUE(rv$dep_has_run)) {
        return(shiny::div(
          "Click RUN DEP to create comparison results.",
          style = "padding:20px;color:#6c757d;"
        ))
      }
      if (!isTRUE(rv$dep_ready)) {
        return(shiny::div(
          "No valid DEP result is available.",
          style = "padding:20px;color:#6c757d;"
        ))
      }

      comparisons <- rv$compare_data
      tabs <- lapply(seq_len(nrow(comparisons)), function(i) {
        g1 <- as.character(comparisons$Group1[[i]])
        g2 <- as.character(comparisons$Group2[[i]])
        key <- paste0(g1, "_vs_", g2)
        if (is.null(rv$dep_results[[key]])) return(NULL)

        table_id <- paste0("dep_table_", i)
        volcano_id <- paste0("volcano_plot_", i)
        heatmap_id <- paste0("heatmap_", i)
        bar_id <- paste0("bar_dep_", i)
        show_id <- paste0("show_volcano_", i)
        volcano_download_id <- paste0("download_volcano_", i)
        heatmap_download_id <- paste0("download_heatmap_", i)
        bar_download_id <- paste0("download_bar_", i)

        output[[table_id]] <- DT::renderDT({
          DT::datatable(
            rv$dep_results[[key]],
            rownames = FALSE,
            extensions = "Buttons",
            options = list(
              scrollX = TRUE,
              pageLength = 10,
              dom = "Bfrtip",
              buttons = c("copy", "csv", "excel")
            )
          )
        })

        output[[volcano_id]] <- shiny::renderPlot({
          baseline <- rv$volcano_baseline[[show_id]] %||% 0
          shiny::validate(
            shiny::need(
              (input[[show_id]] %||% 0) > baseline,
              "Volcano plot is not rendered by default because large result sets can be slow. Click SHOW VOLCANO."
            )
          )
          print(.protvis_dep_volcano_plot(
            rv$dep_results[[key]],
            rv$dep_params,
            up = input[[paste0("color_up_", i)]] %||% "#d62728",
            down = input[[paste0("color_down_", i)]] %||% "#1f77b4",
            ns = input[[paste0("color_ns_", i)]] %||% "#7f7f7f"
          ))
        })

        output[[volcano_download_id]] <- shiny::downloadHandler(
          filename = function() paste0("Volcano_", key, ".pdf"),
          content = function(file) {
            plot <- .protvis_dep_volcano_plot(
              rv$dep_results[[key]],
              rv$dep_params,
              up = input[[paste0("color_up_", i)]] %||% "#d62728",
              down = input[[paste0("color_down_", i)]] %||% "#1f77b4",
              ns = input[[paste0("color_ns_", i)]] %||% "#7f7f7f"
            )
            ggplot2::ggsave(
              file, plot = plot, device = "pdf",
              width = input[[paste0("volcano_width_", i)]] %||% 8,
              height = input[[paste0("volcano_height_", i)]] %||% 6,
              units = "in"
            )
          }
        )

        output[[heatmap_id]] <- shiny::renderPlot({
          result <- rv$dep_results[[key]]
          sig <- result[result$regulation != "Not significant", , drop = FALSE]
          shiny::validate(
            shiny::need(nrow(sig) > 1L, "No significant proteins for the heatmap.")
          )
          metric <- rv$dep_params$p_metric
          ord <- order(sig[[metric]], -abs(sig$logFC), na.last = NA)
          top_n <- as.integer(input[[paste0("heatmap_top_", i)]] %||% 100L)
          sig <- sig[head(ord, top_n), , drop = FALSE]
          samples <- c(
            as.character(rv$sample_info$sample_id[
              as.character(rv$sample_info$group) == g1
            ]),
            as.character(rv$sample_info$sample_id[
              as.character(rv$sample_info$group) == g2
            ])
          )
          samples <- intersect(samples, colnames(rv$normalized_matrix))
          mat <- rv$normalized_matrix[
            intersect(sig$ID, rownames(rv$normalized_matrix)),
            samples,
            drop = FALSE
          ]
          mat <- t(scale(t(mat)))
          mat[!is.finite(mat)] <- 0
          pheatmap::pheatmap(
            mat,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_colnames = isTRUE(
              input[[paste0("heatmap_show_colnames_", i)]]
            ),
            show_rownames = FALSE,
            main = paste(g1, "vs", g2)
          )
        })

        output[[heatmap_download_id]] <- shiny::downloadHandler(
          filename = function() paste0("Heatmap_", key, ".pdf"),
          content = function(file) {
            grDevices::pdf(
              file,
              width = input[[paste0("heatmap_width_", i)]] %||% 8,
              height = input[[paste0("heatmap_height_", i)]] %||% 6
            )
            result <- rv$dep_results[[key]]
            sig <- result[result$regulation != "Not significant", , drop = FALSE]
            metric <- rv$dep_params$p_metric
            ord <- order(sig[[metric]], -abs(sig$logFC), na.last = NA)
            top_n <- as.integer(input[[paste0("heatmap_top_", i)]] %||% 100L)
            sig <- sig[head(ord, top_n), , drop = FALSE]
            samples <- c(
              as.character(rv$sample_info$sample_id[
                as.character(rv$sample_info$group) == g1
              ]),
              as.character(rv$sample_info$sample_id[
                as.character(rv$sample_info$group) == g2
              ])
            )
            samples <- intersect(samples, colnames(rv$normalized_matrix))
            mat <- rv$normalized_matrix[
              intersect(sig$ID, rownames(rv$normalized_matrix)),
              samples,
              drop = FALSE
            ]
            mat <- t(scale(t(mat)))
            mat[!is.finite(mat)] <- 0
            pheatmap::pheatmap(
              mat,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              show_colnames = isTRUE(
                input[[paste0("heatmap_show_colnames_", i)]]
              ),
              show_rownames = FALSE,
              main = paste(g1, "vs", g2)
            )
            grDevices::dev.off()
          }
        )

        output[[bar_id]] <- shiny::renderPlot({
          print(.protvis_dep_count_plot(
            rv$dep_results[[key]],
            title = .protvis_dep_stage_label(g1),
            up = input[[paste0("bar_up_", i)]] %||% "#FA8072",
            ns = input[[paste0("bar_ns_", i)]] %||% "#B3B3B3",
            down = input[[paste0("bar_down_", i)]] %||% "#90EE90"
          ))
        })

        output[[bar_download_id]] <- shiny::downloadHandler(
          filename = function() paste0("DEP_count_", key, ".pdf"),
          content = function(file) {
            plot <- .protvis_dep_count_plot(
              rv$dep_results[[key]],
              title = .protvis_dep_stage_label(g1),
              up = input[[paste0("bar_up_", i)]] %||% "#FA8072",
              ns = input[[paste0("bar_ns_", i)]] %||% "#B3B3B3",
              down = input[[paste0("bar_down_", i)]] %||% "#90EE90"
            )
            ggplot2::ggsave(
              file, plot = plot, device = "pdf",
              width = input[[paste0("bar_width_", i)]] %||% 7,
              height = input[[paste0("bar_height_", i)]] %||% 4.5,
              units = "in"
            )
          }
        )

        bslib::nav_panel(
          paste(g1, "vs", g2),
          bslib::layout_column_wrap(
            width = 1/2,
            gap = "1rem",
            bslib::card(
              height = "560px",
              bslib::card_header(paste("DEP table -", g1, "vs", g2)),
              bslib::card_body(DT::DTOutput(ns(table_id)))
            ),
            bslib::card(
              height = "560px",
              bslib::card_header(paste("Volcano plot -", g1, "vs", g2)),
              bslib::card_body(
                bslib::layout_sidebar(
                  sidebar = bslib::sidebar(
                    width = 240,
                    shiny::actionButton(
                      ns(show_id), "SHOW VOLCANO",
                      class = "btn btn-primary fw-bold w-100"
                    ),
                    shiny::tags$small(
                      "The volcano is intentionally lazy-rendered for large DEP result sets.",
                      style = "color:#6c757d;"
                    ),
                    colourpicker::colourInput(
                      ns(paste0("color_up_", i)),
                      "Upregulated", "#d62728"
                    ),
                    colourpicker::colourInput(
                      ns(paste0("color_down_", i)),
                      "Downregulated", "#1f77b4"
                    ),
                    colourpicker::colourInput(
                      ns(paste0("color_ns_", i)),
                      "Not significant", "#7f7f7f"
                    ),
                    shiny::numericInput(
                      ns(paste0("volcano_width_", i)),
                      "PDF width", 8, min = 4, max = 20
                    ),
                    shiny::numericInput(
                      ns(paste0("volcano_height_", i)),
                      "PDF height", 6, min = 4, max = 20
                    ),
                    shiny::downloadButton(
                      ns(volcano_download_id), "DOWNLOAD VOLCANO"
                    )
                  ),
                  shiny::plotOutput(ns(volcano_id), height = "430px")
                )
              )
            ),
            bslib::card(
              height = "560px",
              bslib::card_header(paste("Heatmap -", g1, "vs", g2)),
              bslib::card_body(
                bslib::layout_sidebar(
                  sidebar = bslib::sidebar(
                    width = 240,
                    shiny::numericInput(
                      ns(paste0("heatmap_top_", i)),
                      "Top significant proteins", 100,
                      min = 10, max = 500, step = 10
                    ),
                    shiny::checkboxInput(
                      ns(paste0("heatmap_show_colnames_", i)),
                      "Show sample names", TRUE
                    ),
                    shiny::numericInput(
                      ns(paste0("heatmap_width_", i)),
                      "PDF width", 8, min = 4, max = 20
                    ),
                    shiny::numericInput(
                      ns(paste0("heatmap_height_", i)),
                      "PDF height", 6, min = 4, max = 20
                    ),
                    shiny::downloadButton(
                      ns(heatmap_download_id), "DOWNLOAD HEATMAP"
                    )
                  ),
                  shiny::plotOutput(ns(heatmap_id), height = "430px")
                )
              )
            ),
            bslib::card(
              height = "560px",
              bslib::card_header(paste("DEP count -", g1, "vs", g2)),
              bslib::card_body(
                bslib::layout_sidebar(
                  sidebar = bslib::sidebar(
                    width = 240,
                    colourpicker::colourInput(
                      ns(paste0("bar_up_", i)), "Up", "#FA8072"
                    ),
                    colourpicker::colourInput(
                      ns(paste0("bar_ns_", i)),
                      "Not significant", "#B3B3B3"
                    ),
                    colourpicker::colourInput(
                      ns(paste0("bar_down_", i)), "Down", "#90EE90"
                    ),
                    shiny::numericInput(
                      ns(paste0("bar_width_", i)),
                      "PDF width", 7, min = 4, max = 20
                    ),
                    shiny::numericInput(
                      ns(paste0("bar_height_", i)),
                      "PDF height", 4.5, min = 3, max = 20
                    ),
                    shiny::downloadButton(
                      ns(bar_download_id), "DOWNLOAD DEP COUNT"
                    )
                  ),
                  shiny::plotOutput(ns(bar_id), height = "430px")
                )
              )
            )
          )
        )
      })

      tabs <- Filter(Negate(is.null), tabs)
      do.call(
        bslib::navset_card_tab,
        c(list(full_screen = TRUE), tabs)
      )
    })
  })
}
