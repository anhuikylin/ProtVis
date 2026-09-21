#' Differential protein analysis UI
#'
#' The default workflow is Recommended DEP: observed Step4 log2 intensities,
#' sample-wise median centering, no imputation, detection filtering, and
#' limma empirical-Bayes statistics with BH FDR control. The historical
#' Step6 workflow is retained as an explicit Archived reproduction mode.
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
              "Recommended DEP is the default. Archived reproduction remains available for historical comparability.",
              style = "color:#6c757d;"
            ),
            shiny::selectInput(
              ns("dep_mode"), "DEP mode",
              choices = c(
                "Recommended DEP" = "recommended",
                "Archived reproduction" = "archived"
              ),
              selected = "recommended"
            ),
            shiny::numericInput(
              ns("dep_logfc"), "|log2FC| threshold",
              value = 1, min = 0, max = 10, step = 0.1
            ),
            shiny::numericInput(
              ns("dep_fdr"), "FDR threshold",
              value = 0.05, min = 0, max = 1, step = 0.01
            ),
            shiny::selectInput(
              ns("dep_p_metric"), "Significance criterion",
              choices = c(
                "Adjusted P-value / FDR (recommended)" = "adj.P.Val",
                "Raw P-value" = "P.Value"
              ),
              selected = "adj.P.Val"
            ),
            shiny::selectInput(
              ns("dep_adjust_method"), "Multiple-testing adjustment",
              choices = c(
                "BH (recommended)" = "BH",
                "BY" = "BY",
                "Bonferroni" = "bonferroni",
                "Holm" = "holm",
                "None" = "none"
              ),
              selected = "BH"
            ),
            shiny::selectInput(
              ns("dep_sort_by"), "Result sorting",
              choices = c(
                "P-value (recommended)" = "P",
                "logFC" = "logFC",
                "B-statistic" = "B",
                "None" = "none"
              ),
              selected = "P"
            ),
            shiny::conditionalPanel(
              condition = paste0(
                "input['", ns("dep_mode"), "'] == 'recommended'"
              ),
              shiny::tags$div(
                class = "alert alert-info py-2 px-3",
                shiny::tags$small(
                  "Input: Step4 observed log2 matrix. No KNN imputation, no row-wise positive shift, and no zero-to-one replacement."
                )
              ),
              shiny::selectInput(
                ns("dep_test_method"), "Statistical test",
                choices = c(
                  "limma robust eBayes (recommended)" = "ebayes_robust",
                  "limma treat (minimum effect-size test)" = "treat"
                ),
                selected = "ebayes_robust"
              ),
              shiny::checkboxInput(
                ns("dep_center_samples"),
                "Median-center samples before DEP",
                value = TRUE
              )
            ),
            shiny::conditionalPanel(
              condition = paste0(
                "input['", ns("dep_mode"), "'] == 'archived'"
              ),
              shiny::tags$div(
                class = "alert alert-secondary py-2 px-3",
                shiny::tags$small(
                  "Uses the historical Step6 normalized matrix and archived limma settings."
                )
              ),
              shiny::checkboxInput(
                ns("dep_matrix_shift"),
                "Apply historical x + abs(min(x)) shift",
                value = TRUE
              )
            ),
            shiny::selectInput(
              ns("dep_protein_universe"), "Detection filter",
              choices = c(
                "Detected in >=2 replicates in both groups (recommended)" = "both_genotypes",
                "Detected in >=2 replicates in either group" = "either_genotype",
                "Detected in any comparison sample (archived)" = "archived_any_detected",
                "All proteins" = "all"
              ),
              selected = "both_genotypes"
            ),
            shiny::conditionalPanel(
              condition = paste0(
                "input['", ns("dep_protein_universe"),
                "'] == 'both_genotypes' || input['",
                ns("dep_protein_universe"), "'] == 'either_genotype'"
              ),
              shiny::numericInput(
                ns("dep_min_detected"),
                "Minimum detected replicates per group",
                value = 2, min = 1, max = 10, step = 1
              )
            ),
            shiny::selectInput(
              ns("dep_volcano_p_metric"), "Volcano y-axis",
              choices = c(
                "Raw P-value (recommended display)" = "P.Value",
                "Adjusted P-value / FDR" = "adj.P.Val"
              ),
              selected = "P.Value"
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
            "DEP Input Data",
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
                shiny::tags$hr(),
                shiny::tags$strong("Presence/absence evidence"),
                colourpicker::colourInput(
                  ns("summary_group1_only"),
                  "Group1 only", value = "#E76F51"
                ),
                colourpicker::colourInput(
                  ns("summary_group2_only"),
                  "Group2 only", value = "#2A9D8F"
                ),
                shiny::numericInput(
                  ns("summary_width"), "PDF width (inch)",
                  value = 12, min = 6, max = 24
                ),
                shiny::numericInput(
                  ns("summary_height"), "PDF height (inch)",
                  value = 5.5, min = 4, max = 20
                ),
                shiny::downloadButton(
                  ns("download_dep_summary"), "DOWNLOAD SUMMARY PDF"
                ),
                shiny::downloadButton(
                  ns("download_dep_counts"), "DOWNLOAD QUANTITATIVE CSV"
                ),
                shiny::downloadButton(
                  ns("download_evidence_counts"),
                  "DOWNLOAD PARALLEL EVIDENCE CSV"
                )
              ),
              shiny::div(
                class = "pv-dep-summary-main",
                style = paste0(
                  "width:100%;max-width:none;min-width:0;",
                  "height:520px;margin:8px 0 0 0;",
                  "padding:8px 18px 0 18px;",
                  "box-sizing:border-box;",
                  "align-self:stretch;justify-self:stretch;"
                ),
                shiny::plotOutput(
                  ns("dep_summary_plot"),
                  height = "500px",
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

.protvis_dep_recommended_defaults <- function() {
  list(
    mode = "recommended",
    logfc = 1,
    fdr = 0.05,
    p_metric = "adj.P.Val",
    volcano_p_metric = "P.Value",
    adjust_method = "BH",
    sort_by = "P",
    test_method = "ebayes_robust",
    center_samples = TRUE,
    matrix_shift = FALSE,
    protein_universe = "both_genotypes",
    min_detected = 2L,
    matrix_source = "Step4_data_transformed"
  )
}

.protvis_dep_prepare_recommended_matrix <- function(
    matrix, center_samples = TRUE) {
  x <- as.matrix(matrix)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- NA_real_

  if (!nrow(x) || !ncol(x)) {
    stop("Recommended DEP input matrix is empty.", call. = FALSE)
  }

  if (isTRUE(center_samples)) {
    medians <- apply(x, 2L, stats::median, na.rm = TRUE)
    if (any(!is.finite(medians))) {
      stop(
        "At least one DEP sample contains no finite observed values.",
        call. = FALSE
      )
    }
    x <- sweep(x, 2L, medians, FUN = "-")
  }
  x
}

.protvis_dep_presence_absence <- function(
    observed_matrix, group1_samples, group2_samples, min_detected = 2L) {
  if (is.null(observed_matrix)) return(data.frame())
  x <- as.matrix(observed_matrix)
  storage.mode(x) <- "numeric"
  required <- c(group1_samples, group2_samples)
  if (!all(required %in% colnames(x))) return(data.frame())

  min_detected <- max(1L, as.integer(min_detected))
  d1 <- rowSums(is.finite(x[, group1_samples, drop = FALSE]))
  d2 <- rowSums(is.finite(x[, group2_samples, drop = FALSE]))
  keep <- (d1 >= min_detected & d2 == 0L) |
    (d2 >= min_detected & d1 == 0L)
  if (!any(keep)) return(data.frame())

  row_median <- function(mat) {
    apply(mat, 1L, function(v) {
      v <- v[is.finite(v)]
      if (length(v)) stats::median(v) else NA_real_
    })
  }

  g1_median <- row_median(x[, group1_samples, drop = FALSE])
  g2_median <- row_median(x[, group2_samples, drop = FALSE])

  data.frame(
    ID = rownames(x)[keep],
    Group1_detected = d1[keep],
    Group2_detected = d2[keep],
    Group1_median_observed = g1_median[keep],
    Group2_median_observed = g2_median[keep],
    Detection_difference = d1[keep] - d2[keep],
    Pattern = ifelse(
      d1[keep] >= min_detected & d2[keep] == 0L,
      "Detected in Group1 only",
      "Detected in Group2 only"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.protvis_dep_run_limma_recommended <- function(
    matrix, group1_samples, group2_samples, group1, group2,
    adjust_method = "BH", sort_by = "P",
    test_method = "ebayes_robust", lfc = 1) {
  x <- as.matrix(matrix[, c(group1_samples, group2_samples), drop = FALSE])
  storage.mode(x) <- "numeric"

  groups <- factor(
    c(rep("G1", length(group1_samples)), rep("G2", length(group2_samples))),
    levels = c("G1", "G2")
  )
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- c("G1", "G2")
  rownames(design) <- colnames(x)

  fit <- limma::lmFit(x, design)
  contrast <- limma::makeContrasts(G1 - G2, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)

  fc_threshold_tested <- FALSE
  if (identical(test_method, "treat")) {
    fit <- limma::treat(
      fit,
      lfc = as.numeric(lfc),
      trend = TRUE,
      robust = TRUE
    )
    out <- limma::topTreat(
      fit,
      coef = 1,
      n = Inf,
      adjust.method = adjust_method,
      sort.by = "none"
    )
    fc_threshold_tested <- TRUE
  } else {
    fit <- limma::eBayes(
      fit,
      trend = TRUE,
      robust = TRUE
    )
    out <- limma::topTable(
      fit,
      coef = 1,
      n = Inf,
      adjust.method = adjust_method,
      sort.by = "none"
    )
  }

  if (identical(sort_by, "P") && "P.Value" %in% names(out)) {
    out <- out[order(out$P.Value, na.last = TRUE), , drop = FALSE]
  } else if (identical(sort_by, "logFC") && "logFC" %in% names(out)) {
    out <- out[order(out$logFC, decreasing = TRUE, na.last = TRUE), ,
               drop = FALSE]
  } else if (identical(sort_by, "B") && "B" %in% names(out)) {
    out <- out[order(out$B, decreasing = TRUE, na.last = TRUE), ,
               drop = FALSE]
  }

  out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
  out <- out[stats::complete.cases(out[, intersect(
    c("logFC", "P.Value", "adj.P.Val"), names(out)
  ), drop = FALSE]), , drop = FALSE]
  out <- tibble::rownames_to_column(out, "ID")
  out$FC <- 2 ^ out$logFC
  out$Group1 <- group1
  out$Group2 <- group2
  out$fc_threshold_tested <- fc_threshold_tested
  out
}

.protvis_dep_archived_defaults <- function() {
  list(
    mode = "archived",
    logfc = 1,
    fdr = 0.05,
    p_metric = "adj.P.Val",
    volcano_p_metric = "adj.P.Val",
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
                                  p_metric = "adj.P.Val",
                                  fc_threshold_tested = FALSE) {
  result <- as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
  if (!p_metric %in% names(result)) {
    stop("Significance column not found: ", p_metric, call. = FALSE)
  }
  p <- suppressWarnings(as.numeric(result[[p_metric]]))
  lfc <- suppressWarnings(as.numeric(result$logFC))
  if (isTRUE(fc_threshold_tested)) {
    result$regulation <- ifelse(
      is.finite(p) & is.finite(lfc) & p < cutoff & lfc > 0,
      "Upregulated",
      ifelse(
        is.finite(p) & is.finite(lfc) & p < cutoff & lfc < 0,
        "Downregulated",
        "Not significant"
      )
    )
  } else {
    result$regulation <- ifelse(
      is.finite(p) & is.finite(lfc) & p < cutoff & lfc > logfc,
      "Upregulated",
      ifelse(
        is.finite(p) & is.finite(lfc) & p < cutoff & lfc < -logfc,
        "Downregulated",
        "Not significant"
      )
    )
  }
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
  if (identical(mode, "either_genotype")) {
    return(rownames(mat)[g1 | g2])
  }
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

.protvis_dep_evidence_summary <- function(
    results, presence_absence, comparison_order = names(results)) {
  if (!length(comparison_order)) return(data.frame())

  rows <- lapply(comparison_order, function(key) {
    parts <- strsplit(key, "_vs_", fixed = TRUE)[[1L]]
    g1 <- if (length(parts) >= 1L) parts[[1L]] else key
    g2 <- if (length(parts) >= 2L) parts[[2L]] else NA_character_
    stage <- .protvis_dep_stage_label(g1)

    quantitative <- results[[key]]
    q_counts <- c(
      "Upregulated" = 0L,
      "Downregulated" = 0L
    )
    if (!is.null(quantitative) && nrow(quantitative)) {
      q_tab <- table(factor(
        quantitative$regulation,
        levels = c("Upregulated", "Downregulated")
      ))
      q_counts[names(q_tab)] <- as.integer(q_tab)
    }

    presence <- presence_absence[[key]]
    p_counts <- c(
      "Detected in Group1 only" = 0L,
      "Detected in Group2 only" = 0L
    )
    if (!is.null(presence) && nrow(presence) &&
        "Pattern" %in% names(presence)) {
      p_tab <- table(factor(
        presence$Pattern,
        levels = names(p_counts)
      ))
      p_counts[names(p_tab)] <- as.integer(p_tab)
    }

    rbind(
      data.frame(
        Comparison = key,
        Group1 = g1,
        Group2 = g2,
        Stage = stage,
        Evidence = "Quantitative DEP",
        Direction = names(q_counts),
        Protein_number = as.integer(q_counts),
        stringsAsFactors = FALSE
      ),
      data.frame(
        Comparison = key,
        Group1 = g1,
        Group2 = g2,
        Stage = stage,
        Evidence = "Presence/absence",
        Direction = names(p_counts),
        Protein_number = as.integer(p_counts),
        stringsAsFactors = FALSE
      )
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}

.protvis_dep_presence_summary_plot <- function(
    counts,
    group1_only = "#E76F51",
    group2_only = "#2A9D8F") {
  presence <- counts[counts$Evidence == "Presence/absence", , drop = FALSE]
  if (!nrow(presence)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text", x = 0, y = 0,
          label = "No presence/absence evidence."
        )
    )
  }

  stage_order <- unique(presence$Stage)
  presence$Stage <- factor(
    presence$Stage,
    levels = rev(stage_order)
  )
  presence$Direction <- factor(
    presence$Direction,
    levels = c(
      "Detected in Group1 only",
      "Detected in Group2 only"
    ),
    labels = c("Group1 only", "Group2 only")
  )

  ggplot2::ggplot(
    presence,
    ggplot2::aes(
      x = Stage,
      y = Protein_number,
      fill = Direction
    )
  ) +
    ggplot2::geom_col(
      colour = "black",
      linewidth = 0.3,
      width = 0.72,
      position = ggplot2::position_dodge(width = 0.78)
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "Group1 only" = group1_only,
        "Group2 only" = group2_only
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Evidence 2 · Presence/absence",
      subtitle = "Proteins detected in one group but absent from the other",
      x = NULL,
      y = "Protein number",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 11
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, size = 8.5, colour = "#6c757d"
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

.protvis_dep_parallel_evidence_plot <- function(
    quantitative_counts,
    evidence_counts,
    up = "#FA8072",
    ns = "#B3B3B3",
    down = "#90EE90",
    group1_only = "#E76F51",
    group2_only = "#2A9D8F") {
  quantitative <- .protvis_dep_summary_plot(
    quantitative_counts,
    up = up,
    ns = ns,
    down = down
  ) +
    ggplot2::labs(
      title = "Evidence 1 · Quantitative DEP",
      subtitle = "Reliable quantitative proteins tested by limma"
    ) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, size = 8.5, colour = "#6c757d"
      )
    )

  presence <- .protvis_dep_presence_summary_plot(
    evidence_counts,
    group1_only = group1_only,
    group2_only = group2_only
  )

  patchwork::wrap_plots(
    quantitative,
    presence,
    ncol = 2,
    widths = c(1.08, 0.92)
  ) +
    patchwork::plot_annotation(
      title = "Two parallel differential-protein evidence streams",
      subtitle = paste0(
        "Quantitative abundance differences and presence/absence evidence ",
        "are reported separately rather than forcing missing proteins into ",
        "the same statistical model."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 13
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5, size = 9, colour = "#6c757d"
        )
      )
    )
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
  p_metric <- params$volcano_p_metric %||% params$p_metric
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
    {
      if (identical(p_metric, params$p_metric)) {
        ggplot2::geom_hline(
          yintercept = -log10(params$fdr),
          linetype = "dashed", colour = "black"
        )
      } else {
        NULL
      }
    } +
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
      subtitle = paste0(
        "Point classification: ", params$p_metric, " < ", params$fdr,
        if (identical(params$test_method %||% "", "treat")) {
          paste0("; limma treat tests |log2FC| > ", params$logfc)
        } else {
          paste0("; |log2FC| > ", params$logfc)
        }
      ),
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
      dep_analysis_matrix = NULL,
      compare_data = NULL,
      dep_results = list(),
      presence_absence = list(),
      dep_summary = data.frame(),
      evidence_summary = data.frame(),
      dep_params = .protvis_dep_recommended_defaults(),
      load_success = FALSE,
      dep_ready = FALSE,
      dep_has_run = FALSE,
      volcano_baseline = list()
    )

    reset_dep <- function() {
      rv$dep_results <- list()
      rv$presence_absence <- list()
      rv$dep_analysis_matrix <- NULL
      rv$dep_summary <- data.frame()
      rv$evidence_summary <- data.frame()
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

    shiny::observeEvent(input$dep_mode, {
      if (identical(input$dep_mode, "archived")) {
        defaults <- .protvis_dep_archived_defaults()
        shiny::updateSelectInput(
          session, "dep_sort_by", selected = defaults$sort_by
        )
        shiny::updateSelectInput(
          session, "dep_protein_universe",
          selected = defaults$protein_universe
        )
        shiny::updateSelectInput(
          session, "dep_volcano_p_metric",
          selected = defaults$volcano_p_metric
        )
        shiny::updateNumericInput(
          session, "dep_logfc", value = defaults$logfc
        )
        shiny::updateNumericInput(
          session, "dep_fdr", value = defaults$fdr
        )
      } else {
        defaults <- .protvis_dep_recommended_defaults()
        shiny::updateSelectInput(
          session, "dep_sort_by", selected = defaults$sort_by
        )
        shiny::updateSelectInput(
          session, "dep_protein_universe",
          selected = defaults$protein_universe
        )
        shiny::updateSelectInput(
          session, "dep_volcano_p_metric",
          selected = defaults$volcano_p_metric
        )
        shiny::updateSelectInput(
          session, "dep_test_method", selected = defaults$test_method
        )
        shiny::updateCheckboxInput(
          session, "dep_center_samples", value = defaults$center_samples
        )
        shiny::updateNumericInput(
          session, "dep_logfc", value = defaults$logfc
        )
        shiny::updateNumericInput(
          session, "dep_fdr", value = defaults$fdr
        )
      }
      reset_dep()
    }, ignoreInit = TRUE)

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
            paste0(
              "✅ Step6 loaded. Step4 was not found, so Recommended DEP is ",
              "unavailable until the Transformation step is completed."
            )
          } else {
            paste0(
              "✅ Step4 observed log2 matrix and Step6 normalized matrix loaded. ",
              "Recommended DEP is ready."
            )
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
        shiny::br(),
        shiny::tags$small(
          if (is.null(rv$pre_knn_matrix)) {
            "Recommended DEP input (Step4 observed log2): unavailable"
          } else {
            paste0(
              "Recommended DEP input: Step4 observed log2 (",
              nrow(rv$pre_knn_matrix), " proteins × ",
              ncol(rv$pre_knn_matrix), " samples)"
            )
          }
        )
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
      mode <- input$dep_mode %||% "recommended"
      matrix <- if (identical(mode, "archived")) {
        rv$normalized_matrix
      } else {
        shiny::req(rv$pre_knn_matrix)
        .protvis_dep_prepare_recommended_matrix(
          rv$pre_knn_matrix,
          center_samples = isTRUE(input$dep_center_samples)
        )
      }
      shiny::req(matrix)
      df <- data.frame(
        Protein_ID = rownames(matrix),
        matrix,
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
        rv$sample_info,
        rv$compare_data
      )

      mode <- input$dep_mode %||% "recommended"
      if (identical(mode, "recommended") && is.null(rv$pre_knn_matrix)) {
        shiny::showNotification(
          paste0(
            "Recommended DEP requires Step4_data_transformed.rda ",
            "(observed log2 values before imputation). Run Transformation first."
          ),
          type = "error",
          duration = 8
        )
        return(invisible(NULL))
      }
      if (identical(mode, "archived") && is.null(rv$normalized_matrix)) {
        shiny::showNotification(
          "Archived reproduction requires the Step6 normalized matrix.",
          type = "error",
          duration = 8
        )
        return(invisible(NULL))
      }

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

      if (identical(mode, "recommended")) {
        params <- .protvis_dep_recommended_defaults()
        params$test_method <- input$dep_test_method %||% params$test_method
        params$center_samples <- isTRUE(input$dep_center_samples)
      } else {
        params <- .protvis_dep_archived_defaults()
        params$matrix_shift <- isTRUE(input$dep_matrix_shift)
      }
      params$logfc <- as.numeric(input$dep_logfc %||% params$logfc)
      params$fdr <- as.numeric(input$dep_fdr %||% params$fdr)
      params$p_metric <- input$dep_p_metric %||% params$p_metric
      params$volcano_p_metric <-
        input$dep_volcano_p_metric %||% params$volcano_p_metric
      params$adjust_method <-
        input$dep_adjust_method %||% params$adjust_method
      params$sort_by <- input$dep_sort_by %||% params$sort_by
      params$protein_universe <-
        input$dep_protein_universe %||% params$protein_universe
      params$min_detected <-
        as.integer(input$dep_min_detected %||% params$min_detected)

      rv$dep_params <- params
      rv$dep_results <- list()
      rv$presence_absence <- list()
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

      analysis_matrix <- if (identical(mode, "recommended")) {
        .protvis_dep_prepare_recommended_matrix(
          rv$pre_knn_matrix,
          center_samples = params$center_samples
        )
      } else {
        rv$normalized_matrix
      }
      rv$dep_analysis_matrix <- analysis_matrix

      progress_message <- if (identical(mode, "recommended")) {
        "Running Recommended DEP"
      } else {
        "Running archived-compatible limma"
      }

      shinyWidgets::updateProgressBar(
        session = session,
        id = "dep_progress",
        value = 0,
        total = 100
      )

      shiny::withProgress(message = progress_message, value = 0, {
        for (i in seq_len(nrow(comparisons))) {
          g1 <- as.character(comparisons$Group1[[i]])
          g2 <- as.character(comparisons$Group2[[i]])
          progress_value <- (i - 1) / nrow(comparisons)
          shiny::setProgress(
            progress_value,
            detail = paste(g1, "vs", g2)
          )
          shinyWidgets::updateProgressBar(
            session = session,
            id = "dep_progress",
            value = round(progress_value * 100),
            total = 100
          )

          s1 <- as.character(info$sample_id[as.character(info$group) == g1])
          s2 <- as.character(info$sample_id[as.character(info$group) == g2])
          s1 <- s1[s1 %in% colnames(analysis_matrix)]
          s2 <- s2[s2 %in% colnames(analysis_matrix)]
          if (!length(s1) || !length(s2)) next

          key <- paste0(g1, "_vs_", g2)

          if (identical(mode, "recommended")) {
            presence <- .protvis_dep_presence_absence(
              rv$pre_knn_matrix,
              s1,
              s2,
              min_detected = params$min_detected
            )
            if (nrow(presence)) {
              presence$Group1 <- g1
              presence$Group2 <- g2
            }
            rv$presence_absence[[key]] <- presence
          } else {
            rv$presence_absence[[key]] <- data.frame()
          }

          matrix_use <- analysis_matrix
          effective_universe <- "all"
          if (!identical(params$protein_universe, "all")) {
            detection_matrix <- rv$pre_knn_matrix
            if (!is.null(detection_matrix)) {
              ids <- .protvis_dep_shared_ids(
                detection_matrix,
                s1,
                s2,
                mode = params$protein_universe,
                min_detected = params$min_detected
              )
              if (!is.null(ids)) {
                ids <- intersect(ids, rownames(matrix_use))
                matrix_use <- matrix_use[ids, , drop = FALSE]
                effective_universe <- params$protein_universe
              }
            }
          }

          if (!nrow(matrix_use)) next

          if (identical(mode, "recommended")) {
            result <- .protvis_dep_run_limma_recommended(
              matrix_use,
              s1, s2, g1, g2,
              adjust_method = params$adjust_method,
              sort_by = params$sort_by,
              test_method = params$test_method,
              lfc = params$logfc
            )
            fc_tested <- if (nrow(result)) {
              isTRUE(result$fc_threshold_tested[[1L]])
            } else {
              FALSE
            }
          } else {
            result <- .protvis_dep_run_limma_archived(
              matrix_use,
              s1, s2, g1, g2,
              adjust_method = params$adjust_method,
              sort_by = params$sort_by,
              matrix_shift = params$matrix_shift
            )
            fc_tested <- FALSE
          }

          if (!nrow(result)) next
          result <- .protvis_dep_classify(
            result,
            logfc = params$logfc,
            cutoff = params$fdr,
            p_metric = params$p_metric,
            fc_threshold_tested = fc_tested
          )
          result$analysis_mode <- mode
          result$protein_universe <- effective_universe
          result$matrix_source <- if (identical(mode, "recommended")) {
            "Step4_data_transformed"
          } else {
            "Step6_data_normalization"
          }
          result$test_method <- if (identical(mode, "recommended")) {
            params$test_method
          } else {
            "archived_eBayes"
          }

          rv$dep_results[[key]] <- result
          rv$volcano_baseline[[paste0("show_volcano_", i)]] <-
            input[[paste0("show_volcano_", i)]] %||% 0
        }
        shiny::setProgress(1)
        shinyWidgets::updateProgressBar(
          session = session,
          id = "dep_progress",
          value = 100,
          total = 100
        )
      })

      comparison_order <- paste0(
        comparisons$Group1, "_vs_", comparisons$Group2
      )
      rv$dep_summary <- .protvis_dep_summary_table(
        rv$dep_results, comparison_order
      )
      rv$evidence_summary <- .protvis_dep_evidence_summary(
        rv$dep_results,
        rv$presence_absence,
        comparison_order
      )
      rv$dep_ready <- length(rv$dep_results) > 0

      if (inherits(shared_state$dataset, "ProtVis_dataset") && rv$dep_ready) {
        dataset <- shared_state$dataset
        dataset$analysis_results$DEP <- list(
          results = rv$dep_results,
          presence_absence = rv$presence_absence,
          summary = rv$dep_summary,
          evidence_summary = rv$evidence_summary,
          parameters = rv$dep_params,
          comparisons = comparisons,
          input_matrix = rv$dep_analysis_matrix,
          provenance = list(
            mode = mode,
            matrix_source = if (identical(mode, "recommended")) {
              "Step4_data_transformed.rda"
            } else {
              "Step6_data_normalization.rda"
            },
            imputation_for_statistics = identical(mode, "archived"),
            sample_median_centering = if (identical(mode, "recommended")) {
              params$center_samples
            } else {
              NA
            },
            row_wise_positive_shift = if (identical(mode, "recommended")) {
              FALSE
            } else {
              params$matrix_shift
            },
            zero_to_one_for_statistics = identical(mode, "archived"),
            detection_filter = params$protein_universe,
            minimum_detected_replicates = params$min_detected,
            test_method = if (identical(mode, "recommended")) {
              params$test_method
            } else {
              "archived_eBayes"
            },
            adjust_method = params$adjust_method,
            significance_metric = params$p_metric,
            fdr_threshold = params$fdr,
            logfc_threshold = params$logfc
          )
        )
        dataset <- .protvis_append_process(
          dataset,
          if (identical(mode, "recommended")) {
            "differential_analysis_recommended"
          } else {
            "differential_analysis_archived"
          },
          status = "success",
          parameters = dataset$analysis_results$DEP$provenance
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
          paste0(
            "✅ ",
            if (identical(mode, "recommended")) {
              "Recommended DEP"
            } else {
              "Archived DEP"
            },
            " completed. Volcano plots remain unloaded until SHOW VOLCANO is clicked."
          )
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
      print(.protvis_dep_parallel_evidence_plot(
        rv$dep_summary,
        rv$evidence_summary,
        up = input$summary_up %||% "#FA8072",
        ns = input$summary_ns %||% "#B3B3B3",
        down = input$summary_down %||% "#90EE90",
        group1_only = input$summary_group1_only %||% "#E76F51",
        group2_only = input$summary_group2_only %||% "#2A9D8F"
      ))
    })

    output$download_dep_summary <- shiny::downloadHandler(
      filename = function() paste0("DEP_summary_", Sys.Date(), ".pdf"),
      content = function(file) {
        shiny::req(nrow(rv$dep_summary) > 0)
        plot <- .protvis_dep_parallel_evidence_plot(
          rv$dep_summary,
          rv$evidence_summary,
          up = input$summary_up %||% "#FA8072",
          ns = input$summary_ns %||% "#B3B3B3",
          down = input$summary_down %||% "#90EE90",
          group1_only = input$summary_group1_only %||% "#E76F51",
          group2_only = input$summary_group2_only %||% "#2A9D8F"
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

    output$download_evidence_counts <- shiny::downloadHandler(
      filename = function() {
        paste0("DEP_parallel_evidence_counts_", Sys.Date(), ".csv")
      },
      content = function(file) {
        shiny::req(nrow(rv$evidence_summary) > 0)
        utils::write.csv(
          rv$evidence_summary,
          file,
          row.names = FALSE
        )
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
        presence_id <- paste0("presence_absence_", i)
        show_id <- paste0("show_volcano_", i)
        volcano_download_id <- paste0("download_volcano_", i)
        heatmap_download_id <- paste0("download_heatmap_", i)
        bar_download_id <- paste0("download_bar_", i)
        presence_download_id <- paste0("download_presence_absence_", i)

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

        output[[presence_id]] <- DT::renderDT({
          presence <- rv$presence_absence[[key]] %||% data.frame()
          if (!nrow(presence)) {
            presence <- data.frame(
              Message = if (identical(rv$dep_params$mode, "recommended")) {
                "No presence/absence candidate met the current detection rule."
              } else {
                "Presence/absence candidates are reported only in Recommended DEP."
              },
              stringsAsFactors = FALSE
            )
          }
          DT::datatable(
            presence,
            rownames = FALSE,
            options = list(
              scrollX = TRUE,
              pageLength = 10
            )
          )
        })

        output[[presence_download_id]] <- shiny::downloadHandler(
          filename = function() {
            paste0("Presence_absence_", key, ".csv")
          },
          content = function(file) {
            presence <- rv$presence_absence[[key]] %||% data.frame()
            utils::write.csv(presence, file, row.names = FALSE)
          }
        )

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
          heatmap_source <- rv$dep_analysis_matrix %||% rv$normalized_matrix
          samples <- intersect(samples, colnames(heatmap_source))
          mat <- heatmap_source[
            intersect(sig$ID, rownames(heatmap_source)),
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
            heatmap_source <- rv$dep_analysis_matrix %||% rv$normalized_matrix
            samples <- intersect(samples, colnames(heatmap_source))
            mat <- heatmap_source[
              intersect(sig$ID, rownames(heatmap_source)),
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
              bslib::card_header(paste("Evidence 1 · Quantitative DEP table -", g1, "vs", g2)),
              bslib::card_body(DT::DTOutput(ns(table_id)))
            ),
            bslib::card(
              height = "560px",
              bslib::card_header(paste("Evidence 1 · Quantitative volcano -", g1, "vs", g2)),
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
              bslib::card_header(paste("Evidence 1 · Quantitative heatmap -", g1, "vs", g2)),
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
              bslib::card_header(paste("Evidence 1 · Quantitative count -", g1, "vs", g2)),
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
            ),
            bslib::card(
              height = "560px",
              bslib::card_header(
                paste("Evidence 2 · Presence/absence candidates -", g1, "vs", g2)
              ),
              bslib::card_body(
                shiny::tags$small(
                  paste0(
                    "Recommended DEP does not impute proteins that are absent ",
                    "from one group. Such candidates are reported separately."
                  ),
                  style = "color:#6c757d;"
                ),
                shiny::br(),
                shiny::br(),
                shiny::downloadButton(
                  ns(presence_download_id),
                  "DOWNLOAD CANDIDATES"
                ),
                shiny::br(),
                shiny::br(),
                DT::DTOutput(ns(presence_id))
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
