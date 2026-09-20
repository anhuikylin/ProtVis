#' Overview UI Module
#' Creates the user interface for the overview analysis module
#' @param id Character string specifying the namespace id
#' @return A Shiny UI tagList containing the overview analysis interface
#' @import shiny
#' @import bslib
#' @importFrom colourpicker colourInput
#' @name overview_ui
#' @export
#'
overview_ui <- function(id) {
  ns <- NS(id)
  shiny::tagList(
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 300,
        shiny::actionButton(ns("load_data"), "LOAD DATA", class = "btn btn-light fw-bold"),
        shiny::uiOutput(ns("load_status_panel")),
        bslib::accordion(
          bslib::accordion_panel(
            title = "Correlation",
            icon = correlation_icon,
            shiny::selectInput(
              inputId = ns("cor_method"),
              label = "Correlation Method:",
              choices = c("Pearson", "Spearman", "Kendall"),
              selected = "Pearson"
            ),
            colourpicker::colourInput(
              ns("cor_high_color"),
              "High Color",
              value = "purple"
            ),
            colourpicker::colourInput(
              ns("cor_mid_color"),
              "middle Color",
              value = "black"
            ),
            colourpicker::colourInput(
              ns("cor_low_color"),
              "Low Color",
              value = "yellow"
            ),
            shiny::numericInput(ns("cor_color_min"), "Set Min Value", value = -1, step = 0.1),
            shiny::numericInput(ns("cor_color_max"), "Set Max Value", value = 1, step = 0.1),
            shiny::checkboxInput(ns("cor_cluster_rows"), "Cluster rows", TRUE),
            shiny::checkboxInput(ns("cor_cluster_columns"), "Cluster columns", TRUE),
            shiny::checkboxInput(ns("cor_show_numbers"), "Show correlation values", TRUE),
            shiny::checkboxInput(ns("cor_show_column_names"), "Show sample names", FALSE),
            shiny::actionButton(ns("run_correlation"), "Run Correlation"),
            shiny::numericInput(ns("cor_plot_width"), "Download Plot Width (inches)", value = 10),
            shiny::numericInput(ns("cor_plot_height"), "Download Plot Height (inches)", value = 7),
            shiny::downloadButton(ns("cor_download_pdf"), "Download PDF"),
            shiny::downloadButton(
              ns("cor_download_table"),
              "Download Correlation Table"
            )
          ),
          bslib::accordion_panel(
            title = "Expression pattern",
            icon = expression_pattern_icon,
            shiny::sliderInput(
              inputId = ns("exp_top_n"),
              label = "Top N Features:",
              min = 50,
              max = 2000,
              value = 500,
              step = 50
            ),
            shiny::checkboxInput(
              inputId = ns("exp_scale"),
              label = "Scale Data",
              value = TRUE
            ),
            shiny::selectInput(
              ns("exp_scale_method"),
              "Scaling direction",
              choices = c("By protein (row)" = "row",
                          "By sample (column)" = "column",
                          "No scaling" = "none"),
              selected = "row"
            ),
            shiny::checkboxInput(ns("exp_cluster_rows"), "Cluster samples", TRUE),
            shiny::checkboxInput(ns("exp_cluster_columns"), "Cluster proteins", TRUE),
            shiny::checkboxInput(ns("exp_show_feature_names"), "Show protein names", FALSE),
            colourpicker::colourInput(
              ns("exp_high_color"),
              "High Color",
              value = "purple"
            ),
            colourpicker::colourInput(
              ns("exp_mid_color"),
              "middle Color",
              value = "black"
            ),
            colourpicker::colourInput(
              ns("exp_low_color"),
              "Low Color",
              value = "yellow"
            ),
            shiny::numericInput(ns("exp_color_min"), "Set Min Value", value = -1, step = 0.1),
            shiny::numericInput(ns("exp_color_max"), "Set Max Value", value = 1, step = 0.1),
            shiny::actionButton(ns("run_expression"), "Run Expression"),
            shiny::numericInput(ns("exp_plot_width"), "Download Plot Width (inches)", value = 10),
            shiny::numericInput(ns("exp_plot_height"), "Download Plot Height (inches)", value = 7),
            shiny::downloadButton(ns("exp_download_pdf"), "Download PDF")
          ),
          bslib::accordion_panel(
            title = "Dimensionality Reduction",
            icon = dimensionality_reduction_icon,
            shiny::helpText(
              "Reproduces the archived Figure 3 UMAP: transformed pre-KNN ",
              "intensity matrix, five early developmental groups, seed 10086, ",
              "group colours, species shapes, and four biological-region ellipses."
            ),
            shiny::actionButton(ns("DR_analyse"), "Run UMAP"),
            shiny::numericInput(
              ns("dr_plot_width"),
              "Download Plot Width (inches)",
              value = 10
            ),
            shiny::numericInput(
              ns("dr_plot_height"),
              "Download Plot Height (inches)",
              value = 7
            ),
            shiny::downloadButton(ns("dr_download_pdf"), "Download PDF")
          ),
          bslib::accordion_panel(
            title = "Proteomics QC",
            icon = bsicons::bs_icon("clipboard-pulse"),
            shiny::selectInput(
              ns("qc_download_plot_type"),
              "QC figure to download",
              choices = c(
                "Normalized intensity density" = "density",
                "Protein coefficient of variation" = "cv"
              ),
              selected = "density"
            ),
            shiny::numericInput(
              ns("qc_plot_width"),
              "Download Plot Width (inches)",
              value = 8
            ),
            shiny::numericInput(
              ns("qc_plot_height"),
              "Download Plot Height (inches)",
              value = 6
            ),
            shiny::downloadButton(ns("qc_download_pdf"), "Download QC PDF"),
            shiny::downloadButton(
              ns("qc_download_matrix"),
              "Download Normalized Matrix"
            )
          )

        )
      ),
      bslib::page_fluid(
        bslib::layout_column_wrap(
          width = 1/2,
          gap = "1rem",
          bslib::card(
            height = "520px",
            bslib::card_header("Correlation"),
            bslib::card_body(
              shiny::plotOutput(ns("cor_res"), height = "430px")
            )
          ),
          bslib::card(
            height = "520px",
            bslib::card_header("Expression pattern"),
            bslib::card_body(
              shiny::plotOutput(ns("expression_pattern"), height = "430px")
            )
          ),
          bslib::card(
            height = "520px",
            bslib::card_header("Normalized intensity density"),
            bslib::card_body(
              shiny::plotOutput(ns("qc_density_plot"), height = "430px")
            )
          ),
          bslib::card(
            height = "520px",
            bslib::card_header("Dimensionality reduction analyse"),
            bslib::card_body(
              shiny::plotOutput(ns("DR_Reproduction"), height = "430px")
            )
          )
        ),
        shiny::div(style = "height: 1rem;"),
        bslib::layout_column_wrap(
          width = 1,
          gap = "1rem",
          bslib::card(
            height = "440px",
            bslib::card_header("Protein coefficient of variation"),
            bslib::card_body(
              shiny::plotOutput(ns("qc_cv_plot"), height = "350px")
            )
          )
        )
      )
    )
  )
}

#' Overview Server Module
#' Server-side logic for the overview analysis module
#' @param id Character string specifying the namespace id
#' @param shared_state Reactive values shared across modules
#' @return A module server function that handles the overview analysis logic
#' @import shiny
#' @importFrom dplyr left_join mutate select case_when
#' @importFrom stringr str_split str_remove_all
#' @importFrom ComplexHeatmap Heatmap rowAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar grid.text
#' @importFrom grDevices pdf dev.off
#' @importFrom matrixStats rowVars
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom vegan metaMDS
#' @importFrom ggsci scale_color_lancet scale_fill_lancet
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_bw labs geom_col geom_boxplot geom_density geom_histogram geom_text theme_minimal theme element_text
#' @name overview_server
#' @export
#'
utils::globalVariables(c(
  "tissue", "tissue2", "species", "Type", "Species",
  "V1", "V2", "SampleType", "UMAP1", "UMAP2", "Group", "Region"
))

overview_server <- function(id, shared_state) {
  shiny::moduleServer(id, function(input, output, session) {

    standardize_overview_matrix <- function(data) {
      matrix <- base::as.matrix(data)
      storage.mode(matrix) <- "numeric"
      matrix[!is.finite(matrix)] <- NA_real_
      if (base::ncol(matrix) == 0L) return(matrix)
      for (j in base::seq_len(base::ncol(matrix))) {
        observed <- matrix[, j]
        center <- if (base::any(is.finite(observed))) {
          stats::median(observed[is.finite(observed)])
        } else {
          0
        }
        matrix[, j] <- observed - center
      }
      matrix
    }

    impute_overview_matrix <- function(data) {
      matrix <- base::as.matrix(data)
      storage.mode(matrix) <- "numeric"
      matrix[!is.finite(matrix)] <- NA_real_
      for (j in base::seq_len(base::ncol(matrix))) {
        missing <- is.na(matrix[, j])
        if (base::any(missing)) {
          replacement <- if (base::any(!missing)) {
            stats::median(matrix[!missing, j])
          } else {
            0
          }
          matrix[missing, j] <- replacement
        }
      }
      matrix
    }

    rv <- shiny::reactiveValues(
      sample_info = NULL,
      load_success = FALSE,
      transformed_matrix = NULL,
      normalized_matrix = NULL,
      cor_results = NULL,
      exp_results = NULL
    )

    shiny::observeEvent(input$load_data, {
      tryCatch({
        step4 <- NULL
        step6 <- NULL
        if (!base::is.null(shared_state$workdir)) {
          step4_path <- base::file.path(
            shared_state$workdir, "Step4_data_transformed.rda"
          )
          step6_path <- base::file.path(
            shared_state$workdir, "Step6_data_normalization.rda"
          )
          step4 <- .protvis_load_stage_dataset(
            step4_path, expression_names = "transformed"
          )
          step6 <- .protvis_load_stage_dataset(
            step6_path, expression_names = "normalized_data"
          )
        }

        if (!base::is.null(step4) && !base::is.null(step6)) {
          transformed_mat <- base::as.matrix(step4$expression_data)
          normalized_mat <- base::as.matrix(step6$expression_data)
          storage.mode(transformed_mat) <- "numeric"
          storage.mode(normalized_mat) <- "numeric"

          rv$sample_info <- step6$sample_info
          rv$transformed_matrix <- transformed_mat
          rv$normalized_matrix <- normalized_mat
          shared_state$dataset <- step6
          notice <- paste0(
            "✅ Transformed (Step4) and normalized (Step6) ",
            "ProtVis_dataset stages loaded."
          )
        } else if (inherits(shared_state$dataset, "ProtVis_dataset")) {
          matrix <- base::as.matrix(shared_state$dataset$expression_data)
          storage.mode(matrix) <- "numeric"
          rv$sample_info <- shared_state$dataset$sample_info
          rv$transformed_matrix <- NULL
          rv$normalized_matrix <- matrix
          notice <- paste0(
            "⚠️ Step4/Step6 snapshots were unavailable. QC uses the current ",
            "ProtVis_dataset; archived UMAP reproduction requires ",
            "Step4_data_transformed.rda."
          )
        } else {
          stop(
            "Run transformation and normalization before loading the overview."
          )
        }

        rv$cor_results <- NULL
        rv$exp_results <- NULL
        rv$load_success <- TRUE
        shiny::showNotification(notice, type = "message")
      }, error = function(e) {
        shiny::showNotification(
          base::paste("Error loading data:", e$message),
          type = "error"
        )
        rv$load_success <- FALSE
        rv$transformed_matrix <- NULL
        rv$normalized_matrix <- NULL
        rv$cor_results <- NULL
        rv$exp_results <- NULL
      })
    })

    output$load_status_panel <- shiny::renderUI({
      if (isTRUE(rv$load_success)) {
        items <- list(
          shiny::span("✅ Data loaded successfully", style = "color: green;"),
          shiny::br()
        )
        if (!base::is.null(rv$transformed_matrix)) {
          items <- base::c(
            items,
            list(base::paste(
              "Transformed data:",
              base::nrow(rv$transformed_matrix), "proteins,",
              base::ncol(rv$transformed_matrix), "samples"
            ), shiny::br())
          )
        }
        items <- base::c(
          items,
          list(base::paste(
            "Normalized data:",
            base::nrow(rv$normalized_matrix), "proteins,",
            base::ncol(rv$normalized_matrix), "samples"
          ))
        )
        do.call(shiny::div, items)
      } else {
        shiny::span("❌ Data not loaded", style = "color: red;")
      }
    })

    shiny::observeEvent(input$run_correlation, {
      shiny::req(rv$normalized_matrix)

      shiny::withProgress(message = "Calculating correlations...", value = 0.5, {
        matrix <- base::as.matrix(rv$normalized_matrix)
        storage.mode(matrix) <- "numeric"
        if (base::ncol(matrix) < 2L) {
          shiny::showNotification(
            "At least two samples are required for correlation analysis.",
            type = "error"
          )
          rv$cor_results <- NULL
          return()
        }
        result <- tryCatch(
          stats::cor(
            matrix,
            method = base::tolower(input$cor_method),
            use = "pairwise.complete.obs"
          ),
          error = function(e) {
            shiny::showNotification(
              paste("Correlation failed:", conditionMessage(e)),
              type = "error"
            )
            NULL
          }
        )
        if (!base::is.null(result)) {
          # Constant or entirely missing samples have undefined correlation.
          # Preserve that information in the result; the heatmap converts it
          # to a neutral display value without modifying the dataset.
          diag(result) <- 1
        }
        rv$cor_results <- result
        if (!base::is.null(result)) {
          .protvis_record_shared_run(
            shared_state,
            module = "overview_correlation",
            method = base::tolower(input$cor_method),
            category = "qc",
            parameters = list(
              method = base::tolower(input$cor_method),
              cluster_rows = isTRUE(input$cor_cluster_rows),
              cluster_columns = isTRUE(input$cor_cluster_columns)
            ),
            matrices = list(correlation_matrix = result),
            plot_data = list(
              correlation_matrix = as.data.frame(result)
            ),
            plot_config = list(
              color_min = input$cor_color_min,
              color_max = input$cor_color_max,
              low_color = input$cor_low_color,
              mid_color = input$cor_mid_color,
              high_color = input$cor_high_color
            )
          )
        }
        shiny::incProgress(1, detail = "Done")
      })
    })

    make_metadata_annotation <- function(matrix) {
      metadata_share <- base::data.frame(
        sample_id = base::colnames(matrix), stringsAsFactors = FALSE
      )
      info <- rv$sample_info
      if (base::is.null(info) || !is.data.frame(info)) {
        info <- base::data.frame(sample_id = character(),
                                 stringsAsFactors = FALSE)
      }
      info_index <- match(metadata_share$sample_id, info$sample_id)
      if ("maxquant_id" %in% base::colnames(info)) {
        fallback_index <- match(metadata_share$sample_id, info$maxquant_id)
        info_index[is.na(info_index)] <- fallback_index[is.na(info_index)]
      }
      for (column in base::setdiff(base::colnames(info), "sample_id")) {
        metadata_share[[column]] <- info[[column]][info_index]
      }
      tissue_values <- if ("tissue2" %in% names(metadata_share)) {
        as.character(metadata_share$tissue2)
      } else if ("tissue" %in% names(metadata_share)) {
        as.character(metadata_share$tissue)
      } else {
        rep(NA_character_, nrow(metadata_share))
      }
      tissue_lower <- tolower(tissue_values)
      tissue_values[grepl("root|below[ ._-]*ground|underground", tissue_lower)] <- "Root"
      tissue_values[grepl("leaf|shoot|stem|above[ ._-]*ground|aerial", tissue_lower)] <- "Shoot"
      sample_lower <- tolower(metadata_share$sample_id)
      fallback_tissue <- ifelse(
        grepl("root|below[ ._-]*ground|underground", sample_lower), "Root",
        ifelse(grepl("leaf|shoot|stem|above[ ._-]*ground|aerial", sample_lower),
               "Shoot", NA_character_)
      )
      channel <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", metadata_share$sample_id)))
      fallback_tissue[is.na(fallback_tissue) & !is.na(channel) & channel <= 3L] <- "Root"
      fallback_tissue[is.na(fallback_tissue) & !is.na(channel) & channel >= 4L] <- "Shoot"
      tissue_values[is.na(tissue_values) | !nzchar(tissue_values) |
                      tissue_values == "NA" | tissue_values == "All samples"] <-
        fallback_tissue[is.na(tissue_values) | !nzchar(tissue_values) |
                         tissue_values == "NA" | tissue_values == "All samples"]
      tissue_values[is.na(tissue_values) | !nzchar(tissue_values)] <- "All samples"
      metadata_share$tissue2 <- tissue_values
      metadata_share$species <- if ("species" %in% names(metadata_share)) {
        as.character(metadata_share$species)
      } else {
        ifelse(grepl("B73", metadata_share$sample_id, ignore.case = TRUE),
               "Zea mays ssp. mays",
               ifelse(grepl("Y12", metadata_share$sample_id,
                            ignore.case = TRUE),
                      "Zea mays ssp. mexicana", "All samples"))
      }
      metadata_share$species[is.na(metadata_share$species) |
                               !nzchar(metadata_share$species)] <- "All samples"
      metadata_share$group <- .protvis_sample_group_values(
        rv$sample_info,
        metadata_share$sample_id,
        mode = "triplicate"
      )
      group_colors <- .protvis_group_palette(metadata_share$group)
      ComplexHeatmap::rowAnnotation(
        Group = base::as.matrix(metadata_share["group"]),
        Tissue = base::as.matrix(metadata_share["tissue2"]),
        Species = base::as.matrix(metadata_share["species"]),
        col = base::list(
          Group = group_colors,
          Tissue = c("Shoot" = "#65a30d", "Root" = "#c2410c",
                     "Leaf" = "#65a30d", "Pulvinus" = "#a16207",
                     "Stem" = "#166534",
                     "Shoot.tip" = "#2563eb", "All samples" = "#94a3b8"),
          Species = c("Zea mays ssp. mays" = "#f59e0b",
                      "Zea mays ssp. mexicana" = "#84cc16",
                      "All samples" = "#94a3b8")
        ),
        show_legend = c(
          Group = FALSE,
          Tissue = TRUE,
          Species = TRUE
        ),
        annotation_name_gp = grid::gpar(
          fontsize = 7, fontfamily = "sans", fontface = "plain"
        ),
        annotation_legend_param = base::list(
          Tissue = base::list(
            title_gp = grid::gpar(
              fontsize = 7, fontfamily = "sans", fontface = "plain"
            ),
            labels_gp = grid::gpar(
              fontsize = 7, fontfamily = "sans", fontface = "plain"
            )
          ),
          Species = base::list(
            title_gp = grid::gpar(
              fontsize = 7, fontfamily = "sans", fontface = "plain"
            ),
            labels_gp = grid::gpar(
              fontsize = 7, fontfamily = "sans", fontface = "plain"
            )
          )
        )
      )
    }

    cor_heatmap <- shiny::reactive({
      shiny::req(isTRUE(rv$load_success))
      shiny::req(!base::is.null(rv$cor_results))
      shiny::req(!base::is.null(rv$sample_info))

      ha <- make_metadata_annotation(rv$normalized_matrix)

      min_break <- input$cor_color_min
      max_break <- input$cor_color_max
      shiny::validate(shiny::need(
        is.finite(min_break) && is.finite(max_break) && min_break < max_break,
        "Correlation color limits must be finite and min < max."
      ))
      mid_break <- (min_break + max_break) / 2
      heatmap_matrix <- rv$cor_results
      heatmap_matrix[!is.finite(heatmap_matrix)] <- 0
      sample_groups <- .protvis_sample_group_values(
        rv$sample_info,
        base::rownames(heatmap_matrix),
        mode = "triplicate"
      )
      sample_groups <- base::factor(
        sample_groups,
        levels = base::unique(sample_groups)
      )

      ComplexHeatmap::Heatmap(
        heatmap_matrix,
        right_annotation = ha,
        row_split = sample_groups,
        column_split = sample_groups,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        row_gap = grid::unit(1.2, "mm"),
        column_gap = grid::unit(1.2, "mm"),
        row_title = NULL,
        column_title = NULL,
        cluster_rows = isTRUE(input$cor_cluster_rows),
        cluster_columns = isTRUE(input$cor_cluster_columns),
        show_row_names = TRUE,
        show_column_names = isTRUE(input$cor_show_column_names),
        row_names_gp = grid::gpar(
          fontsize = 7, fontfamily = "sans", fontface = "plain"
        ),
        column_names_gp = grid::gpar(
          fontsize = 7, fontfamily = "sans", fontface = "plain"
        ),
        border = "black",
        na_col = "#d1d5db",
        name = "r",
        col = circlize::colorRamp2(
          breaks = c(min_break, mid_break, max_break),
          colors = c(
            input$cor_low_color,
            input$cor_mid_color,
            input$cor_high_color
          )
        ),
        heatmap_legend_param = base::list(
          title_gp = grid::gpar(
            fontsize = 7, fontfamily = "sans", fontface = "plain"
          ),
          labels_gp = grid::gpar(
            fontsize = 7, fontfamily = "sans", fontface = "plain"
          )
        ),
        cell_fun = if (isTRUE(input$cor_show_numbers)) function(j, i, x, y, width, height, fill) {
          grid::grid.text(
            label = if (is.finite(rv$cor_results[i, j])) {
              base::round(rv$cor_results[i, j], 2)
            } else {
              "NA"
            },
            x = x,
            y = y,
            gp = grid::gpar(
              fontsize = 6, col = "white",
              fontfamily = "sans", fontface = "plain"
            )
          )
        } else NULL
      )
    })

    output$cor_res <- shiny::renderPlot({
      shiny::validate(
        shiny::need(isTRUE(rv$load_success), "")
      )
      shiny::validate(
        shiny::need(!base::is.null(rv$cor_results), "")
      )

      tryCatch({
        ht <- cor_heatmap()
        shiny::req(!base::is.null(ht))
        ComplexHeatmap::draw(ht)
      }, error = function(e) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste("Correlation unavailable:",
                                       conditionMessage(e)), cex = 0.9)
      })
    })

    output$cor_download_pdf <- shiny::downloadHandler(
      filename = function() {
        base::paste0("correlation_heatmap_", base::Sys.Date(), ".pdf")
      },
      content = function(file) {
        shiny::req(isTRUE(rv$load_success))
        shiny::req(!base::is.null(rv$cor_results))

        grDevices::pdf(
          file,
          width = input$cor_plot_width,
          height = input$cor_plot_height
        )
        ComplexHeatmap::draw(cor_heatmap())
        grDevices::dev.off()
      }
    )

    # Download the computed sample-by-sample correlation matrix without
    # rendering another table in the UI. Keep undefined pairwise correlations
    # as NA in the exported file rather than replacing them with heatmap display
    # values.
    output$cor_download_table <- shiny::downloadHandler(
      filename = function() {
        method <- base::tolower(input$cor_method %||% "pearson")
        base::paste0(
          "correlation_matrix_", method, "_", base::Sys.Date(), ".csv"
        )
      },
      content = function(file) {
        shiny::req(isTRUE(rv$load_success))
        shiny::req(!base::is.null(rv$cor_results))

        result <- base::as.data.frame(
          rv$cor_results,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        result <- base::data.frame(
          Sample = base::rownames(rv$cor_results),
          result,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        utils::write.csv(
          result,
          file = file,
          row.names = FALSE,
          na = "NA"
        )
      }
    )

    shiny::observeEvent(input$run_expression, {
      shiny::req(rv$normalized_matrix)

      shiny::withProgress(message = "Analyzing expression patterns...", value = 0.5, {
        mat <- rv$normalized_matrix
        mat <- base::as.matrix(mat)
        storage.mode(mat) <- "numeric"
        row_vars <- matrixStats::rowVars(mat, na.rm = TRUE)
        row_vars[!is.finite(row_vars)] <- -Inf

        top_n <- base::min(input$exp_top_n, base::nrow(mat))
        top_idx <- base::order(row_vars, decreasing = TRUE)[base::seq_len(top_n)]
        mat <- mat[top_idx, , drop = FALSE]

        if (isTRUE(input$exp_scale) && input$exp_scale_method != "none") {
          if (identical(input$exp_scale_method, "column")) {
            mat <- base::scale(mat)
          } else {
            mat <- base::t(base::scale(base::t(mat)))
          }
          mat[!is.finite(mat)] <- NA_real_
        }

        rv$exp_results <- base::as.data.frame(
          mat, stringsAsFactors = FALSE, check.names = FALSE
        )
        .protvis_record_shared_run(
          shared_state,
          module = "overview_expression_pattern",
          method = if (isTRUE(input$exp_scale)) {
            paste0("top_variance_", input$exp_scale_method)
          } else {
            "top_variance_unscaled"
          },
          category = "qc",
          parameters = list(
            top_n = top_n,
            scale = isTRUE(input$exp_scale),
            scale_method = input$exp_scale_method,
            cluster_rows = isTRUE(input$exp_cluster_rows),
            cluster_columns = isTRUE(input$exp_cluster_columns)
          ),
          matrices = list(expression_pattern = mat),
          tables = list(
            expression_pattern = data.frame(
              protein_id = rownames(mat),
              mat,
              check.names = FALSE,
              stringsAsFactors = FALSE
            )
          ),
          plot_data = list(expression_pattern = as.data.frame(mat))
        )
        shiny::incProgress(1, detail = "Done")
      })
    })

    exp_heatmap <- shiny::reactive({
      shiny::req(isTRUE(rv$load_success))
      shiny::req(!base::is.null(rv$exp_results))
      shiny::req(!base::is.null(rv$sample_info))

      ha <- make_metadata_annotation(rv$normalized_matrix)

      min_break <- input$exp_color_min
      max_break <- input$exp_color_max
      shiny::validate(shiny::need(
        is.finite(min_break) && is.finite(max_break) && min_break < max_break,
        "Expression color limits must be finite and min < max."
      ))
      mid_break <- (min_break + max_break) / 2
      heatmap_matrix <- base::t(base::as.matrix(rv$exp_results))
      cluster_matrix <- heatmap_matrix
      for (i in base::seq_len(base::nrow(cluster_matrix))) {
        missing <- !is.finite(cluster_matrix[i, ])
        if (base::any(missing)) {
          replacement <- stats::median(cluster_matrix[i, !missing], na.rm = TRUE)
          if (!is.finite(replacement)) replacement <- 0
          cluster_matrix[i, missing] <- replacement
        }
      }
      # Rows are samples and are split into biological triplicate groups.
      # ComplexHeatmap does not allow a categorical row_split together with a
      # precomputed row dendrogram. Use a logical cluster_rows flag instead:
      # this clusters samples within each triplicate slice while preserving the
      # group order via cluster_row_slices = FALSE.
      cluster_sample_rows <- isTRUE(input$exp_cluster_rows) &&
        base::nrow(cluster_matrix) > 1L
      column_dend <- if (isTRUE(input$exp_cluster_columns) &&
                         base::ncol(cluster_matrix) > 1L) {
        stats::hclust(stats::dist(base::t(cluster_matrix)))
      } else FALSE

      sample_groups <- .protvis_sample_group_values(
        rv$sample_info,
        base::rownames(heatmap_matrix),
        mode = "triplicate"
      )
      sample_groups <- base::factor(
        sample_groups,
        levels = base::unique(sample_groups)
      )

      ComplexHeatmap::Heatmap(
        heatmap_matrix,
        right_annotation = ha,
        row_split = sample_groups,
        cluster_row_slices = FALSE,
        row_gap = grid::unit(1.2, "mm"),
        row_title = NULL,
        cluster_rows = cluster_sample_rows,
        cluster_columns = column_dend,
        show_row_names = TRUE,
        show_column_names = isTRUE(input$exp_show_feature_names),
        row_names_gp = grid::gpar(
          fontsize = 7, fontfamily = "sans", fontface = "plain"
        ),
        column_names_gp = grid::gpar(
          fontsize = 7, fontfamily = "sans", fontface = "plain"
        ),
        border = "black",
        na_col = "#d1d5db",
        name = ifelse(isTRUE(input$exp_scale) &&
                        input$exp_scale_method != "none",
                      "Z-score", "Intensity"),
        col = circlize::colorRamp2(
          breaks = c(min_break, mid_break, max_break),
          colors = c(
            input$exp_low_color,
            input$exp_mid_color,
            input$exp_high_color
          )
        ),
        heatmap_legend_param = base::list(
          title_gp = grid::gpar(
            fontsize = 7, fontfamily = "sans", fontface = "plain"
          ),
          labels_gp = grid::gpar(
            fontsize = 7, fontfamily = "sans", fontface = "plain"
          )
        )
      )
    })

    output$expression_pattern <- shiny::renderPlot({
      shiny::validate(
        shiny::need(isTRUE(rv$load_success), "")
      )
      shiny::validate(
        shiny::need(!base::is.null(rv$exp_results), "")
      )

      tryCatch({
        ht <- exp_heatmap()
        shiny::req(!base::is.null(ht))
        ComplexHeatmap::draw(ht)
      }, error = function(e) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste("Expression pattern unavailable:",
                                       conditionMessage(e)), cex = 0.9)
      })
    })

    output$exp_download_pdf <- shiny::downloadHandler(
      filename = function() {
        base::paste0("expression_pattern_heatmap_", base::Sys.Date(), ".pdf")
      },
      content = function(file) {
        shiny::req(isTRUE(rv$load_success))
        shiny::req(!base::is.null(rv$exp_results))

        grDevices::pdf(
          file,
          width = input$exp_plot_width,
          height = input$exp_plot_height
        )
        ComplexHeatmap::draw(exp_heatmap())
        grDevices::dev.off()
      }
    )

    DR_results <- shiny::reactiveValues(
      reproduction = NULL
    )

    archived_umap_groups <- c(
      "Root_VE",
      "Root_V1.V2",
      "Root_V4",
      "Leaf_VE.V1.V2",
      "Leaf_V4.V6.V8"
    )

    run_archived_umap <- function(data) {
      if (base::is.null(data)) {
        stop(
          "Step4 transformed data are required for the archived UMAP reproduction.",
          call. = FALSE
        )
      }

      matrix <- base::as.matrix(data)
      storage.mode(matrix) <- "numeric"

      # Match the archived Figure 3 code: select columns containing VE, V2 or
      # V4, then sort sample names before UMAP.
      selected <- base::grepl("VE|V2|V4", base::colnames(matrix))
      sample_names <- base::sort(base::colnames(matrix)[selected])
      if (base::length(sample_names) < 3L) {
        stop(
          "Archived UMAP reproduction requires the five early developmental groups.",
          call. = FALSE
        )
      }
      matrix <- matrix[, sample_names, drop = FALSE]

      # Match dplyr::filter(rowSums(.) > 0) from the archived script. Rows with
      # non-finite sums are discarded just as filter() discards NA conditions.
      keep <- base::rowSums(matrix) > 0
      keep[base::is.na(keep)] <- FALSE
      matrix <- matrix[keep, , drop = FALSE]
      if (base::nrow(matrix) < 2L) {
        stop("Too few complete positive features remain for UMAP.", call. = FALSE)
      }

      base::set.seed(10086)
      fit <- umap::umap(base::t(matrix))
      df <- base::as.data.frame(
        fit$layout[, 1:2, drop = FALSE],
        stringsAsFactors = FALSE
      )
      base::colnames(df) <- c("UMAP1", "UMAP2")
      df$sample <- base::rownames(df)
      if (base::is.null(df$sample) || base::any(!base::nzchar(df$sample))) {
        df$sample <- sample_names
      }
      df$Genotype <- base::sub("_.*$", "", df$sample)
      df$Group <- base::sub(
        "^(B73|Y12)_", "",
        base::sub("_[123]$", "", df$sample)
      )
      df$Species <- base::ifelse(
        df$Genotype == "B73",
        "Zea mays ssp. mays",
        "Zea mays ssp. mexicana"
      )
      df$Group <- base::factor(df$Group, levels = archived_umap_groups)
      df$Region <- base::ifelse(
        df$Group %in% c("Root_VE", "Root_V1.V2"),
        "Root VE-V2",
        base::ifelse(
          df$Group == "Root_V4",
          "Root V4",
          base::ifelse(
            df$Group == "Leaf_VE.V1.V2",
            "Leaf VE-V2",
            base::ifelse(
              df$Group == "Leaf_V4.V6.V8",
              "Leaf V4-V8",
              NA_character_
            )
          )
        )
      )
      df$Region <- base::factor(
        df$Region,
        levels = c("Root VE-V2", "Root V4", "Leaf VE-V2", "Leaf V4-V8")
      )

      df <- df[!base::is.na(df$Group) & !base::is.na(df$Region), , drop = FALSE]
      if (base::nrow(df) < 3L) {
        stop("No archived early-development sample groups were identified.",
             call. = FALSE)
      }
      df
    }

    plot_archived_umap <- function(df) {
      region_centres <- df |>
        dplyr::group_by(Region) |>
        dplyr::summarise(
          UMAP1 = base::mean(UMAP1),
          UMAP2 = base::mean(UMAP2),
          .groups = "drop"
        )

      ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2)) +
        ggplot2::stat_ellipse(
          ggplot2::aes(group = Region, fill = Region),
          geom = "polygon",
          type = "norm",
          alpha = 0.18,
          colour = "black",
          linewidth = 0.35,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          ggplot2::aes(colour = Group, shape = Species),
          size = 2.3,
          alpha = 0.85
        ) +
        ggplot2::geom_label(
          data = region_centres,
          ggplot2::aes(
            x = UMAP1, y = UMAP2, label = Region
          ),
          inherit.aes = FALSE,
          size = 2.6,
          label.size = NA,
          fill = scales::alpha("white", 0.65)
        ) +
        ggsci::scale_color_lancet() +
        ggplot2::scale_fill_manual(
          values = c(
            "Root VE-V2" = "#F4A6A1",
            "Root V4" = "#E7C570",
            "Leaf VE-V2" = "#9CD6B0",
            "Leaf V4-V8" = "#9CD6B0"
          )
        ) +
        ggplot2::labs(
          title = "UMAP analysis",
          x = "UMAP 1",
          y = "UMAP 2",
          colour = "Group",
          shape = NULL
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5, face = "bold"
          ),
          panel.border = ggplot2::element_rect(
            colour = "black", linewidth = 0.8
          )
        )
    }

    shiny::observeEvent(input$DR_analyse, {
      if (!isTRUE(rv$load_success)) {
        shiny::showNotification(
          "Load data before dimensionality reduction.",
          type = "warning"
        )
        return(invisible(NULL))
      }
      if (base::is.null(rv$transformed_matrix)) {
        shiny::showNotification(
          paste0(
            "Archived UMAP reproduction requires Step4_data_transformed.rda. ",
            "Run the Transformation step first."
          ),
          type = "error",
          duration = 8
        )
        return(invisible(NULL))
      }

      shiny::withProgress(
        message = "Reproducing archived UMAP...",
        value = 0.5,
        {
          result <- tryCatch(
            run_archived_umap(rv$transformed_matrix),
            error = function(e) {
              shiny::showNotification(
                base::paste(
                  "UMAP reproduction failed:",
                  base::conditionMessage(e)
                ),
                type = "error",
                duration = 8
              )
              NULL
            }
          )
          DR_results$reproduction <- result
          shiny::incProgress(1, detail = "Done")
        }
      )

      if (!base::is.null(DR_results$reproduction)) {
        .protvis_record_shared_run(
          shared_state,
          module = "dimensionality_reduction",
          method = "UMAP_archived_Figure3",
          category = "dimensionality_reduction",
          parameters = list(
            source_stage = "Step4_data_transformed",
            seed = 10086L,
            groups = archived_umap_groups,
            ellipse_regions = c(
              "Root VE-V2", "Root V4", "Leaf VE-V2", "Leaf V4-V8"
            )
          ),
          tables = list(
            archived_umap = DR_results$reproduction
          ),
          plot_data = list(
            archived_umap = DR_results$reproduction
          )
        )
      }
    })

    output$DR_Reproduction <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          !base::is.null(DR_results$reproduction),
          "Run UMAP to display the archived Figure 3 reproduction."
        )
      )
      print(plot_archived_umap(DR_results$reproduction))
    })

    output$dr_download_pdf <- shiny::downloadHandler(
      filename = function() {
        base::paste0(
          "UMAP_archived_Figure3_",
          base::Sys.Date(),
          ".pdf"
        )
      },
      content = function(file) {
        shiny::req(DR_results$reproduction)
        grDevices::pdf(
          file,
          width = input$dr_plot_width,
          height = input$dr_plot_height
        )
        print(plot_archived_umap(DR_results$reproduction))
        grDevices::dev.off()
      }
    )

    qc_matrix <- shiny::reactive({
      shiny::req(isTRUE(rv$load_success))
      shiny::req(rv$normalized_matrix)
      base::as.matrix(rv$normalized_matrix)
    })

    qc_long_intensity <- shiny::reactive({
      mat <- qc_matrix()
      long <- base::data.frame(
        Sample = rep(base::colnames(mat), each = base::nrow(mat)),
        Intensity = as.vector(mat),
        stringsAsFactors = FALSE
      )
      long[is.finite(long$Intensity), , drop = FALSE]
    })

    qc_sample_total_plot <- shiny::reactive({
      mat <- qc_matrix()
      total_df <- base::data.frame(
        Sample = base::colnames(mat),
        TotalIntensity = base::colSums(mat, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      ggplot2::ggplot(total_df, ggplot2::aes(x = Sample, y = TotalIntensity)) +
        ggplot2::geom_col(fill = "#2563eb") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Sample total normalized intensity", x = NULL, y = "Total intensity")
    })

    qc_missing_rate_plot <- shiny::reactive({
      mat <- qc_matrix()
      missing_df <- base::data.frame(
        Sample = base::colnames(mat),
        MissingRate = base::colMeans(base::is.na(mat)),
        stringsAsFactors = FALSE
      )
      ggplot2::ggplot(missing_df, ggplot2::aes(x = Sample, y = MissingRate)) +
        ggplot2::geom_col(fill = "#dc2626") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Missing value rate by sample", x = NULL, y = "Missing rate")
    })

    qc_boxplot <- shiny::reactive({
      ggplot2::ggplot(qc_long_intensity(), ggplot2::aes(x = Sample, y = Intensity)) +
        ggplot2::geom_boxplot(fill = "#38bdf8", outlier.size = 0.6,
                              na.rm = TRUE) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Normalized intensity distribution", x = NULL, y = "Intensity")
    })

    qc_density_plot <- shiny::reactive({
      ggplot2::ggplot(
        qc_long_intensity(),
        ggplot2::aes(x = Intensity, color = Sample)
      ) +
        ggplot2::geom_density(na.rm = TRUE, linewidth = 0.65) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
          legend.position = "none"
        ) +
        ggplot2::labs(
          title = "Normalized intensity density",
          x = "Intensity",
          y = "Density"
        )
    })

    qc_pca_plot <- shiny::reactive({
      mat <- qc_matrix()
      row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
      keep <- is.finite(row_sds) & row_sds > 0
      top_n <- base::min(input$qc_top_n, base::sum(keep))
      if (top_n > 0 && base::sum(keep) > top_n) {
        top_idx <- base::order(row_sds, decreasing = TRUE)[base::seq_len(top_n)]
        keep <- base::seq_along(row_sds) %in% top_idx
      }
      pca_mat <- base::t(mat[keep, , drop = FALSE])
      pca_mat[base::is.na(pca_mat)] <- 0
      shiny::validate(shiny::need(base::nrow(pca_mat) >= 2 && base::ncol(pca_mat) >= 2, "Need at least two samples and two variable features for PCA."))
      pca <- stats::prcomp(pca_mat, center = TRUE, scale. = TRUE)
      var_exp <- base::round(100 * (pca$sdev^2 / base::sum(pca$sdev^2))[1:2], 1)
      pca_df <- base::data.frame(Sample = base::rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
      ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample)) +
        ggplot2::geom_point(size = 3, color = "#7c3aed") +
        ggplot2::geom_text(vjust = -0.7, size = 3) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = "PCA of normalized proteomics samples", x = base::paste0("PC1 (", var_exp[1], "%)"), y = base::paste0("PC2 (", var_exp[2], "%)"))
    })

    qc_cv_plot <- shiny::reactive({
      mat <- qc_matrix()
      row_mean <- base::rowMeans(mat, na.rm = TRUE)
      row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
      cv_df <- base::data.frame(CV = row_sd / base::abs(row_mean), stringsAsFactors = FALSE)
      cv_df <- cv_df[base::is.finite(cv_df$CV), , drop = FALSE]
      ggplot2::ggplot(cv_df, ggplot2::aes(x = CV)) +
        ggplot2::geom_histogram(bins = 50, fill = "#22c55e", color = "white") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = "Protein coefficient of variation", x = "CV", y = "Protein count")
    })

    qc_plot_by_type <- function(type) {
      switch(
        type,
        density = qc_density_plot(),
        cv = qc_cv_plot(),
        qc_density_plot()
      )
    }

    output$qc_summary <- shiny::renderPrint({
      mat <- qc_matrix()
      cat("Proteomics QC summary\n")
      cat("Proteins/features:", base::nrow(mat), "\n")
      cat("Samples:", base::ncol(mat), "\n")
      cat("Overall missing rate:", base::round(base::mean(base::is.na(mat)), 4), "\n")
      cat("Median sample intensity range:", base::paste(base::round(base::range(apply(mat, 2, stats::median, na.rm = TRUE)), 4), collapse = " - "), "\n")
    })

    safe_qc_plot <- function(plot_function) {
      tryCatch({
        print(plot_function())
      }, error = function(e) {
        graphics::plot.new()
        graphics::text(
          0.5, 0.5, paste("QC plot unavailable:", conditionMessage(e)),
          cex = 0.85
        )
      })
    }

    output$qc_sample_total_plot <- shiny::renderPlot(
      safe_qc_plot(qc_sample_total_plot), height = 280
    )
    output$qc_missing_rate_plot <- shiny::renderPlot(
      safe_qc_plot(qc_missing_rate_plot), height = 180
    )
    output$qc_boxplot <- shiny::renderPlot(
      safe_qc_plot(qc_boxplot), height = 260
    )
    output$qc_density_plot <- shiny::renderPlot(
      safe_qc_plot(qc_density_plot), height = 430
    )
    output$qc_pca_plot <- shiny::renderPlot(
      safe_qc_plot(qc_pca_plot), height = 220
    )
    output$qc_cv_plot <- shiny::renderPlot(
      safe_qc_plot(qc_cv_plot), height = 350
    )

    output$qc_download_pdf <- shiny::downloadHandler(
      filename = function() {
        base::paste0("overview_proteomics_qc_", input$qc_download_plot_type, "_", base::Sys.Date(), ".pdf")
      },
      content = function(file) {
        grDevices::pdf(file, width = input$qc_plot_width, height = input$qc_plot_height)
        print(qc_plot_by_type(input$qc_download_plot_type))
        grDevices::dev.off()
      }
    )

    output$qc_download_matrix <- shiny::downloadHandler(
      filename = function() {
        base::paste0("overview_normalized_matrix_", base::Sys.Date(), ".csv")
      },
      content = function(file) {
        mat <- qc_matrix()
        out <- base::data.frame(ID = base::rownames(mat), mat, check.names = FALSE)
        utils::write.csv(out, file, row.names = FALSE)
      }
    )
  })
}
