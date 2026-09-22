#' @import shiny
#' @import bslib
#' @importFrom colourpicker colourInput
#' @importFrom bsicons bs_icon
#' @name DEG_ui
#' @title DEG Analysis UI
#' @description This function creates the user interface for DEG (Differential Expression Analysis) that includes file upload options, PCA plot settings, and volcano plot settings.
#' @param id A unique ID for the Shiny module, used to create input/output bindings.
#' @export
#'
DEG_ui <- function(id) {
  ns <- NS(id)

  shiny::tagList(
    shiny::tags$style(
      shiny::HTML("
        .pv-transcriptome-gsea-status {
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
          align-items: center;
          margin-bottom: 10px;
        }
        .pv-transcriptome-gsea-badge {
          display: inline-flex;
          align-items: center;
          gap: 5px;
          border: 1px solid #d9e7f0;
          border-radius: 999px;
          background: #f7fbfd;
          padding: 4px 9px;
          color: #52697c;
          font-size: 0.76rem;
        }
        .pv-transcriptome-gsea-badge strong {
          color: #17394f;
        }
        .pv-transcriptome-gsea-card .card-body {
          min-height: 620px;
        }
        .pv-transcriptome-gsea-note {
          color: #62798b;
          font-size: 0.78rem;
          line-height: 1.45;
        }
        .pv-transcriptome-gsea-reference {
          border: 1px solid #dce8f0;
          border-radius: 10px;
          padding: 12px 14px;
          background: #f8fbfd;
        }
      ")
    ),

    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 300,

        shiny::fileInput(
          ns("data_count_file"),
          "Upload Gene Expression Data",
          accept = c(".xlsx")
        ),
        shiny::helpText(
          "The first column name in the expression file should be 'GeneID'."
        ),

        shiny::fileInput(
          ns("group_file"),
          "Upload Group Information",
          accept = c(".xlsx")
        ),
        shiny::helpText(
          "The first column in the group file should be 'Sample', ",
          "and the second column should be 'Group'."
        ),

        shiny::actionButton(
          ns("generate_plot"),
          "Run Analysis",
          class = "btn-primary"
        ),
        shiny::hr(),

        bslib::accordion(
          bslib::accordion_panel(
            title = "PCA Settings",
            icon = pca_icon,

            shiny::selectInput(
              ns("pca_colby"),
              "Color by:",
              choices = c("None" = "none"),
              selected = "none"
            ),
            shiny::selectInput(
              ns("pca_shapeby"),
              "Shape by:",
              choices = c("None" = "none"),
              selected = "none"
            ),
            shiny::selectInput(
              ns("pca_pointsize"),
              "Point Size:",
              choices = c(
                "Small" = 2,
                "Medium" = 3,
                "Large" = 4
              ),
              selected = 3
            ),
            shiny::uiOutput(ns("group_colors_ui")),
            colourpicker::colourInput(
              ns("pca_base_color"),
              "Base Color (when no grouping)",
              value = "#2E86AB"
            ),
            shiny::checkboxInput(
              ns("pca_show_labels"),
              "Show Sample Labels",
              value = FALSE
            ),
            shiny::checkboxInput(
              ns("pca_encircle"),
              "Encircle Groups",
              value = TRUE
            ),
            shiny::checkboxInput(
              ns("pca_show_ellipse"),
              "Show Confidence Ellipse",
              value = TRUE
            ),
            shiny::numericInput(
              ns("pca_ellipse_alpha"),
              "Ellipse Transparency",
              value = 0.2,
              min = 0,
              max = 1,
              step = 0.1
            ),
            shiny::numericInput(
              ns("pca_legend_size"),
              "Legend Text Size",
              value = 12,
              min = 8,
              max = 20,
              step = 1
            ),
            shiny::hr(),
            shiny::numericInput(
              ns("download_width_pca"),
              "Width of PCA Plot (inches)",
              value = 8,
              min = 3,
              max = 20
            ),
            shiny::numericInput(
              ns("download_height_pca"),
              "Height of PCA Plot (inches)",
              value = 7,
              min = 3,
              max = 20
            ),
            shiny::downloadButton(
              ns("download_pca"),
              "Download PCA Plot PDF",
              class = "btn-sm"
            ),
            shiny::downloadButton(
              ns("download_pca_data"),
              "Download PCA Data",
              class = "btn-sm"
            )
          ),

          bslib::accordion_panel(
            title = "Volcano Plot Settings",
            icon = volcano_icon,

            colourpicker::colourInput(
              ns("color_up"),
              "Color for Up-regulated",
              value = "salmon"
            ),
            colourpicker::colourInput(
              ns("color_down"),
              "Color for Down-regulated",
              value = "lightblue"
            ),
            colourpicker::colourInput(
              ns("color_not_sig"),
              "Color for Not Significant",
              value = "grey"
            ),
            shiny::numericInput(
              ns("volcano_point_size"),
              "Point Size",
              value = 2,
              min = 1,
              max = 5,
              step = 0.5
            ),
            shiny::sliderInput(
              ns("volcano_alpha"),
              "Point Transparency",
              min = 0.1,
              max = 1,
              value = 0.7,
              step = 0.1
            ),
            shiny::checkboxInput(
              ns("volcano_show_grid"),
              "Show Grid",
              value = FALSE
            ),
            shiny::hr(),
            shiny::numericInput(
              ns("download_width_voc"),
              "Width of Volcano Plot (inches)",
              value = 8,
              min = 3,
              max = 20
            ),
            shiny::numericInput(
              ns("download_height_voc"),
              "Height of Volcano Plot (inches)",
              value = 7,
              min = 3,
              max = 20
            ),
            shiny::downloadButton(
              ns("download_pdf"),
              "Download Volcano Plot PDF",
              class = "btn-sm"
            ),
            shiny::downloadButton(
              ns("download_deg_data"),
              "Download DEG Data",
              class = "btn-sm"
            )
          ),

          bslib::accordion_panel(
            title = "GSEA Settings",
            icon = shiny::icon("chart-line"),

            shiny::div(
              class = "pv-transcriptome-gsea-note",
              "KEGG gene sets are bundled from the maize-teosinte ",
              "Enrichmentdb2 annotation supplied with the built-in RNA-seq workflow. ",
              "Expression and group tables are always supplied by the user."
            ),
            shiny::hr(),

            shiny::selectInput(
              ns("gsea_rank_metric"),
              "Ranking metric",
              choices = c(
                "log2 fold change (RNAseq.R default)" =
                  "log2FoldChange",
                "DESeq2 Wald statistic" = "stat",
                "Signed -log10(P-value)" = "signed_log10_p"
              ),
              selected = "log2FoldChange"
            ),

            shiny::numericInput(
              ns("gsea_min_size"),
              "Minimum gene-set size",
              value = 5,
              min = 1,
              step = 1
            ),
            shiny::numericInput(
              ns("gsea_max_size"),
              "Maximum gene-set size",
              value = 500,
              min = 5,
              step = 10
            ),
            shiny::numericInput(
              ns("gsea_pvalue_cutoff"),
              "GSEA p-value cutoff",
              value = 1,
              min = 0,
              max = 1,
              step = 0.05
            ),
            shiny::numericInput(
              ns("gsea_top_n"),
              "Pathways in dotplot",
              value = 10,
              min = 1,
              max = 50,
              step = 1
            ),

            shiny::uiOutput(
              ns("gsea_pathway_ui")
            ),

            shiny::actionButton(
              ns("run_gsea"),
              "Run / Refresh GSEA",
              icon = shiny::icon("play"),
              class = "btn-outline-primary w-100"
            ),

            shiny::hr(),

            shiny::downloadButton(
              ns("download_gsea_dotplot"),
              "Download GSEA Dotplot PDF",
              class = "btn-sm"
            ),
            shiny::downloadButton(
              ns("download_gsea_curve"),
              "Download GSEA Curve PDF",
              class = "btn-sm"
            ),
            shiny::downloadButton(
              ns("download_gsea_results"),
              "Download GSEA Results",
              class = "btn-sm"
            ),
            shiny::downloadButton(
              ns("download_gsea_ranks"),
              "Download Ranked List",
              class = "btn-sm"
            )
          )
        )
      ),

      shiny::div(
        bslib::layout_columns(
          col_widths = c(6, 6),

          bslib::card(
            height = "800px",
            bslib::card_header(
              "PCA Analysis",
              icon = shiny::icon("chart-pie")
            ),
            bslib::card_body(
              shiny::tabsetPanel(
                type = "tabs",
                shiny::tabPanel(
                  "Plot",
                  shiny::plotOutput(
                    ns("pca_plot"),
                    height = "650px"
                  )
                ),
                shiny::tabPanel(
                  "PCA Data",
                  shiny::div(
                    style = "margin-bottom: 10px;",
                    shiny::downloadButton(
                      ns("download_pca_table"),
                      "Download as CSV",
                      class = "btn-sm btn-success",
                      style = "float: right;"
                    )
                  ),
                  DT::DTOutput(
                    ns("pca_data_table"),
                    height = "600px"
                  )
                )
              )
            )
          ),

          bslib::card(
            height = "800px",
            bslib::card_header(
              "Volcano Plot",
              icon = shiny::icon("fire")
            ),
            bslib::card_body(
              shiny::tabsetPanel(
                type = "tabs",
                shiny::tabPanel(
                  "Plot",
                  shiny::plotOutput(
                    ns("voc_plot"),
                    height = "650px"
                  )
                ),
                shiny::tabPanel(
                  "DEG Results",
                  shiny::div(
                    style = "margin-bottom: 10px;",
                    shiny::downloadButton(
                      ns("download_degs"),
                      "Download as CSV",
                      class = "btn-sm btn-success",
                      style = "float: right;"
                    )
                  ),
                  DT::DTOutput(
                    ns("deg_table"),
                    height = "600px"
                  )
                ),
                shiny::tabPanel(
                  "Statistics",
                  bslib::card(
                    bslib::card_header(
                      "DEG Summary Statistics"
                    ),
                    shiny::tableOutput(ns("deg_stats"))
                  ),
                  bslib::card(
                    bslib::card_header("Top DEGs"),
                    DT::DTOutput(
                      ns("top_degs_table"),
                      height = "300px"
                    )
                  )
                )
              )
            )
          )
        ),

        bslib::card(
          class = "pv-transcriptome-gsea-card",
          full_screen = TRUE,
          bslib::card_header(
            "KEGG Gene Set Enrichment Analysis",
            icon = shiny::icon("chart-line")
          ),
          bslib::card_body(
            shiny::uiOutput(ns("gsea_status")),

            shiny::tabsetPanel(
              type = "tabs",

              shiny::tabPanel(
                "Dotplot",
                shiny::plotOutput(
                  ns("gsea_dotplot"),
                  height = "560px"
                )
              ),

              shiny::tabPanel(
                "Enrichment Curve",
                shiny::plotOutput(
                  ns("gsea_curve"),
                  height = "620px"
                )
              ),

              shiny::tabPanel(
                "GSEA Results",
                DT::DTOutput(
                  ns("gsea_table"),
                  height = "560px"
                )
              ),

              shiny::tabPanel(
                "Ranked List",
                DT::DTOutput(
                  ns("gsea_rank_table"),
                  height = "560px"
                )
              ),

              shiny::tabPanel(
                "Built-in reference",
                shiny::uiOutput(ns("gsea_reference_ui"))
              )
            )
          )
        )
      )
    )
  )
}



.protvis_deg_gsea_locate <- function(file_name) {
  candidates <- c(
    system.file(
      "extdata",
      "transcriptome_gsea",
      file_name,
      package = "ProtVis"
    ),
    file.path(
      "inst",
      "extdata",
      "transcriptome_gsea",
      file_name
    ),
    file.path(
      getwd(),
      "inst",
      "extdata",
      "transcriptome_gsea",
      file_name
    )
  )
  candidates <- candidates[
    nzchar(candidates) & file.exists(candidates)
  ]
  if (!length(candidates)) {
    stop(
      "The bundled transcriptome GSEA asset is unavailable: ",
      file_name,
      call. = FALSE
    )
  }
  candidates[[1L]]
}


.protvis_deg_gsea_normalize_term <- function(x) {
  x <- trimws(as.character(x))
  numeric_only <- grepl("^[0-9]+$", x)
  short_numeric <- numeric_only & nchar(x) < 5L
  x[short_numeric] <- sprintf(
    "%05d",
    as.integer(x[short_numeric])
  )
  x
}


.protvis_deg_load_gsea_background <- function() {
  candidates <- c(
    system.file(
      "extdata",
      "maize_teosinte_KEGG_background.tsv.xz",
      package = "ProtVis"
    ),
    file.path(
      "inst",
      "extdata",
      "maize_teosinte_KEGG_background.tsv.xz"
    ),
    file.path(
      getwd(),
      "inst",
      "extdata",
      "maize_teosinte_KEGG_background.tsv.xz"
    )
  )
  candidates <- candidates[
    nzchar(candidates) & file.exists(candidates)
  ]
  if (!length(candidates)) {
    stop(
      "The built-in maize KEGG background is unavailable.",
      call. = FALSE
    )
  }

  bg <- utils::read.delim(
    candidates[[1L]],
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("TERM", "GENE", "NAME") %in% names(bg))) {
    stop(
      "The built-in KEGG background must contain TERM, GENE and NAME.",
      call. = FALSE
    )
  }

  bg <- bg[, c("TERM", "GENE", "NAME"), drop = FALSE]
  bg$TERM <- .protvis_deg_gsea_normalize_term(bg$TERM)
  bg$GENE <- trimws(as.character(bg$GENE))
  bg$NAME <- trimws(as.character(bg$NAME))
  bg <- bg[
    nzchar(bg$TERM) &
      nzchar(bg$GENE) &
      nzchar(bg$NAME),
    ,
    drop = FALSE
  ]
  bg <- unique(bg)

  t2g <- unique(bg[, c("TERM", "GENE"), drop = FALSE])
  t2n <- bg[
    !duplicated(bg$TERM),
    c("TERM", "NAME"),
    drop = FALSE
  ]

  list(
    TERM2GENE = t2g,
    TERM2NAME = t2n,
    background = bg
  )
}


.protvis_deg_load_gsea_reference <- function() {
  ref <- utils::read.delim(
    .protvis_deg_gsea_locate(
      "RNAseq_GSEA_reference_00940.tsv"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ref$ID <- .protvis_deg_gsea_normalize_term(ref$ID)
  ref
}


.protvis_deg_prepare_gsea_rank <- function(
    res_tbl,
    metric = c(
      "log2FoldChange",
      "stat",
      "signed_log10_p"
    )) {
  metric <- match.arg(metric)
  df <- as.data.frame(
    res_tbl,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("GeneID", "log2FoldChange") %in% names(df))) {
    stop(
      "DESeq2 results must contain GeneID and log2FoldChange.",
      call. = FALSE
    )
  }

  score <- switch(
    metric,
    log2FoldChange =
      suppressWarnings(as.numeric(df$log2FoldChange)),
    stat = {
      if (!"stat" %in% names(df)) {
        stop(
          "DESeq2 Wald statistic is unavailable.",
          call. = FALSE
        )
      }
      suppressWarnings(as.numeric(df$stat))
    },
    signed_log10_p = {
      if (!"pvalue" %in% names(df)) {
        stop(
          "Raw P-values are unavailable.",
          call. = FALSE
        )
      }
      lfc <- suppressWarnings(
        as.numeric(df$log2FoldChange)
      )
      pv <- suppressWarnings(
        as.numeric(df$pvalue)
      )
      sign(lfc) *
        -log10(pmax(pv, .Machine$double.xmin))
    }
  )

  ids <- trimws(as.character(df$GeneID))
  keep <- !is.na(ids) & nzchar(ids) & is.finite(score)
  ids <- ids[keep]
  score <- score[keep]

  if (anyDuplicated(ids)) {
    ord_abs <- order(
      abs(score),
      decreasing = TRUE,
      na.last = NA
    )
    ids <- ids[ord_abs]
    score <- score[ord_abs]
    unique_keep <- !duplicated(ids)
    ids <- ids[unique_keep]
    score <- score[unique_keep]
  }

  ord <- order(score, decreasing = TRUE, na.last = NA)
  score <- score[ord]
  ids <- ids[ord]
  names(score) <- ids

  score
}


.protvis_deg_default_gsea_pathway <- function(
    result_table,
    current = NULL,
    preferred = "00940") {
  df <- as.data.frame(
    result_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!nrow(df) || !"ID" %in% names(df)) {
    return("")
  }

  ids <- as.character(df$ID)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) {
    return("")
  }

  current <- as.character(current %||% "")
  if (length(current) && nzchar(current[[1L]]) &&
      current[[1L]] %in% ids) {
    return(current[[1L]])
  }

  preferred <- as.character(preferred %||% "")
  if (length(preferred) && nzchar(preferred[[1L]]) &&
      preferred[[1L]] %in% ids) {
    return(preferred[[1L]])
  }

  ids[[1L]]
}


.protvis_deg_run_gsea <- function(
    res_tbl,
    metric = "log2FoldChange",
    min_size = 5L,
    max_size = 500L,
    pvalue_cutoff = 1) {
  background <- .protvis_deg_load_gsea_background()
  ranks <- .protvis_deg_prepare_gsea_rank(
    res_tbl,
    metric = metric
  )

  overlap <- intersect(
    names(ranks),
    unique(background$TERM2GENE$GENE)
  )
  if (length(overlap) < 5L) {
    stop(
      "Too few ranked genes overlap the built-in KEGG background.",
      call. = FALSE
    )
  }

  # Important for historical reproduction: RNAseq.R passed the complete
  # DESeq2 ranked list to clusterProfiler::GSEA. Unannotated genes remain in
  # the ranked universe and only pathway membership comes from TERM2GENE.
  ranks <- sort(ranks, decreasing = TRUE)

  result <- clusterProfiler::GSEA(
    geneList = ranks,
    minGSSize = as.integer(min_size),
    maxGSSize = as.integer(max_size),
    pvalueCutoff = as.numeric(pvalue_cutoff),
    TERM2GENE = background$TERM2GENE,
    TERM2NAME = background$TERM2NAME,
    verbose = FALSE
  )

  list(
    result = result,
    ranks = ranks,
    result_table = as.data.frame(result),
    background = background
  )
}


#' @title DEG Analysis Server Logic
#' @description This function contains the server-side logic for performing DEG (Differential Expression Analysis), including PCA and volcano plot generation, and DEG result calculations.
#' @param id A unique ID for the Shiny module, used to create input/output bindings.
#' @import shiny
#' @importFrom readxl read_xlsx
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr mutate across everything case_when arrange desc filter select where
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom DT renderDT datatable formatStyle styleEqual
#' @importFrom grDevices pdf dev.off
#' @importFrom colourpicker updateColourInput
#' @name DEG_server
#' @export

utils::globalVariables(c("padj", "log2FoldChange", "regular",
                         "GeneID","baseMean","lfcSE","pvalue","Regulation"))

DEG_server <- function(id, shared_state = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # State management
    analysis_ready <- shiny::reactiveVal(FALSE)
    # Reactive loading data
    expression_data <- shiny::reactive({
      shiny::req(input$data_count_file)
      df <- readxl::read_xlsx(input$data_count_file$datapath)
      shiny::validate(
        shiny::need("GeneID" %in% base::colnames(df), "Error: Expression file must contain 'GeneID' column"),
        shiny::need(base::ncol(df) > 1, "Error: Expression file must contain sample columns")
      )
      return(df)
    })
    # Load grouping information
    group_data <- shiny::reactive({
      shiny::req(input$group_file)
      df <- readxl::read_xlsx(input$group_file$datapath)
      shiny::validate(
        shiny::need("Sample" %in% base::colnames(df), "Error: Group file must contain 'Sample' column"),
        shiny::need("Group" %in% base::colnames(df), "Error: Group file must contain 'Group' column")
      )
      return(df)
    })
    # Gets the column names of grouped data (used for PCA color and shape selection)
    group_columns <- shiny::reactive({
      shiny::req(group_data())
      cols <- base::colnames(group_data())
      # Exclude Sample column
      cols <- cols[cols != "Sample"]
      return(cols)
    })
    # Get different groups of the current grouping variable.
    selected_groups <- shiny::reactive({
      shiny::req(group_data(), input$pca_colby)
      if (input$pca_colby != "none") {
        groups <- base::unique(group_data()[[input$pca_colby]])
        return(base::sort(base::as.character(groups)))  # 确保是字符型并排序
      }
      return(NULL)
    })
    # Observe the change of packet data and update PCA setting options.
    shiny::observeEvent(group_data(), {
      cols <- group_columns()
      if(base::length(cols) > 0) {
        # Update color selection
        shiny::updateSelectInput(session, "pca_colby",
                          choices = c("None" = "none", cols),
                          selected = "Group")
        # Update shape selection
        shiny::updateSelectInput(session, "pca_shapeby",
                          choices = c("None" = "none", cols),
                          selected = "none")
      }
    })
    # Observe the change of grouping column selection
    shiny::observeEvent(input$pca_colby, {
      if (input$pca_colby != "none" && !is.null(group_data())) {
        # Clear the color input that may exist before.
        shiny::removeUI(
          selector = paste0("#", ns("group_colors_title")),
          immediate = TRUE
        )
      }
    })
    # Generate dynamic color selector
    output$group_colors_ui <- shiny::renderUI({
      groups <- selected_groups()
      if (is.null(groups) || input$pca_colby == "none") {
        return(NULL)  # If no grouping is selected or the grouping is "none", the color selector is not displayed.
      }
      # Generate a set of beautiful default colors
      default_colors <- c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
        "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA"
      )
      # Create color selectors for each group.
      color_pickers <- base::lapply(base::seq_along(groups), function(i) {
        group <- groups[i]
        default_color <- default_colors[(i-1) %% base::length(default_colors) + 1]
        # Create a unique ID for each group.
        group_id <- base::gsub("[^A-Za-z0-9]", "_", group)
        shiny::tagList(
          colourpicker::colourInput(
            ns(paste0("color_", group_id)),
            label = paste("Color for:", group),
            value = default_color
          )
        )
      })
      # Add a reset button
      reset_button <- shiny::actionButton(
        ns("reset_colors"),
        "Reset Colors to Default",
        icon = shiny::icon("refresh"),
        class = "btn-sm btn-outline-secondary"
      )
      shiny::tagList(
        shiny::h5("Customize Group Colors:", id = ns("group_colors_title")),
        shiny::br(),
        color_pickers,
        shiny::br(),
        reset_button
      )
    })
    # Handle color reset button
    shiny::observeEvent(input$reset_colors, {
      groups <- selected_groups()
      if (!is.null(groups)) {
        default_colors <- c(
          "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
          "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
          "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
          "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA"
        )
        for (i in base::seq_along(groups)) {
          group <- groups[i]
          default_color <- default_colors[(i-1) %% base::length(default_colors) + 1]
          group_id <- base::gsub("[^A-Za-z0-9]", "_", group)
          colourpicker::updateColourInput(
            session,
            base::paste0("color_", group_id),
            value = default_color
          )
        }
      }
    })
    # Gets the color selected by the user.
    get_group_colors <- shiny::reactive({
      groups <- selected_groups()
      if (is.null(groups) || input$pca_colby == "none") {
        return(NULL)
      }
      colors <- base::character(0)
      for (group in groups) {
        group_id <- base::gsub("[^A-Za-z0-9]", "_", group)
        color_input <- base::paste0("color_", group_id)
        if (!is.null(input[[color_input]])) {
          colors <- c(colors, input[[color_input]])
        } else {
          # If the color is not set, the default color is used.
          default_colors <- c(
            "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
            "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5"
          )
          default_color <- default_colors[(which(groups == group) - 1) %% base::length(default_colors) + 1]
          colors <- c(colors, default_color)
        }
      }
      base::names(colors) <- groups
      return(colors)
    })
    # Create sample information
    sample_info <- shiny::reactive({
      shiny::req(group_data())
      col_data <- group_data() %>%
        tibble::column_to_rownames("Sample")
      return(col_data)
    })
    # Operational analysis
    shiny::observeEvent(input$generate_plot, {
      # Verification data
      shiny::validate(
        shiny::need(!is.null(expression_data()), "Please upload expression data"),
        shiny::need(!is.null(sample_info()), "Please upload group information"),
        shiny::need(nrow(expression_data()) > 0, "Expression data is empty"),
        shiny::need(nrow(sample_info()) > 0, "Group information is empty")
      )
      # Check whether the sample names match.
      expr_samples <- base::colnames(expression_data())[-1]  # Exclude GeneID column
      group_samples <- base::rownames(sample_info())
      shiny::validate(
        shiny::need(all(expr_samples %in% group_samples),
                    base::paste("Error: Sample names in expression data do not match group data.\n",
                   "Expression samples:", base::paste(expr_samples, collapse = ", "), "\n",
                   "Group samples:", base::paste(group_samples, collapse = ", ")))
      )
      analysis_ready(TRUE)
    })
    # Perform PCA analysis
    pca_result <- shiny::reactive({
      shiny::req(analysis_ready(), expression_data(), sample_info())
      # Extract expression matrix
      expr_mat <- expression_data() %>%
        tibble::column_to_rownames("GeneID") %>%
        base::as.matrix()
      # Ensure that the sample order is consistent.
      expr_mat <- expr_mat[, base::rownames(sample_info()), drop = FALSE]
      # Run PCA
      pca <- PCAtools::pca(expr_mat, metadata = sample_info(), removeVar = 0.1)
      return(pca)
    })
    # Draw PCA diagram
    pca_plot_obj <- shiny::reactive({
      shiny::req(pca_result())
      # Get color and shape settings
      colby <- input$pca_colby
      shapeby <- input$pca_shapeby
      show_labels <- input$pca_show_labels
      encircle <- input$pca_encircle
      show_ellipse <- input$pca_show_ellipse
      ellipse_alpha <- input$pca_ellipse_alpha
      point_size <- base::as.numeric(input$pca_pointsize)
      base_color <- input$pca_base_color
      legend_size <- input$pca_legend_size
      # Basic PCA diagram setting
      pca_args <- base::list(
        pca_result(),
        x = "PC1",
        y = "PC2",
        legendPosition = "right",
        legendLabSize = legend_size,
        legendIconSize = 6,
        pointSize = point_size,
        title = "PCA Plot",
        subtitle = "Principal Component Analysis"
      )
      # Set color
      if (colby != "none") {
        pca_args$colby <- colby
        # Gets the user-defined color.
        group_colors <- get_group_colors()
        if (base::length(group_colors) > 0) {
          pca_args$colkey <- group_colors
        }
      } else {
        pca_args$colby <- NULL
        pca_args$colkey <- base_color
      }
      # Set the shape
      if (shapeby != "none") {
        pca_args$shape <- shapeby
      } else {
        pca_args$shape <- NULL
      }
      # Set sample label
      if (show_labels) {
        pca_args$lab <- base::rownames(pca_result()$metadata)
      } else {
        pca_args$lab <- NULL
      }
      # Set ellipse
      if (encircle && colby != "none" && show_ellipse) {
        pca_args$encircle <- TRUE
        pca_args$encircleFill <- TRUE
        pca_args$encircleAlpha <- ellipse_alpha
        pca_args$encircleLineSize <- 1
      } else {
        pca_args$encircle <- FALSE
      }
      # draw a graph
      pca_plot <- base::do.call(PCAtools::biplot, pca_args)
      return(pca_plot)
    })
    # Prepare PCA data table (show pca_result$rotated)
    pca_rotated_data <- shiny::reactive({
      shiny::req(pca_result())
      # Get the rotated coordinates.
      rotated_data <- base::as.data.frame(pca_result()$rotated)
      rotated_data <- rotated_data[, 1:base::min(10, base::ncol(rotated_data))]  # Only the top 10 principal components are displayed.
      # Add sample name
      rotated_data <- base::cbind(
        Sample = base::rownames(rotated_data),
        rotated_data
      )
      # Add grouping information
      if (!is.null(sample_info())) {
        rotated_data <- base::cbind(
          rotated_data,
          sample_info()
        )
      }
      return(rotated_data)
    })
    # Perform differential expression analysis.
    deseq_results <- shiny::eventReactive(input$generate_plot, {
      shiny::req(expression_data(), sample_info())
      shiny::withProgress(message = 'Running DESeq2 analysis...', value = 0.3, {
        # Prepare counting matrix
        count_mat <- expression_data() %>%
          tibble::column_to_rownames("GeneID") %>%
          dplyr::mutate(dplyr::across(dplyr::everything(), ceiling)) %>%
          base::as.matrix()
        # Ensure that the sample order is consistent.
        count_mat <- count_mat[, base::rownames(sample_info()), drop = FALSE]
        # Create a DESeq2 object
        shiny::incProgress(0.2, detail = "Creating DESeq2 object...")
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = count_mat,
          colData = sample_info(),
          design = ~ Group
        )
        # run DESeq2
        shiny::incProgress(0.3, detail = "Running DESeq2...")
        dds <- DESeq2::DESeq(dds)
        # Get results
        shiny::incProgress(0.2, detail = "Extracting results...")
        res <- DESeq2::results(dds, contrast = c("Group", "B73", "Y12"))
        # Collate results
        res_tbl <- res %>%
          base::as.data.frame() %>%
          tibble::rownames_to_column("GeneID") %>%
          dplyr::mutate(
            regular = dplyr::case_when(
              padj < 0.05 & log2FoldChange > 1 ~ "up",
              padj < 0.05 & log2FoldChange < -1 ~ "down",
              TRUE ~ "not sig"
            ),
            significant = base::ifelse(padj < 0.05 & base::abs(log2FoldChange) > 1, "yes", "no"),
            Regulation = dplyr::case_when(
              regular == "up" ~ "Up-regulated",
              regular == "down" ~ "Down-regulated",
              TRUE ~ "Not significant"
            )
          ) %>%
          dplyr::arrange(padj, dplyr::desc(base::abs(log2FoldChange)))
        return(res_tbl)
      })
    })
    # KEGG GSEA: reproduces the archived RNAseq.R method by
    # ranking the complete DESeq2 result with log2FoldChange by default.
    gsea_bundle <- shiny::eventReactive(
      {
        input$generate_plot
        input$run_gsea
      },
      {
        shiny::req(deseq_results())

        tryCatch(
          .protvis_deg_run_gsea(
            deseq_results(),
            metric =
              input$gsea_rank_metric %||% "log2FoldChange",
            min_size =
              input$gsea_min_size %||% 5L,
            max_size =
              input$gsea_max_size %||% 500L,
            pvalue_cutoff =
              input$gsea_pvalue_cutoff %||% 1
          ),
          error = function(e) {
            list(
              result = NULL,
              ranks = numeric(),
              result_table = data.frame(),
              background = NULL,
              error = conditionMessage(e)
            )
          }
        )
      },
      ignoreInit = TRUE
    )

    output$gsea_pathway_ui <- shiny::renderUI({
      bundle <- gsea_bundle()

      if (is.null(bundle) ||
          !is.null(bundle$error) ||
          is.null(bundle$result) ||
          !nrow(bundle$result_table)) {
        return(
          shiny::selectInput(
            ns("gsea_pathway"),
            "Enrichment curve pathway",
            choices = c(
              "Run analysis first" = ""
            ),
            selected = ""
          )
        )
      }

      df <- bundle$result_table
      ids <- as.character(df$ID)
      labels <- paste0(
        df$Description,
        " [",
        ids,
        "]"
      )
      choices <- stats::setNames(
        ids,
        labels
      )

      selected <- .protvis_deg_default_gsea_pathway(
        df,
        current = isolate(input$gsea_pathway),
        preferred = "00940"
      )

      shiny::selectizeInput(
        ns("gsea_pathway"),
        "Enrichment curve pathway",
        choices = choices,
        selected = selected,
        options = list(
          placeholder = "Search KEGG pathway..."
        )
      )
    })

    selected_gsea_pathway <- shiny::reactive({
      bundle <- gsea_bundle()
      shiny::req(
        bundle,
        is.null(bundle$error),
        nrow(bundle$result_table)
      )

      .protvis_deg_default_gsea_pathway(
        bundle$result_table,
        current = input$gsea_pathway,
        preferred = "00940"
      )
    })


    # Get DEG statistics

    # Get DEG statistics
    deg_stats <- shiny::reactive({
      shiny::req(deseq_results())
      res_tbl <- deseq_results()
      stats <- base::list(
        total_genes = base::nrow(res_tbl),
        up_regulated = base::sum(res_tbl$regular == "up", na.rm = TRUE),
        down_regulated = base::sum(res_tbl$regular == "down", na.rm = TRUE),
        significant = base::sum(res_tbl$regular %in% c("up", "down"), na.rm = TRUE),
        percent_sig = base::round(base::sum(res_tbl$regular %in% c("up", "down"), na.rm = TRUE) / base::nrow(res_tbl) * 100, 2)
      )
      return(stats)
    })
    # obtain top DEGs
    top_degs <- shiny::reactive({
      shiny::req(deseq_results())
      res_tbl <- deseq_results()
      # Obtaining significantly differentially expressed genes
      sig_genes <- res_tbl %>%
        dplyr::filter(regular %in% c("up", "down")) %>%
        dplyr::arrange(padj, dplyr::desc(base::abs(log2FoldChange))) %>%
        utils::head(20)
      return(sig_genes)
    })
    # Persist PCA + DESeq2 results as one immutable schema-v4 run for
    # every explicit analysis request.
    shiny::observeEvent(input$generate_plot, {
      shiny::req(analysis_ready())
      res_tbl <- deseq_results()
      pca_tbl <- pca_rotated_data()
      stats <- deg_stats()
      top_tbl <- top_degs()
      gsea_now <- gsea_bundle()
      gsea_ok <- is.null(gsea_now$error) &&
        !is.null(gsea_now$result) &&
        nrow(gsea_now$result_table)

      .protvis_record_shared_run(
        shared_state,
        module = "deg_deseq2",
        method = "DESeq2 + KEGG GSEA",
        category = "differential_analysis",
        parameters = list(
          contrast = c("Group", "B73", "Y12"),
          padj_threshold = 0.05,
          log2fc_threshold = 1,
          pca_colby = input$pca_colby,
          pca_shapeby = input$pca_shapeby,
          gsea_rank_metric =
            input$gsea_rank_metric %||% "log2FoldChange",
          gsea_min_size =
            input$gsea_min_size %||% 5L,
          gsea_max_size =
            input$gsea_max_size %||% 500L,
          gsea_pvalue_cutoff =
            input$gsea_pvalue_cutoff %||% 1,
          gsea_background =
            "bundled Enrichmentdb2 KEGG annotation"
        ),
        tables = c(
          list(
            differential_expression =
              as.data.frame(res_tbl),
            top_DEGs =
              as.data.frame(top_tbl),
            PCA_scores =
              as.data.frame(pca_tbl)
          ),
          if (gsea_ok) {
            list(
              KEGG_GSEA =
                as.data.frame(gsea_now$result_table),
              GSEA_ranked_list =
                data.frame(
                  GeneID = names(gsea_now$ranks),
                  Rank_metric =
                    as.numeric(gsea_now$ranks),
                  stringsAsFactors = FALSE
                )
            )
          } else {
            list()
          }
        ),
        statistics = c(
          stats,
          list(
            gsea_pathways = if (gsea_ok) {
              nrow(gsea_now$result_table)
            } else {
              0L
            }
          )
        ),
        plot_data = c(
          list(
            volcano = as.data.frame(res_tbl),
            PCA = as.data.frame(pca_tbl)
          ),
          if (gsea_ok) {
            list(
              GSEA = as.data.frame(
                gsea_now$result_table
              )
            )
          } else {
            list()
          }
        ),
        plot_config = list(
          volcano_colors = c(
            up = input$color_up,
            down = input$color_down,
            not_significant = input$color_not_sig
          ),
          point_size = input$volcano_point_size,
          alpha = input$volcano_alpha,
          gsea_top_n =
            input$gsea_top_n %||% 10L,
          gsea_reference_pathway =
            "00940 Phenylpropanoid biosynthesis"
        )
      )
    }, ignoreInit = TRUE, priority = -10)

    # Draw a volcano map
    voc_plot_obj <- shiny::reactive({
      shiny::req(deseq_results())
      res_tbl <- deseq_results()
      # Calculation statistics are used for subheadings
      stats <- deg_stats()
      # Create a volcano map
      p <- ggplot2::ggplot(res_tbl, ggplot2::aes(x = log2FoldChange, y = -log10(padj))) +
        ggplot2::geom_point(ggplot2::aes(color = regular),
                            size = input$volcano_point_size,
                            alpha = input$volcano_alpha) +
        ggplot2::scale_color_manual(
          values = c(
            "up" = input$color_up,
            "down" = input$color_down,
            "not sig" = input$color_not_sig
          ),
          name = "Expression"
        ) +
        ggplot2::geom_hline(
          yintercept = -log10(0.05),
          linetype = "dashed",
          color = "black",
          alpha = 0.5
        ) +
        ggplot2::geom_vline(
          xintercept = c(-1, 1),
          linetype = "dashed",
          color = "black",
          alpha = 0.5
        ) +
        ggplot2::labs(
          title = "Volcano Plot",
          subtitle = paste(
            "Up-regulated:", stats$up_regulated,
            "| Down-regulated:", stats$down_regulated,
            "| Total significant:", stats$significant,
            base::paste0("(", stats$percent_sig, "%)")
          ),
          x = "log2(Fold Change)",
          y = "-log10(Adjusted p-value)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 16, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 12, color = "gray50"),
          axis.title = ggplot2::element_text(size = 12),
          legend.position = "right",
          panel.grid = if(input$volcano_show_grid) ggplot2::element_line(color = "gray90") else ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5)
        ) +
        ggplot2::coord_cartesian(ylim = c(0, base::max(-log10(res_tbl$padj[base::is.finite(-log10(res_tbl$padj))]), na.rm = TRUE) * 1.1))
      return(p)
    })
    # Rendering PCA diagram
    output$pca_plot <- shiny::renderPlot({
      shiny::req(pca_plot_obj())
      pca_plot_obj()
    })
    # Render a volcano map
    output$voc_plot <- shiny::renderPlot({
      shiny::req(voc_plot_obj())
      voc_plot_obj()
    })
    output$gsea_status <- shiny::renderUI({
      bundle <- gsea_bundle()

      if (is.null(bundle)) {
        return(
          shiny::div(
            class = "pv-transcriptome-gsea-note",
            "Run the transcriptome analysis to calculate KEGG GSEA."
          )
        )
      }

      if (!is.null(bundle$error)) {
        return(
          shiny::div(
            class = "text-danger small",
            paste0("GSEA failed: ", bundle$error)
          )
        )
      }

      table <- bundle$result_table
      pathway <- selected_gsea_pathway()
      selected <- table[
        as.character(table$ID) == pathway,
        ,
        drop = FALSE
      ]

      shiny::div(
        class = "pv-transcriptome-gsea-status",

        shiny::span(
          class = "pv-transcriptome-gsea-badge",
          shiny::span("Ranked genes"),
          shiny::tags$strong(
            format(
              length(bundle$ranks),
              big.mark = ","
            )
          )
        ),

        shiny::span(
          class = "pv-transcriptome-gsea-badge",
          shiny::span("Pathways"),
          shiny::tags$strong(nrow(table))
        ),

        shiny::span(
          class = "pv-transcriptome-gsea-badge",
          shiny::span("Ranking"),
          shiny::tags$strong(
            input$gsea_rank_metric %||%
              "log2FoldChange"
          )
        ),

        if (nrow(selected)) {
          shiny::span(
            class = "pv-transcriptome-gsea-badge",
            shiny::span("Selected NES"),
            shiny::tags$strong(
              signif(selected$NES[[1L]], 4)
            )
          )
        },

        shiny::span(
          class = "pv-transcriptome-gsea-badge",
          shiny::span("Pathway"),
          shiny::tags$strong(
            if (nzchar(pathway)) pathway else "—"
          )
        ),

        shiny::span(
          class = "pv-transcriptome-gsea-badge",
          shiny::span("Background"),
          shiny::tags$strong("Built-in maize KEGG")
        )
      )
    })

    gsea_dotplot_obj <- shiny::reactive({
      bundle <- gsea_bundle()
      shiny::req(
        bundle,
        is.null(bundle$error),
        bundle$result,
        nrow(bundle$result_table)
      )

      top_n <- as.integer(
        input$gsea_top_n %||% 10L
      )

      tryCatch(
        {
          enrichplot::dotplot(
            bundle$result,
            showCategory = top_n,
            split = ".sign"
          ) +
            ggplot2::facet_grid(. ~ .sign) +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::theme(
              axis.text.y =
                ggplot2::element_text(
                  colour = "black"
                ),
              axis.text.x =
                ggplot2::element_text(
                  colour = "black"
                )
            )
        },
        error = function(e) {
          # Some enrichplot versions handle the internal .sign split
          # differently. Fall back to the standard GSEA dotplot rather than
          # leaving an empty panel.
          enrichplot::dotplot(
            bundle$result,
            showCategory = top_n
          ) +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::theme(
              axis.text.y =
                ggplot2::element_text(
                  colour = "black"
                ),
              axis.text.x =
                ggplot2::element_text(
                  colour = "black"
                )
            )
        }
      )
    })

    gsea_curve_obj <- shiny::reactive({
      bundle <- gsea_bundle()
      shiny::req(
        bundle,
        is.null(bundle$error),
        bundle$result,
        nrow(bundle$result_table)
      )

      pathway <- selected_gsea_pathway()
      shiny::req(nzchar(pathway))

      df <- bundle$result_table
      gene_set_index <- match(
        pathway,
        as.character(df$ID)
      )
      shiny::validate(
        shiny::need(
          !is.na(gene_set_index),
          "The selected pathway is not present in the current GSEA result."
        )
      )

      title <- df$Description[[gene_set_index]]

      # Match 02.MaizeTeosintePro/01.src/RNAseq.R: the archived script used
      # gseaplot(..., by = "all", geneSetID = <result row index>).
      enrichplot::gseaplot(
        bundle$result,
        by = "all",
        geneSetID = as.integer(gene_set_index),
        title = title
      )
    })


    output$gsea_dotplot <- shiny::renderPlot({
      bundle <- gsea_bundle()
      shiny::validate(
        shiny::need(
          !is.null(bundle) &&
            is.null(bundle$error) &&
            nrow(bundle$result_table),
          if (!is.null(bundle$error)) {
            paste0("GSEA failed: ", bundle$error)
          } else {
            "Run GSEA to display the pathway dotplot."
          }
        )
      )
      print(gsea_dotplot_obj())
    }, res = 100)

    output$gsea_curve <- shiny::renderPlot({
      bundle <- gsea_bundle()
      shiny::validate(
        shiny::need(
          !is.null(bundle) &&
            is.null(bundle$error) &&
            nrow(bundle$result_table),
          if (!is.null(bundle$error)) {
            paste0("GSEA failed: ", bundle$error)
          } else {
            "Run GSEA to display the enrichment curve."
          }
        )
      )
      print(gsea_curve_obj())
    }, res = 100)

    output$gsea_table <- DT::renderDT({
      bundle <- gsea_bundle()
      shiny::req(
        bundle,
        is.null(bundle$error),
        nrow(bundle$result_table)
      )

      table <- bundle$result_table
      keep <- intersect(
        c(
          "ID",
          "Description",
          "setSize",
          "enrichmentScore",
          "NES",
          "pvalue",
          "p.adjust",
          "qvalues",
          "rank",
          "leading_edge",
          "core_enrichment"
        ),
        names(table)
      )
      table <- table[, keep, drop = FALSE]

      DT::datatable(
        table,
        extensions = c("Buttons", "Scroller"),
        options = list(
          pageLength = 10,
          dom = "Bfrtip",
          buttons = c(
            "copy",
            "csv",
            "excel"
          ),
          scrollX = TRUE,
          scrollY = 480,
          scroller = TRUE
        ),
        rownames = FALSE,
        class = "display compact"
      )
    })

    output$gsea_rank_table <- DT::renderDT({
      bundle <- gsea_bundle()
      shiny::req(
        bundle,
        is.null(bundle$error),
        length(bundle$ranks)
      )

      ranks <- data.frame(
        Rank = seq_along(bundle$ranks),
        GeneID = names(bundle$ranks),
        Rank_metric = as.numeric(bundle$ranks),
        stringsAsFactors = FALSE
      )

      DT::datatable(
        ranks,
        extensions = c("Buttons", "Scroller"),
        options = list(
          pageLength = 10,
          dom = "Bfrtip",
          buttons = c("copy", "csv"),
          scrollY = 480,
          scroller = TRUE
        ),
        rownames = FALSE,
        class = "display compact"
      )
    })

    output$gsea_reference_ui <- shiny::renderUI({
      ref <- .protvis_deg_load_gsea_reference()
      row <- ref[
        ref$ID == "00940",
        ,
        drop = FALSE
      ]

      shiny::div(
        class = "pv-transcriptome-gsea-reference",

        shiny::h5(
          "Built-in RNA-seq workflow"
        ),
        shiny::p(
          class = "pv-transcriptome-gsea-note",
          "This is a regression/reference record extracted from ",
          "02.MaizeTeosintePro/03.progress/04.RNAseq/KEGGenrich.xlsx. ",
          "It is not used as the active GSEA result. Uploading the original ",
          "expression and group tables should follow this workflow."
        ),

        shiny::tags$table(
          class = "table table-sm",
          shiny::tags$tbody(
            shiny::tags$tr(
              shiny::tags$th("Pathway"),
              shiny::tags$td(
                paste0(
                  row$Description[[1L]],
                  " [",
                  row$ID[[1L]],
                  "]"
                )
              )
            ),
            shiny::tags$tr(
              shiny::tags$th("Set size"),
              shiny::tags$td(row$setSize[[1L]])
            ),
            shiny::tags$tr(
              shiny::tags$th("ES"),
              shiny::tags$td(
                signif(
                  row$enrichmentScore[[1L]],
                  5
                )
              )
            ),
            shiny::tags$tr(
              shiny::tags$th("NES"),
              shiny::tags$td(
                signif(row$NES[[1L]], 5)
              )
            ),
            shiny::tags$tr(
              shiny::tags$th("P value"),
              shiny::tags$td(
                signif(row$pvalue[[1L]], 5)
              )
            ),
            shiny::tags$tr(
              shiny::tags$th("Adjusted P"),
              shiny::tags$td(
                signif(row$p.adjust[[1L]], 5)
              )
            ),
            shiny::tags$tr(
              shiny::tags$th("Rank"),
              shiny::tags$td(row$rank[[1L]])
            )
          )
        )
      )
    })

    # Render PCA data table (show pca_result$rotated)
    output$pca_data_table <- DT::renderDT({
      shiny::req(pca_rotated_data())
      DT::datatable(
        pca_rotated_data(),
        extensions = c('Buttons', 'Scroller'),
        options = list(
          pageLength = 10,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
          scrollX = TRUE,
          scrollY = 550,
          scroller = TRUE
        ),
        rownames = FALSE,
        class = 'display compact'
      )
    })
    # Render DEG result table (display res_tbl)
    output$deg_table <- DT::renderDT({
      shiny::req(deseq_results())
      res_tbl <- deseq_results() %>%
        dplyr::select(GeneID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, Regulation) %>%
        dplyr::mutate(
          dplyr::across(dplyr::where(is.numeric), ~ base::round(., 4)),
          padj = base::format(padj, scientific = TRUE, digits = 3)
        )
      DT::datatable(
        res_tbl,
        extensions = c('Buttons', 'Scroller'),
        options = list(
          pageLength = 10,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
          scrollX = TRUE,
          scrollY = 550,
          scroller = TRUE
        ),
        rownames = FALSE,
        class = 'display compact'
      ) %>%
        DT::formatStyle(
          'Regulation',
          backgroundColor = DT::styleEqual(
            c('Up-regulated', 'Down-regulated', 'Not significant'),
            c('#FFCCCC', '#CCE5FF', '#F2F2F2')
          )
        )
    })
    # Render DEG statistics table
    output$deg_stats <- shiny::renderTable({
      shiny::req(deg_stats())
      stats <- deg_stats()
      base::data.frame(
        Statistic = c("Total Genes", "Up-regulated", "Down-regulated",
                      "Total Significant", "Percentage Significant"),
        Value = c(
          stats$total_genes,
          paste(stats$up_regulated, "genes"),
          paste(stats$down_regulated, "genes"),
          paste(stats$significant, "genes"),
          paste(stats$percent_sig, "%")
        )
      )
    }, align = 'lr')
    # Render Top DEGs table
    output$top_degs_table <- DT::renderDT({
      shiny::req(top_degs())
      top_genes <- top_degs() %>%
        dplyr::select(GeneID, log2FoldChange, padj, Regulation) %>%
        dplyr::mutate(
          log2FoldChange = base::round(log2FoldChange, 3),
          padj = base::format(padj, scientific = TRUE, digits = 3)
        )
      DT::datatable(
        top_genes,
        extensions = c('Buttons', 'Scroller'),
        options = list(
          pageLength = 5,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
          scrollX = TRUE
        ),
        rownames = FALSE,
        class = 'display compact'
      ) %>%
        DT::formatStyle(
          'Regulation',
          backgroundColor = DT::styleEqual(
            c('Up-regulated', 'Down-regulated', 'Not significant'),
            c('#FFCCCC', '#CCE5FF', '#F2F2F2')
          )
        )
    })
    # Download PCA diagram
    output$download_pca <- shiny::downloadHandler(
      filename = function() {
        base::paste("PCA_plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        shiny::req(pca_plot_obj())
        grDevices::pdf(file, width = input$download_width_pca, height = input$download_height_pca)
        print(pca_plot_obj())
        grDevices::dev.off()
      }
    )
    # Download volcano map
    output$download_pdf <- shiny::downloadHandler(
      filename = function() {
        base::paste("volcano_plot_", base::Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        shiny::req(voc_plot_obj())
        grDevices::pdf(file, width = input$download_width_voc, height = input$download_height_voc)
        print(voc_plot_obj())
        grDevices::dev.off()
      }
    )
    output$download_pca_table <- shiny::downloadHandler(
      filename = function() {
        base::paste("pca_rotated_data_", base::Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        shiny::req(pca_rotated_data())
        utils::write.csv(pca_rotated_data(), file, row.names = FALSE)
      }
    )
    # Download PCA data (from the sidebar button)
    output$download_pca_data <- shiny::downloadHandler(
      filename = function() {
        base::paste("pca_rotated_data_", base::Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        shiny::req(pca_rotated_data())
        utils::write.csv(pca_rotated_data(), file, row.names = FALSE)
      }
    )
    # Download DEG data (res_tbl)-from the sidebar button
    output$download_deg_data <- shiny::downloadHandler(
      filename = function() {
        base::paste("deg_analysis_results_", base::Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        shiny::req(deseq_results())
        utils::write.csv(deseq_results(), file, row.names = FALSE)
      }
    )
    output$download_gsea_results <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "transcriptome_KEGG_GSEA_",
          base::Sys.Date(),
          ".csv"
        )
      },
      content = function(file) {
        bundle <- gsea_bundle()
        shiny::req(
          bundle,
          is.null(bundle$error),
          nrow(bundle$result_table)
        )
        utils::write.csv(
          bundle$result_table,
          file,
          row.names = FALSE
        )
      }
    )

    output$download_gsea_ranks <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "transcriptome_GSEA_ranked_list_",
          base::Sys.Date(),
          ".csv"
        )
      },
      content = function(file) {
        bundle <- gsea_bundle()
        shiny::req(
          bundle,
          is.null(bundle$error),
          length(bundle$ranks)
        )
        out <- data.frame(
          Rank = seq_along(bundle$ranks),
          GeneID = names(bundle$ranks),
          Rank_metric = as.numeric(bundle$ranks),
          stringsAsFactors = FALSE
        )
        utils::write.csv(
          out,
          file,
          row.names = FALSE
        )
      }
    )

    output$download_gsea_dotplot <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "transcriptome_GSEA_dotplot_",
          base::Sys.Date(),
          ".pdf"
        )
      },
      content = function(file) {
        shiny::req(gsea_dotplot_obj())
        grDevices::pdf(
          file,
          width = 10,
          height = 7
        )
        print(gsea_dotplot_obj())
        grDevices::dev.off()
      }
    )

    output$download_gsea_curve <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "transcriptome_GSEA_curve_",
          base::Sys.Date(),
          ".pdf"
        )
      },
      content = function(file) {
        shiny::req(gsea_curve_obj())
        grDevices::pdf(
          file,
          width = 9,
          height = 7
        )
        print(gsea_curve_obj())
        grDevices::dev.off()
      }
    )

    # Download DEG data (res_tbl)-from the button in the table
    output$download_degs <- shiny::downloadHandler(
      filename = function() {
        base::paste("deg_results_", base::Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        shiny::req(deseq_results())
        utils::write.csv(deseq_results(), file, row.names = FALSE)
      }
    )
  })
}
