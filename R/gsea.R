#' GSEA module UI
#' @import shiny
#' @import bslib
#' @name gsea_ui
#' @export
#'
gsea_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_sidebar(
    title = "GSEA Analysis",
    sidebar = bslib::sidebar(
      shiny::selectInput(
        ns("analysis_mode"),
        "Analysis mode",
        choices = c(
          "Standard GSEA" = "standard",
          "B73-Y12 Figure reproduction" = "archived"
        ),
        selected = "standard"
      ),

      shiny::conditionalPanel(
        condition = "input.analysis_mode === 'standard'",
        ns = ns,
        shiny::fileInput(
          ns("expr_file"),
          "Upload Expression Matrix",
          accept = ".txt"
        ),
        shiny::fileInput(
          ns("group_file"),
          "Upload Group Info",
          accept = ".txt"
        ),
        shiny::fileInput(
          ns("ko_file"),
          "Upload KO Pathway Info",
          accept = ".txt"
        ),
        shiny::numericInput(
          ns("minGSSize"),
          "minGSSize",
          value = 1,
          min = 1,
          step = 1
        ),
        shiny::numericInput(
          ns("maxGSSize"),
          "maxGSSize",
          value = 5000,
          min = 10,
          step = 10
        ),
        shiny::numericInput(
          ns("plot_top_x"),
          "Show Top X Pathways",
          value = 20,
          min = 1,
          step = 1
        ),
        shiny::checkboxInput(
          ns("sig_only"),
          "Show only significant (p.adjust < 0.05)",
          value = FALSE
        ),
        shiny::radioButtons(
          ns("sort_by"),
          "Dotplot Sort By",
          choices = c("NES", "p.adjust"),
          inline = TRUE
        ),
        shiny::selectInput(
          ns("pathway"),
          "Select a Pathway",
          choices = NULL
        )
      ),

      shiny::conditionalPanel(
        condition = "input.analysis_mode === 'archived'",
        ns = ns,
        shiny::selectInput(
          ns("archive_stage"),
          "Archived comparison",
          choices = c(
            "Root VE: B73 vs Y12" = "Root_VE"
          ),
          selected = "Root_VE"
        ),
        shiny::selectInput(
          ns("archive_pathway"),
          "Pathway",
          choices = c(
            "Phenylpropanoid biosynthesis" =
              "Phenylpropanoid biosynthesis"
          ),
          selected = "Phenylpropanoid biosynthesis"
        ),
        shiny::div(
          class = "alert alert-info py-2 small",
          shiny::tags$b("Figure-compatible workflow"),
          shiny::tags$br(),
          "Historical retained-protein rank list; B73 - Y12 log2FC ",
          "ranking; weighted GSEA (p = 1); fgseaMultilevel; ",
          "seed = 20260920; minSize = 5; maxSize = 500."
        )
      ),

      shiny::actionButton(
        ns("run"),
        "Run Analysis",
        class = "btn-primary w-100 mt-2",
        icon = shiny::icon("play")
      ),
      shiny::hr(),

      shiny::downloadButton(
        ns("download_csv"),
        "Download Results CSV",
        class = "btn-success w-100"
      ),

      shiny::conditionalPanel(
        condition = "input.analysis_mode === 'standard'",
        ns = ns,
        shiny::downloadButton(
          ns("download_dotplot"),
          "Download Dotplot PDF",
          class = "btn-success w-100 mt-1"
        )
      ),

      shiny::downloadButton(
        ns("download_curve"),
        "Download Curve PDF",
        class = "btn-success w-100 mt-1"
      )
    ),

    shiny::conditionalPanel(
      condition = "input.analysis_mode === 'standard'",
      ns = ns,
      bslib::card(
        bslib::card_header("GSEA Results Table"),
        bslib::card_body(
          DT::DTOutput(ns("gsea_table"))
        )
      ),
      bslib::layout_columns(
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Top Pathways Dotplot"),
          bslib::card_body(
            shiny::plotOutput(ns("dotplot"), height = "600px")
          )
        ),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Single Pathway Enrichment Curve"),
          bslib::card_body(
            shiny::plotOutput(ns("gsea_plot"), height = "600px")
          )
        )
      )
    ),

    shiny::conditionalPanel(
      condition = "input.analysis_mode === 'archived'",
      ns = ns,
      bslib::card(
        bslib::card_header(
          "Archived B73-Y12 GSEA reproduction"
        ),
        bslib::card_body(
          shiny::uiOutput(ns("archive_status"))
        )
      ),
      bslib::layout_columns(
        col_widths = c(4, 8),
        bslib::card(
          bslib::card_header("GSEA statistics"),
          bslib::card_body(
            DT::DTOutput(ns("archive_table"))
          )
        ),
        bslib::card(
          bslib::card_header("Reproduction metadata"),
          bslib::card_body(
            DT::DTOutput(ns("archive_metadata"))
          )
        )
      ),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header(
          "Phenylpropanoid biosynthesis"
        ),
        bslib::card_body(
          shiny::plotOutput(
            ns("archive_curve"),
            height = "650px"
          )
        )
      )
    )
  )
}


.protvis_gsea_base36 <- function(x) {
  alphabet <- c(as.character(0:9), letters)
  decode_one <- function(value) {
    value <- as.character(value)
    sign <- 1L
    if (startsWith(value, "-")) {
      sign <- -1L
      value <- substring(value, 2L)
    }
    chars <- strsplit(tolower(value), "", fixed = TRUE)[[1L]]
    digits <- match(chars, alphabet) - 1L
    if (!length(digits) || anyNA(digits)) {
      stop("Invalid base36 value in archived GSEA data.", call. = FALSE)
    }
    powers <- rev(seq_along(digits) - 1L)
    sign * sum(digits * (36 ^ powers))
  }
  vapply(x, decode_one, numeric(1))
}


.protvis_gsea_archive_path <- function() {
  candidates <- c(
    system.file(
      "extdata", "figure4_archive",
      "Root_VE_phenylpropanoid.pvg",
      package = "ProtVis"
    ),
    file.path(
      "inst", "extdata", "figure4_archive",
      "Root_VE_phenylpropanoid.pvg"
    ),
    file.path(
      getwd(), "inst", "extdata", "figure4_archive",
      "Root_VE_phenylpropanoid.pvg"
    )
  )
  candidates <- candidates[
    nzchar(candidates) & file.exists(candidates)
  ]
  if (!length(candidates)) {
    stop(
      "The bundled Root_VE GSEA reproduction archive is unavailable.",
      call. = FALSE
    )
  }
  candidates[[1L]]
}


.protvis_gsea_load_archived_rootve <- function() {
  lines <- readLines(
    .protvis_gsea_archive_path(),
    warn = FALSE
  )
  get_field <- function(name) {
    prefix <- paste0(name, "=")
    hit <- lines[startsWith(lines, prefix)]
    if (!length(hit)) {
      stop(
        "Missing field in archived GSEA data: ",
        name,
        call. = FALSE
      )
    }
    substring(hit[[1L]], nchar(prefix) + 1L)
  }

  scale <- as.numeric(get_field("scale"))
  n_expected <- as.integer(get_field("n"))
  deltas <- strsplit(
    get_field("rank_deltas"),
    ",",
    fixed = TRUE
  )[[1L]]
  delta_values <- .protvis_gsea_base36(deltas)

  ranks_int <- numeric(length(delta_values))
  ranks_int[[1L]] <- delta_values[[1L]]
  if (length(delta_values) > 1L) {
    for (i in 2:length(delta_values)) {
      ranks_int[[i]] <-
        ranks_int[[i - 1L]] - delta_values[[i]]
    }
  }
  rank_metric <- ranks_int / scale

  if (length(rank_metric) != n_expected) {
    stop(
      "Archived GSEA ranked-list length failed validation.",
      call. = FALSE
    )
  }

  hit_tokens <- strsplit(
    get_field("hits"),
    ",",
    fixed = TRUE
  )[[1L]]
  hit_parts <- strsplit(
    hit_tokens,
    ":",
    fixed = TRUE
  )
  position_delta <- .protvis_gsea_base36(
    vapply(hit_parts, `[[`, character(1), 1L)
  )
  hit_position <- cumsum(position_delta)
  hit_id_token <- vapply(
    hit_parts,
    `[[`,
    character(1),
    2L
  )

  decode_id <- function(token) {
    prefix <- substring(token, 1L, 1L)
    suffix <- substring(token, 2L)
    if (identical(prefix, "P")) {
      return(
        paste0(
          "PZ00001a",
          sprintf(
            "%06d",
            as.integer(.protvis_gsea_base36(suffix))
          )
        )
      )
    }
    if (identical(prefix, "Z")) {
      return(
        paste0(
          "Zm00001d",
          sprintf(
            "%06d",
            as.integer(.protvis_gsea_base36(suffix))
          )
        )
      )
    }
    if (identical(prefix, "X")) {
      return(suffix)
    }
    stop("Unknown archived protein-ID token.", call. = FALSE)
  }

  hit_id <- vapply(
    hit_id_token,
    decode_id,
    character(1)
  )
  hit_mask <- rep(FALSE, length(rank_metric))
  hit_mask[hit_position] <- TRUE

  list(
    stage = get_field("stage"),
    pathway = get_field("pathway"),
    term = get_field("term"),
    contrast = get_field("contrast"),
    rank_metric_name = get_field("rank_metric"),
    rank_metric = rank_metric,
    hit_position = as.integer(hit_position),
    hit_id = hit_id,
    hit_mask = hit_mask
  )
}


.protvis_gsea_prepare_archived <- function() {
  archive <- .protvis_gsea_load_archived_rootve()
  n_all <- length(archive$rank_metric)
  n_hit <- sum(archive$hit_mask)

  weights <- abs(archive$rank_metric)
  hit_weight <- sum(weights[archive$hit_mask])
  if (!is.finite(hit_weight) || hit_weight <= 0 ||
      n_all <= n_hit) {
    stop(
      "Archived GSEA data cannot produce a weighted running score.",
      call. = FALSE
    )
  }

  increment <- ifelse(
    archive$hit_mask,
    weights / hit_weight,
    -1 / (n_all - n_hit)
  )
  running_es <- cumsum(increment)
  extreme_position <- which.max(abs(running_es))

  rank_tbl <- data.frame(
    position = seq_len(n_all),
    rank_metric = archive$rank_metric,
    in_pathway = archive$hit_mask,
    running_ES = running_es,
    stringsAsFactors = FALSE
  )
  rank_tbl$Protein_ID <- NA_character_
  rank_tbl$Protein_ID[archive$hit_position] <-
    archive$hit_id

  list(
    archive = archive,
    rank_tbl = rank_tbl,
    extreme_position = as.integer(extreme_position),
    extreme_es = running_es[[extreme_position]]
  )
}


.protvis_gsea_run_archived <- function(seed = 20260920L) {
  prepared <- .protvis_gsea_prepare_archived()
  archive <- prepared$archive

  ranked_ids <- sprintf(
    "rank_%05d",
    seq_along(archive$rank_metric)
  )
  stats <- archive$rank_metric
  names(stats) <- ranked_ids
  pathway_ids <- ranked_ids[archive$hit_mask]

  set.seed(as.integer(seed))
  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = stats::setNames(
      list(pathway_ids),
      archive$pathway
    ),
    stats = stats,
    minSize = 5,
    maxSize = 500,
    eps = 0
  )
  fgsea_result <- as.data.frame(fgsea_result)
  if ("leadingEdge" %in% names(fgsea_result)) {
    fgsea_result$leadingEdge <- vapply(
      fgsea_result$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }

  metadata <- data.frame(
    Stage = archive$stage,
    Contrast = "B73 - Y12",
    Pathway = archive$pathway,
    KEGG_term = archive$term,
    Rank_metric = "historical B73 - Y12 log2FC",
    Ranked_proteins = length(archive$rank_metric),
    Pathway_members = sum(archive$hit_mask),
    Weighted_ES_extreme = prepared$extreme_es,
    Extreme_position = prepared$extreme_position,
    Seed = as.integer(seed),
    minSize = 5L,
    maxSize = 500L,
    stringsAsFactors = FALSE
  )

  c(
    prepared,
    list(
      fgsea_result = fgsea_result,
      metadata = metadata
    )
  )
}


.protvis_gsea_archived_plot <- function(result) {
  rank_tbl <- result$rank_tbl
  archive <- result$archive
  extreme_position <- result$extreme_position

  base_theme <-
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey45",
        linewidth = 0.45
      ),
      plot.margin = ggplot2::margin(4, 7, 4, 7)
    )

  rank_plot <-
    ggplot2::ggplot(
      rank_tbl,
      ggplot2::aes(x = position)
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey65",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = position,
        y = 0,
        yend = rank_metric
      ),
      colour = "grey20",
      linewidth = 0.22
    ) +
    ggplot2::labs(
      title = archive$pathway,
      x = NULL,
      y = "Ranked List Metric"
    ) +
    base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 17
      )
    )

  es_range <- range(rank_tbl$running_ES, 0)
  tick_height <- diff(es_range) * 0.07

  es_plot <-
    ggplot2::ggplot(
      rank_tbl,
      ggplot2::aes(
        x = position,
        y = running_ES
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 0.45
    ) +
    ggplot2::geom_vline(
      xintercept = extreme_position,
      colour = "#F28E8E",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    ggplot2::geom_segment(
      data = rank_tbl[rank_tbl$in_pathway, , drop = FALSE],
      ggplot2::aes(
        x = position,
        xend = position,
        y = 0,
        yend = -tick_height
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.32
    ) +
    ggplot2::geom_line(
      colour = "#00D800",
      linewidth = 0.85
    ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = es_range[[1L]] + diff(es_range) * 0.08,
      hjust = 0,
      label = "Zea mays ssp. mays",
      colour = "red",
      size = 5.2,
      fontface = "italic"
    ) +
    ggplot2::annotate(
      "text",
      x = nrow(rank_tbl),
      y = es_range[[2L]] - diff(es_range) * 0.08,
      hjust = 1,
      label = "Zea mays ssp. mexicana",
      colour = "blue",
      size = 5.2,
      fontface = "italic"
    ) +
    ggplot2::labs(
      x = "Position in the Ranked List of Genes",
      y = "Running Enrichment Score"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0.01, 0.01)
      )
    ) +
    base_theme

  patchwork::wrap_plots(
    rank_plot,
    es_plot,
    ncol = 1,
    heights = c(1, 1)
  )
}


#' GSEA module server
#' @import shiny
#' @name gsea_server
#' @export
#'
gsea_server <- function(id, shared_state = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    gsea_res_val <- shiny::reactiveVal()
    top_df_val <- shiny::reactiveVal()
    archived_val <- shiny::reactiveVal()

    shiny::observeEvent(input$run, {
      mode <- input$analysis_mode %||% "standard"

      if (identical(mode, "archived")) {
        result <- tryCatch(
          .protvis_gsea_run_archived(),
          error = function(e) e
        )
        if (inherits(result, "error")) {
          archived_val(NULL)
          shiny::showNotification(
            conditionMessage(result),
            type = "error",
            duration = 8
          )
          return()
        }

        archived_val(result)

        .protvis_record_shared_run(
          shared_state,
          module = "gsea",
          method = "archived_B73_Y12_weighted_GSEA",
          category = "enrichment",
          parameters = list(
            stage = result$archive$stage,
            contrast = "B73 - Y12",
            pathway = result$archive$pathway,
            rank_metric = "historical B73 - Y12 log2FC",
            weighted_p = 1,
            seed = 20260920L,
            minSize = 5L,
            maxSize = 500L
          ),
          tables = list(
            gsea_results = result$fgsea_result,
            metadata = result$metadata,
            ranked_background = result$rank_tbl
          ),
          plot_data = list(
            enrichment_curve = result$rank_tbl
          ),
          plot_config = list(
            extreme_position = result$extreme_position,
            left_label = "Zea mays ssp. mays",
            right_label = "Zea mays ssp. mexicana"
          )
        )
        return()
      }

      shiny::req(
        input$expr_file,
        input$group_file,
        input$ko_file
      )

      expr <- data.table::fread(
        input$expr_file$datapath
      )
      group <- data.table::fread(
        input$group_file$datapath
      )
      ko <- data.table::fread(
        input$ko_file$datapath
      )

      mat <- base::as.matrix(expr[, -1])
      base::rownames(mat) <- expr$ID
      group_df <- base::data.frame(
        row.names = group$Sample,
        condition = group$Group
      )

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = mat,
        colData = group_df,
        design = ~ condition
      )
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::results(dds)
      res <- res[!is.na(res$log2FoldChange), ]

      gene_ranks <- res$log2FoldChange
      base::names(gene_ranks) <- base::rownames(res)
      gene_ranks <- base::sort(
        gene_ranks,
        decreasing = TRUE
      )

      term2gene <- data.table::rbindlist(
        base::lapply(
          seq_len(base::nrow(ko)),
          function(i) {
            genes <- base::unlist(
              base::strsplit(
                ko[i, 3][[1]],
                ";"
              )
            )
            data.table::data.table(
              term = ko[i, 1][[1]],
              gene = genes
            )
          }
        )
      )

      gsea_res <- clusterProfiler::GSEA(
        geneList = gene_ranks,
        TERM2GENE = term2gene,
        pvalueCutoff = 1,
        minGSSize = input$minGSSize,
        maxGSSize = input$maxGSSize
      )
      gsea_res_val(gsea_res)

      df <- gsea_res@result
      if (isTRUE(input$sig_only)) {
        df <- df[df$p.adjust < 0.05, , drop = FALSE]
      }
      if (identical(input$sort_by, "NES")) {
        df <- df[order(-abs(df$NES)), , drop = FALSE]
      } else {
        df <- df[order(df$p.adjust), , drop = FALSE]
      }
      top_df <- df[
        seq_len(min(nrow(df), input$plot_top_x)),
        ,
        drop = FALSE
      ]
      top_df_val(top_df)

      if (nrow(top_df) > 0) {
        shiny::updateSelectInput(
          session,
          "pathway",
          choices = top_df$Description
        )
      }

      .protvis_record_shared_run(
        shared_state,
        module = "gsea",
        method = "clusterProfiler_GSEA",
        category = "enrichment",
        parameters = list(
          minGSSize = input$minGSSize,
          maxGSSize = input$maxGSSize,
          significant_only = isTRUE(input$sig_only),
          sort_by = input$sort_by,
          top_n = input$plot_top_x
        ),
        tables = list(
          gsea_results =
            as.data.frame(gsea_res@result),
          top_pathways =
            as.data.frame(top_df),
          gene_ranks = data.frame(
            gene = names(gene_ranks),
            rank = as.numeric(gene_ranks),
            stringsAsFactors = FALSE
          )
        ),
        plot_data = list(
          top_pathways = as.data.frame(top_df)
        ),
        plot_config = list(
          selected_pathway =
            input$pathway %||% NA_character_
        )
      )
    })

    output$gsea_table <- DT::renderDT({
      top_df_val()
    })

    output$dotplot <- shiny::renderPlot({
      gsea_res <- gsea_res_val()
      top_df <- top_df_val()
      shiny::req(gsea_res, top_df)
      if (nrow(top_df) > 0) {
        sub_res <- gsea_res
        sub_res@result <- top_df
        print(
          clusterProfiler::dotplot(
            sub_res,
            showCategory = nrow(top_df)
          )
        )
      }
    })

    output$gsea_plot <- shiny::renderPlot({
      gsea_res <- gsea_res_val()
      shiny::req(gsea_res, input$pathway)
      print(
        GseaVis::gseaNb(
          object = gsea_res,
          geneSetID = input$pathway,
          subPlot = 3
        )
      )
    })

    output$archive_status <- shiny::renderUI({
      result <- archived_val()
      if (is.null(result)) {
        return(
          shiny::span(
            "Select the archived reproduction mode and click Run Analysis.",
            class = "text-muted"
          )
        )
      }
      fg <- result$fgsea_result
      shiny::tagList(
        shiny::span(
          class = "text-success fw-semibold",
          "✓ Archived Root_VE rank list reproduced"
        ),
        shiny::tags$br(),
        shiny::span(
          paste0(
            "Ranked proteins: ",
            format(
              nrow(result$rank_tbl),
              big.mark = ","
            ),
            " · pathway hits: ",
            sum(result$rank_tbl$in_pathway),
            " · extreme position: ",
            result$extreme_position,
            " · ES: ",
            signif(result$extreme_es, 4),
            if (nrow(fg)) {
              paste0(
                " · NES: ",
                signif(fg$NES[[1L]], 4),
                " · padj: ",
                signif(fg$padj[[1L]], 4)
              )
            } else {
              ""
            }
          ),
          class = "text-muted"
        )
      )
    })

    output$archive_table <- DT::renderDT({
      result <- archived_val()
      shiny::req(result)
      DT::datatable(
        result$fgsea_result,
        rownames = FALSE,
        options = list(
          dom = "t",
          scrollX = TRUE
        )
      )
    })

    output$archive_metadata <- DT::renderDT({
      result <- archived_val()
      shiny::req(result)
      meta <- data.frame(
        Parameter = names(result$metadata),
        Value = vapply(
          result$metadata,
          function(x) as.character(x[[1L]]),
          character(1)
        ),
        stringsAsFactors = FALSE
      )
      DT::datatable(
        meta,
        rownames = FALSE,
        options = list(
          dom = "t",
          pageLength = nrow(meta)
        )
      )
    })

    output$archive_curve <- shiny::renderPlot({
      result <- archived_val()
      shiny::req(result)
      print(.protvis_gsea_archived_plot(result))
    }, res = 110)

    output$download_csv <- shiny::downloadHandler(
      filename = function() {
        if (identical(
          input$analysis_mode,
          "archived"
        )) {
          return(
            paste0(
              "GSEA_Root_VE_phenylpropanoid_",
              base::Sys.Date(),
              ".csv"
            )
          )
        }
        paste0(
          "GSEA_result_",
          base::Sys.Date(),
          ".csv"
        )
      },
      content = function(file) {
        if (identical(
          input$analysis_mode,
          "archived"
        )) {
          result <- archived_val()
          shiny::req(result)
          utils::write.csv(
            result$fgsea_result,
            file,
            row.names = FALSE
          )
          return()
        }

        res <- gsea_res_val()
        shiny::req(res)
        utils::write.csv(
          res@result,
          file,
          row.names = FALSE
        )
      }
    )

    output$download_dotplot <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "dotplot_",
          base::Sys.Date(),
          ".pdf"
        )
      },
      content = function(file) {
        gsea_res <- gsea_res_val()
        top_df <- top_df_val()
        shiny::req(gsea_res, top_df)
        sub_res <- gsea_res
        sub_res@result <- top_df
        grDevices::pdf(
          file,
          width = 8,
          height = 6
        )
        print(
          clusterProfiler::dotplot(
            sub_res,
            showCategory = nrow(top_df)
          )
        )
        grDevices::dev.off()
      }
    )

    output$download_curve <- shiny::downloadHandler(
      filename = function() {
        if (identical(
          input$analysis_mode,
          "archived"
        )) {
          return(
            paste0(
              "GSEA_Root_VE_phenylpropanoid_",
              base::Sys.Date(),
              ".pdf"
            )
          )
        }
        paste0(
          "gsea_curve_",
          base::Sys.Date(),
          ".pdf"
        )
      },
      content = function(file) {
        if (identical(
          input$analysis_mode,
          "archived"
        )) {
          result <- archived_val()
          shiny::req(result)
          grDevices::pdf(
            file,
            width = 8.2,
            height = 5
          )
          print(
            .protvis_gsea_archived_plot(result)
          )
          grDevices::dev.off()
          return()
        }

        gsea_res <- gsea_res_val()
        shiny::req(
          gsea_res,
          input$pathway
        )
        grDevices::pdf(
          file,
          width = 8,
          height = 6
        )
        print(
          GseaVis::gseaNb(
            object = gsea_res,
            geneSetID = input$pathway,
            subPlot = 3
          )
        )
        grDevices::dev.off()
      }
    )
  })
}
