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
        shiny::actionButton(
          ns("load_previous_dep"),
          "LOAD PREVIOUS DEP RESULTS",
          icon = shiny::icon("folder-open"),
          class = "btn-outline-primary w-100"
        ),
        shiny::uiOutput(ns("dep_load_status")),
        shiny::hr(),
        shiny::tags$h6("DEP comparisons"),
        shiny::uiOutput(ns("dep_comparisons_ui")),
        shiny::div(
          class = "d-flex gap-2 mb-2",
          shiny::actionLink(
            ns("dep_select_all"),
            "Select all"
          ),
          shiny::tags$span("·"),
          shiny::actionLink(
            ns("dep_clear_all"),
            "Clear"
          )
        ),
        shiny::selectInput(
          ns("focus_comparison"),
          "Focus comparison",
          choices = NULL
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
        shiny::selectInput(
          ns("archive_rank_metric"),
          "Ranking metric",
          choices = c(
            "Moderated t statistic (script default)" = "t",
            "log2 fold change" = "logFC",
            "Signed -log10(P-value)" = "signed_log10_p"
          ),
          selected = "t"
        ),
        shiny::radioButtons(
          ns("archive_filter_mode"),
          "Protein inclusion",
          choices = c(
            "All tested proteins (GSEA default)" = "all",
            "Custom thresholds" = "custom"
          ),
          selected = "all"
        ),
        shiny::conditionalPanel(
          condition = "input.archive_filter_mode === 'custom'",
          ns = ns,
          shiny::selectInput(
            ns("archive_p_metric"),
            "P-value filter",
            choices = c(
              "None" = "none",
              "Raw P-value" = "P.Value",
              "BH-adjusted P-value" = "adj.P.Val"
            ),
            selected = "none"
          ),
          shiny::numericInput(
            ns("archive_p_cutoff"),
            "P-value cutoff",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),
          shiny::numericInput(
            ns("archive_abs_logfc"),
            "Minimum |log2FC|",
            value = 0,
            min = 0,
            step = 0.1
          )
        ),
        shiny::numericInput(
          ns("archive_seed"),
          "Random seed",
          value = 20260920,
          min = 1,
          max = .Machine$integer.max,
          step = 1
        ),
        shiny::div(
          class = "alert alert-info py-2 small",
          shiny::tags$b("Figure-reproduction defaults"),
          shiny::tags$br(),
          "Load all available DEP comparisons, then focus on Root_VE. ",
          "The original script uses the complete Root_VE tested-protein ",
          "result, moderated t ranking, no significance pre-filter, ",
          "weighted GSEA (p = 1), fgseaMultilevel, minSize = 5, ",
          "maxSize = 500 and eps = 0."
        )
      ),

      shiny::actionButton(
        ns("run"),      shiny::actionButton(
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
          "Previous DEP → GSEA"
        ),
        bslib::card_body(
          shiny::uiOutput(ns("archive_status"))
        )
      ),
      bslib::layout_columns(
        col_widths = c(7, 5),
        bslib::card(
          bslib::card_header(
            "Selected-comparison GSEA summary"
          ),
          bslib::card_body(
            DT::DTOutput(ns("archive_table"))
          )
        ),
        bslib::card(
          bslib::card_header(
            "Focus comparison metadata"
          ),
          bslib::card_body(
            DT::DTOutput(ns("archive_metadata"))
          )
        )
      ),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header(
          shiny::uiOutput(ns("archive_curve_title"))
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



.protvis_gsea_dep_id_column <- function(x) {
  candidates <- c("ID", "Protein_ID", "protein_id", "Protein", "Gene")
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) {
    stop("DEP result does not contain a protein ID column.", call. = FALSE)
  }
  hit[[1L]]
}


.protvis_gsea_rootve_comparison <- function(comparisons) {
  comparisons <- as.character(comparisons)
  if (!length(comparisons)) return(NA_character_)
  exact <- comparisons[
    comparisons == "B73_Root_VE_vs_Y12_Root_VE"
  ]
  if (length(exact)) return(exact[[1L]])
  root <- comparisons[
    grepl("B73_Root_VE", comparisons, fixed = TRUE) &
      grepl("Y12_Root_VE", comparisons, fixed = TRUE)
  ]
  if (length(root)) return(root[[1L]])
  comparisons[[1L]]
}


.protvis_gsea_dep_compatibility <- function(result) {
  result <- as.data.frame(
    result,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  get_unique <- function(column) {
    if (!column %in% names(result)) return(character())
    out <- unique(as.character(result[[column]]))
    out[!is.na(out) & nzchar(out)]
  }
  mode <- get_unique("analysis_mode")
  universe <- get_unique("protein_universe")
  matrix_source <- get_unique("matrix_source")
  test_method <- get_unique("test_method")

  exact <- length(mode) &&
    all(mode == "archived") &&
    length(universe) &&
    all(universe == "archived_any_detected") &&
    length(matrix_source) &&
    all(matrix_source == "Step6_data_normalization") &&
    length(test_method) &&
    all(test_method == "archived_eBayes")

  list(
    exact = isTRUE(exact),
    label = if (isTRUE(exact)) {
      "Archived DEP compatible with the historical limma workflow"
    } else {
      paste0(
        "Loaded DEP is not the archived historical workflow; ",
        "GSEA is valid for the loaded result but exact figure ",
        "reproduction is not guaranteed."
      )
    }
  )
}


.protvis_gsea_prepare_dep_rank <- function(
    dep_result,
    rank_metric = c("t", "logFC", "signed_log10_p"),
    filter_mode = c("all", "custom"),
    p_metric = c("none", "P.Value", "adj.P.Val"),
    p_cutoff = 0.05,
    abs_logfc = 0) {
  rank_metric <- match.arg(rank_metric)
  filter_mode <- match.arg(filter_mode)
  p_metric <- match.arg(p_metric)

  df <- as.data.frame(
    dep_result,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  id_col <- .protvis_gsea_dep_id_column(df)
  ids <- trimws(as.character(df[[id_col]]))

  required <- switch(
    rank_metric,
    t = "t",
    logFC = "logFC",
    signed_log10_p = c("logFC", "P.Value")
  )
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(
      "Ranking metric requires missing DEP column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (identical(filter_mode, "custom")) {
    if (!"logFC" %in% names(df)) {
      stop("Custom filtering requires logFC.", call. = FALSE)
    }
    keep <- is.finite(suppressWarnings(as.numeric(df$logFC))) &
      abs(suppressWarnings(as.numeric(df$logFC))) >=
        as.numeric(abs_logfc)

    if (!identical(p_metric, "none")) {
      if (!p_metric %in% names(df)) {
        stop(
          "Custom filtering requires ", p_metric, ".",
          call. = FALSE
        )
      }
      pv <- suppressWarnings(as.numeric(df[[p_metric]]))
      keep <- keep & is.finite(pv) & pv <= as.numeric(p_cutoff)
    }
    df <- df[keep, , drop = FALSE]
    ids <- ids[keep]
  }

  metric <- switch(
    rank_metric,
    t = suppressWarnings(as.numeric(df$t)),
    logFC = suppressWarnings(as.numeric(df$logFC)),
    signed_log10_p = {
      lfc <- suppressWarnings(as.numeric(df$logFC))
      pv <- suppressWarnings(as.numeric(df$P.Value))
      sign(lfc) * -log10(pmax(pv, .Machine$double.xmin))
    }
  )

  keep <- !is.na(ids) & nzchar(ids) & is.finite(metric)
  df <- df[keep, , drop = FALSE]
  ids <- ids[keep]
  metric <- metric[keep]

  # DEP tables are complete comparison results. Keep one row per ID, then
  # perform the GSEA-specific ranking. This mirrors the source script's
  # distinct(Protein_ID) -> arrange(desc(t)) behavior.
  duplicate <- duplicated(ids)
  if (any(duplicate)) {
    df <- df[!duplicate, , drop = FALSE]
    ids <- ids[!duplicate]
    metric <- metric[!duplicate]
  }

  ord <- order(metric, decreasing = TRUE, na.last = NA)
  df <- df[ord, , drop = FALSE]
  ids <- ids[ord]
  metric <- metric[ord]

  rank_vector <- metric
  names(rank_vector) <- ids

  list(
    rank_vector = rank_vector,
    table = df,
    n_after_filter = length(rank_vector)
  )
}


.protvis_gsea_kegg_pathway <- function(
    kegg_background,
    rank_vector,
    pathway_name = "Phenylpropanoid biosynthesis") {
  bg <- as.data.frame(
    kegg_background,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("TERM", "GENE", "NAME") %in% names(bg))) {
    stop(
      "Built-in KEGG background must contain TERM, GENE and NAME.",
      call. = FALSE
    )
  }

  bg$TERM <- trimws(as.character(bg$TERM))
  bg$GENE <- trimws(as.character(bg$GENE))
  bg$NAME <- trimws(as.character(bg$NAME))
  bg <- bg[
    nzchar(bg$TERM) & nzchar(bg$GENE) & nzchar(bg$NAME),
    ,
    drop = FALSE
  ]

  pathway_rows <- bg[
    bg$TERM %in% c("map00940", "00940") |
      grepl(
        pathway_name,
        bg$NAME,
        ignore.case = TRUE,
        fixed = TRUE
      ),
    ,
    drop = FALSE
  ]
  if (!nrow(pathway_rows)) {
    stop(
      "Phenylpropanoid biosynthesis (map00940) was not found ",
      "in the built-in KEGG background.",
      call. = FALSE
    )
  }

  term <- pathway_rows$TERM[[1L]]
  canonical_name <- pathway_rows$NAME[[1L]]
  annotated_ids <- unique(bg$GENE)
  background_ids <- intersect(
    names(rank_vector),
    annotated_ids
  )
  rank_vector <- sort(
    rank_vector[background_ids],
    decreasing = TRUE
  )
  pathway_genes <- intersect(
    unique(bg$GENE[bg$TERM == term]),
    names(rank_vector)
  )

  if (length(pathway_genes) < 5L) {
    stop(
      "Fewer than five pathway proteins overlap the tested ",
      "GSEA background; check the DEP source and protein IDs.",
      call. = FALSE
    )
  }

  list(
    term = term,
    name = canonical_name,
    rank_vector = rank_vector,
    pathway_genes = pathway_genes,
    background_size = length(rank_vector)
  )
}


.protvis_gsea_run_dep_comparison <- function(
    dep_result,
    comparison,
    kegg_background,
    pathway_name = "Phenylpropanoid biosynthesis",
    rank_metric = "t",
    filter_mode = "all",
    p_metric = "none",
    p_cutoff = 0.05,
    abs_logfc = 0,
    seed = 20260920L) {
  ranked <- .protvis_gsea_prepare_dep_rank(
    dep_result,
    rank_metric = rank_metric,
    filter_mode = filter_mode,
    p_metric = p_metric,
    p_cutoff = p_cutoff,
    abs_logfc = abs_logfc
  )

  pathway <- .protvis_gsea_kegg_pathway(
    kegg_background,
    ranked$rank_vector,
    pathway_name = pathway_name
  )
  rank_vector <- pathway$rank_vector
  pathway_genes <- pathway$pathway_genes

  set.seed(as.integer(seed))
  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = stats::setNames(
      list(pathway_genes),
      pathway$name
    ),
    stats = rank_vector,
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

  hit <- names(rank_vector) %in% pathway_genes
  n_all <- length(rank_vector)
  n_hit <- sum(hit)
  weights <- abs(as.numeric(rank_vector))
  hit_weight <- sum(weights[hit])
  if (!is.finite(hit_weight) || hit_weight <= 0 ||
      n_all <= n_hit) {
    stop(
      "The selected DEP result cannot produce a weighted GSEA curve.",
      call. = FALSE
    )
  }

  increment <- ifelse(
    hit,
    weights / hit_weight,
    -1 / (n_all - n_hit)
  )
  running_es <- cumsum(increment)
  extreme_position <- which.max(abs(running_es))

  rank_tbl <- data.frame(
    position = seq_len(n_all),
    Protein_ID = names(rank_vector),
    rank_metric = as.numeric(rank_vector),
    in_pathway = hit,
    running_ES = running_es,
    stringsAsFactors = FALSE
  )

  result_df <- as.data.frame(
    dep_result,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  group1 <- if ("Group1" %in% names(result_df)) {
    unique(as.character(result_df$Group1))
  } else {
    character()
  }
  group2 <- if ("Group2" %in% names(result_df)) {
    unique(as.character(result_df$Group2))
  } else {
    character()
  }
  group1 <- group1[!is.na(group1) & nzchar(group1)]
  group2 <- group2[!is.na(group2) & nzchar(group2)]
  group1 <- if (length(group1)) group1[[1L]] else "Group1"
  group2 <- if (length(group2)) group2[[1L]] else "Group2"

  compatibility <- .protvis_gsea_dep_compatibility(dep_result)
  metadata <- data.frame(
    Comparison = comparison,
    Group1 = group1,
    Group2 = group2,
    Pathway = pathway$name,
    KEGG_term = pathway$term,
    Rank_metric = rank_metric,
    Protein_filter = filter_mode,
    P_filter = p_metric,
    P_cutoff = if (identical(p_metric, "none")) {
      NA_real_
    } else {
      as.numeric(p_cutoff)
    },
    Min_abs_log2FC = as.numeric(abs_logfc),
    Ranked_proteins = length(rank_vector),
    Pathway_members = n_hit,
    Extreme_position = as.integer(extreme_position),
    Extreme_ES = running_es[[extreme_position]],
    Seed = as.integer(seed),
    Historical_DEP_compatible = compatibility$exact,
    stringsAsFactors = FALSE
  )

  list(
    archive = list(
      stage = .protvis_dep_stage_label(group1, group2),
      pathway = pathway$name,
      term = pathway$term,
      contrast = paste(group1, "-", group2),
      rank_metric_name = rank_metric
    ),
    rank_tbl = rank_tbl,
    extreme_position = as.integer(extreme_position),
    extreme_es = running_es[[extreme_position]],
    fgsea_result = fgsea_result,
    metadata = metadata,
    compatibility = compatibility
  )
}


.protvis_gsea_combined_summary <- function(results) {
  if (!length(results)) return(data.frame())
  rows <- lapply(names(results), function(comparison) {
    item <- results[[comparison]]
    fg <- item$fgsea_result
    meta <- item$metadata
    if (is.null(fg) || !nrow(fg)) {
      return(data.frame(
        Comparison = comparison,
        Pathway = meta$Pathway[[1L]],
        ES = NA_real_,
        NES = NA_real_,
        pval = NA_real_,
        padj = NA_real_,
        size = meta$Pathway_members[[1L]],
        Ranked_proteins = meta$Ranked_proteins[[1L]],
        Historical_DEP_compatible =
          meta$Historical_DEP_compatible[[1L]],
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      Comparison = comparison,
      Pathway = as.character(fg$pathway[[1L]]),
      ES = as.numeric(fg$ES[[1L]]),
      NES = as.numeric(fg$NES[[1L]]),
      pval = as.numeric(fg$pval[[1L]]),
      padj = as.numeric(fg$padj[[1L]]),
      size = as.integer(fg$size[[1L]]),
      Ranked_proteins = meta$Ranked_proteins[[1L]],
      Historical_DEP_compatible =
        meta$Historical_DEP_compatible[[1L]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' GSEA module server#' GSEA module server
#' @import shiny
#' @name gsea_server
#' @export
#'
gsea_server <- function(id, shared_state = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    gsea_res_val <- shiny::reactiveVal()
    top_df_val <- shiny::reactiveVal()
    archived_val <- shiny::reactiveVal()
    previous_dep_val <- shiny::reactiveVal(
      list(
        results = list(),
        provenance = list(),
        source = NULL
      )
    )

    load_previous_dep_results <- function() {
      results <- shared_state$dep_results %||% list()
      dataset <- shared_state$dataset
      provenance <- list()
      source <- NULL

      if (length(results)) {
        source <- "Current session DEP results"
        if (inherits(dataset, "ProtVis_dataset")) {
          provenance <-
            dataset$analysis_results$DEP$provenance %||% list()
        }
      }

      if (!length(results) &&
          inherits(dataset, "ProtVis_dataset")) {
        dep <- dataset$analysis_results$DEP %||% list()
        results <- dep$results %||% list()
        provenance <- dep$provenance %||% list()
        if (length(results)) {
          source <- "Current ProtVis_dataset"
        }
      }

      if (!length(results)) {
        workdir <- as.character(shared_state$workdir %||% "")
        step7 <- if (nzchar(workdir)) {
          file.path(
            workdir,
            "Step7_differential_analysis.rda"
          )
        } else {
          ""
        }
        if (nzchar(step7) && file.exists(step7)) {
          dataset7 <- .protvis_load_stage_dataset(step7)
          if (!is.null(dataset7)) {
            dep <- dataset7$analysis_results$DEP %||% list()
            results <- dep$results %||% list()
            provenance <- dep$provenance %||% list()
            if (length(results)) {
              source <- "Step7_differential_analysis.rda"
            }
          }
        }
      }

      if (!length(results)) {
        stop(
          "No previous DEP results were found. Run Differential ",
          "Analysis first or load a project containing ",
          "Step7_differential_analysis.rda.",
          call. = FALSE
        )
      }

      list(
        results = results,
        provenance = provenance,
        source = source %||% "Previous DEP results"
      )
    }

    shiny::observeEvent(input$load_previous_dep, {
      bundle <- tryCatch(
        load_previous_dep_results(),
        error = function(e) e
      )
      if (inherits(bundle, "error")) {
        previous_dep_val(list(
          results = list(),
          provenance = list(),
          source = NULL
        ))
        archived_val(NULL)
        shiny::showNotification(
          conditionMessage(bundle),
          type = "error",
          duration = 8
        )
        return()
      }

      previous_dep_val(bundle)
      comparisons <- names(bundle$results)
      rootve <- .protvis_gsea_rootve_comparison(comparisons)

      shiny::updateSelectInput(
        session,
        "focus_comparison",
        choices = comparisons,
        selected = rootve
      )

      archived_val(NULL)
      shiny::showNotification(
        paste0(
          "Loaded ", length(comparisons),
          " DEP comparison(s). Root_VE is selected as the focus when available."
        ),
        type = "message",
        duration = 5
      )
    }, ignoreInit = TRUE)

    output$dep_load_status <- shiny::renderUI({
      bundle <- previous_dep_val()
      results <- bundle$results %||% list()
      if (!length(results)) {
        return(
          shiny::div(
            class = "text-muted small mt-2",
            "No previous DEP results loaded."
          )
        )
      }

      rootve <- .protvis_gsea_rootve_comparison(names(results))
      root_status <- .protvis_gsea_dep_compatibility(
        results[[rootve]]
      )
      shiny::div(
        class = "small mt-2",
        shiny::span(
          class = "text-success fw-semibold",
          paste0(
            "✓ ", length(results),
            " comparison(s) loaded"
          )
        ),
        shiny::tags$br(),
        shiny::span(
          bundle$source %||% "Previous DEP results",
          class = "text-muted"
        ),
        shiny::tags$br(),
        shiny::span(
          root_status$label,
          class = if (isTRUE(root_status$exact)) {
            "text-success"
          } else {
            "text-warning"
          }
        )
      )
    })

    output$dep_comparisons_ui <- shiny::renderUI({
      comparisons <- names(
        previous_dep_val()$results %||% list()
      )
      if (!length(comparisons)) {
        return(
          shiny::helpText(
            "Load previous DEP results to list comparisons."
          )
        )
      }
      selected <- isolate(input$dep_comparisons)
      selected <- intersect(
        selected %||% comparisons,
        comparisons
      )
      if (!length(selected)) selected <- comparisons

      shiny::checkboxGroupInput(
        session$ns("dep_comparisons"),
        label = NULL,
        choices = comparisons,
        selected = selected,
        inline = FALSE
      )
    })

    shiny::observeEvent(input$dep_select_all, {
      comparisons <- names(
        previous_dep_val()$results %||% list()
      )
      shiny::updateCheckboxGroupInput(
        session,
        "dep_comparisons",
        choices = comparisons,
        selected = comparisons,
        inline = FALSE
      )
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$dep_clear_all, {
      comparisons <- names(
        previous_dep_val()$results %||% list()
      )
      shiny::updateCheckboxGroupInput(
        session,
        "dep_comparisons",
        choices = comparisons,
        selected = character(),
        inline = FALSE
      )
    }, ignoreInit = TRUE)

    shiny::observe({
      comparisons <- names(
        previous_dep_val()$results %||% list()
      )
      if (!length(comparisons)) return()

      selected <- intersect(
        input$dep_comparisons %||% comparisons,
        comparisons
      )
      if (!length(selected)) selected <- comparisons

      current <- input$focus_comparison %||% ""
      rootve <- .protvis_gsea_rootve_comparison(selected)
      focus <- if (current %in% selected) {
        current
      } else {
        rootve
      }

      shiny::updateSelectInput(
        session,
        "focus_comparison",
        choices = selected,
        selected = focus
      )
    })

    shiny::observeEvent(input$run, {
      mode <- input$analysis_mode %||% "standard"

      if (identical(mode, "archived")) {
        bundle <- previous_dep_val()
        dep_results <- bundle$results %||% list()

        if (!length(dep_results)) {
          shiny::showNotification(
            paste0(
              "Load previous DEP results first. Multiple comparisons ",
              "can be loaded; Root_VE will be used as the default focus."
            ),
            type = "error",
            duration = 7
          )
          return()
        }

        selected <- intersect(
          input$dep_comparisons %||% character(),
          names(dep_results)
        )
        if (!length(selected)) {
          shiny::showNotification(
            "Select at least one DEP comparison.",
            type = "error",
            duration = 6
          )
          return()
        }

        focus <- input$focus_comparison %||%
          .protvis_gsea_rootve_comparison(selected)
        if (!focus %in% selected) {
          focus <- .protvis_gsea_rootve_comparison(selected)
        }

        built_in <- tryCatch(
          .protvis_load_builtin_enrichment_background(),
          error = function(e) e
        )
        if (inherits(built_in, "error")) {
          shiny::showNotification(
            conditionMessage(built_in),
            type = "error",
            duration = 8
          )
          return()
        }

        seed <- suppressWarnings(
          as.integer(input$archive_seed %||% 20260920L)
        )
        if (!is.finite(seed) || seed < 1L) {
          seed <- 20260920L
        }

        results <- list()
        errors <- character()
        for (comparison in selected) {
          item <- tryCatch(
            .protvis_gsea_run_dep_comparison(
              dep_result = dep_results[[comparison]],
              comparison = comparison,
              kegg_background =
                built_in$KEGG_background,
              pathway_name =
                input$archive_pathway %||%
                  "Phenylpropanoid biosynthesis",
              rank_metric =
                input$archive_rank_metric %||% "t",
              filter_mode =
                input$archive_filter_mode %||% "all",
              p_metric =
                input$archive_p_metric %||% "none",
              p_cutoff =
                input$archive_p_cutoff %||% 0.05,
              abs_logfc =
                input$archive_abs_logfc %||% 0,
              seed = seed
            ),
            error = function(e) e
          )
          if (inherits(item, "error")) {
            errors[[comparison]] <- conditionMessage(item)
          } else {
            results[[comparison]] <- item
          }
        }

        if (!length(results)) {
          archived_val(NULL)
          shiny::showNotification(
            paste(
              "GSEA failed for all selected comparisons:",
              paste(
                paste0(names(errors), ": ", errors),
                collapse = " | "
              )
            ),
            type = "error",
            duration = 10
          )
          return()
        }

        if (!focus %in% names(results)) {
          focus <- .protvis_gsea_rootve_comparison(
            names(results)
          )
        }

        archived_val(list(
          results = results,
          summary = .protvis_gsea_combined_summary(results),
          focus = focus,
          source = bundle$source,
          errors = errors
        ))

        focus_result <- results[[focus]]
        .protvis_record_shared_run(
          shared_state,
          module = "gsea",
          method = "previous_DEP_weighted_GSEA",
          category = "enrichment",
          parameters = list(
            comparisons = names(results),
            focus_comparison = focus,
            pathway = focus_result$archive$pathway,
            rank_metric =
              input$archive_rank_metric %||% "t",
            filter_mode =
              input$archive_filter_mode %||% "all",
            p_metric =
              input$archive_p_metric %||% "none",
            p_cutoff =
              input$archive_p_cutoff %||% 0.05,
            minimum_abs_log2FC =
              input$archive_abs_logfc %||% 0,
            weighted_p = 1,
            seed = seed,
            minSize = 5L,
            maxSize = 500L,
            eps = 0
          ),
          tables = list(
            gsea_summary =
              .protvis_gsea_combined_summary(results),
            focus_gsea_results =
              focus_result$fgsea_result,
            focus_metadata =
              focus_result$metadata,
            focus_ranked_background =
              focus_result$rank_tbl
          ),
          plot_data = list(
            focus_enrichment_curve =
              focus_result$rank_tbl
          ),
          plot_config = list(
            focus_comparison = focus,
            extreme_position =
              focus_result$extreme_position,
            left_label = "Zea mays ssp. mays",
            right_label = "Zea mays ssp. mexicana"
          )
        )

        if (length(errors)) {
          shiny::showNotification(
            paste0(
              "Completed ", length(results),
              " comparison(s); ",
              length(errors),
              " comparison(s) failed. See the summary/status."
            ),
            type = "warning",
            duration = 7
          )
        }
        return()
      }

      shiny::req(
        input$expr_file,      shiny::req(
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
      run <- archived_val()
      bundle <- previous_dep_val()

      if (is.null(run)) {
        loaded <- length(bundle$results %||% list())
        return(
          shiny::span(
            if (loaded) {
              paste0(
                "Loaded ", loaded,
                " comparison(s). Select comparisons and click Run Analysis. ",
                "Root_VE is the default focus when available."
              )
            } else {
              paste0(
                "Load previous DEP results first. ",
                "The module can load multiple comparisons; ",
                "the final figure can focus on Root_VE."
              )
            },
            class = "text-muted"
          )
        )
      }

      focus <- run$focus
      result <- run$results[[focus]]
      fg <- result$fgsea_result
      compatible <- isTRUE(
        result$metadata$Historical_DEP_compatible[[1L]]
      )

      shiny::tagList(
        shiny::span(
          class = "text-success fw-semibold",
          paste0(
            "✓ ", length(run$results),
            " comparison(s) analysed · focus: ",
            focus
          )
        ),
        shiny::tags$br(),
        shiny::span(
          paste0(
            "Rank metric: ",
            result$metadata$Rank_metric[[1L]],
            " · ranked proteins: ",
            format(
              nrow(result$rank_tbl),
              big.mark = ","
            ),
            " · pathway hits: ",
            sum(result$rank_tbl$in_pathway),
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
        ),
        shiny::tags$br(),
        shiny::span(
          result$compatibility$label,
          class = if (compatible) {
            "text-success"
          } else {
            "text-warning"
          }
        ),
        if (length(run$errors)) {
          shiny::tagList(
            shiny::tags$br(),
            shiny::span(
              paste0(
                "Failed comparisons: ",
                paste(names(run$errors), collapse = ", ")
              ),
              class = "text-warning"
            )
          )
        }
      )
    })

    output$archive_table <- DT::renderDT({
      run <- archived_val()
      shiny::req(run)
      DT::datatable(
        run$summary,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          scrollX = TRUE,
          pageLength = max(5L, nrow(run$summary)),
          dom = "Bfrtip",
          buttons = c("copy", "csv")
        )
      )
    })

    output$archive_metadata <- DT::renderDT({
      run <- archived_val()
      shiny::req(run)
      result <- run$results[[run$focus]]
      meta <- result$metadata
      display <- data.frame(
        Parameter = names(meta),
        Value = vapply(
          meta,
          function(x) as.character(x[[1L]]),
          character(1)
        ),
        stringsAsFactors = FALSE
      )
      DT::datatable(
        display,
        rownames = FALSE,
        options = list(
          dom = "t",
          pageLength = nrow(display),
          scrollX = TRUE
        )
      )
    })

    output$archive_curve_title <- shiny::renderUI({
      run <- archived_val()
      if (is.null(run)) {
        return("Phenylpropanoid biosynthesis · focus comparison")
      }
      result <- run$results[[run$focus]]
      shiny::tags$span(
        paste0(
          result$archive$pathway,
          " · ",
          run$focus
        )
      )
    })

    output$archive_curve <- shiny::renderPlot({
      run <- archived_val()
      shiny::req(run)
      result <- run$results[[run$focus]]
      print(.protvis_gsea_archived_plot(result))
    }, res = 110)

    output$download_csv <- shiny::downloadHandler(    output$download_csv <- shiny::downloadHandler(
      filename = function() {
        if (identical(
          input$analysis_mode,
          "archived"
        )) {
          return(
            paste0(
              "GSEA_previous_DEP_summary_",
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
          run <- archived_val()
          shiny::req(run)
          utils::write.csv(
            run$summary,
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
              "GSEA_",
              gsub(
                "[^A-Za-z0-9._-]+",
                "_",
                archived_val()$focus %||% "Root_VE"
              ),
              "_phenylpropanoid_",
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
          run <- archived_val()
          shiny::req(run)
          result <- run$results[[run$focus]]
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
