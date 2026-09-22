
#' @importFrom circlize colorRamp2 circos.clear circos.genomicInitialize circos.trackPlotRegion
#'   get.cell.meta.data circos.text get.all.sector.index circos.axis circos.genomicTrack
#'   circos.genomicRect circos.genomicText
#' @importFrom ComplexHeatmap draw Legend
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid unit gpar grid.lines
#' @importFrom utils head write.csv read.csv read.table
#' @importFrom grDevices pdf dev.off
#' @importFrom shiny NS moduleServer reactiveValues reactive req observeEvent renderUI
#'   uiOutput fileInput actionButton textAreaInput conditionalPanel selectInput
#'   checkboxGroupInput numericInput sliderInput downloadButton plotOutput renderPlot
#'   renderText showNotification textOutput tagList span div hr
#' @import bslib
#' @importFrom shinyWidgets switchInput
#' @importFrom colourpicker colourInput
#' @importFrom dplyr filter pull select mutate across everything
#' @importFrom tidyr separate
#' @importFrom stringr str_trim
#' @importFrom readxl read_excel excel_sheets
#' @importFrom tools file_ext
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom clusterProfiler enricher
#' @importFrom ggplot2 ggplot aes geom_col geom_point coord_flip theme_bw labs
#'   theme element_text scale_size_continuous
#' @importFrom graphics plot text par
#' @title plot_go_circos
#' @name plot_go_circos
#' @export
plot_go_circos <- function(go_data, top_n = 15, output_pdf = NULL) {
  if (base::is.null(go_data) || !base::is.data.frame(go_data) || base::nrow(go_data) == 0) {
    message("No enrichment data available to plot.")
    return(NULL)
  }

  required_cols <- c("ID", "Description", "pvalue", "BgRatio", "GeneRatio")
  if (!base::all(required_cols %in% base::colnames(go_data))) {
    message("Input data does not contain required columns.")
    return(NULL)
  }

  data <- go_data[base::order(go_data$pvalue), , drop = FALSE]
  datasig <- data[data$pvalue < 0.05, , drop = FALSE]
  data <- utils::head(datasig, top_n)

  if (base::nrow(data) == 0) {
    message("No significant GO/KEGG terms to plot.")
    return(NULL)
  }

  BgGene <- base::as.numeric(base::sapply(base::strsplit(base::as.character(data$BgRatio), "/"), `[`, 1))
  Gene <- base::as.numeric(base::sapply(base::strsplit(base::as.character(data$GeneRatio), "/"), `[`, 1))
  ratio <- Gene / BgGene
  logpvalue <- -base::log10(data$pvalue)

  logpvalue.col <- RColorBrewer::brewer.pal(n = 8, name = "Reds")

  circlize::colorRamp2(
    breaks = c(0, 2, 4, 6, 8, 10, 15, 20),
    colors = logpvalue.col
  )

  BgGene.col <- logpvalue.col[
    base::cut(
      base::pmin(logpvalue, 20),
      breaks = c(0, 2, 4, 6, 8, 10, 15, 20, Inf),
      include.lowest = TRUE,
      labels = FALSE
    )
  ]

  df_circos <- base::data.frame(
    GO = data$ID,
    start = 1,
    end = base::max(BgGene),
    stringsAsFactors = FALSE
  )
  base::rownames(df_circos) <- df_circos$GO

  bed2 <- base::data.frame(
    GO = data$ID,
    start = 1,
    end = BgGene,
    label = BgGene,
    col = BgGene.col,
    stringsAsFactors = FALSE
  )

  bed3 <- base::data.frame(
    GO = data$ID,
    start = 1,
    end = Gene,
    label = Gene,
    stringsAsFactors = FALSE
  )

  bed4 <- base::data.frame(
    GO = data$ID,
    start = 1,
    end = base::max(BgGene),
    ratio = ratio / base::max(ratio) * 9.5,
    col = "#00AFBB",
    stringsAsFactors = FALSE
  )

  if (!base::is.null(output_pdf)) {
    grDevices::pdf(output_pdf, width = 10, height = 6)
  }

  circlize::circos.clear()
  circlize::circos.genomicInitialize(df_circos, plotType = "none")

  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      desc <- data[data$ID == sector.index, "Description"]
      desc <- base::paste(base::strwrap(desc, width = 20), collapse = "\n")
      circlize::circos.text(
        base::mean(xlim), base::mean(ylim),
        desc,
        cex = 0.6,
        facing = "bending.inside",
        niceFacing = TRUE
      )
    },
    track.height = 0.12,
    bg.border = NA,
    bg.col = "grey95"
  )

  for (si in circlize::get.all.sector.index()) {
    circlize::circos.axis(
      h = "top",
      labels.cex = 0.5,
      sector.index = si,
      track.index = 1,
      major.at = base::seq(0, base::max(BgGene), by = 100),
      labels.facing = "clockwise"
    )
  }

  circlize::circos.genomicTrack(
    bed2,
    ylim = c(0, 1),
    track.height = 0.10,
    bg.border = "white",
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(
        region, value,
        ytop = 1, ybottom = 0,
        col = value$col, border = NA, ...
      )
      circlize::circos.genomicText(
        region, value,
        y = 0.4,
        labels = value$label,
        adj = 0,
        cex = 0.6,
        ...
      )
    }
  )

  circlize::circos.genomicTrack(
    bed3,
    ylim = c(0, 1),
    track.height = 0.10,
    bg.border = "white",
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(
        region, value,
        ytop = 1, ybottom = 0,
        col = "#BA55D3", border = NA, ...
      )
      circlize::circos.genomicText(
        region, value,
        y = 0.4,
        labels = value$label,
        adj = 0,
        cex = 0.6,
        ...
      )
    }
  )

  circlize::circos.genomicTrack(
    bed4,
    ylim = c(0, 10),
    track.height = 0.35,
    bg.border = "white",
    bg.col = "grey90",
    panel.fun = function(region, value, ...) {
      cell.xlim <- circlize::get.cell.meta.data("cell.xlim")
      cell.ylim <- circlize::get.cell.meta.data("cell.ylim")
      for (j in 1:9) {
        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
        grid::grid.lines(
          cell.xlim,
          c(y, y),
          gp = grid::gpar(col = "#FFFFFF", lwd = 0.3)
        )
      }
      circlize::circos.genomicRect(
        region, value,
        ytop = value$ratio,
        ybottom = 0,
        col = value$col,
        border = NA,
        ...
      )
    }
  )

  circlize::circos.clear()

  circle_size <- grid::unit(1, "snpc")

  ComplexHeatmap::draw(
    ComplexHeatmap::Legend(
      labels = c("Number of Genes", "Number of Select", "Rich Factor(0-1)"),
      type = "points",
      pch = c(15, 15, 17),
      legend_gp = grid::gpar(col = c("pink", "#BA55D3", "#00AFBB")),
      title = "",
      nrow = 3,
      size = grid::unit(3, "mm")
    ),
    x = circle_size * 0.83,
    y = circle_size * 0.5,
    just = "center"
  )

  ComplexHeatmap::draw(
    ComplexHeatmap::Legend(
      labels = c("(0,2]", "(2,4]", "(4,6]", "(6,8]", "(8,10]", "(10,15]", "(15,20]", ">=20"),
      type = "points",
      pch = 16,
      legend_gp = grid::gpar(col = logpvalue.col),
      title = "-log10(Pvalue)",
      title_position = "topcenter",
      grid_height = grid::unit(5, "mm"),
      grid_width = grid::unit(5, "mm"),
      size = grid::unit(3, "mm")
    ),
    x = circle_size * 1.4,
    y = circle_size * 0.5,
    just = "left"
  )

  if (!base::is.null(output_pdf)) {
    grDevices::dev.off()
  }

  message("GO/KEGG Circos plot finished!")
}


#' @title plot_enrichment_bar
#' @name plot_enrichment_bar
#' @keywords internal
plot_enrichment_bar <- function(enrich_df, top_n = 10, fill_color = "#2c7bb6", title = NULL) {
  if (base::is.null(enrich_df) || !base::is.data.frame(enrich_df) || base::nrow(enrich_df) == 0) {
    return(NULL)
  }

  df <- enrich_df
  if ("p.adjust" %in% base::colnames(df)) {
    df <- df[base::order(df$p.adjust, df$pvalue), , drop = FALSE]
  } else {
    df <- df[base::order(df$pvalue), , drop = FALSE]
  }

  df <- utils::head(df, top_n)
  df$Description <- base::factor(df$Description, levels = base::rev(df$Description))

  ggplot2::ggplot(df, ggplot2::aes(x = Description, y = Count)) +
    ggplot2::geom_col(fill = fill_color) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Count"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}


#' @title plot_enrichment_dot
#' @name plot_enrichment_dot
#' @keywords internal
plot_enrichment_dot <- function(enrich_df, top_n = 10, point_color = "#2c7bb6", title = NULL) {
  if (base::is.null(enrich_df) || !base::is.data.frame(enrich_df) || base::nrow(enrich_df) == 0) {
    return(NULL)
  }

  df <- enrich_df
  if ("p.adjust" %in% base::colnames(df)) {
    df <- df[base::order(df$p.adjust, df$pvalue), , drop = FALSE]
  } else {
    df <- df[base::order(df$pvalue), , drop = FALSE]
  }

  df <- utils::head(df, top_n)

  df$GeneRatio_num <- base::vapply(
    base::strsplit(base::as.character(df$GeneRatio), "/"),
    function(x) {
      base::as.numeric(x[1]) / base::as.numeric(x[2])
    },
    numeric(1)
  )

  df$Description <- base::factor(df$Description, levels = base::rev(df$Description))

  ggplot2::ggplot(df, ggplot2::aes(x = GeneRatio_num, y = Description, size = Count)) +
    ggplot2::geom_point(color = point_color) +
    ggplot2::theme_bw() +
    ggplot2::scale_size_continuous(range = c(3, 8)) +
    ggplot2::labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      size = "Count"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}


# Directional KEGG reproduction is calculated from the current archived DEP
# results and the annotation background. No precomputed Figure 3 enrichment
# result table or static plot is bundled in ProtVis.

.protvis_load_builtin_enrichment_background <- function() {
  read_table <- function(file_name) {
    utils::read.delim(
      .protvis_data_file(
        "backgrounds", "maize_teosinte", file_name
      ), sep = "\t", header = TRUE, quote = "",
      comment.char = "", stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  list(
    GO_background = read_table("maize_teosinte_GO_background.tsv.xz"),
    KEGG_background = read_table("maize_teosinte_KEGG_background.tsv.xz")
  )
}

.protvis_read_enrichment_background <- function(file) {
  sheets <- readxl::excel_sheets(file)
  standard <- c("GO_background", "KEGG_background")
  legacy <- c("t2g.go", "t2n.go", "t2g.kegg", "t2n.kegg")

  if (all(standard %in% sheets)) {
    out <- lapply(standard, function(sheet) {
      as.data.frame(readxl::read_excel(file, sheet = sheet),
                    stringsAsFactors = FALSE)
    })
    names(out) <- standard
  } else if (all(legacy %in% sheets)) {
    combine_legacy <- function(t2g_sheet, t2n_sheet) {
      t2g <- as.data.frame(readxl::read_excel(file, sheet = t2g_sheet),
                           stringsAsFactors = FALSE)
      t2n <- as.data.frame(readxl::read_excel(file, sheet = t2n_sheet),
                           stringsAsFactors = FALSE)
      if (!all(c("TERM", "GENE") %in% names(t2g)) ||
          !all(c("TERM", "NAME") %in% names(t2n))) {
        stop("Legacy annotation sheets must contain TERM/GENE and TERM/NAME.",
             call. = FALSE)
      }
      t2n <- t2n[!duplicated(t2n$TERM), c("TERM", "NAME"), drop = FALSE]
      merge(t2g[, c("TERM", "GENE"), drop = FALSE], t2n,
            by = "TERM", all.x = TRUE, sort = FALSE)
    }
    out <- list(
      GO_background = combine_legacy("t2g.go", "t2n.go"),
      KEGG_background = combine_legacy("t2g.kegg", "t2n.kegg")
    )
  } else {
    stop(
      "Background workbook needs GO_background/KEGG_background sheets or ",
      "legacy t2g.go/t2n.go/t2g.kegg/t2n.kegg sheets.",
      call. = FALSE
    )
  }

  for (name in names(out)) {
    table <- out[[name]]
    if (!all(c("TERM", "GENE", "NAME") %in% names(table))) {
      stop(name, " must contain TERM, GENE, and NAME columns.", call. = FALSE)
    }
    table <- table[, c("TERM", "GENE", "NAME"), drop = FALSE]
    table[] <- lapply(table, function(x) trimws(as.character(x)))
    table <- table[
      !is.na(table$TERM) & nzchar(table$TERM) &
        !is.na(table$GENE) & nzchar(table$GENE) &
        !is.na(table$NAME) & nzchar(table$NAME),
      , drop = FALSE
    ]
    out[[name]] <- unique(table)
  }
  out
}

.protvis_directional_expand_ids <- function(x) {
  ids <- unlist(
    strsplit(as.character(x), "[,;|]", perl = TRUE),
    use.names = FALSE
  )
  ids <- trimws(ids)
  ids <- sub("^CON__", "", ids)
  unique(ids[!is.na(ids) & nzchar(ids)])
}


.protvis_directional_background <- function(
    kegg_background, archived = FALSE) {
  background <- as.data.frame(
    kegg_background,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("TERM", "GENE", "NAME") %in% names(background))) {
    stop(
      "KEGG background must contain TERM, GENE, and NAME columns.",
      call. = FALSE
    )
  }

  background <- background[, c("TERM", "GENE", "NAME"), drop = FALSE]
  background[] <- lapply(
    background,
    function(x) trimws(as.character(x))
  )
  background <- background[
    !is.na(background$TERM) & nzchar(background$TERM) &
      !is.na(background$GENE) & nzchar(background$GENE) &
      !is.na(background$NAME) & nzchar(background$NAME),
    ,
    drop = FALSE
  ]

  if (isTRUE(archived)) {
    # Exact Enrichmentdb2 preparation used by the archived Figure 3 script:
    #   t2n: distinct TERM/NAME, then keep the first TERM for each NAME.
    #   t2g: original TERM/GENE rows restricted to retained t2n TERM values.
    t2n <- unique(background[, c("TERM", "NAME"), drop = FALSE])
    t2n <- t2n[!duplicated(t2n$NAME), , drop = FALSE]

    t2g <- unique(background[, c("TERM", "GENE"), drop = FALSE])
    t2g <- t2g[t2g$TERM %in% t2n$TERM, , drop = FALSE]

    return(list(
      TERM2GENE = unique(t2g),
      TERM2NAME = unique(t2n)
    ))
  }

  # Standard ProtVis ORA accepts protein-group strings and expands them to
  # individual identifiers before testing.
  gene_rows <- lapply(seq_len(nrow(background)), function(i) {
    genes <- .protvis_directional_expand_ids(background$GENE[[i]])
    if (!length(genes)) return(NULL)
    data.frame(
      TERM = rep(background$TERM[[i]], length(genes)),
      GENE = genes,
      stringsAsFactors = FALSE
    )
  })
  gene_rows <- gene_rows[!vapply(gene_rows, is.null, logical(1))]
  t2g <- if (length(gene_rows)) {
    unique(do.call(rbind, gene_rows))
  } else {
    data.frame(TERM = character(), GENE = character())
  }

  t2n <- unique(background[, c("TERM", "NAME"), drop = FALSE])
  t2n <- t2n[!duplicated(t2n$TERM), , drop = FALSE]

  list(
    TERM2GENE = unique(t2g),
    TERM2NAME = unique(t2n)
  )
}


.protvis_directional_figure3_spec <- function() {
  stages <- c(
    "Root_VE",
    "Root_V1.V2",
    "Root_V4",
    "Leaf_VE.V1.V2",
    "Leaf_V4.V6.V8"
  )
  data.frame(
    Stage = stages,
    Comparison = paste0(
      "B73_", stages, "_vs_Y12_", stages
    ),
    stringsAsFactors = FALSE
  )
}


.protvis_base36_to_integer <- function(x) {
  alphabet <- c(as.character(0:9), letters)
  decode_one <- function(value) {
    chars <- strsplit(tolower(as.character(value)), "", fixed = TRUE)[[1L]]
    digits <- match(chars, alphabet) - 1L
    if (!length(digits) || anyNA(digits)) {
      stop("Invalid value in bundled enrichment data.", call. = FALSE)
    }
    powers <- rev(seq_along(digits) - 1L)
    as.integer(sum(digits * (36 ^ powers)))
  }
  vapply(x, decode_one, integer(1))
}


.protvis_directional_load_figure3_gene_lists <- function(
    comparisons = NULL) {
  path <- .protvis_data_file("gsea", "gene_membership.b36")

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) {
    stop("Bundled DEP membership file is empty.", call. = FALSE)
  }

  decode_line <- function(line) {
    prefix <- substr(line, 1L, 1L)
    payload <- sub("^[PM]=", "", line)
    tokens <- strsplit(payload, ",", fixed = TRUE)[[1L]]
    pieces <- strsplit(tokens, ".", fixed = TRUE)
    if (any(lengths(pieces) != 2L)) {
    stop("Bundled DEP membership file is malformed.", call. = FALSE)
    }
    delta <- .protvis_base36_to_integer(
      vapply(pieces, `[[`, character(1), 1L)
    )
    mask <- .protvis_base36_to_integer(
      vapply(pieces, `[[`, character(1), 2L)
    )
    suffix <- cumsum(delta)
    id <- if (identical(prefix, "M")) {
      paste0("Zm00001d", sprintf("%06d", suffix))
    } else if (identical(prefix, "P")) {
      paste0("PZ00001a", sprintf("%06d", suffix))
    } else {
      stop("Unknown protein-ID prefix in bundled enrichment data.", call. = FALSE)
    }
    data.frame(
      ID = id,
      mask = as.integer(mask),
      stringsAsFactors = FALSE
    )
  }

  membership <- do.call(rbind, lapply(lines, decode_line))
  membership <- unique(membership)

  spec <- .protvis_directional_figure3_spec()
  all_up <- setNames(vector("list", nrow(spec)), spec$Stage)
  all_down <- setNames(vector("list", nrow(spec)), spec$Stage)

  for (i in seq_len(nrow(spec))) {
    up_bit <- bitwShiftL(1L, i - 1L)
    down_bit <- bitwShiftL(1L, i + 4L)
    all_up[[i]] <- unique(membership$ID[
      bitwAnd(membership$mask, up_bit) != 0L
    ])
    all_down[[i]] <- unique(membership$ID[
      bitwAnd(membership$mask, down_bit) != 0L
    ])
  }

  expected_up <- c(2287L, 2167L, 2442L, 2291L, 2283L)
  expected_down <- c(1315L, 1477L, 1350L, 1425L, 1419L)
  if (!identical(as.integer(lengths(all_up)), expected_up) ||
      !identical(as.integer(lengths(all_down)), expected_down)) {
    stop(
      "Bundled DEP membership failed its internal count check.",
      call. = FALSE
    )
  }

  if (is.null(comparisons)) {
    comparisons <- spec$Comparison
  }
  comparisons <- intersect(
    as.character(comparisons),
    spec$Comparison
  )
  if (!length(comparisons)) {
    stop("Select at least one comparison.", call. = FALSE)
  }

  selected <- spec$Comparison %in% comparisons
  spec_selected <- spec[selected, , drop = FALSE]
  stages <- spec_selected$Stage

  list(
    group1 = all_up[stages],
    group2 = all_down[stages],
    cluster_to_comparison = setNames(
      spec_selected$Comparison,
      spec_selected$Stage
    ),
    group1_names = paste0("B73_", stages),
    group2_names = paste0("Y12_", stages),
    source = paste0(
      "Frozen significant DEP lists from the original ",
      "02.MaizeTeosintePro/03.progress/03.dep outputs"
    ),
    counts = data.frame(
      Stage = spec$Stage,
      B73_higher = as.integer(lengths(all_up)),
      Y12_higher = as.integer(lengths(all_down)),
      stringsAsFactors = FALSE
    )
  )
}


.protvis_directional_archived_compatibility <- function(
    dep_results, comparisons = NULL) {
  if (!is.list(dep_results) || !length(dep_results)) {
    return(list(ok = FALSE, message = "Run DEP first."))
  }

  if (!is.null(comparisons)) {
    comparisons <- intersect(
      as.character(comparisons),
      names(dep_results)
    )
    dep_results <- dep_results[comparisons]
  }

  if (!length(dep_results)) {
    return(list(
      ok = FALSE,
      message = "Select at least one DEP comparison."
    ))
  }

  required <- c(
    "logFC", "P.Value", "Group1", "Group2",
    "analysis_mode", "protein_universe",
    "matrix_source", "test_method"
  )

  for (comparison in names(dep_results)) {
    result <- as.data.frame(
      dep_results[[comparison]],
      stringsAsFactors = FALSE
    )
    id_col <- c("ID", "protein_id", "Protein", "Gene")[
      c("ID", "protein_id", "Protein", "Gene") %in% names(result)
    ][1L]

    if (is.na(id_col)) {
      return(list(
        ok = FALSE,
        message = paste0(
          comparison, ": protein ID column is unavailable."
        )
      ))
    }

    missing <- setdiff(required, names(result))
    if (length(missing)) {
      return(list(
        ok = FALSE,
        message = paste0(
          comparison,
          ": built-in DEP workflow provenance is incomplete (",
          paste(missing, collapse = ", "),
          ")."
        )
      ))
    }

    mode <- unique(as.character(result$analysis_mode))
    mode <- mode[!is.na(mode) & nzchar(mode)]
    universe <- unique(as.character(result$protein_universe))
    universe <- universe[!is.na(universe) & nzchar(universe)]
    matrix_source <- unique(as.character(result$matrix_source))
    matrix_source <- matrix_source[
      !is.na(matrix_source) & nzchar(matrix_source)
    ]
    test_method <- unique(as.character(result$test_method))
    test_method <- test_method[
      !is.na(test_method) & nzchar(test_method)
    ]

    if (!length(mode) || any(mode != "archived") ||
        !length(universe) ||
        any(universe != "archived_any_detected") ||
        !length(matrix_source) ||
        any(matrix_source != "Step6_data_normalization") ||
        !length(test_method) ||
        any(test_method != "archived_eBayes")) {
      return(list(
        ok = FALSE,
        message = paste0(
          "The built-in KEGG workflow requires DEP > Built-in DEP workflow, ",
          "Detection filter = 'Detected in any comparison sample ",
          "(built-in workflow), and the predefined Step6 limma workflow."
        )
      ))
    }
  }

  list(
    ok = TRUE,
    message = paste0(
      "✓ Built-in DEP workflow compatible · Step6 limma · ",
      "|log2FC| > 1 · BH < 0.05"
    )
  )
}


.protvis_directional_fix_figure3_sample_map <- function(matrix) {
  x <- as.matrix(matrix)
  swap <- c("B73_Root_VE_3", "B73_Root_V1.V2_2")
  if (!all(swap %in% colnames(x))) {
    return(x)
  }

  tmp <- x[, swap[[1L]]]
  x[, swap[[1L]]] <- x[, swap[[2L]]]
  x[, swap[[2L]]] <- tmp
  attr(x, "figure3_sample_map_corrected") <- TRUE
  x
}


.protvis_directional_figure3_expected_retained <- function() {
  c(
    B73_Root_VE_vs_Y12_Root_VE = 11049L,
    B73_Root_V1.V2_vs_Y12_Root_V1.V2 = 11049L,
    B73_Root_V4_vs_Y12_Root_V4 = 11044L,
    B73_Leaf_VE.V1.V2_vs_Y12_Leaf_VE.V1.V2 = 11046L,
    B73_Leaf_V4.V6.V8_vs_Y12_Leaf_V4.V6.V8 = 11036L
  )
}


.protvis_directional_rebuild_archived_dep <- function(
    dep_results, comparisons,
    normalized_matrix, detection_matrix,
    sample_info) {
  if (!is.list(dep_results) || !length(dep_results)) {
    stop(
      "Current DEP results are required to define the selected comparisons.",
      call. = FALSE
    )
  }

  comparisons <- intersect(
    as.character(comparisons),
    names(dep_results)
  )
  if (!length(comparisons)) {
    stop("Select at least one DEP comparison.", call. = FALSE)
  }

  normalized_matrix <- .protvis_directional_fix_figure3_sample_map(
    normalized_matrix
  )
  detection_matrix <- .protvis_directional_fix_figure3_sample_map(
    detection_matrix
  )
  storage.mode(normalized_matrix) <- "numeric"
  storage.mode(detection_matrix) <- "numeric"

  retained_audit <- integer()

  sample_info <- as.data.frame(
    sample_info,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("sample_id", "group") %in% names(sample_info))) {
    stop(
      "The built-in DEP workflow requires sample_info columns: sample_id and group.",
      call. = FALSE
    )
  }

  if (is.null(rownames(normalized_matrix)) ||
      is.null(rownames(detection_matrix))) {
    stop(
      "The built-in DEP workflow requires protein IDs as matrix row names.",
      call. = FALSE
    )
  }

  rebuilt <- list()

  first_label <- function(x, fallback) {
    x <- unique(trimws(as.character(x)))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[[1L]] else fallback
  }

  for (comparison in comparisons) {
    source <- as.data.frame(
      dep_results[[comparison]],
      stringsAsFactors = FALSE
    )

    group1 <- first_label(
      source$Group1 %||% character(),
      NA_character_
    )
    group2 <- first_label(
      source$Group2 %||% character(),
      NA_character_
    )

    if (is.na(group1) || is.na(group2)) {
      parts <- strsplit(
        comparison,
        "_vs_",
        fixed = TRUE
      )[[1L]]
      if (length(parts) == 2L) {
        group1 <- parts[[1L]]
        group2 <- parts[[2L]]
      }
    }

    if (is.na(group1) || is.na(group2) ||
        !nzchar(group1) || !nzchar(group2)) {
      stop(
        "Could not recover Group1/Group2 for comparison: ",
        comparison,
        call. = FALSE
      )
    }

    group1_samples <- as.character(
      sample_info$sample_id[
        as.character(sample_info$group) == group1
      ]
    )
    group2_samples <- as.character(
      sample_info$sample_id[
        as.character(sample_info$group) == group2
      ]
    )

    group1_samples <- intersect(
      group1_samples,
      colnames(normalized_matrix)
    )
    group2_samples <- intersect(
      group2_samples,
      colnames(normalized_matrix)
    )

    if (!length(group1_samples) || !length(group2_samples)) {
      stop(
        "The built-in DEP workflow could not find samples for ",
        group1, " vs ", group2, ".",
        call. = FALSE
      )
    }

    required_samples <- c(
      group1_samples,
      group2_samples
    )
    if (!all(required_samples %in% colnames(detection_matrix))) {
      stop(
        "Step4 detection matrix is missing samples required for ",
        comparison, ".",
        call. = FALSE
      )
    }

    keep_ids <- .protvis_dep_shared_ids(
      detection_matrix,
      group1_samples,
      group2_samples,
      mode = "archived_any_detected"
    )
    keep_ids <- intersect(
      keep_ids,
      rownames(normalized_matrix)
    )
    if (!length(keep_ids)) {
      stop(
        "No proteins passed the reference-workflow detection filter for ",
        comparison, ".",
        call. = FALSE
      )
    }

    retained_audit[[comparison]] <- length(keep_ids)

    result <- .protvis_dep_run_limma_archived(
      normalized_matrix[
        keep_ids,
        ,
        drop = FALSE
      ],
      group1_samples,
      group2_samples,
      group1,
      group2,
      adjust_method = "BH",
      sort_by = "logFC",
      matrix_shift = TRUE
    )

    result <- .protvis_dep_classify(
      result,
      logfc = 1,
      cutoff = 0.05,
      p_metric = "adj.P.Val",
      fc_threshold_tested = FALSE
    )
    result$analysis_mode <- "archived"
    result$protein_universe <- "archived_any_detected"
    result$matrix_source <- "Step6_data_normalization"
    result$test_method <- "archived_eBayes"
    result$reproduction_source <-
      "Directional KEGG auto-rebuild from Step4 + Step6"

    rebuilt[[comparison]] <- result
  }

  attr(rebuilt, "retained_audit") <- retained_audit
  attr(rebuilt, "figure3_sample_map_corrected") <- TRUE
  rebuilt
}


.protvis_directional_archived_lists <- function(
    dep_results, comparisons) {
  compatibility <- .protvis_directional_archived_compatibility(
    dep_results,
    comparisons
  )
  if (!isTRUE(compatibility$ok)) {
    stop(compatibility$message, call. = FALSE)
  }

  group1_lists <- list()
  group2_lists <- list()
  cluster_to_comparison <- character()
  group1_names <- character()
  group2_names <- character()

  for (comparison in comparisons) {
    result <- as.data.frame(
      dep_results[[comparison]],
      stringsAsFactors = FALSE
    )
    id_col <- c("ID", "protein_id", "Protein", "Gene")[
      c("ID", "protein_id", "Protein", "Gene") %in% names(result)
    ][1L]

    group1 <- unique(trimws(as.character(result$Group1)))
    group1 <- group1[!is.na(group1) & nzchar(group1)]
    group2 <- unique(trimws(as.character(result$Group2)))
    group2 <- group2[!is.na(group2) & nzchar(group2)]
    group1 <- if (length(group1)) group1[[1L]] else "Group1"
    group2 <- if (length(group2)) group2[[1L]] else "Group2"

    cluster <- .protvis_dep_stage_label(group1)
    if (cluster %in% names(group1_lists)) {
      cluster <- comparison
    }

    # Historical Figure 3 used Protein_ID exactly as stored in the archived
    # matrix. Do not expand semicolon-delimited protein-group strings here.
    ids <- trimws(as.character(result[[id_col]]))
    logfc <- suppressWarnings(as.numeric(result$logFC))
    raw_p <- suppressWarnings(as.numeric(result$P.Value))

    # Recreate adj.P.Val from the complete comparison. This fixes the exact
    # historical threshold independently of the current DEP display settings.
    bh <- stats::p.adjust(raw_p, method = "BH")
    finite <- !is.na(ids) & nzchar(ids) &
      is.finite(logfc) & is.finite(bh)

    group1_lists[[cluster]] <- unique(ids[
      finite & bh < 0.05 & logfc > 1
    ])
    group2_lists[[cluster]] <- unique(ids[
      finite & bh < 0.05 & logfc < -1
    ])

    cluster_to_comparison[[cluster]] <- comparison
    group1_names <- c(group1_names, group1)
    group2_names <- c(group2_names, group2)
  }

  group1_lists <- group1_lists[lengths(group1_lists) > 0L]
  group2_lists <- group2_lists[lengths(group2_lists) > 0L]

  list(
    group1 = group1_lists,
    group2 = group2_lists,
    cluster_to_comparison = cluster_to_comparison,
    group1_names = group1_names,
    group2_names = group2_names
  )
}


.protvis_directional_group_root <- function(x, fallback) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(fallback)

  roots <- sub("_.*$", "", x)
  roots <- unique(roots[nzchar(roots)])
  if (length(roots) == 1L) roots[[1L]] else fallback
}


.protvis_directional_archived_title <- function(
    data, direction) {
  meta <- attr(data, "directional_meta") %||% list()
  root <- if (identical(direction, "Upregulated")) {
    meta$group1_root %||% "Group 1"
  } else {
    meta$group2_root %||% "Group 2"
  }

  if (identical(root, "B73")) {
    return(
      expression(italic("Zea mays ssp. mays") ~ "– enriched")
    )
  }
  if (identical(root, "Y12")) {
    return(
      expression(italic("Zea mays ssp. mexicana") ~ "– enriched")
    )
  }

  paste0(root, " – enriched")
}


.protvis_directional_kegg_data <- function(
    dep_results = NULL, kegg_background,
    presence_absence = list(),
    evidence_mode = c("quantitative", "all_retained"),
    comparisons = NULL,
    top_n = 10L,
    p_adjust_cutoff = 0.05,
    method = c(
      "figure3_archive",
      "archived_comparecluster",
      "standard_ora"
    ),
    pvalue_cutoff = 0.05) {
  evidence_mode <- match.arg(evidence_mode)
  method <- match.arg(method)
  if (!is.list(presence_absence)) {
    presence_absence <- list()
  }
  top_n <- max(1L, as.integer(top_n))

  archived_style <- method %in% c(
    "figure3_archive",
    "archived_comparecluster"
  )

  if (identical(method, "figure3_archive")) {
    spec <- .protvis_directional_figure3_spec()
    if (is.null(comparisons)) {
      comparisons <- spec$Comparison
    }
    comparisons <- intersect(
      as.character(comparisons),
      spec$Comparison
    )
    if (!length(comparisons)) {
      return(data.frame())
    }
    lists <- .protvis_directional_load_figure3_gene_lists(
      comparisons
    )
  } else {
    if (!is.list(dep_results) || !length(dep_results)) {
      return(data.frame())
    }
    if (is.null(comparisons)) {
      comparisons <- names(dep_results)
    }
    comparisons <- intersect(
      as.character(comparisons),
      names(dep_results)
    )
    if (!length(comparisons)) {
      return(data.frame())
    }

    if (method %in% c("figure3_archive", "archived_comparecluster")) {
      if (!identical(evidence_mode, "quantitative")) {
        stop(
          "The built-in KEGG workflow uses quantitative DEP only.",
          call. = FALSE
        )
      }
      lists <- .protvis_directional_archived_lists(
        dep_results,
        comparisons
      )
    }
  }

  if (archived_style) {
    background <- .protvis_directional_background(
      kegg_background,
      archived = TRUE
    )

    run_comparecluster <- function(gene_clusters) {
      if (!length(gene_clusters)) {
        return(NULL)
      }

      tryCatch(
        clusterProfiler::compareCluster(
          geneCluster = gene_clusters,
          fun = clusterProfiler::enricher,
          TERM2GENE = background$TERM2GENE,
          TERM2NAME = background$TERM2NAME,
          pvalueCutoff = as.numeric(pvalue_cutoff),
          qvalueCutoff = 1
        ),
        error = function(e) {
          stop(
            "Built-in KEGG workflow failed: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
    }

    objects <- list(
      Upregulated = run_comparecluster(lists$group1),
      Downregulated = run_comparecluster(lists$group2)
    )

    rows <- lapply(names(objects), function(direction) {
      object <- objects[[direction]]
      if (is.null(object)) return(NULL)

      out <- tryCatch(
        as.data.frame(object),
        error = function(e) NULL
      )
      if (is.null(out) || !nrow(out)) return(NULL)

      out$Comparison <- unname(
        lists$cluster_to_comparison[
          as.character(out$Cluster)
        ]
      )
      out$Direction <- direction
      out$Direction_label <- if (
        identical(direction, "Upregulated")
      ) {
        paste0(
          .protvis_directional_group_root(
            lists$group1_names,
            "Group 1"
          ),
          " higher"
        )
      } else {
        paste0(
          .protvis_directional_group_root(
            lists$group2_names,
            "Group 2"
          ),
          " higher"
        )
      }
      out$Evidence <- if (identical(
        method,
        "figure3_archive"
      )) {
        "Built-in DEP lists"
      } else {
        "Built-in DEP workflow"
      }
      out$Analysis_mode <- if (identical(
        method,
        "figure3_archive"
      )) {
        "Built-in KEGG workflow"
      } else {
        "Built-in DEP workflow"
      }
      out$Universe <-
        "KEGG annotation background (clusterProfiler default)"
      out
    })

    rows <- rows[!vapply(rows, is.null, logical(1))]
    out <- if (length(rows)) {
      do.call(rbind, rows)
    } else {
      data.frame()
    }

    attr(out, "analysis_method") <- method
    attr(out, "comparecluster_objects") <- objects
    attr(out, "plot_top_n") <- top_n
    attr(out, "directional_meta") <- list(
      group1_root = .protvis_directional_group_root(
        lists$group1_names,
        "Group 1"
      ),
      group2_root = .protvis_directional_group_root(
        lists$group2_names,
        "Group 2"
      ),
      pvalue_cutoff = as.numeric(pvalue_cutoff),
      qvalue_cutoff = 1,
      universe = "annotation background",
      source = if (identical(method, "figure3_archive")) {
        lists$source
      } else {
        "Built-in DEP workflow"
      }
    )

    return(out)
  }

  dep_results <- dep_results[comparisons]
  background <- .protvis_directional_background(
    kegg_background,
    archived = FALSE
  )
  t2g <- background$TERM2GENE
  t2n <- background$TERM2NAME

  p_adjust_cutoff <- as.numeric(p_adjust_cutoff)
  rows <- lapply(names(dep_results), function(comparison) {
    result <- as.data.frame(
      dep_results[[comparison]],
      stringsAsFactors = FALSE
    )
    id_col <- c("ID", "protein_id", "Protein", "Gene")[
      c("ID", "protein_id", "Protein", "Gene") %in% names(result)
    ][1L]
    if (is.na(id_col) ||
        !"regulation" %in% names(result)) {
      return(NULL)
    }

    tested <- .protvis_directional_expand_ids(
      result[[id_col]]
    )
    presence <- presence_absence[[comparison]] %||% data.frame()
    presence_id_col <- c("ID", "protein_id", "Protein", "Gene")[
      c("ID", "protein_id", "Protein", "Gene") %in% names(presence)
    ][1L]
    presence_ids <- if (!is.na(presence_id_col)) {
      .protvis_directional_expand_ids(
        presence[[presence_id_col]]
      )
    } else {
      character()
    }

    if (identical(evidence_mode, "all_retained")) {
      tested <- unique(c(tested, presence_ids))
    }
    if (!length(tested)) {
      return(NULL)
    }

    first_label <- function(x, fallback) {
      x <- unique(trimws(as.character(x)))
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x)) x[[1L]] else fallback
    }

    labels <- c(
      Upregulated = first_label(
        result$Group1 %||% character(),
        "Group 1"
      ),
      Downregulated = first_label(
        result$Group2 %||% character(),
        "Group 2"
      )
    )

    lapply(names(labels), function(direction) {
      genes <- .protvis_directional_expand_ids(
        result[[id_col]][
          as.character(result$regulation) == direction
        ]
      )

      if (identical(evidence_mode, "all_retained") &&
          nrow(presence) &&
          "Pattern" %in% names(presence) &&
          !is.na(presence_id_col)) {
        pattern <- if (
          identical(direction, "Upregulated")
        ) {
          "Detected in Group1 only"
        } else {
          "Detected in Group2 only"
        }

        genes <- unique(c(
          genes,
          .protvis_directional_expand_ids(
            presence[[presence_id_col]][
              presence$Pattern == pattern
            ]
          )
        ))
      }
      if (!length(genes)) {
        return(NULL)
      }

      enriched <- tryCatch(
        clusterProfiler::enricher(
          gene = genes,
          universe = tested,
          TERM2GENE = t2g,
          TERM2NAME = t2n,
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          minGSSize = 1,
          maxGSSize = Inf
        ),
        error = function(e) NULL
      )
      out <- tryCatch(
        as.data.frame(enriched@result),
        error = function(e) NULL
      )
      if (is.null(out) || !nrow(out)) {
        return(NULL)
      }

      out <- out[
        is.finite(out$p.adjust) &
          out$p.adjust <= p_adjust_cutoff,
        ,
        drop = FALSE
      ]
      if (!nrow(out)) {
        return(NULL)
      }

      out <- out[
        order(
          out$p.adjust,
          out$pvalue,
          -out$Count
        ),
        ,
        drop = FALSE
      ]
      out <- utils::head(out, top_n)
      out$Comparison <- comparison
      out$Cluster <- .protvis_dep_stage_label(
        labels[[direction]]
      )
      out$Direction <- direction
      out$Direction_label <- paste0(
        labels[[direction]],
        " higher"
      )
      out$Evidence <- if (
        identical(evidence_mode, "all_retained")
      ) {
        "All retained proteins"
      } else {
        "Evidence 1 · Quantitative DEP"
      }
      out$Analysis_mode <- "Standard ORA"
      out$Universe <- "Comparison-specific tested proteins"
      out
    })
  })

  rows <- unlist(rows, recursive = FALSE)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  out <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame()
  }
  attr(out, "analysis_method") <- method
  attr(out, "plot_top_n") <- top_n
  out
}


.protvis_directional_plot_colour <- function(
    value, fallback) {
  value <- as.character(value %||% "")
  if (length(value) != 1L || !nzchar(value)) {
    return(fallback)
  }
  ok <- tryCatch({
    grDevices::col2rgb(value)
    TRUE
  }, error = function(e) FALSE)
  if (isTRUE(ok)) value else fallback
}


.protvis_directional_kegg_panel <- function(
    data, direction,
    show_legend = TRUE,
    top_n = NULL,
    colour_low = "#2C7FB8",
    colour_high = "#F7FBFF") {
  colour_low <- .protvis_directional_plot_colour(
    colour_low, "#2C7FB8"
  )
  colour_high <- .protvis_directional_plot_colour(
    colour_high, "#F7FBFF"
  )
  method <- attr(data, "analysis_method") %||% "standard_ora"

  if (method %in% c("figure3_archive", "archived_comparecluster")) {
    objects <- attr(data, "comparecluster_objects") %||% list()
    object <- objects[[direction]]

    if (is.null(object)) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate(
            "text",
            x = 0,
            y = 0,
            label = paste(
              "No significant",
              tolower(direction),
              "KEGG pathways."
            ),
            colour = "#64748B",
            size = 4
          )
      )
    }

    if (is.null(top_n)) {
      top_n <- attr(data, "plot_top_n") %||% 10L
    }

    plot <- enrichplot::dotplot(
      object,
      showCategory = max(1L, as.integer(top_n)),
      size = "count",
      color = "pvalue",
      label_format = 100
    ) +
      ggplot2::scale_colour_gradient(
        low = colour_low,
        high = colour_high,
        name = "pvalue"
      ) +
      ggplot2::labs(
        title = .protvis_directional_archived_title(
          data,
          direction
        ),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          hjust = 1,
          colour = "black"
        ),
        axis.text.y = ggplot2::element_text(
          colour = "black"
        ),
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          colour = "black",
          size = 12
        ),
        panel.border = ggplot2::element_rect(
          colour = "black",
          linewidth = 0.8
        ),
        legend.position = if (isTRUE(show_legend)) {
          "right"
        } else {
          "none"
        }
      )

    return(plot)
  }

  df <- data[
    as.character(data$Direction) == direction,
    ,
    drop = FALSE
  ]
  if (!nrow(df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = paste(
            "No significant",
            tolower(direction),
            "KEGG pathways."
          ),
          colour = "#64748B",
          size = 4
        ) +
        ggplot2::labs(
          title = if (identical(direction, "Upregulated")) {
            "↑ Upregulated proteins"
          } else {
            "↓ Downregulated proteins"
          }
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0,
            face = "bold",
            colour = "#24384D",
            size = 12
          )
        )
    )
  }

  x_levels <- unique(df$Cluster)
  y_levels <- unique(
    df$Description[
      order(df$p.adjust, df$pvalue)
    ]
  )
  df$Cluster <- factor(df$Cluster, levels = x_levels)
  df$Description <- factor(
    df$Description,
    levels = rev(y_levels)
  )

  title <- if (identical(direction, "Upregulated")) {
    "↑ Upregulated proteins"
  } else {
    "↓ Downregulated proteins"
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Cluster,
      y = Description,
      size = Count,
      fill = p.adjust
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "#355268",
      stroke = 0.4,
      alpha = 0.95
    ) +
    ggplot2::scale_fill_gradient(
      low = colour_low,
      high = colour_high,
      name = "BH-adjusted P"
    ) +
    ggplot2::scale_size_continuous(
      range = c(3.5, 10.5),
      name = "Protein count"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        order = 1,
        reverse = TRUE
      ),
      size = ggplot2::guide_legend(
        order = 2,
        override.aes = list(fill = "white")
      )
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1,
        colour = "#334155"
      ),
      axis.text.y = ggplot2::element_text(
        colour = "#334155"
      ),
      plot.title = ggplot2::element_text(
        hjust = 0,
        face = "bold",
        colour = "#24384D",
        size = 12
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "#D6E4EE",
        fill = NA,
        linewidth = 0.4
      ),
      legend.position = if (isTRUE(show_legend)) {
        "bottom"
      } else {
        "none"
      },
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(10, 12, 8, 10)
    )
}


.protvis_plot_directional_kegg <- function(
    data, ncol = 2L, top_n = NULL,
    colour_low = "#2C7FB8",
    colour_high = "#F7FBFF") {
  method <- attr(data, "analysis_method") %||% "standard_ora"

  if (method %in% c("figure3_archive", "archived_comparecluster")) {
    return(
      patchwork::wrap_plots(
        .protvis_directional_kegg_panel(
          data,
          "Upregulated",
          top_n = top_n,
          colour_low = colour_low,
          colour_high = colour_high
        ),
        .protvis_directional_kegg_panel(
          data,
          "Downregulated",
          top_n = top_n,
          colour_low = colour_low,
          colour_high = colour_high
        ),
        ncol = ncol,
        widths = c(1.05, 0.95)
      ) +
        patchwork::plot_annotation(
          tag_levels = list(c("C", "D")),
          theme = ggplot2::theme(
            plot.tag = ggplot2::element_text(
              size = 14,
              face = "bold",
              colour = "black"
            )
          )
        )
    )
  }

  evidence <- unique(as.character(data$Evidence))
  evidence <- evidence[
    !is.na(evidence) & nzchar(evidence)
  ]

  plot <- patchwork::wrap_plots(
    .protvis_directional_kegg_panel(
      data,
      "Upregulated",
      colour_low = colour_low,
      colour_high = colour_high
    ),
    .protvis_directional_kegg_panel(
      data,
      "Downregulated",
      colour_low = colour_low,
      colour_high = colour_high
    ),
    ncol = ncol
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Directional KEGG enrichment across DEP comparisons",
      subtitle = paste0(
        if (length(evidence)) {
          evidence[[1L]]
        } else {
          "Selected differential evidence"
        },
        ": each comparison uses its retained tested proteins ",
        "as the enrichment universe."
      )
    )

  plot & ggplot2::theme(legend.position = "bottom")
}


#' Enrichment Analysis Module UI
#'
#' This function creates the user interface for the enrichment analysis module.
#'
#' @param id The namespace identifier for the module
#' @return A Shiny UI tagList containing the enrichment analysis interface
#' @import shiny
#' @import bslib
#' @importFrom shinyWidgets switchInput
#' @importFrom colourpicker colourInput
#' @name enrichment_analysis_ui
#' @export
enrichment_analysis_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    bslib::navset_card_tab(
      id = ns("enrichment_main_tabs"),

      bslib::nav_panel(
        "Enrichment workflow",
        bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 300,

        shiny::actionButton(
          ns("load_data"),
          "LOAD DATA",
          class = "btn btn-light fw-bold"
        ),

        shiny::uiOutput(ns("load_status_panel")),
        shiny::uiOutput(ns("compare_select_ui")),

        shiny::div(
          style = "margin-bottom: 15px;",
          shiny::fileInput(
            ns("enrichment_analysis_file"),
            "Upload Enrichment Analysis File (.xlsx)",
            accept = c(".xlsx"),
            buttonLabel = "Browse..."
          ),
          shiny::tags$small(
            "Maximum upload size: 2 GB (configured globally for ProtVis).",
            class = "text-muted"
          ),
          shiny::actionButton(
            ns("check_file"),
            "Check File",
            class = "btn btn-success fw-bold mb-2"
          )
        ),

        shiny::hr(),

        shiny::tags$small(
          "The genelist requires an ID column (.xlsx or .csv).",
          style = "color: #6c757d"
        ),

        shinyWidgets::switchInput(
          inputId = ns("input_mode"),
          label = "Input Manually",
          value = TRUE,
          onLabel = "Upload",
          offLabel = "Paste",
          width = "100%"
        ),

        shiny::conditionalPanel(
          condition = base::paste0("input['", ns("input_mode"), "'] == true"),
          shiny::tags$small("Upload Genelist", style = "color: #6c757d"),
          shiny::fileInput(
            inputId = ns("genelist_file"),
            label = NULL,
            multiple = FALSE,
            accept = c(".csv", ".xlsx")
          )
        ),

        shiny::conditionalPanel(
          condition = base::paste0("input['", ns("input_mode"), "'] == false"),
          shiny::div(
            shiny::tags$small("Paste Genelist", style = "color: #6c757d"),
            shiny::textAreaInput(
              inputId = ns("paste_data"),
              label = NULL,
              placeholder = "Copy and paste Excel data here.",
              rows = 5
            ),
            shiny::actionButton(
              ns("apply_paste"),
              "Apply paste data",
              class = "btn btn-light fw-bold"
            )
          )
        ),

        bslib::accordion(
          bslib::accordion_panel(
            title = "Enrichment analysis",
            icon = enrichment_bubble_icon,

            shiny::selectInput(
              inputId = ns("species"),
              label = "Select taxonomic group:",
              choices = c(
                "Plant" = "Plant",
                "Animals" = "Animals",
                "Bacteria" = "Bacteria",
                "Fungi" = "Fungi",
                "Eukaryotes" = "Eukaryotes",
                "Hsa" = "Hsa"
              ),
              selected = "Plant"
            ),

            shiny::checkboxGroupInput(
              inputId = ns("choices"),
              label = "Please select the analysis content:",
              choices = c(
                "GO" = "go_analysis",
                "KEGG" = "kegg_analysis"
              ),
              selected = c("go_analysis", "kegg_analysis")
            ),

            shiny::actionButton(
              ns("run_enrichment_analysis"),
              "Analysis"
            ),
            shiny::textOutput(ns("analysis_message"))
          )
        )
      ),

      bslib::card(
        bslib::card_header("File Check Result"),
        bslib::card_body(
          shiny::textOutput(ns("file_check_result"))
        )
      ),

      bslib::layout_column_wrap(
        width = 1 / 2,
        bslib::card(
          bslib::card_header("Enrichment Analysis File"),
          bslib::card_body(
            shiny::tags$small(
              "Preview of the validated GO/KEGG background tables (first 100 rows per sheet).",
              class = "text-muted"
            ),
            shiny::uiOutput(ns("background_preview_ui"))
          )
        ),
        bslib::card(
          bslib::card_header("Genelist (ID column required)"),
          bslib::card_body(
            shiny::tags$small(
              "Uploaded gene/protein IDs are shown here; otherwise the DEP gene list is shown.",
              class = "text-muted"
            ),
            DT::DTOutput(ns("genelist_preview"))
          )
        )
      ),

      bslib::layout_column_wrap(
        width = 1 / 2,
        height = 600,

        bslib::card(
          height = "800px",
          bslib::card_header("GO Enrichment Analysis"),
          bslib::card_body(
            shiny::tabsetPanel(
              id = ns("go_tabs"),
              type = "tabs",

              shiny::tabPanel(
                "Visualization",
                bslib::layout_sidebar(
                  sidebar = bslib::sidebar(
                    width = 250,
                    position = "left",
                    open = "open",

                    shiny::selectInput(
                      ns("go_plot_type"),
                      "Select plot type:",
                      choices = c(
                        "Bar plot" = "bar",
                        "Dot plot" = "dot",
                        "Circle plot" = "circle"
                      ),
                      selected = "bar"
                    ),

                    shiny::sliderInput(
                      ns("go_top_n"),
                      "Top N terms:",
                      min = 5,
                      max = 20,
                      value = 10
                    ),

                    colourpicker::colourInput(
                      ns("go_color"),
                      "Select color:",
                      value = "#2c7bb6"
                    ),

                    shiny::numericInput(
                      ns("go_width"),
                      "Plot width (inch)",
                      value = 8,
                      min = 4,
                      max = 20
                    ),

                    shiny::numericInput(
                      ns("go_height"),
                      "Plot height (inch)",
                      value = 6,
                      min = 4,
                      max = 20
                    ),

                    shiny::downloadButton(
                      ns("download_go_plot"),
                      "Download Plot (PDF)"
                    ),

                    shiny::downloadButton(
                      ns("download_go_table"),
                      "Download Table (CSV)"
                    )
                  ),

                  bslib::card_body(
                    shiny::plotOutput(ns("go_plot"))
                  )
                )
              ),

              shiny::tabPanel(
                "Result Table",
                DT::DTOutput(ns("go_res_table"))
              )
            )
          )
        ),

        bslib::card(
          height = "800px",
          bslib::card_header("KEGG Enrichment Analysis"),
          bslib::card_body(
            shiny::tabsetPanel(
              id = ns("kegg_tabs"),
              type = "tabs",

              shiny::tabPanel(
                "Visualization",
                bslib::layout_sidebar(
                  sidebar = bslib::sidebar(
                    width = 250,
                    position = "left",
                    open = "open",

                    shiny::selectInput(
                      ns("kegg_plot_type"),
                      "Select plot type:",
                      choices = c(
                        "Bar plot" = "bar",
                        "Dot plot" = "dot",
                        "Circle plot" = "circle"
                      ),
                      selected = "bar"
                    ),

                    shiny::sliderInput(
                      ns("kegg_top_n"),
                      "Top N pathways:",
                      min = 5,
                      max = 20,
                      value = 10
                    ),

                    colourpicker::colourInput(
                      ns("kegg_color"),
                      "Select color:",
                      value = "#d7191c"
                    ),

                    shiny::numericInput(
                      ns("kegg_width"),
                      "Plot width (inch)",
                      value = 8,
                      min = 4,
                      max = 20
                    ),

                    shiny::numericInput(
                      ns("kegg_height"),
                      "Plot height (inch)",
                      value = 6,
                      min = 4,
                      max = 20
                    ),

                    shiny::downloadButton(
                      ns("download_kegg_plot"),
                      "Download Plot (PDF)"
                    ),

                    shiny::downloadButton(
                      ns("download_kegg_table"),
                      "Download Table (CSV)"
                    )
                  ),

                  bslib::card_body(
                    shiny::plotOutput(ns("kegg_plot"))
                  )
                )
              ),

              shiny::tabPanel(
                "Result Table",
                DT::DTOutput(ns("kegg_res_table"))
              )
            )
          )
        )
      )
        )
      ),

      bslib::nav_panel(
        "Directional KEGG",
        shiny::div(
          class = "protvis-directional-kegg",
          shiny::tags$style(shiny::HTML("
            .protvis-directional-kegg {
              --pv-blue: #2095CF;
              --pv-blue-dark: #147FB8;
              --pv-blue-soft: #F3F9FD;
              --pv-border: #D6E4EE;
              --pv-text: #24384D;
              --pv-muted: #64748B;
            }
            .protvis-directional-kegg .directional-kegg-shell {
              border: 1px solid var(--pv-border);
              border-radius: 14px;
              box-shadow: none;
              overflow: hidden;
            }
            .protvis-directional-kegg .directional-kegg-title {
              font-size: 1.05rem;
              font-weight: 700;
              color: var(--pv-text);
              margin: 0;
            }
            .protvis-directional-kegg .directional-kegg-intro {
              color: var(--pv-muted);
              margin: 4px 0 0;
              font-size: 0.92rem;
              font-weight: 400;
            }
            .protvis-directional-kegg .directional-kegg-config-grid {
              display: grid;
              grid-template-columns: minmax(230px, 0.9fr) minmax(420px, 1.9fr) minmax(235px, 0.85fr);
              gap: 16px;
              align-items: stretch;
              margin-bottom: 14px;
            }
            .protvis-directional-kegg .directional-config-card {
              min-width: 0;
              background: #FFFFFF;
              border: 1px solid var(--pv-border);
              border-radius: 12px;
              padding: 16px 18px;
            }
            .protvis-directional-kegg .directional-card-heading {
              display: flex;
              align-items: flex-start;
              justify-content: space-between;
              gap: 12px;
              margin-bottom: 10px;
            }
            .protvis-directional-kegg .directional-card-heading h5 {
              color: var(--pv-text);
              font-size: 0.96rem;
              font-weight: 700;
              margin: 0;
            }
            .protvis-directional-kegg .directional-card-heading p {
              color: var(--pv-muted);
              font-size: 0.78rem;
              line-height: 1.35;
              margin: 3px 0 0;
            }
            .protvis-directional-kegg .directional-selection-count {
              flex: 0 0 auto;
              color: var(--pv-blue-dark);
              background: var(--pv-blue-soft);
              border-radius: 999px;
              padding: 3px 9px;
              font-size: 0.76rem;
              font-weight: 700;
              white-space: nowrap;
            }
            .protvis-directional-kegg .directional-evidence-card .form-group,
            .protvis-directional-kegg .directional-settings-card .form-group {
              margin-bottom: 10px;
            }
            .protvis-directional-kegg .directional-evidence-card .shiny-options-group {
              margin-top: 2px;
            }
            .protvis-directional-kegg .directional-card-divider {
              border-top: 1px solid #E7EEF4;
              margin: 12px 0;
            }
            .protvis-directional-kegg .directional-background-row {
              display: flex;
              align-items: flex-end;
              justify-content: space-between;
              gap: 10px;
            }
            .protvis-directional-kegg .directional-background-label {
              color: var(--pv-text);
              font-size: 0.82rem;
              font-weight: 700;
              margin-bottom: 4px;
            }
            .protvis-directional-kegg .directional-state {
              font-size: 0.78rem;
              line-height: 1.3;
            }
            .protvis-directional-kegg .directional-state-loaded {
              color: #2E7D55;
            }
            .protvis-directional-kegg .directional-state-waiting {
              color: var(--pv-muted);
            }
            .protvis-directional-kegg .directional-fixed-settings {
              color: var(--pv-muted);
              background: #F8FBFD;
              border: 1px solid #E7EEF4;
              border-radius: 8px;
              padding: 8px 10px;
              font-size: 0.78rem;
              line-height: 1.4;
              margin-top: 6px;
            }
            .protvis-directional-kegg .directional-fixed-settings .directional-state {
              display: block;
              margin-top: 5px;
            }
            .protvis-directional-kegg .directional-comparison-scroll {
              max-height: 156px;
              overflow-y: auto;
              padding-right: 8px;
              scrollbar-width: thin;
            }
            .protvis-directional-kegg .directional-comparison-scroll .checkbox {
              display: block;
              margin: 0 0 7px;
              line-height: 1.25;
            }
            .protvis-directional-kegg .directional-comparison-scroll .form-group {
              margin-bottom: 0;
            }
            .protvis-directional-kegg .directional-comparison-actions {
              display: flex;
              gap: 8px;
              margin-top: 10px;
              padding-top: 10px;
              border-top: 1px solid #E7EEF4;
            }
            .protvis-directional-kegg .directional-settings-grid {
              display: grid;
              grid-template-columns: 1fr;
              gap: 2px;
            }
            .protvis-directional-kegg .directional-colour-block {
              margin-top: 8px;
              padding-top: 10px;
              border-top: 1px solid #E7EEF4;
            }
            .protvis-directional-kegg .directional-colour-title {
              font-size: 0.78rem;
              font-weight: 700;
              color: #24384D;
              margin-bottom: 6px;
            }
            .protvis-directional-kegg .directional-colour-controls {
              display: grid;
              grid-template-columns: repeat(2, minmax(0, 1fr));
              gap: 8px;
              align-items: end;
            }
            .protvis-directional-kegg .directional-colour-controls .form-group {
              margin-bottom: 4px;
            }
            .protvis-directional-kegg .directional-colour-controls label {
              font-size: 0.74rem;
              color: var(--pv-muted);
              margin-bottom: 3px;
            }
            .protvis-directional-kegg .directional-colour-reset {
              width: 100%;
              margin-top: 4px;
            }
            .protvis-directional-kegg .directional-run-button {
              width: 100%;
              min-height: 42px;
              margin-top: 6px;
              font-weight: 700;
              border-radius: 8px;
            }
            .protvis-directional-kegg .directional-result-toolbar {
              display: flex;
              align-items: center;
              justify-content: space-between;
              gap: 14px;
              flex-wrap: wrap;
              background: var(--pv-blue-soft);
              border: 1px solid var(--pv-border);
              border-radius: 10px;
              padding: 9px 12px;
              margin-bottom: 12px;
            }
            .protvis-directional-kegg .directional-result-status {
              display: flex;
              align-items: center;
              gap: 8px;
              min-width: 240px;
              color: var(--pv-muted);
              font-size: 0.82rem;
            }
            .protvis-directional-kegg .directional-result-status i {
              color: var(--pv-blue-dark);
            }
            .protvis-directional-kegg .directional-export-cluster {
              display: flex;
              align-items: center;
              gap: 7px;
              margin-left: auto;
            }
            .protvis-directional-kegg .directional-export-label {
              color: var(--pv-muted);
              font-size: 0.78rem;
              font-weight: 700;
              margin-right: 2px;
            }
            .protvis-directional-kegg .directional-export-cluster .shiny-download-link {
              min-width: 64px;
            }
            .protvis-directional-kegg .directional-results-tabs > .tab-content {
              padding-top: 12px;
            }
            .protvis-directional-kegg .directional-plot-grid {
              display: grid;
              grid-template-columns: repeat(2, minmax(0, 1fr));
              gap: 18px;
              align-items: stretch;
            }
            .protvis-directional-kegg .directional-plot-card {
              min-width: 0;
              background: #FFFFFF;
              border: 1px solid var(--pv-border);
              border-radius: 12px;
              padding: 8px 10px 2px;
            }
            .protvis-directional-kegg .directional-plot-card .shiny-plot-output {
              width: 100% !important;
            }
            @media (max-width: 1399.98px) {
              .protvis-directional-kegg .directional-plot-grid {
                grid-template-columns: 1fr;
              }
            }
            @media (max-width: 1199.98px) {
              .protvis-directional-kegg .directional-kegg-config-grid {
                grid-template-columns: minmax(240px, 0.9fr) minmax(360px, 1.4fr);
              }
              .protvis-directional-kegg .directional-settings-card {
                grid-column: 1 / -1;
              }
              .protvis-directional-kegg .directional-settings-grid {
                grid-template-columns: minmax(140px, 1fr) minmax(140px, 1fr) minmax(200px, 1.2fr);
                gap: 12px;
                align-items: end;
              }
              .protvis-directional-kegg .directional-run-button {
                margin-bottom: 10px;
              }
            }
            @media (max-width: 767.98px) {
              .protvis-directional-kegg .directional-kegg-config-grid,
              .protvis-directional-kegg .directional-settings-grid {
                grid-template-columns: 1fr;
              }
              .protvis-directional-kegg .directional-settings-card {
                grid-column: auto;
              }
              .protvis-directional-kegg .directional-result-toolbar {
                align-items: flex-start;
              }
              .protvis-directional-kegg .directional-export-cluster {
                width: 100%;
                margin-left: 0;
                flex-wrap: wrap;
              }
            }
          ")),
          bslib::card(
            class = "directional-kegg-shell",
            bslib::card_header(
              shiny::div(
                shiny::tags$h4(
                  "Directional KEGG enrichment across DEP comparisons",
                  class = "directional-kegg-title"
                ),
                shiny::tags$p(
                  "Run the built-in KEGG workflow or a custom ORA workflow across selected DEP contrasts.",
                  class = "directional-kegg-intro"
                )
              )
            ),
            bslib::card_body(
              shiny::div(
                class = "directional-kegg-config-grid",

                shiny::div(
                  class = "directional-config-card directional-evidence-card",
                  shiny::div(
                    class = "directional-card-heading",
                    shiny::div(
                      shiny::tags$h5("Analysis mode"),
                      shiny::tags$p(
                        "Choose a built-in workflow with predefined settings or a custom ProtVis ORA workflow."
                      )
                    )
                  ),
                  shiny::selectInput(
                    ns("directional_method"),
                    label = NULL,
                    choices = c(
                      "Built-in KEGG workflow" =
                        "figure3_archive",
                      "Standard ORA (tested universe + BH)" =
                        "standard_ora"
                    ),
                    selected = "figure3_archive"
                  ),
                  shiny::conditionalPanel(
                    condition = paste0(
                      "input['",
                      ns("directional_method"),
                      "'] == 'figure3_archive'"
                    ),
                    shiny::div(
                      class = "directional-fixed-settings",
                      shiny::tags$small(
                        "Uses bundled DEP protein lists and predefined settings; KEGG is recalculated when you run the analysis."
                      ),
                      shiny::uiOutput(
                        ns("directional_dep_compatibility")
                      )
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = paste0(
                      "input['",
                      ns("directional_method"),
                      "'] == 'standard_ora'"
                    ),
                    shiny::radioButtons(
                      ns("directional_evidence"),
                      label = "Evidence source",
                      choices = c(
                        "Quantitative DEP" = "quantitative",
                        "All retained proteins" = "all_retained"
                      ),
                      selected = "quantitative",
                      inline = FALSE
                    )
                  ),
                  shiny::div(
                    class = "directional-card-divider"
                  ),
                  shiny::div(
                    class = "directional-background-row",
                    shiny::div(
                      shiny::div(
                        "Background",
                        class = "directional-background-label"
                      ),
                      shiny::uiOutput(
                        ns("directional_background_status")
                      )
                    ),
                    shiny::actionButton(
                      ns("load_maize_teosinte_background"),
                      "Load built-in",
                      icon = shiny::icon("database"),
                      class = "btn-outline-primary btn-sm"
                    )
                  )
                ),

                shiny::div(
                  class = "directional-config-card directional-comparisons-card",
                  shiny::div(
                    class = "directional-card-heading",
                    shiny::div(
                      shiny::tags$h5("DEP comparisons"),
                      shiny::tags$p("Select one or more contrasts; long lists stay inside this panel.")
                    ),
                    shiny::span(
                      shiny::textOutput(ns("directional_selected_count"), inline = TRUE),
                      class = "directional-selection-count"
                    )
                  ),
                  shiny::div(
                    class = "directional-comparison-scroll",
                    shiny::uiOutput(ns("directional_comparisons_ui"))
                  ),
                  shiny::div(
                    class = "directional-comparison-actions",
                    shiny::actionButton(
                      ns("directional_select_all"), "Select all",
                      class = "btn-outline-secondary btn-sm"
                    ),
                    shiny::actionButton(
                      ns("directional_clear_all"), "Clear",
                      class = "btn-outline-secondary btn-sm"
                    )
                  )
                ),

                shiny::div(
                  class = "directional-config-card directional-settings-card",
                  shiny::div(
                    class = "directional-card-heading",
                    shiny::div(
                      shiny::tags$h5("KEGG settings"),
                      shiny::tags$p(
                        "Reporting depth is shared; statistical filtering follows the selected mode."
                      )
                    )
                  ),
                  shiny::div(
                    class = "directional-settings-grid",
                    shiny::numericInput(
                      ns("directional_top_n"),
                      "Top pathways",
                      10,
                      min = 1,
                      max = 20
                    ),
                    shiny::conditionalPanel(
                      condition = paste0(
                        "input['",
                        ns("directional_method"),
                        "'] == 'standard_ora'"
                      ),
                      shiny::numericInput(
                        ns("directional_p_adjust"),
                        "BH FDR cutoff",
                        0.05,
                        min = 0,
                        max = 1,
                        step = 0.01
                      )
                    ),
                    shiny::conditionalPanel(
                      condition = paste0(
                        "input['",
                        ns("directional_method"),
                        "'] == 'figure3_archive'"
                      ),
                      shiny::div(
                        class = "directional-fixed-settings",
                        shiny::tags$small(
                          "Built-in settings: bundled DEP lists; p-value cutoff = 0.05; q-value cutoff = 1; raw p-value statistic; bundled Enrichmentdb2."
                        )
                      )
                    ),
                    shiny::div(
                      class = "directional-colour-block",
                      shiny::div(
                        "Colour gradient",
                        class = "directional-colour-title"
                      ),
                      shiny::div(
                        class = "directional-colour-controls",
                        colourpicker::colourInput(
                          ns("directional_colour_low"),
                          "Low value",
                          value = "#2C7FB8",
                          showColour = "background",
                          allowTransparent = FALSE
                        ),
                        colourpicker::colourInput(
                          ns("directional_colour_high"),
                          "High value",
                          value = "#F7FBFF",
                          showColour = "background",
                          allowTransparent = FALSE
                        )
                      ),
                      shiny::actionButton(
                        ns("reset_directional_colours"),
                        "RESET BLUE",
                        icon = shiny::icon("rotate-left"),
                        class = paste(
                          "btn-outline-secondary btn-sm",
                          "directional-colour-reset"
                        )
                      )
                    ),
                    shiny::actionButton(
                      ns("run_directional_kegg"),
                      "RUN KEGG",
                      icon = shiny::icon("play"),
                      class = "btn-primary directional-run-button"
                    )
                  )
                )
              ),

              shiny::div(
                class = "directional-result-toolbar",
                shiny::div(
                  class = "directional-result-status",
                  shiny::icon("circle-info"),
                  shiny::textOutput(ns("directional_kegg_status"), inline = TRUE)
                ),
                shiny::div(
                  class = "directional-export-cluster",
                  shiny::span("Export", class = "directional-export-label"),
                  shiny::downloadButton(
                    ns("download_directional_kegg_pdf"), "PDF",
                    class = "btn-outline-secondary btn-sm"
                  ),
                  shiny::downloadButton(
                    ns("download_directional_kegg_svg"), "SVG",
                    class = "btn-outline-secondary btn-sm"
                  ),
                  shiny::downloadButton(
                    ns("download_directional_kegg_data"), "CSV",
                    class = "btn-outline-secondary btn-sm"
                  )
                )
              ),

              shiny::div(
                class = "directional-results-tabs",
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "Figure",
                    shiny::div(
                      class = "directional-plot-grid",
                      shiny::div(
                        class = "directional-plot-card",
                        shiny::plotOutput(
                          ns("directional_kegg_up_plot"),
                          height = "470px"
                        )
                      ),
                      shiny::div(
                        class = "directional-plot-card",
                        shiny::plotOutput(
                          ns("directional_kegg_down_plot"),
                          height = "470px"
                        )
                      )
                    )
                  ),
                  shiny::tabPanel(
                    "Result table",
                    DT::DTOutput(ns("directional_kegg_table"))
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

utils::globalVariables(c(
  "regulation", "V3", "Pathway_ID", "TERM", "GENE", "NAME",
  "Panel", "Cluster_label", "Cluster", "Description", "Count", "pvalue",
  "Direction"
))

# Store the complete, validated enrichment workbook in the project object.
# The two sheets remain separate so a saved ProtVis_dataset can be reopened
# and used for GO/KEGG analysis without re-uploading the workbook.
.protvis_add_enrichment_background <- function(dataset, background, file_name) {
  dataset <- as_protvis_dataset(dataset)
  required <- c("GO_background", "KEGG_background")
  if (!is.list(background) || !all(required %in% names(background)) ||
      !all(vapply(background[required], is.data.frame, logical(1)))) {
    stop("background must contain GO_background and KEGG_background data.frames.",
         call. = FALSE)
  }

  dataset <- .protvis_new_analysis_dataset(
    dataset,
    "enrichment_background",
    parameters = list(
      file_name = as.character(file_name),
      worksheets = required
    )
  )
  dataset$other_files$enrichment_background <- list(
    file_name = as.character(file_name),
    uploaded_at = as.character(Sys.time()),
    sheets = background[required]
  )
  dataset <- .protvis_append_process(
    dataset,
    "enrichment_background",
    status = "success",
    parameters = dataset$process_info$parameters$enrichment_background,
    message = paste0("Stored enrichment background workbook: ", file_name)
  )
  validate_protvis_dataset(dataset)
  dataset
}


#' Enrichment Analysis Module Server
#'
#' This function provides the server-side logic for the enrichment analysis module.
#'
#' @param id The namespace identifier for the module
#' @param shared_state A reactive values list for sharing state between modules
#' @return A module server function
#' @import shiny
#' @name enrichment_analysis_server
#' @export
enrichment_analysis_server <- function(id, shared_state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    rv <- shiny::reactiveValues(
      sample_info = NULL,
      load_success = FALSE,
      normalized_matrix = NULL,
      compare_data = NULL,
      dep_results = base::list(),
      file_check_msg = NULL,
      background_data = NULL,
      go_res = NULL,
      kegg_res = NULL,
      directional_kegg = NULL,
      directional_kegg_message = NULL,
      analysis_message = NULL,
      pasted_genelist = NULL
    )

    set_background <- function(background, file_name, message) {
      rv$background_data <- background
      shared_state$pending_enrichment_background <- background
      shared_state$pending_enrichment_background_name <- file_name
      if (inherits(shared_state$dataset, "ProtVis_dataset")) {
        saved <- tryCatch({
          dataset <- .protvis_add_enrichment_background(
            shared_state$dataset, background, file_name
          )
          if (!base::is.null(shared_state$workdir) &&
              base::dir.exists(shared_state$workdir)) {
            dataset <- protvis_auto_export_dataset(
              dataset, directory = shared_state$workdir
            )
          }
          .protvis_ui_sync_state(dataset, shared_state)
          TRUE
        }, error = function(e) {
          rv$file_check_msg <- paste("❌ Background valid but could not be saved:", e$message)
          FALSE
        })
        if (!saved) return(FALSE)
        rv$file_check_msg <- paste0(message, " Saved to ProtVis_dataset.")
      } else {
        rv$file_check_msg <- paste0(
          message,
          " It will be added to ProtVis_dataset when Project init is completed."
        )
      }
      TRUE
    }

    shiny::observeEvent(input$load_maize_teosinte_background, {
      background <- tryCatch(
        .protvis_load_builtin_enrichment_background(),
        error = function(e) e
      )
      if (inherits(background, "error")) {
        shiny::showNotification(conditionMessage(background), type = "error")
        return()
      }
      set_background(
        background,
        "maize_teosinte_background.tsv.xz",
        "✅ Built-in maize-teosinte GO/KEGG background loaded."
      )
      shiny::showNotification("Maize-teosinte enrichment background loaded.", type = "message")
    }, ignoreInit = TRUE)

    get_result_df <- function(enrich_obj) {
      if (base::is.null(enrich_obj)) {
        return(NULL)
      }

      res <- tryCatch(
        enrich_obj@result,
        error = function(e) NULL
      )

      if (!base::is.null(res) && base::nrow(res) > 0) {
        return(base::as.data.frame(res))
      }

      NULL
    }

    normalize_dep_results <- function(dep_obj) {
      if (base::is.null(dep_obj) || !base::is.list(dep_obj) ||
          base::length(dep_obj) == 0L) {
        return(NULL)
      }
      dep_obj <- dep_obj[!vapply(dep_obj, is.null, logical(1))]
      if (base::length(dep_obj) == 0L) return(NULL)
      dep_obj
    }

    directional_dep_bundle <- function() {
      results <- normalize_dep_results(shared_state$dep_results)
      presence <- shared_state$presence_absence %||% list()
      if (base::is.null(results) && base::length(rv$dep_results) > 0L) {
        results <- rv$dep_results
      }
      if (base::is.null(results) && inherits(shared_state$dataset, "ProtVis_dataset")) {
        stored_dep <- shared_state$dataset@analysis_results$DEP
        if (base::is.list(stored_dep) && base::is.list(stored_dep$results)) {
          results <- normalize_dep_results(stored_dep$results)
          presence <- stored_dep$presence_absence %||% list()
        }
      }
      list(results = results, presence = presence)
    }

    directional_available_comparisons <- function() {
      if (identical(
        input$directional_method %||% "figure3_archive",
        "figure3_archive"
      )) {
        return(
          .protvis_directional_figure3_spec()$Comparison
        )
      }
      bundle <- directional_dep_bundle()
      names(bundle$results %||% list())
    }

    directional_archived_stage_inputs <- function(load = FALSE) {
      workdir <- as.character(shared_state$workdir %||% "")
      if (!nzchar(workdir) || !dir.exists(workdir)) {
        return(list(
          ok = FALSE,
          message = paste0(
            "Working directory is unavailable; the built-in DEP ",
            "Step4/Step6 inputs cannot be reconstructed."
          )
        ))
      }

      step4_path <- file.path(
        workdir,
        "Step4_data_transformed.rda"
      )
      step6_path <- file.path(
        workdir,
        "Step6_data_normalization.rda"
      )

      missing <- c(
        if (!file.exists(step4_path)) {
          "Step4_data_transformed.rda"
        },
        if (!file.exists(step6_path)) {
          "Step6_data_normalization.rda"
        }
      )

      if (length(missing)) {
        return(list(
          ok = FALSE,
          message = paste0(
            "The built-in KEGG workflow needs ",
            paste(missing, collapse = " + "),
            " in the project working directory."
          )
        ))
      }

      if (!isTRUE(load)) {
        return(list(
          ok = TRUE,
          message = paste0(
            "✓ Step4 + Step6 available · built-in DEP workflow will ",
            "be rebuilt automatically"
          ),
          step4_path = step4_path,
          step6_path = step6_path
        ))
      }

      step4 <- .protvis_load_stage_dataset(
        step4_path,
        expression_names = "transformed"
      )
      step6 <- .protvis_load_stage_dataset(
        step6_path,
        expression_names = "normalized_data"
      )

      list(
        ok = TRUE,
        message = paste0(
          "✓ Step4 + Step6 loaded · built-in DEP workflow prepared ",
          "inside Directional KEGG"
        ),
        step4 = step4,
        step6 = step6,
        step4_path = step4_path,
        step6_path = step6_path
      )
    }

    normalise_term2gene <- function(background) {
      background <- base::as.data.frame(background, stringsAsFactors = FALSE)
      if (!base::all(c("TERM", "GENE") %in% base::names(background))) {
        return(base::data.frame(TERM = character(), GENE = character()))
      }
      pieces <- base::lapply(base::seq_len(base::nrow(background)), function(i) {
        genes <- unlist(base::strsplit(base::as.character(background$GENE[[i]]),
                                       "[,;|]", perl = TRUE), use.names = FALSE)
        genes <- trimws(genes)
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (!length(genes)) return(NULL)
        base::data.frame(TERM = as.character(background$TERM[[i]]),
                         GENE = genes, stringsAsFactors = FALSE)
      })
      pieces <- pieces[!vapply(pieces, is.null, logical(1))]
      if (!length(pieces)) return(base::data.frame(TERM = character(), GENE = character()))
      unique(do.call(rbind, pieces))
    }

    shiny::observeEvent(input$load_data, {
      dep_obj <- normalize_dep_results(shared_state$dep_results)
      source_label <- "current DEP results"

      # ProtVis_dataset is the canonical source for saved analyses.
      if (base::is.null(dep_obj) && inherits(shared_state$dataset, "ProtVis_dataset")) {
        stored_dep <- shared_state$dataset@analysis_results$DEP
        if (base::is.list(stored_dep) && base::is.list(stored_dep$results)) {
          dep_obj <- normalize_dep_results(stored_dep$results)
          if (!base::is.null(dep_obj)) source_label <- "ProtVis_dataset DEP"
        }
        stored <- shared_state$dataset@analysis_results$differential_analysis
        if (base::is.null(dep_obj) && base::is.list(stored) && base::is.list(stored$comparisons)) {
          dep_obj <- normalize_dep_results(stored$comparisons)
          if (!base::is.null(dep_obj)) source_label <- "ProtVis_dataset"
        }
        if (base::is.null(dep_obj) && base::is.list(stored) && base::is.data.frame(stored$table)) {
          tab <- stored$table
          id_col <- if ("protein_id" %in% names(tab)) "protein_id" else if ("ID" %in% names(tab)) "ID" else NULL
          if (!base::is.null(id_col)) {
            logfc_col <- if ("log2FC" %in% names(tab)) "log2FC" else if ("logFC" %in% names(tab)) "logFC" else NULL
            p_col <- if ("p_value" %in% names(tab)) "p_value" else if ("P.Value" %in% names(tab)) "P.Value" else NULL
            if (!base::is.null(logfc_col) && !base::is.null(p_col)) {
              tab$ID <- base::as.character(tab[[id_col]])
              tab$logFC <- base::as.numeric(tab[[logfc_col]])
              tab$P.Value <- base::as.numeric(tab[[p_col]])
              significant <- if ("significant" %in% names(tab)) {
                !is.na(tab$significant) & base::as.logical(tab$significant)
              } else {
                rep(TRUE, base::nrow(tab))
              }
              tab$regulation <- ifelse(
                significant & tab$logFC > 0, "Upregulated",
                ifelse(significant & tab$logFC < 0, "Downregulated", "Not significant")
              )
              comparison <- base::paste(stored$group1 %||% "Group1", "vs", stored$group2 %||% "Group2")
              dep_obj <- base::list()
              dep_obj[[comparison]] <- tab
              source_label <- "ProtVis_dataset"
            }
          }
        }
      }

      # Load the result written in the current working directory.  Keep the
      # old Step7 filename as a read-only migration fallback.
      if (base::is.null(dep_obj) && !base::is.null(shared_state$workdir)) {
        for (rda_name in c("differential_analysis.rda", "Step7_DEP_result.rda")) {
          if (!base::is.null(dep_obj)) break
          rda_path <- base::file.path(shared_state$workdir, rda_name)
          if (!base::file.exists(rda_path)) next
          e <- base::new.env()
          loaded <- tryCatch({ base::load(rda_path, envir = e); TRUE }, error = function(e) FALSE)
          if (isTRUE(loaded) && base::exists("dep_results2", envir = e, inherits = FALSE)) {
            dep_obj <- normalize_dep_results(base::get("dep_results2", envir = e, inherits = FALSE))
            source_label <- rda_name
          }
        }
      }

      if (!base::is.null(dep_obj)) {
        rv$dep_results <- dep_obj
        rv$load_success <- TRUE

        dep_names <- base::names(rv$dep_results)

        if (!base::is.null(dep_names) && base::length(dep_names) > 0) {
          rv$compare_data <- rv$dep_results[[dep_names[1]]]
        } else {
          rv$compare_data <- NULL
        }

        shiny::showNotification(base::paste("Data loaded successfully from", source_label, "."), type = "message")
      } else {
        rv$dep_results <- base::list()
        rv$compare_data <- NULL
        rv$load_success <- FALSE

        shiny::showNotification(
          paste0(
            "No DEP results are available for the standard enrichment workflow. ",
            "Run DEP first, or upload a gene list for an independent analysis."
          ),
          type = "message"
        )
      }
    })

    output$load_status_panel <- shiny::renderUI({
      if (base::isTRUE(rv$load_success)) {
        shiny::span("✅ Data loaded", style = "color: green;")
      } else {
        shiny::span("❌ Data not loaded", style = "color: red;")
      }
    })

    output$directional_comparisons_ui <- shiny::renderUI({
      choices <- directional_available_comparisons()
      if (!length(choices)) {
        return(
          shiny::helpText(
            "Run DEP first; available comparisons will appear here."
          )
        )
      }
      selected <- isolate(input$directional_comparisons)
      selected <- intersect(selected %||% choices, choices)
      if (!length(selected)) selected <- choices
      shiny::checkboxGroupInput(
        ns("directional_comparisons"),
        label = NULL,
        choices = choices,
        selected = selected,
        inline = FALSE
      )
    })

    output$directional_selected_count <- shiny::renderText({
      choices <- directional_available_comparisons()
      selected <- intersect(
        input$directional_comparisons %||% character(),
        choices
      )
      paste0(length(selected), " selected")
    })

    output$directional_background_status <- shiny::renderUI({
      if (identical(
        input$directional_method %||% "figure3_archive",
        "figure3_archive"
      )) {
        return(
          shiny::span(
            "✓ Fixed · bundled Enrichmentdb2",
            class = "directional-state directional-state-loaded"
          )
        )
      }

      background <- rv$background_data
      if (base::is.null(background) ||
          base::is.null(background$KEGG_background)) {
        return(
          shiny::span(
            "○ Not loaded",
            class = "directional-state directional-state-waiting"
          )
        )
      }

      table <- base::as.data.frame(background$KEGG_background)
      n_ids <- if ("GENE" %in% base::names(table)) {
        base::length(base::unique(
          base::as.character(table$GENE[
            !base::is.na(table$GENE) & base::nzchar(table$GENE)
          ])
        ))
      } else {
        base::nrow(table)
      }

      shiny::span(
        paste0(
          "✓ Loaded · ",
          base::format(n_ids, big.mark = ","),
          " annotated IDs"
        ),
        class = "directional-state directional-state-loaded"
      )
    })

    output$directional_dep_compatibility <- shiny::renderUI({
      if (identical(
        input$directional_method %||% "figure3_archive",
        "figure3_archive"
      )) {
        lists <- tryCatch(
          .protvis_directional_load_figure3_gene_lists(),
          error = function(e) e
        )
        if (inherits(lists, "error")) {
          return(
            shiny::span(
              conditionMessage(lists),
              class = "directional-state directional-state-waiting"
            )
          )
        }
        return(
          shiny::span(
            paste0(
              "✓ Built-in DEP lists verified · B73-higher ",
              sum(lists$counts$B73_higher),
              " memberships · Y12-higher ",
              sum(lists$counts$Y12_higher),
              " memberships"
            ),
            class = "directional-state directional-state-loaded"
          )
        )
      }

      shiny::span(
        "Standard ORA uses the currently loaded DEP results.",
        class = "directional-state directional-state-loaded"
      )
    })

    shiny::observeEvent(input$directional_method, {
      if (identical(
        input$directional_method,
        "figure3_archive"
      )) {
        shiny::updateRadioButtons(
          session,
          "directional_evidence",
          selected = "quantitative"
        )
        shiny::updateNumericInput(
          session,
          "directional_top_n",
          value = 10
        )
      }

      choices <- directional_available_comparisons()
      shiny::updateCheckboxGroupInput(
        session,
        "directional_comparisons",
        choices = choices,
        selected = choices,
        inline = FALSE
      )

      rv$directional_kegg <- NULL
      rv$directional_kegg_message <- if (identical(
        input$directional_method,
        "figure3_archive"
      )) {
        paste0(
          "Built-in KEGG workflow selected. ",
          "Bundled DEP lists and Enrichmentdb2 ",
          "will be used; KEGG is recalculated when RUN KEGG is clicked."
        )
      } else {
        paste0(
          "Standard ORA selected. ",
          "Choose evidence and comparisons, then run KEGG."
        )
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$reset_directional_colours, {
      colourpicker::updateColourInput(
        session,
        "directional_colour_low",
        value = "#2C7FB8"
      )
      colourpicker::updateColourInput(
        session,
        "directional_colour_high",
        value = "#F7FBFF"
      )
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$directional_select_all, {
      choices <- directional_available_comparisons()
      shiny::updateCheckboxGroupInput(
        session,
        "directional_comparisons",
        choices = choices,
        selected = choices,
        inline = FALSE
      )
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$directional_clear_all, {
      choices <- directional_available_comparisons()
      shiny::updateCheckboxGroupInput(
        session,
        "directional_comparisons",
        choices = choices,
        selected = character(),
        inline = FALSE
      )
    }, ignoreInit = TRUE)

    output$background_preview_ui <- shiny::renderUI({
      if (base::is.null(rv$background_data)) {
        return(shiny::tags$p("Upload and check an Enrichment Analysis File first.", class = "text-muted"))
      }
      shiny::tabsetPanel(
        id = ns("background_tabs"),
        shiny::tabPanel("GO_background", DT::DTOutput(ns("go_background_preview"))),
        shiny::tabPanel("KEGG_background", DT::DTOutput(ns("kegg_background_preview")))
      )
    })

    output$go_background_preview <- DT::renderDT({
      table <- utils::head(base::as.data.frame(rv$background_data$GO_background), 100L)
      DT::datatable(table, options = base::list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })

    output$kegg_background_preview <- DT::renderDT({
      table <- utils::head(base::as.data.frame(rv$background_data$KEGG_background), 100L)
      DT::datatable(table, options = base::list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })

    output$background_preview <- DT::renderDT({
      background <- rv$background_data
      if (base::is.null(background)) {
        return(DT::datatable(
          base::data.frame(Message = "Upload and check an Enrichment Analysis File first."),
          options = base::list(dom = "t"), rownames = FALSE
        ))
      }
      pieces <- base::lapply(names(background), function(sheet) {
        table <- base::as.data.frame(background[[sheet]], stringsAsFactors = FALSE)
        table <- utils::head(table, 100L)
        table$Sheet <- sheet
        table[, c("Sheet", setdiff(names(table), "Sheet")), drop = FALSE]
      })
      DT::datatable(
        do.call(rbind, pieces),
        options = base::list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$genelist_preview <- DT::renderDT({
      ids <- tryCatch(genelist(), error = function(e) character())
      if (base::length(ids) == 0L) {
        return(DT::datatable(
          base::data.frame(Message = "No valid IDs available."),
          options = base::list(dom = "t"), rownames = FALSE
        ))
      }
      DT::datatable(
        base::data.frame(ID = ids, stringsAsFactors = FALSE),
        options = base::list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$compare_select_ui <- shiny::renderUI({
      if (!base::isTRUE(rv$load_success) || base::length(rv$dep_results) == 0) {
        return(NULL)
      }

      choices <- base::names(rv$dep_results)

      if (base::is.null(choices) || base::length(choices) == 0) {
        choices <- base::paste0("Comparison_", base::seq_along(rv$dep_results))
      }

      shiny::selectInput(
        ns("dep_compare"),
        label = "Select DEP comparison",
        choices = choices,
        selected = choices[1]
      )
    })

    shiny::observeEvent(input$dep_compare, {
      shiny::req(rv$dep_results)

      dep_names <- base::names(rv$dep_results)

      if (base::is.null(dep_names) || base::length(dep_names) == 0) {
        return()
      }

      if (!base::is.null(input$dep_compare) && input$dep_compare %in% dep_names) {
        rv$compare_data <- rv$dep_results[[input$dep_compare]]
      }
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$apply_paste, {
      if (base::is.null(input$paste_data) || base::nchar(base::trimws(input$paste_data)) == 0) {
        rv$pasted_genelist <- NULL
        shiny::showNotification("Paste data is empty.", type = "warning")
        return()
      }

      df <- tryCatch(
        utils::read.table(
          text = input$paste_data,
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE,
          check.names = FALSE,
          quote = "",
          comment.char = ""
        ),
        error = function(e) NULL
      )

      if (base::is.null(df) || !"ID" %in% base::colnames(df)) {
        rv$pasted_genelist <- NULL
        shiny::showNotification("Paste data must contain an ID column.", type = "error")
        return()
      }

      ids <- base::as.character(df$ID)
      ids <- ids[!base::is.na(ids) & base::nzchar(ids)]

      if (base::length(ids) == 0) {
        rv$pasted_genelist <- NULL
        shiny::showNotification("No valid IDs found in pasted data.", type = "error")
        return()
      }

      rv$pasted_genelist <- base::unique(ids)
      shiny::showNotification("Paste data applied successfully.", type = "message")
    })

    extract_genelist_ids <- function(df) {
      if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0L) {
        return(character())
      }
      names_clean <- tolower(gsub("[^a-z0-9]", "", base::names(df)))
      candidates <- c("id", "geneid", "proteinid", "gene", "protein")
      index <- match(candidates, names_clean, nomatch = 0L)
      index <- index[index > 0L][1L]
      if (is.na(index) && base::ncol(df) == 1L) index <- 1L
      if (is.na(index)) return(character())
      genes <- trimws(base::as.character(df[[index]]))
      genes <- unique(unlist(base::strsplit(genes, "[,;|]", perl = TRUE), use.names = FALSE))
      genes[!is.na(genes) & nzchar(genes)]
    }

    genelist <- shiny::reactive({
      # An uploaded file takes precedence over the automatic DEP list.
      if (!base::is.null(input$genelist_file)) {
        ext <- base::tolower(tools::file_ext(input$genelist_file$name))

        df <- tryCatch(
          {
            if (ext == "csv") {
              utils::read.csv(
                input$genelist_file$datapath,
                stringsAsFactors = FALSE,
                check.names = FALSE
              )
            } else if (ext == "xlsx") {
              readxl::read_excel(input$genelist_file$datapath)
            } else {
              NULL
            }
          },
          error = function(e) NULL
        )

        genes <- extract_genelist_ids(df)
        if (base::length(genes) > 0) return(genes)
      }

      # Pasted IDs are used only while the manual/paste mode is active.
      if (isFALSE(input$input_mode) &&
          !base::is.null(rv$pasted_genelist) &&
          base::length(rv$pasted_genelist) > 0) {
        return(rv$pasted_genelist)
      }

      # No uploaded (or pasted) list: use the selected DEP comparison.  Read
      # shared_state directly as well, so Analysis works even when the user
      # has not opened the separate LOAD DATA panel in this module.
      dep_data <- rv$compare_data
      if ((base::is.null(dep_data) || !base::is.data.frame(dep_data)) &&
          base::is.list(shared_state$dep_results) &&
          base::length(shared_state$dep_results) > 0L) {
        dep_data <- shared_state$dep_results[[1L]]
      }
      if (base::is.data.frame(dep_data) && base::nrow(dep_data) > 0L) {
        id_candidates <- c("ID", "protein_id", "Protein", "Gene", "gene")
        id_col <- id_candidates[id_candidates %in% base::colnames(dep_data)][1L]
        if (!base::is.na(id_col)) {
          keep <- rep(TRUE, base::nrow(dep_data))
          if ("regulation" %in% base::colnames(dep_data)) {
            keep <- dep_data$regulation != "Not significant"
          } else if ("significant" %in% base::colnames(dep_data)) {
            keep <- !is.na(dep_data$significant) & base::as.logical(dep_data$significant)
          }
          genes <- base::unique(trimws(base::as.character(dep_data[[id_col]][keep])))
          genes <- genes[!base::is.na(genes) & nzchar(genes)]
          if (base::length(genes) > 0) return(genes)
        }
      }

      NULL
    })

    shiny::observeEvent(input$check_file, {
      shiny::req(input$enrichment_analysis_file)

      ext <- base::tolower(tools::file_ext(input$enrichment_analysis_file$name))
      if (ext != "xlsx") {
        rv$file_check_msg <- "❌ Background file must be an .xlsx file."
        rv$background_data <- NULL
        return()
      }

      file <- input$enrichment_analysis_file$datapath

      background <- tryCatch(
        .protvis_read_enrichment_background(file),
        error = function(e) e
      )
      if (inherits(background, "error")) {
        rv$file_check_msg <- paste("❌", conditionMessage(background))
        rv$background_data <- NULL
        return()
      }
      set_background(
        background,
        input$enrichment_analysis_file$name,
        "✅ Background file valid."
      )
    })

    output$file_check_result <- shiny::renderText({
      rv$file_check_msg
    })

    selected_kegg_background <- shiny::reactive({
      shiny::req(input$species)
      shiny::req(rv$background_data)

      background_data <- switch(
        input$species,
        "Plant" = ProtVisDatabase::Plant_KEGG_Background,
        "Animals" = ProtVisDatabase::Animals_KEGG_Background,
        "Bacteria" = ProtVisDatabase::Bacteria_KEGG_Background,
        "Fungi" = ProtVisDatabase::Fungi_KEGG_Background,
        "Eukaryotes" = ProtVisDatabase::Eukaryotes_KEGG_Background,
        "Hsa" = ProtVisDatabase::hsa_KEGG_Background,
        NULL
      )

      shiny::req(background_data)

      bg_df <- tryCatch(
        base::as.data.frame(background_data, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      shiny::req(bg_df)

      if (!"V3" %in% base::colnames(bg_df)) {
        return(rv$background_data$KEGG_background)
      }

      map_id <- bg_df %>%
        tidyr::separate(
          col = V3,
          into = c("Pathway_ID", "Pathway_Name"),
          sep = " ",
          extra = "merge",
          fill = "right"
        ) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), stringr::str_trim)) %>%
        dplyr::pull(Pathway_ID)

      map_id <- base::unique(map_id)
      map_id <- map_id[!base::is.na(map_id) & base::nzchar(map_id)]
      map_id <- base::paste0("map", map_id)

      filtered_bg <- rv$background_data$KEGG_background %>%
        dplyr::filter(TERM %in% map_id)

      if (base::nrow(filtered_bg) == 0) {
        return(rv$background_data$KEGG_background)
      }

      filtered_bg
    })

    shiny::observeEvent(input$run_enrichment_analysis, {
      genes <- unique(trimws(unlist(base::strsplit(base::as.character(genelist()),
                                                   "[,;|]", perl = TRUE), use.names = FALSE)))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      if (base::length(genes) == 0) {
        shiny::showNotification("No valid genelist found.", type = "error")
        return()
      }

      if (base::is.null(rv$background_data)) {
        shiny::showNotification("Please check and load the background file first.", type = "error")
        return()
      }

      rv$go_res <- NULL
      rv$kegg_res <- NULL
      rv$analysis_message <- NULL

      if ("go_analysis" %in% input$choices) {
        t2g.go <- normalise_term2gene(rv$background_data$GO_background)

        t2n.go <- rv$background_data$GO_background %>%
          dplyr::select(TERM, NAME)

        rv$go_res <- tryCatch(
          clusterProfiler::enricher(
            gene = genes,
            TERM2GENE = t2g.go,
            TERM2NAME = t2n.go,
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 1,
            maxGSSize = Inf
          ),
          error = function(e) NULL
        )
      }

      if ("kegg_analysis" %in% input$choices) {
        filtered_bg <- selected_kegg_background()

        t2g.kegg <- normalise_term2gene(filtered_bg)

        t2n.kegg <- filtered_bg %>%
          dplyr::select(TERM, NAME)

        rv$kegg_res <- tryCatch(
          clusterProfiler::enricher(
            gene = genes,
            TERM2GENE = t2g.kegg,
            TERM2NAME = t2n.kegg,
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 1,
            maxGSSize = Inf
          ),
          error = function(e) NULL
        )
      }

      available <- c(
        GO = base::nrow(get_result_df(rv$go_res) %||% data.frame()),
        KEGG = base::nrow(get_result_df(rv$kegg_res) %||% data.frame())
      )
      if (!base::any(available > 0L)) {
        rv$analysis_message <- paste0(
          "No enriched terms were found. Checked genes: ", base::length(genes),
          ". Verify that the ID column uses the same identifiers as the background workbook."
        )
      } else {
        rv$analysis_message <- paste0(
          "Enrichment completed: ",
          if (available[["GO"]] > 0L) paste0("GO ", available[["GO"]], " terms") else "GO 0 terms",
          "; ",
          if (available[["KEGG"]] > 0L) paste0("KEGG ", available[["KEGG"]], " pathways") else "KEGG 0 pathways",
          "."
        )
      }

      shiny::showNotification("Enrichment analysis completed.", type = "message")
    })

    output$analysis_message <- shiny::renderText({
      rv$analysis_message %||% ""
    })

    shiny::observeEvent(input$run_directional_kegg, {
      method <- input$directional_method %||%
        "figure3_archive"
      choices <- directional_available_comparisons()
      selected_comparisons <- intersect(
        input$directional_comparisons %||% character(),
        choices
      )

      if (!length(selected_comparisons)) {
        rv$directional_kegg <- NULL
        rv$directional_kegg_message <-
          "Select at least one DEP comparison."
        shiny::showNotification(
          rv$directional_kegg_message,
          type = "error"
        )
        return()
      }

      if (identical(method, "figure3_archive")) {
        bundled <- tryCatch(
          .protvis_load_builtin_enrichment_background(),
          error = function(e) e
        )
        if (inherits(bundled, "error")) {
          rv$directional_kegg <- NULL
          rv$directional_kegg_message <-
            conditionMessage(bundled)
          shiny::showNotification(
            rv$directional_kegg_message,
            type = "error"
          )
          return()
        }

        result <- tryCatch(
          .protvis_directional_kegg_data(
            dep_results = NULL,
            kegg_background = bundled$KEGG_background,
            evidence_mode = "quantitative",
            comparisons = selected_comparisons,
            top_n = input$directional_top_n %||% 10L,
            method = "figure3_archive",
            pvalue_cutoff = 0.05
          ),
          error = function(e) e
        )
      } else {
        if (base::is.null(rv$background_data)) {
          rv$directional_kegg <- NULL
          rv$directional_kegg_message <-
            "Load or validate a background workbook first."
          shiny::showNotification(
            rv$directional_kegg_message,
            type = "error"
          )
          return()
        }

        bundle <- directional_dep_bundle()
        dep_obj <- bundle$results
        if (base::is.null(dep_obj)) {
          rv$directional_kegg <- NULL
          rv$directional_kegg_message <-
            "Run DEP first, then open this panel."
          shiny::showNotification(
            rv$directional_kegg_message,
            type = "error"
          )
          return()
        }

        evidence <- input$directional_evidence %||%
          "quantitative"
        result <- tryCatch(
          .protvis_directional_kegg_data(
            dep_results = dep_obj,
            kegg_background =
              rv$background_data$KEGG_background,
            presence_absence = bundle$presence,
            evidence_mode = evidence,
            comparisons = selected_comparisons,
            top_n = input$directional_top_n %||% 10L,
            p_adjust_cutoff =
              input$directional_p_adjust %||% 0.05,
            method = "standard_ora",
            pvalue_cutoff = 0.05
          ),
          error = function(e) e
        )
      }

      if (inherits(result, "error")) {
        rv$directional_kegg <- NULL
        rv$directional_kegg_message <-
          conditionMessage(result)
        shiny::showNotification(
          rv$directional_kegg_message,
          type = "error",
          duration = 9
        )
        return()
      }

      rv$directional_kegg <- result

      if (!nrow(result)) {
        rv$directional_kegg_message <- if (identical(
          method,
          "figure3_archive"
        )) {
          paste0(
            "No pathways passed the built-in KEGG ",
            "compareCluster workflow (pvalueCutoff = 0.05; ",
            "qvalueCutoff = 1)."
          )
        } else {
          paste0(
            "No directional KEGG pathways passed BH ≤ ",
            input$directional_p_adjust %||% 0.05,
            ". Comparison-specific tested protein universes ",
            "were applied."
          )
        }
      } else if (identical(
        method,
        "figure3_archive"
      )) {
        rv$directional_kegg_message <- paste0(
          "Built-in KEGG workflow calculated from bundled DEP lists: ",
          length(unique(result$Description)),
          " pathways across ",
          length(unique(result$Cluster)),
          " enriched comparison clusters · ",
          "compareCluster p ≤ 0.05 / q ≤ 1 · ",
          "raw pvalue colour · bundled Enrichmentdb2."
        )
      } else {
        rv$directional_kegg_message <- paste0(
          "Standard directional KEGG completed: ",
          nrow(result),
          " pathways retained across ",
          length(unique(result$Comparison)),
          " DEP comparisons."
        )
      }
    }, ignoreInit = TRUE)

    output$directional_kegg_status <- shiny::renderText({
      rv$directional_kegg_message %||%
        paste0(
          "Select one or more comparisons, then run KEGG."
        )
    })

    render_directional_panel <- function(direction) {
      data <- rv$directional_kegg
      if (base::is.null(data) || !nrow(data)) {
        return(
          ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::annotate(
              "text",
              x = 0,
              y = 0,
              label = rv$directional_kegg_message %||%
                "No directional KEGG result available.",
              colour = "#64748B",
              size = 4
            )
        )
      }

      .protvis_directional_kegg_panel(
        data,
        direction,
        top_n = input$directional_top_n %||% 10L,
        colour_low =
          input$directional_colour_low %||% "#2C7FB8",
        colour_high =
          input$directional_colour_high %||% "#F7FBFF"
      )
    }

    output$directional_kegg_up_plot <- shiny::renderPlot({
      print(render_directional_panel("Upregulated"))
    }, res = 96)

    output$directional_kegg_down_plot <- shiny::renderPlot({
      print(render_directional_panel("Downregulated"))
    }, res = 96)

    output$directional_kegg_plot <- shiny::renderPlot({
      data <- rv$directional_kegg
      if (base::is.null(data) || !nrow(data)) {
        plot(
          0, 0,
          type = "n",
          axes = FALSE,
          xlab = "",
          ylab = ""
        )
        text(
          0, 0,
          rv$directional_kegg_message %||%
            "No directional KEGG result available."
        )
        return(invisible(NULL))
      }

      print(
        .protvis_plot_directional_kegg(
          data,
          top_n = input$directional_top_n %||% 10L,
          colour_low =
            input$directional_colour_low %||% "#2C7FB8",
          colour_high =
            input$directional_colour_high %||% "#F7FBFF"
        )
      )
    }, res = 96)

    output$directional_kegg_table <- DT::renderDT({
      data <- rv$directional_kegg
      if (base::is.null(data) || !nrow(data)) {
        return(DT::datatable(
          base::data.frame(Message = rv$directional_kegg_message %||% "No result available."),
          options = base::list(dom = "t"), rownames = FALSE
        ))
      }
      DT::datatable(data, options = base::list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
    })

    output$download_directional_kegg_pdf <- shiny::downloadHandler(
      filename = function() paste0("directional_KEGG_", Sys.Date(), ".pdf"),
      content = function(file) {
        shiny::req(rv$directional_kegg)
        grDevices::pdf(file, width = 14, height = 7)
        print(.protvis_plot_directional_kegg(
          rv$directional_kegg,
          top_n = input$directional_top_n %||% 10L,
          colour_low =
            input$directional_colour_low %||% "#2C7FB8",
          colour_high =
            input$directional_colour_high %||% "#F7FBFF"
        ))
        grDevices::dev.off()
      }
    )

    output$download_directional_kegg_svg <- shiny::downloadHandler(
      filename = function() paste0("directional_KEGG_", Sys.Date(), ".svg"),
      content = function(file) {
        shiny::req(rv$directional_kegg)
        grDevices::svg(file, width = 14, height = 7, onefile = TRUE)
        print(.protvis_plot_directional_kegg(
          rv$directional_kegg,
          top_n = input$directional_top_n %||% 10L,
          colour_low =
            input$directional_colour_low %||% "#2C7FB8",
          colour_high =
            input$directional_colour_high %||% "#F7FBFF"
        ))
        grDevices::dev.off()
      }
    )

    output$download_directional_kegg_data <- shiny::downloadHandler(
      filename = function() paste0("directional_KEGG_", Sys.Date(), ".csv"),
      content = function(file) {
        shiny::req(rv$directional_kegg)
        utils::write.csv(rv$directional_kegg, file, row.names = FALSE)
      }
    )

    output$go_plot <- shiny::renderPlot({
      go_df <- get_result_df(rv$go_res)

      if (base::is.null(go_df) || base::nrow(go_df) == 0) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(0, 0, rv$analysis_message %||% "No GO enrichment results.")
        return(invisible(NULL))
      }

      if (input$go_plot_type == "bar") {
        p <- plot_enrichment_bar(
          enrich_df = go_df,
          top_n = input$go_top_n,
          fill_color = input$go_color,
          title = "GO Enrichment"
        )
        print(p)
      } else if (input$go_plot_type == "dot") {
        p <- plot_enrichment_dot(
          enrich_df = go_df,
          top_n = input$go_top_n,
          point_color = input$go_color,
          title = "GO Enrichment"
        )
        print(p)
      } else {
        plot_go_circos(
          go_data = go_df,
          top_n = input$go_top_n
        )
      }
    })

    output$kegg_plot <- shiny::renderPlot({
      kegg_df <- get_result_df(rv$kegg_res)

      if (base::is.null(kegg_df) || base::nrow(kegg_df) == 0) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(0, 0, rv$analysis_message %||% "No KEGG enrichment results.")
        return(invisible(NULL))
      }

      if (input$kegg_plot_type == "bar") {
        p <- plot_enrichment_bar(
          enrich_df = kegg_df,
          top_n = input$kegg_top_n,
          fill_color = input$kegg_color,
          title = "KEGG Enrichment"
        )
        print(p)
      } else if (input$kegg_plot_type == "dot") {
        p <- plot_enrichment_dot(
          enrich_df = kegg_df,
          top_n = input$kegg_top_n,
          point_color = input$kegg_color,
          title = "KEGG Enrichment"
        )
        print(p)
      } else {
        plot_go_circos(
          go_data = kegg_df,
          top_n = input$kegg_top_n
        )
      }
    })

    output$go_res_table <- DT::renderDT({
      go_df <- get_result_df(rv$go_res)

      if (base::is.null(go_df) || base::nrow(go_df) == 0) {
        return(DT::datatable(
          base::data.frame(Message = rv$analysis_message %||% "No GO enrichment results."),
          options = base::list(dom = "t"), rownames = FALSE
        ))
      }

      DT::datatable(
        go_df,
        options = base::list(pageLength = 10, scrollX = TRUE)
      )
    })

    output$kegg_res_table <- DT::renderDT({
      kegg_df <- get_result_df(rv$kegg_res)

      if (base::is.null(kegg_df) || base::nrow(kegg_df) == 0) {
        return(DT::datatable(
          base::data.frame(Message = rv$analysis_message %||% "No KEGG enrichment results."),
          options = base::list(dom = "t"), rownames = FALSE
        ))
      }

      DT::datatable(
        kegg_df,
        options = base::list(pageLength = 10, scrollX = TRUE)
      )
    })

    output$download_go_plot <- shiny::downloadHandler(
      filename = function() {
        base::paste0("GO_enrichment_plot_", base::Sys.Date(), ".pdf")
      },
      content = function(file) {
        go_df <- get_result_df(rv$go_res)
        shiny::req(go_df)

        if (input$go_plot_type == "circle") {
          plot_go_circos(
            go_data = go_df,
            top_n = input$go_top_n,
            output_pdf = file
          )
          return()
        }

        grDevices::pdf(file, width = input$go_width, height = input$go_height)

        if (input$go_plot_type == "bar") {
          p <- plot_enrichment_bar(
            enrich_df = go_df,
            top_n = input$go_top_n,
            fill_color = input$go_color,
            title = "GO Enrichment"
          )
          print(p)
        } else if (input$go_plot_type == "dot") {
          p <- plot_enrichment_dot(
            enrich_df = go_df,
            top_n = input$go_top_n,
            point_color = input$go_color,
            title = "GO Enrichment"
          )
          print(p)
        }

        grDevices::dev.off()
      }
    )

    output$download_go_table <- shiny::downloadHandler(
      filename = function() {
        base::paste0("GO_enrichment_table_", base::Sys.Date(), ".csv")
      },
      content = function(file) {
        go_df <- get_result_df(rv$go_res)
        shiny::req(go_df)
        utils::write.csv(go_df, file, row.names = FALSE)
      }
    )

    output$download_kegg_plot <- shiny::downloadHandler(
      filename = function() {
        base::paste0("KEGG_enrichment_plot_", base::Sys.Date(), ".pdf")
      },
      content = function(file) {
        kegg_df <- get_result_df(rv$kegg_res)
        shiny::req(kegg_df)

        if (input$kegg_plot_type == "circle") {
          plot_go_circos(
            go_data = kegg_df,
            top_n = input$kegg_top_n,
            output_pdf = file
          )
          return()
        }

        grDevices::pdf(file, width = input$kegg_width, height = input$kegg_height)

        if (input$kegg_plot_type == "bar") {
          p <- plot_enrichment_bar(
            enrich_df = kegg_df,
            top_n = input$kegg_top_n,
            fill_color = input$kegg_color,
            title = "KEGG Enrichment"
          )
          print(p)
        } else if (input$kegg_plot_type == "dot") {
          p <- plot_enrichment_dot(
            enrich_df = kegg_df,
            top_n = input$kegg_top_n,
            point_color = input$kegg_color,
            title = "KEGG Enrichment"
          )
          print(p)
        }

        grDevices::dev.off()
      }
    )

    output$download_kegg_table <- shiny::downloadHandler(
      filename = function() {
        base::paste0("KEGG_enrichment_table_", base::Sys.Date(), ".csv")
      },
      content = function(file) {
        kegg_df <- get_result_df(rv$kegg_res)
        shiny::req(kegg_df)
        utils::write.csv(kegg_df, file, row.names = FALSE)
      }
    )
  })
}
