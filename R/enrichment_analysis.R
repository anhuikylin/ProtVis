
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


# Built-in source for the archived maize-teosinte Figure 3C-D KEGG panel.
# The 29 displayed pathway-stage points were reconstructed from the archived
# directional DEP/KEGG workflow and checked against
# 03.Maize_Teosinte_Jul02_2024/04.result/01.Publish_figures/KEGG enrichment.png.
.protvis_maize_teosinte_kegg_data <- function() {
  file_name <- "maize_teosinte_kegg_figure3_cd.csv"
  path <- system.file(
    "extdata",
    file_name,
    package = "ProtVis"
  )

  if (!nzchar(path)) {
    candidates <- c(
      file.path("inst", "extdata", file_name),
      file.path(getwd(), "inst", "extdata", file_name)
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) path <- hit[[1L]]
  }

  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "Built-in maize-teosinte KEGG reproduction data are unavailable.",
      call. = FALSE
    )
  }

  out <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  numeric_cols <- intersect(c("pvalue", "p.adjust", "Count"), names(out))
  out[numeric_cols] <- lapply(out[numeric_cols], as.numeric)
  .protvis_validate_maize_teosinte_kegg_data(out)
  out
}


.protvis_validate_maize_teosinte_kegg_data <- function(data) {
  required <- c(
    "Panel", "Direction", "Cluster", "Cluster_label", "ID",
    "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Built-in maize-teosinte KEGG data are missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(data) != 29L ||
      sum(data$Panel == "C") != 18L ||
      sum(data$Panel == "D") != 11L) {
    stop(
      "Built-in maize-teosinte KEGG data failed the 29-point ",
      "(Panel C = 18; Panel D = 11) integrity check.",
      call. = FALSE
    )
  }

  expected_c <- c(
    "Pentose and glucuronate interconversions",
    "Phenylpropanoid biosynthesis",
    "Galactose metabolism",
    "Lipid biosynthesis proteins",
    "Fatty acid biosynthesis"
  )
  expected_d <- c(
    "Phenylpropanoid biosynthesis",
    "Metabolism of terpenoids and polyketides",
    "Monoterpenoid biosynthesis",
    "Glutathione metabolism",
    "Photosynthesis"
  )
  observed_c <- unique(data$Description[data$Panel == "C"])
  observed_d <- unique(data$Description[data$Panel == "D"])
  if (!setequal(observed_c, expected_c) ||
      !setequal(observed_d, expected_d)) {
    stop(
      "Built-in maize-teosinte KEGG pathway labels failed integrity checks.",
      call. = FALSE
    )
  }

  if (any(!is.finite(data$pvalue)) ||
      any(!is.finite(data$Count)) ||
      any(data$pvalue <= 0) ||
      any(data$Count <= 0) ||
      anyDuplicated(data[c("Panel", "Cluster", "Description")])) {
    stop(
      "Built-in maize-teosinte KEGG numeric/key values failed integrity checks.",
      call. = FALSE
    )
  }

  # Publication panel D has no significant Root_V1.V2 pathway point.
  if (any(data$Panel == "D" & data$Cluster == "Root_V1.V2")) {
    stop(
      "Panel D integrity check failed: Root_V1.V2 should be absent.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.protvis_maize_teosinte_kegg_panel <- function(data, panel = c("C", "D")) {
  panel <- match.arg(panel)
  df <- data[data$Panel == panel, , drop = FALSE]

  if (!nrow(df)) {
    stop("No built-in KEGG data are available for panel ", panel, ".",
         call. = FALSE)
  }

  if (identical(panel, "C")) {
    x_levels <- c(
      "Root_VE", "Root_V1.V2", "Root_V4", "Leaf_VE-V2", "Leaf_V4-V8"
    )
    y_levels <- c(
      "Pentose and glucuronate interconversions",
      "Phenylpropanoid biosynthesis",
      "Galactose metabolism",
      "Lipid biosynthesis proteins",
      "Fatty acid biosynthesis"
    )
    p_breaks <- c(0.0005, 0.0010, 0.0015)
    p_labels <- c("0.0005", "0.0010", "0.0015")
    count_breaks <- c(20, 30, 40)
    title <- expression(italic("Zea mays ssp. mays") ~ "– enriched")
  } else {
    x_levels <- c("Root_VE", "Root_V4", "Leaf_VE-V2", "Leaf_V4-V8")
    y_levels <- c(
      "Phenylpropanoid biosynthesis",
      "Metabolism of terpenoids and polyketides",
      "Monoterpenoid biosynthesis",
      "Glutathione metabolism",
      "Photosynthesis"
    )
    p_breaks <- c(0.0003, 0.0006, 0.0009)
    p_labels <- c("0.0003", "0.0006", "0.0009")
    count_breaks <- c(5, 20, 30)
    title <- expression(italic("Zea mays ssp. mexicana") ~ "– enriched")
  }

  df$Cluster_label <- factor(df$Cluster_label, levels = x_levels)
  df$Description <- factor(df$Description, levels = rev(y_levels))

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Cluster_label,
      y = Description,
      size = Count,
      fill = pvalue
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.45
    ) +
    ggplot2::scale_fill_gradient(
      low = "red",
      high = "blue",
      breaks = p_breaks,
      labels = p_labels,
      name = "pvalue"
    ) +
    ggplot2::scale_size_continuous(
      range = c(3.5, 10.5),
      breaks = count_breaks,
      name = "Count"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        order = 1,
        reverse = TRUE
      ),
      size = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          fill = "white",
          colour = "black"
        )
      )
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 9,
        colour = "black",
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = 9,
        colour = "black"
      ),
      plot.title = ggplot2::element_text(
        size = 12,
        hjust = 0.5,
        colour = "black"
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
        linewidth = 0.8
      ),
      panel.grid.major = ggplot2::element_line(
        colour = "#e7e7e7",
        linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.45,
        colour = "black"
      ),
      legend.text = ggplot2::element_text(size = 9, colour = "black"),
      legend.title = ggplot2::element_text(size = 10, colour = "black"),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
}


plot_maize_teosinte_kegg_reproduction <- function(data = NULL) {
  if (is.null(data)) {
    data <- .protvis_maize_teosinte_kegg_data()
  }

  panel_c <- .protvis_maize_teosinte_kegg_panel(data, "C")
  panel_d <- .protvis_maize_teosinte_kegg_panel(data, "D")

  patchwork::wrap_plots(
    panel_c,
    panel_d,
    ncol = 2,
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
}

.protvis_load_builtin_enrichment_background <- function() {
  locate <- function(file_name) {
    candidates <- c(
      system.file("extdata", file_name, package = "ProtVis"),
      file.path("inst", "extdata", file_name),
      file.path(getwd(), "inst", "extdata", file_name)
    )
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    if (!length(candidates)) {
      stop("The built-in maize-teosinte enrichment background is unavailable.",
           call. = FALSE)
    }
    candidates[[1L]]
  }
  read_table <- function(file_name) {
    utils::read.delim(
      locate(file_name), sep = "\t", header = TRUE, quote = "",
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

.protvis_directional_kegg_data <- function(
    dep_results, kegg_background, presence_absence = list(),
    evidence_mode = c("quantitative", "all_retained"),
    comparisons = NULL, top_n = 5L, p_adjust_cutoff = 0.05) {
  if (!is.list(dep_results) || !length(dep_results)) return(data.frame())
  evidence_mode <- match.arg(evidence_mode)
  if (!is.list(presence_absence)) presence_absence <- list()
  if (!is.null(comparisons)) {
    comparisons <- intersect(as.character(comparisons), names(dep_results))
    dep_results <- dep_results[comparisons]
  }
  background <- as.data.frame(kegg_background, stringsAsFactors = FALSE)
  if (!all(c("TERM", "GENE", "NAME") %in% names(background))) {
    stop("KEGG background must contain TERM, GENE, and NAME columns.",
         call. = FALSE)
  }
  t2g <- unique(background[, c("TERM", "GENE"), drop = FALSE])
  t2g$TERM <- trimws(as.character(t2g$TERM))
  t2g$GENE <- trimws(as.character(t2g$GENE))
  t2g <- t2g[nzchar(t2g$TERM) & nzchar(t2g$GENE), , drop = FALSE]
  t2n <- unique(background[, c("TERM", "NAME"), drop = FALSE])
  t2n$TERM <- trimws(as.character(t2n$TERM))
  t2n$NAME <- trimws(as.character(t2n$NAME))
  t2n <- t2n[!duplicated(t2n$TERM) & nzchar(t2n$NAME), , drop = FALSE]

  # MaxQuant protein groups can contain several identifiers separated by
  # semicolons.  Enrichment must test each mapped identifier, rather than
  # treating the entire protein-group string as a new, unmatched ID.
  expand_ids <- function(x) {
    ids <- unlist(strsplit(as.character(x), "[,;|]", perl = TRUE), use.names = FALSE)
    ids <- trimws(ids)
    ids <- sub("^CON__", "", ids)
    unique(ids[!is.na(ids) & nzchar(ids)])
  }

  top_n <- max(1L, as.integer(top_n))
  p_adjust_cutoff <- as.numeric(p_adjust_cutoff)
  rows <- lapply(names(dep_results), function(comparison) {
    result <- as.data.frame(dep_results[[comparison]], stringsAsFactors = FALSE)
    id_col <- c("ID", "protein_id", "Protein", "Gene")[
      c("ID", "protein_id", "Protein", "Gene") %in% names(result)
    ][1L]
    if (is.na(id_col) || !"regulation" %in% names(result)) return(NULL)
    tested <- expand_ids(result[[id_col]])
    presence <- presence_absence[[comparison]] %||% data.frame()
    presence_id_col <- c("ID", "protein_id", "Protein", "Gene")[
      c("ID", "protein_id", "Protein", "Gene") %in% names(presence)
    ][1L]
    presence_ids <- if (!is.na(presence_id_col)) {
      expand_ids(presence[[presence_id_col]])
    } else {
      character()
    }
    if (identical(evidence_mode, "all_retained")) {
      tested <- unique(c(tested, presence_ids))
    }
    if (!length(tested)) return(NULL)
    first_label <- function(x, fallback) {
      x <- unique(trimws(as.character(x)))
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x)) x[[1L]] else fallback
    }
    labels <- c(
      Upregulated = first_label(result$Group1 %||% character(), "Group 1"),
      Downregulated = first_label(result$Group2 %||% character(), "Group 2")
    )
    lapply(names(labels), function(direction) {
      genes <- expand_ids(
        result[[id_col]][as.character(result$regulation) == direction]
      )
      if (identical(evidence_mode, "all_retained") && nrow(presence) &&
          "Pattern" %in% names(presence) && !is.na(presence_id_col)) {
        pattern <- if (identical(direction, "Upregulated")) {
          "Detected in Group1 only"
        } else {
          "Detected in Group2 only"
        }
        genes <- unique(c(
          genes,
          expand_ids(presence[[presence_id_col]][presence$Pattern == pattern])
        ))
      }
      if (!length(genes)) return(NULL)
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
      out <- tryCatch(as.data.frame(enriched@result), error = function(e) NULL)
      if (is.null(out) || !nrow(out)) return(NULL)
      out <- out[is.finite(out$p.adjust) & out$p.adjust <= p_adjust_cutoff, , drop = FALSE]
      if (!nrow(out)) return(NULL)
      out <- out[order(out$p.adjust, out$pvalue, -out$Count), , drop = FALSE]
      out <- utils::head(out, top_n)
      out$Comparison <- comparison
      out$Cluster <- .protvis_dep_stage_label(labels[[direction]])
      out$Direction <- direction
      out$Direction_label <- paste0(labels[[direction]], " higher")
      out$Evidence <- if (identical(evidence_mode, "all_retained")) {
        "All retained proteins"
      } else {
        "Evidence 1 · Quantitative DEP"
      }
      out
    })
  })
  rows <- unlist(rows, recursive = FALSE)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

.protvis_directional_kegg_panel <- function(data, direction) {
  df <- data[as.character(data$Direction) == direction, , drop = FALSE]
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = paste("No significant", tolower(direction), "KEGG pathways.")))
  }
  x_levels <- unique(df$Cluster)
  y_levels <- unique(df$Description[order(df$p.adjust, df$pvalue)])
  df$Cluster <- factor(df$Cluster, levels = x_levels)
  df$Description <- factor(df$Description, levels = rev(y_levels))
  label <- unique(df$Direction_label)
  title <- if (length(label) == 1L) {
    paste(label, "– enriched")
  } else {
    paste(direction, "proteins – enriched")
  }
  ggplot2::ggplot(df, ggplot2::aes(
    x = Cluster, y = Description, size = Count, fill = pvalue
  )) +
    ggplot2::geom_point(shape = 21, colour = "black", stroke = 0.45) +
    ggplot2::scale_fill_gradient(low = "red", high = "blue", name = "p value") +
    ggplot2::scale_size_continuous(range = c(3.5, 10.5), name = "Count") +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1, reverse = TRUE),
      size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "white"))
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, colour = "black"),
      axis.text.y = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "italic"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.protvis_plot_directional_kegg <- function(data) {
  evidence <- unique(as.character(data$Evidence))
  evidence <- evidence[!is.na(evidence) & nzchar(evidence)]
  patchwork::wrap_plots(
    .protvis_directional_kegg_panel(data, "Upregulated"),
    .protvis_directional_kegg_panel(data, "Downregulated"),
    ncol = 2
  ) + patchwork::plot_annotation(
    title = "Directional KEGG enrichment across DEP comparisons",
    subtitle = paste0(
      if (length(evidence)) evidence[[1L]] else "Selected differential evidence",
      ": each comparison uses its retained tested proteins as the enrichment universe."
    )
  )
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
        bslib::card(
          bslib::card_header("Directional KEGG enrichment across DEP comparisons"),
          bslib::card_body(
            shiny::tags$p(
              "Choose one differential-evidence source and one or more DEP comparisons. Upregulated and downregulated proteins are enriched separately; each comparison uses its own retained protein universe.",
              class = "text-muted"
            ),
            bslib::layout_columns(
              col_widths = c(4, 4, 2, 2),
              shiny::radioButtons(
                ns("directional_evidence"),
                "Differential evidence",
                choices = c(
                  "Evidence 1 · Quantitative DEP" = "quantitative",
                  "All retained proteins" = "all_retained"
                ),
                selected = "quantitative",
                inline = TRUE
              ),
              shiny::uiOutput(ns("directional_comparisons_ui")),
              shiny::numericInput(
                ns("directional_top_n"), "Top pathways", 5,
                min = 1, max = 20
              ),
              shiny::numericInput(
                ns("directional_p_adjust"), "BH FDR cutoff", 0.05,
                min = 0, max = 1, step = 0.01
              )
            ),
            bslib::layout_columns(
              col_widths = c(3, 5, 2, 2),
              shiny::actionButton(
                ns("load_maize_teosinte_background"),
                "LOAD BUILT-IN BACKGROUND",
                class = "btn btn-outline-primary fw-bold w-100"
              ),
              shiny::div(
                class = "pt-2 text-muted",
                shiny::textOutput(ns("directional_kegg_status"))
              ),
              shiny::actionButton(
                ns("run_directional_kegg"), "RUN KEGG",
                class = "btn btn-primary fw-bold w-100"
              ),
              shiny::tags$div(
                class = "d-grid gap-2",
                shiny::downloadButton(
                  ns("download_directional_kegg_pdf"), "PDF",
                  class = "btn btn-outline-secondary w-100"
                ),
                shiny::downloadButton(
                  ns("download_directional_kegg_data"), "CSV",
                  class = "btn btn-outline-secondary w-100 mt-2"
                )
              )
            ),
            shiny::br(), shiny::br(),
            shiny::tabsetPanel(
              shiny::tabPanel("Figure", shiny::plotOutput(ns("directional_kegg_plot"), height = "580px")),
              shiny::tabPanel("Result table", DT::DTOutput(ns("directional_kegg_table")))
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
      bundle <- directional_dep_bundle()
      choices <- names(bundle$results %||% list())
      if (!length(choices)) {
        return(shiny::helpText("Run DEP to select comparisons."))
      }
      selected <- isolate(input$directional_comparisons)
      selected <- intersect(selected %||% choices, choices)
      if (!length(selected)) selected <- choices
      shiny::checkboxGroupInput(
        ns("directional_comparisons"),
        "DEP comparisons",
        choices = choices,
        selected = selected,
        inline = TRUE
      )
    })

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
      if (base::is.null(rv$background_data)) {
        rv$directional_kegg <- NULL
        rv$directional_kegg_message <- "Load or validate a background workbook first."
        shiny::showNotification(rv$directional_kegg_message, type = "error")
        return()
      }

      bundle <- directional_dep_bundle()
      dep_obj <- bundle$results
      selected_comparisons <- intersect(
        input$directional_comparisons %||% character(),
        names(dep_obj %||% list())
      )
      if (base::is.null(dep_obj)) {
        rv$directional_kegg <- NULL
        rv$directional_kegg_message <- "Run DEP first, then open this panel."
        shiny::showNotification(rv$directional_kegg_message, type = "error")
        return()
      }
      if (!length(selected_comparisons)) {
        rv$directional_kegg <- NULL
        rv$directional_kegg_message <- "Select at least one DEP comparison."
        shiny::showNotification(rv$directional_kegg_message, type = "error")
        return()
      }

      result <- tryCatch(
        .protvis_directional_kegg_data(
          dep_obj,
          rv$background_data$KEGG_background,
          presence_absence = bundle$presence,
          evidence_mode = input$directional_evidence %||% "quantitative",
          comparisons = selected_comparisons,
          top_n = input$directional_top_n,
          p_adjust_cutoff = input$directional_p_adjust
        ),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        rv$directional_kegg <- NULL
        rv$directional_kegg_message <- conditionMessage(result)
        shiny::showNotification(rv$directional_kegg_message, type = "error")
        return()
      }
      rv$directional_kegg <- result
      if (!nrow(result)) {
        rv$directional_kegg_message <- paste0(
          "No directional KEGG pathways passed BH ≤ ", input$directional_p_adjust,
          ". The comparison-specific tested protein universes were still applied."
        )
      } else {
        rv$directional_kegg_message <- paste0(
          "Directional KEGG completed: ", nrow(result),
          " pathways retained across ", length(unique(result$Comparison)),
          " DEP comparisons."
        )
      }
    }, ignoreInit = TRUE)

    output$directional_kegg_status <- shiny::renderText({
      rv$directional_kegg_message %||%
        "Load the built-in background, then run this after DEP."
    })

    output$directional_kegg_plot <- shiny::renderPlot({
      data <- rv$directional_kegg
      if (base::is.null(data) || !nrow(data)) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(0, 0, rv$directional_kegg_message %||% "No directional KEGG result available.")
        return(invisible(NULL))
      }
      print(.protvis_plot_directional_kegg(data))
    })

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
        print(.protvis_plot_directional_kegg(rv$directional_kegg))
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
