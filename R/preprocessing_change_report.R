# Reports use finite values: zeros and negative values are valid after
# transformation and normalization, unlike raw-intensity detection counts.
.protvis_matrix_change_report <- function(before, after, method) {
  clean <- function(x) {
    x <- as.data.frame(x, check.names = FALSE)
    if ("ID" %in% names(x)) {
      rownames(x) <- as.character(x$ID)
      x$ID <- NULL
    }
    m <- as.matrix(x)
    suppressWarnings(storage.mode(m) <- "double")
    m
  }
  a <- clean(before); b <- clean(after)
  stats <- function(m) {
    v <- m[is.finite(m)]
    c(nrow(m), ncol(m), sum(is.finite(m)), sum(!is.finite(m)),
      if (length(m)) 100 * mean(!is.finite(m)) else NA_real_,
      if (length(v)) median(v) else NA_real_,
      if (length(v)) min(v) else NA_real_,
      if (length(v)) max(v) else NA_real_)
  }
  summary <- data.frame(Metric = c("Protein rows", "Samples", "Finite values",
    "Missing / non-finite values", "Missing / non-finite (%)", "Median", "Minimum", "Maximum"),
    Before = stats(a), After = stats(b))
  summary$Change <- summary$After - summary$Before
  samples <- intersect(colnames(a), colnames(b))
  rows <- intersect(rownames(a), rownames(b))
  per_sample <- do.call(rbind, lapply(samples, function(sample) {
    old <- a[rows, sample]; new <- b[rows, sample]
    valid_old <- is.finite(old); valid_new <- is.finite(new)
    med <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
    data.frame(Sample = sample, Missing_before = sum(!is.finite(a[, sample])),
      Missing_after = sum(!is.finite(b[, sample])),
      Filled = sum(!valid_old & valid_new), Masked = sum(valid_old & !valid_new),
      Changed_observed = sum(old[valid_old & valid_new] != new[valid_old & valid_new]),
      Median_before = med(a[, sample]), Median_after = med(b[, sample]))
  }))
  list(method = method, summary = summary, per_sample = per_sample,
       created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

.protvis_save_change_report <- function(dataset, stage, before, after, method) {
  report <- .protvis_matrix_change_report(before, after, method)
  key <- paste0(stage, "_report")
  dataset$analysis_results[[key]] <- report
  history_key <- paste0(stage, "_reports")
  dataset$analysis_results[[history_key]] <- c(
    dataset$analysis_results[[history_key]] %||% list(), list(report))
  dataset
}

.protvis_change_report_ui <- function(ns) {
  shiny::tagList(shiny::uiOutput(ns("processing_report_method")),
    DT::DTOutput(ns("processing_report_summary")),
    shiny::tags$h5("Changes by sample"), DT::DTOutput(ns("processing_report_samples")))
}

.protvis_change_report_outputs <- function(output, shared_state, stage) {
  report <- shiny::reactive({
    dataset <- shared_state$dataset
    if (!inherits(dataset, "ProtVis_dataset")) return(NULL)
    dataset$analysis_results[[paste0(stage, "_report")]]
  })
  output$processing_report_method <- shiny::renderUI({
    r <- report()
    if (is.null(r)) return(shiny::tags$p("Complete and save this processing step to generate its report."))
    shiny::tagList(shiny::tags$p(paste("Method:", r$method)),
      shiny::tags$small("Finite values include zero and negative values. Changes in value scale should be interpreted according to the selected method."))
  })
  output$processing_report_summary <- DT::renderDT({
    shiny::req(report())
    DT::formatRound(DT::datatable(report()$summary, rownames = FALSE,
      options = list(dom = "t", scrollX = TRUE)), c("Before", "After", "Change"), digits = 3)
  })
  output$processing_report_samples <- DT::renderDT({
    shiny::req(report())
    DT::datatable(report()$per_sample, rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10))
  })
  invisible(NULL)
}
