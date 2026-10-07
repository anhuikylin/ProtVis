# Explicit scale tracking prevents standardized/raw values being reported as log2FC.
.protvis_set_expression_scale <- function(dataset, scale) {
  dataset$metadata$expression_scale <- scale
  dataset
}

.protvis_require_log2_scale <- function(dataset, input_scale = "auto") {
  explicit <- tolower(as.character(input_scale))
  scale <- dataset$metadata$expression_scale %||% NULL
  if (is.null(scale)) {
    normalization <- dataset$analysis_results$normalization$method %||% ""
    transformation <- .protvis_norm_transformation_method(dataset)
    scale <- if (normalization %in% c("zscore", "z_score", "standardize", "vsn")) {
      "other"
    } else if (tolower(transformation) %in% c("log2", "maxquant_log2", "maxquant_recommended")) {
      "log2"
    } else if (tolower(transformation) %in% c("unknown", "none", "identity", "")) {
      "unknown"
    } else "other"
  }
  # Explicit declarations are permitted for imported matrices with no scale record.
  if (identical(scale, "unknown") && identical(explicit, "log2")) scale <- "log2"
  if (!identical(scale, "log2")) {
    stop(paste0("Differential analysis requires log2 intensities; current scale: ",
      scale, ". Transform raw intensities first, or declare input_scale='log2' ",
      "for a verified imported log2 matrix. Standardized/ln/log10 values cannot be used."),
      call. = FALSE)
  }
  invisible(TRUE)
}
