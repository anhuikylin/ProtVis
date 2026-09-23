# Apply the upload limit as soon as the package namespace is loaded. This is
# required for launchers that do not call run_ProtVis() explicitly.
.onLoad <- function(libname, pkgname) {
  current <- getOption("shiny.maxRequestSize", 0)
  if (!is.numeric(current) || length(current) != 1L || !is.finite(current)) {
    current <- 0
  }
  options(shiny.maxRequestSize = max(current, 2 * 1024^3))
}

# Load feature-specific packages only when the corresponding module is used.
# Keeping these packages out of the package-level import graph makes
# library(ProtVis) quiet and avoids unrelated startup messages/warnings.
.protvis_require_optional <- function(package, feature = package) {
  available <- suppressWarnings(
    suppressPackageStartupMessages(
      requireNamespace(package, quietly = TRUE)
    )
  )

  if (!isTRUE(available)) {
    stop(
      "The optional package '", package,
      "' is required for ", feature,
      ". Reinstall ProtVis with dependencies = TRUE or install the package separately.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
