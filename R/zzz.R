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

# Show a compact package-loading progress indicator when ProtVis is attached
# with library(ProtVis). Heavy feature packages remain lazy-loaded.
.onAttach <- function(libname, pkgname) {
  show_progress <- isTRUE(getOption("protvis.startup.progress", TRUE))

  if (show_progress && interactive()) {
    base::cat("ProtVis loading ")
    pb <- utils::txtProgressBar(
      min = 0,
      max = 100,
      initial = 0,
      style = 3,
      width = 28
    )

    # Keep the bar tied to lightweight attach-time checks rather than adding
    # artificial delays solely for animation.
    utils::setTxtProgressBar(pb, 35)

    ProtVisDatabase::protvis_database_path(
      "extdata",
      "manifest.tsv",
      must_work = TRUE
    )
    utils::setTxtProgressBar(pb, 75)

    version <- as.character(utils::packageVersion(pkgname))
    utils::setTxtProgressBar(pb, 100)
    base::close(pb)

    packageStartupMessage(
      sprintf("ProtVis %s loaded successfully.", version)
    )
  } else if (show_progress) {
    packageStartupMessage(
      sprintf(
        "ProtVis %s loaded successfully.",
        as.character(utils::packageVersion(pkgname))
      )
    )
  }

  invisible()
}

