# Clean installer for ProtVis -------------------------------------------------
# Run this file in a fresh R session before calling library(ProtVis).

if ("ProtVis" %in% loadedNamespaces() ||
    "ProtVisDatabase" %in% loadedNamespaces()) {
  stop(
    "ProtVis or ProtVisDatabase is loaded in this R session. Restart R, do ",
    "not call library(ProtVis), and run this installer again.",
    call. = FALSE
  )
}

protvis_library <- Sys.getenv("PROTVIS_LIBRARY", unset = "")
if (!nzchar(protvis_library)) {
  existing_install <- find.package("ProtVis", quiet = TRUE)
  protvis_library <- if (nzchar(existing_install)) {
    dirname(existing_install)
  } else {
    .libPaths()[[1L]]
  }
}
protvis_library <- normalizePath(
  protvis_library, winslash = "/", mustWork = FALSE
)

if (!dir.exists(protvis_library)) {
  dir.create(protvis_library, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(protvis_library) || file.access(protvis_library, 2L) != 0L) {
  stop("R library is not writable: ", protvis_library, call. = FALSE)
}
.libPaths(unique(c(protvis_library, .libPaths())))
cran_repo <- unname(getOption("repos")["CRAN"])
if (length(cran_repo) != 1L || is.na(cran_repo) || !nzchar(cran_repo) ||
    identical(cran_repo, "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

protvis_target <- file.path(protvis_library, "ProtVis")
protvis_lock <- file.path(protvis_library, "00LOCK-ProtVis")
protvis_database_target <- file.path(protvis_library, "ProtVisDatabase")
protvis_database_lock <- file.path(protvis_library, "00LOCK-ProtVisDatabase")

# Guard the recursive deletion so it can never target the library itself.
if (!identical(basename(protvis_target), "ProtVis") ||
    identical(normalizePath(protvis_target, winslash = "/", mustWork = FALSE),
              protvis_library)) {
  stop("Refusing to clean an unexpected installation path.", call. = FALSE)
}
if (!identical(basename(protvis_database_target), "ProtVisDatabase") ||
    identical(normalizePath(protvis_database_target, winslash = "/", mustWork = FALSE),
              protvis_library)) {
  stop("Refusing to clean an unexpected companion-package path.", call. = FALSE)
}

if (dir.exists(protvis_target)) {
  unlink(protvis_target, recursive = TRUE, force = TRUE)
}
if (dir.exists(protvis_lock)) {
  unlink(protvis_lock, recursive = TRUE, force = TRUE)
}
if (dir.exists(protvis_database_target)) {
  unlink(protvis_database_target, recursive = TRUE, force = TRUE)
}
if (dir.exists(protvis_database_lock)) {
  unlink(protvis_database_lock, recursive = TRUE, force = TRUE)
}
if (dir.exists(protvis_target) || dir.exists(protvis_lock) ||
    dir.exists(protvis_database_target) || dir.exists(protvis_database_lock)) {
  stop(
    "Could not remove the previous ProtVis installation. Close all R/RStudio ",
    "sessions that use either ProtVis package and retry.",
    call. = FALSE
  )
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = protvis_library)
}

remotes::install_github(
  "anhuikylin/ProtVisDatabase",
  lib = protvis_library,
  force = TRUE,
  upgrade = "never",
  dependencies = FALSE
)

remotes::install_github(
  "anhuikylin/ProtVis",
  ref = "dev",
  lib = protvis_library,
  force = TRUE,
  upgrade = "never",
  dependencies = TRUE
)

# A partially restored R library can contain bslib without its cachem runtime
# dependency. Repair that state before the clean-process verification below.
if (!requireNamespace("cachem", quietly = TRUE)) {
  install.packages("cachem", lib = protvis_library, dependencies = TRUE)
}
tryCatch(
  loadNamespace("bslib"),
  error = function(error) {
    install.packages(
      c("cachem", "bslib"),
      lib = protvis_library,
      dependencies = TRUE
    )
    loadNamespace("bslib")
  }
)

database_path <- find.package(
  "ProtVisDatabase", lib.loc = protvis_library, quiet = TRUE
)
if (!nzchar(database_path)) {
  stop("ProtVisDatabase installation did not produce an installed package.",
       call. = FALSE)
}

# Verify every bundled resource against the companion-package manifest before
# ProtVis starts. This catches truncated binary fixtures (especially .xlsx)
# immediately instead of failing later when a Built-in example is loaded.
database_manifest <- ProtVisDatabase::protvis_database_manifest()
database_resources <- file.path(
  database_path, "extdata", database_manifest$path
)
database_exists <- file.exists(database_resources)
database_sizes <- suppressWarnings(file.info(database_resources)$size)
database_expected_sizes <- suppressWarnings(as.numeric(database_manifest$bytes))
database_bad <- !database_exists |
  is.na(database_sizes) |
  is.na(database_expected_sizes) |
  database_sizes != database_expected_sizes

if (any(database_bad)) {
  bad_rows <- which(database_bad)
  details <- paste0(
    database_manifest$path[bad_rows],
    " (expected ", database_expected_sizes[bad_rows],
    " bytes; installed ",
    ifelse(database_exists[bad_rows], database_sizes[bad_rows], "missing"),
    ")"
  )
  stop(
    "ProtVisDatabase installation verification failed for bundled resources:\n",
    paste0(" - ", details, collapse = "\n"),
    "\nRe-run the clean installer so the companion package is downloaded again.",
    call. = FALSE
  )
}

xlsx_resources <- database_resources[
  grepl("\\.xlsx$", database_manifest$path, ignore.case = TRUE)
]
if (length(xlsx_resources)) {
  xlsx_ok <- vapply(
    xlsx_resources,
    function(path) {
      tryCatch(
        {
          listing <- suppressWarnings(utils::unzip(path, list = TRUE))
          is.data.frame(listing) && nrow(listing) > 0L
        },
        error = function(error) FALSE
      )
    },
    logical(1L)
  )
  if (!all(xlsx_ok)) {
    stop(
      "ProtVisDatabase contains an unreadable XLSX resource: ",
      paste(basename(xlsx_resources[!xlsx_ok]), collapse = ", "),
      ". Re-run the clean installer.",
      call. = FALSE
    )
  }
}

installed_path <- find.package(
  "ProtVis", lib.loc = protvis_library, quiet = TRUE
)
if (!nzchar(installed_path)) {
  stop("ProtVis installation did not produce an installed package.",
       call. = FALSE)
}

lazy_databases <- file.path(installed_path, "R", c("ProtVis.rdb", "ProtVis.rdx"))
if (!all(file.exists(lazy_databases))) {
  stop(
    "Installation verification failed: the package lazy-load database is ",
    "incomplete. Please report the R version and installation log.",
    call. = FALSE
  )
}

# Verify the new database in a clean process. Calling body() forces the
# exported launcher promise to be decompressed instead of merely registering
# the namespace, which detects the exact R_decompress1 failure reported by
# interrupted or overwritten installations.
probe_file <- tempfile("protvis_install_probe_", fileext = ".R")
on.exit(unlink(probe_file, force = TRUE), add = TRUE)
writeLines(
  c(
    paste0("protvis_library <- ", deparse(protvis_library)),
    # --vanilla deliberately ignores the user's R_LIBS_USER/Rprofile. Restore
    # the library search path used for installation so transitive dependencies
    # (such as bslib -> cachem) remain visible to this independent process.
    paste0(".libPaths(", paste(capture.output(dput(.libPaths())), collapse = "\n"), ")"),
    "library(ProtVisDatabase, lib.loc = protvis_library)",
    "stopifnot(nzchar(ProtVisDatabase::protvis_database_path('extdata', 'manifest.tsv', must_work = TRUE)))",
    "manifest <- ProtVisDatabase::protvis_database_manifest()",
    "database_path <- find.package('ProtVisDatabase', lib.loc = protvis_library)",
    "resources <- file.path(database_path, 'extdata', manifest$path)",
    "stopifnot(all(file.exists(resources)))",
    "stopifnot(all(file.info(resources)$size == as.numeric(manifest$bytes)))",
    "xlsx <- resources[endsWith(tolower(manifest$path), '.xlsx')]",
    "if (length(xlsx)) stopifnot(all(vapply(xlsx, function(path) tryCatch({ z <- suppressWarnings(utils::unzip(path, list = TRUE)); is.data.frame(z) && nrow(z) > 0L }, error = function(e) FALSE), logical(1L))))",
    "library(ProtVis, lib.loc = protvis_library)",
    "launcher <- getExportedValue('ProtVis', 'run_ProtVis')",
    "stopifnot(is.function(launcher), is.call(body(launcher)) || is.expression(body(launcher)))"
  ),
  probe_file,
  useBytes = TRUE
)
rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
probe_output <- system2(
  rscript,
  c("--vanilla", shQuote(probe_file)),
  stdout = TRUE,
  stderr = TRUE
)
probe_status <- attr(probe_output, "status")
if (is.null(probe_status)) probe_status <- 0L
if (!identical(as.integer(probe_status), 0L)) {
  stop(
    "ProtVis was installed but failed the clean-process load test:\n",
    paste(probe_output, collapse = "\n"),
    call. = FALSE
  )
}

installed_description <- read.dcf(file.path(installed_path, "DESCRIPTION"))
message(
  "ProtVis ", installed_description[1L, "Version"],
  " and its locally installed ProtVisDatabase resources were installed successfully in ", installed_path,
  ". The lazy-load database passed a clean-process verification. Restart R ",
  "before loading the package."
)
, manifest$path, ignore.case = TRUE)]",
    "if (length(xlsx)) stopifnot(all(vapply(xlsx, function(path) tryCatch({ z <- suppressWarnings(utils::unzip(path, list = TRUE)); is.data.frame(z) && nrow(z) > 0L }, error = function(e) FALSE), logical(1L))))",
    "library(ProtVis, lib.loc = protvis_library)",
    "launcher <- getExportedValue('ProtVis', 'run_ProtVis')",
    "stopifnot(is.function(launcher), is.call(body(launcher)) || is.expression(body(launcher)))"
  ),
  probe_file,
  useBytes = TRUE
)
rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
probe_output <- system2(
  rscript,
  c("--vanilla", shQuote(probe_file)),
  stdout = TRUE,
  stderr = TRUE
)
probe_status <- attr(probe_output, "status")
if (is.null(probe_status)) probe_status <- 0L
if (!identical(as.integer(probe_status), 0L)) {
  stop(
    "ProtVis was installed but failed the clean-process load test:\n",
    paste(probe_output, collapse = "\n"),
    call. = FALSE
  )
}

installed_description <- read.dcf(file.path(installed_path, "DESCRIPTION"))
message(
  "ProtVis ", installed_description[1L, "Version"],
  " and its locally installed ProtVisDatabase resources were installed successfully in ", installed_path,
  ". The lazy-load database passed a clean-process verification. Restart R ",
  "before loading the package."
)
