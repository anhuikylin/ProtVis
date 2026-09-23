#' Run the Shiny Application
#' @description Launches the ProtVis Shiny application for proteomics data analysis.
#' @param onStart A function that will be called before the app is actually run.
#' @param options Named list of values passed to `shiny::shinyOptions`.
#' @param enableBookmarking Can be "url", "server", or "disable".
#' @param uiPattern A regular expression used to determine which requests should be handled by the UI.
#' @param progress Logical; show a compact console progress bar while ProtVis initializes.
#' @param ... Arguments to pass to `golem_opts`. See `?golem::get_golem_options` for more details.
#' @return An object that represents the app.
#' @import shiny
#' @importFrom golem with_golem_options
#' @importFrom utils modifyList
#' @name run_ProtVis
#' @export

run_ProtVis <- function(
    onStart = NULL,
    options = list(),
    enableBookmarking = NULL,
    uiPattern = "/",
    progress = getOption("protvis.startup.progress", TRUE),
    ...
) {

  show_progress <- isTRUE(progress) && interactive()
  startup_bar <- NULL

  if (show_progress) {
    base::cat("ProtVis startup ")
    startup_bar <- utils::txtProgressBar(
      min = 0,
      max = 3,
      initial = 0,
      style = 3,
      width = 28
    )
    on.exit(
      if (!is.null(startup_bar)) {
        base::close(startup_bar)
      },
      add = TRUE
    )
  }

  advance_startup <- function(step) {
    if (!is.null(startup_bar)) {
      utils::setTxtProgressBar(startup_bar, step)
    }
    invisible(step)
  }

  # Step 1: validate the local companion-data installation before building UI.
  ProtVisDatabase::protvis_database_path(
    "extdata",
    "manifest.tsv",
    must_work = TRUE
  )
  advance_startup(1)

  # Set a generous default for large proteomics/background workbooks.
  # modifyList ensures that user-defined options are preserved
  default_options <- list(shiny.maxRequestSize = 2 * 1024^3)
  combined_options <- utils::modifyList(default_options, options)
  advance_startup(2)

  app <- shiny::shinyApp(
    ui = app_ui,
    server = app_server,
    onStart = onStart,
    options = combined_options,
    enableBookmarking = enableBookmarking,
    uiPattern = uiPattern
  )

  advance_startup(3)
  if (!is.null(startup_bar)) {
    base::close(startup_bar)
    startup_bar <- NULL
  }

  golem::with_golem_options(
    app = app,
    golem_opts = list(...)
  )
}
