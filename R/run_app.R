#' Run the Shiny Application
#' @description Launches the ProtVis Shiny application for proteomics data analysis.
#' @param onStart A function that will be called before the app is actually run.
#' @param options Named list of values passed to `shiny::shinyOptions`.
#' @param enableBookmarking Can be "url", "server", or "disable".
#' @param uiPattern A regular expression used to determine which requests should be handled by the UI.
#' @param ... Reserved for backward compatibility; currently ignored.
#' @return An object that represents the app.
#' @import shiny
#' @importFrom utils modifyList
#' @name run_ProtVis
#' @export

run_ProtVis <- function(
    onStart = NULL,
    options = list(),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...
) {

  # Validate the local companion-data installation before building UI.
  ProtVisDatabase::protvis_database_path(
    "extdata",
    "manifest.tsv",
    must_work = TRUE
  )

  # Set a generous default for large proteomics/background workbooks.
  # modifyList ensures that user-defined options are preserved
  default_options <- list(shiny.maxRequestSize = 2 * 1024^3)
  combined_options <- utils::modifyList(default_options, options)

  app <- shiny::shinyApp(
    ui = app_ui,
    server = app_server,
    onStart = onStart,
    options = combined_options,
    enableBookmarking = enableBookmarking,
    uiPattern = uiPattern
  )

  invisible(list(...))
  app
}
