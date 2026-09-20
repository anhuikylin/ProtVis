test_that("directory roots contain only accessible unique paths", {
  missing <- file.path(tempdir(), "protvis-directory-that-does-not-exist")
  roots <- ProtVis:::.protvis_sanitize_directory_roots(c(
    Home = tempdir(),
    Duplicate = tempdir(),
    Missing = missing
  ))

  expect_named(roots)
  expect_length(roots, 1L)
  expect_true(all(dir.exists(roots)))
  expect_identical(names(roots), "Home")
})

test_that("the chooser always receives at least one valid root", {
  roots <- ProtVis:::.protvis_directory_roots(
    home = tempdir(),
    r_home = tempdir(),
    system_roots = character()
  )

  expect_length(roots, 1L)
  expect_true(all(dir.exists(roots)))
  expect_identical(names(roots), "Home")
})

test_that("platform roots satisfy the shinyFiles roots contract", {
  roots <- getVolumes_win()

  expect_type(roots, "character")
  expect_named(roots)
  expect_gt(length(roots), 0L)
  expect_true(all(dir.exists(roots)))
  expect_false(anyDuplicated(unname(roots)) > 0L)
})


test_that("directory parser ignores transient shinyDirButton clicks", {
  roots <- c(Home = tempdir())
  called <- FALSE
  parser <- function(roots, selection) {
    called <<- TRUE
    stop("parser should not run for a button click")
  }

  result <- ProtVis:::.protvis_parse_directory_selection(
    roots, 1L, parser = parser
  )

  expect_null(result)
  expect_false(called)
})

test_that("directory parser accepts a confirmed structured selection", {
  roots <- c(Home = tempdir())
  parser <- function(roots, selection) tempdir()
  selection <- list(root = "Home", path = "selected")

  result <- ProtVis:::.protvis_parse_directory_selection(
    roots, selection, parser = parser
  )

  expect_identical(
    result,
    normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  )
})

test_that("directory parser ignores incomplete structured payloads", {
  roots <- c(Home = tempdir())

  expect_null(
    ProtVis:::.protvis_parse_directory_selection(
      roots,
      list(root = NULL, path = NULL),
      parser = function(...) stop("should not parse")
    )
  )
})
