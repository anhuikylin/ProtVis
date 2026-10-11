test_that("replicate filtering preserves observations and only optionally imputes", {
  raw <- data.frame(ID = c("partial", "single", "complete", "absent"),
                    A_1 = c(100, 100, 10, 0),
                    A_2 = c(120, 0, 20, NA),
                    A_3 = c(0, 0, 30, 0),
                    A_4 = c(140, 0, 40, 0))
  filtered <- .protvis_filter_replicates(raw)
  expect_identical(filtered$ID, c("partial", "complete"))
  expect_equal(as.numeric(filtered[1, -1]), c(100, 120, NA, 140))
  filled <- .protvis_filter_replicates(raw, TRUE)
  expect_equal(as.numeric(filled[1, -1]), c(100, 120, 120, 140))
  expect_equal(filtered[2, ], filled[2, ])
})

test_that("replicate filtering masks failed groups without removing other groups", {
  raw <- data.frame(ID = "p", A_1 = 5, A_2 = 0, B_1 = 10, B_2 = 20)
  filtered <- .protvis_filter_replicates(raw, TRUE)
  expect_true(all(is.na(filtered[1, c("A_1", "A_2")])))
  expect_equal(as.numeric(filtered[1, c("B_1", "B_2")]), c(10, 20))
  expect_error(.protvis_filter_replicates(data.frame(ID = "p", A = 1, B = 2)),
               "ending")
  expect_error(.protvis_filter_replicates(data.frame(ID = "p", A_1 = 1, B_1 = 2)),
               "at least two")
})
