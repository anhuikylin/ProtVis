test_that("Protein Workbench domain tracks use labels and adaptive height", {
  domains <- data.frame(
    accession = c("IPR0001", "IPR0002", "IPR0001"),
    name = c(
      "ATPase catalytic domain",
      "A very long transmembrane-associated domain annotation used for display testing",
      "ATPase catalytic domain"
    ),
    start = c(12L, 155L, 420L),
    end = c(120L, 255L, 590L),
    stringsAsFactors = FALSE
  )

  plot_data <- ProtVis:::.protvis_pw_domain_plot_data(domains)
  expect_equal(nrow(plot_data), 3)
  expect_equal(plot_data$label[[1]], "ATPase catalytic domain")
  expect_true(all(nchar(plot_data$label_display) <= 52L))
  expect_equal(length(unique(plot_data$track)), 3)
  expect_gt(
    ProtVis:::.protvis_pw_domain_plot_height(nrow(plot_data)),
    ProtVis:::.protvis_pw_domain_plot_height(0L)
  )
  expect_equal(ProtVis:::.protvis_pw_domain_plot_height(100L), 1600L)
})
