# Shared sample-group helpers for plots and heatmaps ------------------------

.protvis_triplicate_group <- function(sample_ids) {
  sample_ids <- base::as.character(sample_ids)
  base::sub("_[0-9]+$", "", sample_ids)
}

.protvis_normalise_tissue_labels <- function(values) {
  values <- base::as.character(values)
  lower <- base::tolower(values)
  values[base::grepl("root|below[ ._-]*ground|underground", lower)] <- "Root"
  values[base::grepl(
    "leaf|shoot|stem|above[ ._-]*ground|aerial", lower
  )] <- "Shoot"
  values
}

.protvis_sample_metadata_values <- function(sample_info, sample_ids, column) {
  sample_ids <- base::as.character(sample_ids)
  if (base::is.null(sample_info) || !base::is.data.frame(sample_info) ||
      !"sample_id" %in% base::names(sample_info) ||
      !column %in% base::names(sample_info)) {
    return(base::rep(NA_character_, base::length(sample_ids)))
  }
  info <- base::as.data.frame(
    sample_info, stringsAsFactors = FALSE, check.names = FALSE
  )
  index <- base::match(sample_ids, base::as.character(info$sample_id))
  if ("maxquant_id" %in% base::names(info)) {
    fallback <- base::match(sample_ids, base::as.character(info$maxquant_id))
    index[base::is.na(index)] <- fallback[base::is.na(index)]
  }
  base::as.character(info[[column]])[index]
}

.protvis_sample_group_values <- function(sample_info, sample_ids,
                                         mode = "triplicate") {
  sample_ids <- base::as.character(sample_ids)
  mode <- base::tolower(base::as.character(mode %||% "triplicate"))
  triplicate <- .protvis_triplicate_group(sample_ids)

  if (mode %in% c("triplicate", "auto")) {
    return(triplicate)
  }

  if (mode %in% c("experimental_group", "group")) {
    values <- .protvis_sample_metadata_values(sample_info, sample_ids, "group")
    missing <- base::is.na(values) | !base::nzchar(base::trimws(values)) |
      values == "Unassigned"
    values[missing] <- triplicate[missing]
    return(values)
  }

  if (mode == "species") {
    values <- .protvis_sample_metadata_values(sample_info, sample_ids, "species")
    fallback <- base::ifelse(
      base::grepl("^B73_", sample_ids, ignore.case = TRUE),
      "Zea mays ssp. mays",
      base::ifelse(
        base::grepl("^Y12_", sample_ids, ignore.case = TRUE),
        "Zea mays ssp. mexicana",
        base::sub("_.*$", "", sample_ids)
      )
    )
    encoded <- base::grepl("^(B73|Y12)_", sample_ids, ignore.case = TRUE)
    values[encoded] <- fallback[encoded]
    missing <- base::is.na(values) | !base::nzchar(base::trimws(values)) |
      values %in% c("Unassigned", "All samples")
    values[missing] <- fallback[missing]
    return(values)
  }

  if (mode %in% c("tissue", "tissue2")) {
    values <- .protvis_sample_metadata_values(sample_info, sample_ids, mode)
    if (base::all(base::is.na(values)) && mode == "tissue2") {
      values <- .protvis_sample_metadata_values(sample_info, sample_ids, "tissue")
    }
    values <- .protvis_normalise_tissue_labels(values)
    fallback <- base::ifelse(
      base::grepl("root", sample_ids, ignore.case = TRUE),
      "Root",
      base::ifelse(
        base::grepl("leaf|shoot|stem", sample_ids, ignore.case = TRUE),
        "Shoot",
        triplicate
      )
    )
    encoded <- base::grepl(
      "root|leaf|shoot|stem", sample_ids, ignore.case = TRUE
    )
    values[encoded] <- fallback[encoded]
    missing <- base::is.na(values) | !base::nzchar(base::trimws(values)) |
      values %in% c("Unassigned", "All samples")
    values[missing] <- fallback[missing]
    return(values)
  }

  values <- .protvis_sample_metadata_values(sample_info, sample_ids, mode)
  missing <- base::is.na(values) | !base::nzchar(base::trimws(values)) |
    values == "Unassigned"
  values[missing] <- triplicate[missing]
  values
}

.protvis_group_palette <- function(groups) {
  groups <- base::unique(base::as.character(groups))
  groups <- groups[!base::is.na(groups) & base::nzchar(groups)]
  if (!base::length(groups)) return(base::setNames(character(), character()))

  palette <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
    "#76B7B2", "#EDC948", "#B07AA1", "#FF9DA7",
    "#9C755F", "#BAB0AC", "#2F6B99", "#8CD17D",
    "#B6992D", "#499894", "#D37295", "#79706E"
  )
  if (base::length(groups) > base::length(palette)) {
    palette <- grDevices::hcl.colors(base::length(groups), "Dynamic")
  } else {
    palette <- palette[base::seq_along(groups)]
  }
  stats::setNames(palette, groups)
}

.protvis_shape_palette <- function(groups) {
  groups <- base::unique(base::as.character(groups))
  groups <- groups[!base::is.na(groups) & base::nzchar(groups)]
  shapes <- c(16, 17, 15, 3, 8, 4, 7, 9, 10, 12, 13, 14, 0, 1, 2, 5, 6, 11)
  if (base::length(groups) > base::length(shapes)) {
    shapes <- base::rep(shapes, length.out = base::length(groups))
  } else {
    shapes <- shapes[base::seq_along(groups)]
  }
  stats::setNames(shapes, groups)
}

.protvis_group_label <- function(mode) {
  switch(
    base::tolower(base::as.character(mode %||% "")),
    triplicate = "Experimental group",
    group = "Experimental group",
    experimental_group = "Experimental group",
    species = "Species",
    tissue = "Tissue",
    tissue2 = "Tissue",
    batch = "Batch",
    condition = "Condition",
    "Group"
  )
}
