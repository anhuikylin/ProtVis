# Subcellular localization analysis for plant proteins
#
# This module upgrades the former single-service Plant-mPLoc page into a
# fault-tolerant consensus workflow. Plant-mPLoc remains available as one
# evidence source; WoLF PSORT and BUSCA are queried independently and failures
# are reported per provider instead of failing the whole analysis.

.subcellular_provider_urls <- c(
  "Plant-mPLoc" = "http://www.csbio.sjtu.edu.cn/bioinf/plant-multi/",
  "WoLF PSORT" = "https://wolfpsort.hgc.jp/",
  "BUSCA" = "https://busca.biocomp.unibo.it/"
)

.subcellular_location_aliases <- list(
  "Cell membrane" = c("cell membrane", "plasma membrane", "plasma_membrane", "plas"),
  "Cell wall" = c("cell wall", "cell_wall"),
  "Chloroplast" = c("chloroplast", "chloroplastic", "chlo"),
  "Cytoplasm" = c("cytoplasm", "cytoplasmic", "cytosol", "cyto"),
  "Endoplasmic reticulum" = c("endoplasmic reticulum", "endoplasmic_reticulum", "e.r.", "er"),
  "Extracellular" = c("extracellular", "secreted", "extr"),
  "Golgi apparatus" = c("golgi apparatus", "golgi", "golg"),
  "Mitochondrion" = c("mitochondrion", "mitochondrial", "mito"),
  "Nucleus" = c("nucleus", "nuclear", "nucl"),
  "Peroxisome" = c("peroxisome", "peroxisomal", "pero"),
  "Plastid" = c("plastid", "plastidial"),
  "Lysosome/Vacuole" = c("lysosome/vacuole", "lysosome", "lysosomal"),
  "Vacuole" = c("vacuole", "vacuolar", "vacu")
)

.subcellular_normalize_locations <- function(x) {
  if (is.null(x) || !length(x)) return(character())
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x) & !is.na(x)]
  if (!length(x)) return(character())
  out <- character()
  for (value in x) {
    low <- tolower(value)
    matched <- names(.subcellular_location_aliases)[vapply(
      .subcellular_location_aliases,
      function(aliases) any(vapply(aliases, function(a) {
        a <- tolower(a)
        identical(low, a) || (nchar(a) >= 4L && grepl(a, low, fixed = TRUE))
      }, logical(1))),
      logical(1)
    )]
    if (length(matched)) out <- c(out, matched[[1L]])
  }
  unique(out)
}

.subcellular_provider_result <- function(
    source, status = "unavailable", prediction = character(),
    raw_result = "", result_url = NA_character_, error = NULL,
    details = NULL, score_table = NULL) {
  prediction <- .subcellular_normalize_locations(prediction)
  structure(
    list(
      source = source,
      status = status,
      prediction = prediction,
      raw_result = raw_result %||% "",
      result_url = result_url %||% NA_character_,
      error = error,
      details = details,
      score_table = score_table
    ),
    class = "ProtVis_localization_provider"
  )
}

.subcellular_form_field <- function(forms, pattern, types = NULL) {
  candidates <- list()
  for (fi in seq_along(forms)) {
    fields <- forms[[fi]]$fields
    if (!length(fields)) next
    nms <- names(fields)
    for (i in seq_along(fields)) {
      field <- fields[[i]]
      type <- tolower(as.character(field$type %||% ""))
      if (!is.null(types) && !type %in% types) next
      label <- paste(
        nms[[i]] %||% "",
        field$name %||% "",
        field$id %||% "",
        field$placeholder %||% "",
        field$value %||% "",
        collapse = " "
      )
      score <- if (grepl(pattern, label, ignore.case = TRUE, perl = TRUE)) 100L else 0L
      if (type == "textarea") score <- score + 30L
      if (score > 0L) {
        candidates[[length(candidates) + 1L]] <- list(
          form = fi, field = nms[[i]], index = i, score = score
        )
      }
    }
  }
  if (!length(candidates)) return(NULL)
  candidates[[which.max(vapply(candidates, `[[`, numeric(1), "score"))]]
}

.subcellular_form_submit <- function(form) {
  fields <- form$fields
  if (!length(fields)) return(NULL)
  nms <- names(fields)
  idx <- which(vapply(fields, function(x) {
    tolower(as.character(x$type %||% "")) %in% c("submit", "button", "image")
  }, logical(1)))
  if (!length(idx)) return(NULL)
  labels <- vapply(idx, function(i) {
    paste(nms[[i]] %||% "", fields[[i]]$value %||% "", fields[[i]]$name %||% "")
  }, character(1))
  preferred <- which(grepl("start|submit|predict|run", labels, ignore.case = TRUE))
  name <- nms[[idx[if (length(preferred)) preferred[[1L]] else 1L]]]
  if (is.null(name) || !nzchar(name)) NULL else name
}

.subcellular_select_plant_value <- function(field) {
  values <- unique(c(
    as.character(field$values %||% character()),
    as.character(field$value %||% character())
  ))
  values <- values[nzchar(values) & !is.na(values)]
  hit <- values[grepl("plant|viridi", values, ignore.case = TRUE)]
  if (length(hit)) return(hit[[1L]])
  # Common server-side encodings. Used only when field choices are not exposed.
  "plant"
}

.subcellular_session_text <- function(session) {
  body <- tryCatch(rvest::html_element(session, "body"), error = function(e) NULL)
  if (is.null(body) || !length(body)) return("")
  tryCatch(rvest::html_text2(body), error = function(e) "")
}

.subcellular_parse_refresh <- function(refresh) {
  if (is.null(refresh) || !length(refresh) || is.na(refresh[[1L]])) return("")
  refresh <- as.character(refresh[[1L]])
  if (!nzchar(refresh) || !grepl("url[[:space:]]*=", refresh, ignore.case = TRUE)) {
    return("")
  }
  target <- sub(".*url[[:space:]]*=[[:space:]]*", "", refresh, ignore.case = TRUE)
  target <- gsub("^[\"']|[\"']$", "", trimws(target))
  if (!nzchar(target) || grepl("^(javascript|data):", target, ignore.case = TRUE)) "" else target
}

.subcellular_refresh_target <- function(session) {
  headers <- tryCatch(httr::headers(session$response), error = function(e) list())
  header <- headers[which(tolower(names(headers)) == "refresh")]
  refresh <- if (length(header)) header[[1L]] else ""
  if (is.null(refresh) || !length(refresh) || is.na(refresh[[1L]]) ||
      !nzchar(as.character(refresh[[1L]]))) {
    meta <- tryCatch(
      rvest::html_element(session, "meta[http-equiv='refresh']"),
      error = function(e) NULL
    )
    if (!is.null(meta) && length(meta)) {
      refresh <- rvest::html_attr(meta, "content") %||% ""
    }
  }
  .subcellular_parse_refresh(refresh)
}

.subcellular_follow_result <- function(session, timeout = 90, provider = "") {
  started <- Sys.time()
  current <- session
  repeat {
    txt <- .subcellular_session_text(current)
    if (identical(provider, "WoLF PSORT")) {
      if (!is.null(.subcellular_parse_wolfpsort(txt)$scores)) return(current)
    } else if (identical(provider, "BUSCA")) {
      if (grepl("/showresult/", current$url, fixed = TRUE) &&
          !grepl("queued|pending|running|refresh every minute", txt, ignore.case = TRUE)) {
        return(current)
      }
    } else if (grepl("final results|predicted localization", txt, ignore.case = TRUE) &&
               !grepl("processing|please wait|queued|running", txt, ignore.case = TRUE)) {
      return(current)
    }
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (!is.finite(elapsed) || elapsed >= timeout) return(current)

    target <- .subcellular_refresh_target(current)
    if (nzchar(target) && grepl("^https?://", target, ignore.case = TRUE)) {
      host <- sub("^https?://([^/]+).*$", "\\1", target, ignore.case = TRUE)
      allowed <- c("wolfpsort.hgc.jp", "busca.biocomp.unibo.it")
      if (!tolower(host) %in% allowed) target <- ""
    }
    if (!nzchar(target)) {
      Sys.sleep(if (identical(provider, "BUSCA")) 10 else 2)
      target <- current$url
    }
    next_page <- tryCatch(rvest::session_jump_to(current, target), error = function(e) NULL)
    if (!is.null(next_page)) current <- next_page
  }
}

.subcellular_submit_generic <- function(
    url, fasta, provider, timeout = 90, ignore_ssl = FALSE,
    organism_pattern = "organism|tax|kingdom|type") {
  if (!requireNamespace("rvest", quietly = TRUE)) stop("Package 'rvest' is required.")
  if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' is required.")

  opts <- list(
    httr::timeout(timeout),
    httr::user_agent("ProtVis subcellular-localization client")
  )
  if (isTRUE(ignore_ssl)) {
    opts <- c(opts, list(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)))
  }
  ses <- do.call(rvest::session, c(list(url = url), opts))
  forms <- rvest::html_form(ses)
  seq_field <- .subcellular_form_field(
    forms, "seq|fasta|protein|query|input",
    c("textarea", "text", "search")
  )
  if (is.null(seq_field)) stop(provider, " FASTA input field was not found.")
  form <- forms[[seq_field$form]]
  values <- stats::setNames(list(fasta), seq_field$field)

  org <- .subcellular_form_field(
    list(form), organism_pattern,
    c("select", "radio", "text", "hidden")
  )
  if (!is.null(org)) {
    values[[org$field]] <- .subcellular_select_plant_value(form$fields[[org$index]])
  }
  form <- do.call(rvest::html_form_set, c(list(form), values))
  submit <- .subcellular_form_submit(form)
  result <- if (is.null(submit)) {
    rvest::session_submit(ses, form)
  } else {
    rvest::session_submit(ses, form, submit = submit)
  }
  .subcellular_follow_result(result, timeout = timeout, provider = provider)
}

.subcellular_parse_wolfpsort <- function(text) {
  if (is.null(text) || !nzchar(text)) return(list(prediction = character(), scores = NULL))
  clean <- gsub("[\r\n]+", " ", text)
  # WoLF PSORT commonly reports: query chlo: 12, cyto: 8, nucl: 4 ...
  keys <- c(
    chlo = "Chloroplast", cyto = "Cytoplasm", nucl = "Nucleus",
    mito = "Mitochondrion", plas = "Cell membrane", extr = "Extracellular",
    pero = "Peroxisome", vacu = "Vacuole", golg = "Golgi apparatus",
    E.R. = "Endoplasmic reticulum"
  )
  rows <- list()
  for (key in names(keys)) {
    escaped <- gsub("\\.", "\\\\.", key)
    m <- regexec(
      paste0("(?i)(?:^|[^[:alnum:]_])", escaped, "\\s*:\\s*([0-9]+(?:\\.[0-9]+)?)"),
      clean, perl = TRUE
    )
    hit <- regmatches(clean, m)[[1L]]
    if (length(hit) >= 2L) {
      rows[[length(rows) + 1L]] <- data.frame(
        location = unname(keys[[key]]),
        score = as.numeric(hit[[2L]]),
        stringsAsFactors = FALSE
      )
    }
  }
  scores <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(scores) && nrow(scores)) {
    pred <- scores$location[scores$score == max(scores$score, na.rm = TRUE)]
    return(list(prediction = unique(pred), scores = scores[order(-scores$score), , drop = FALSE]))
  }

  final <- sub(".*(?:prediction for|predicted(?: location)?(?: is|:))\\s*", "", clean, ignore.case = TRUE)
  pred <- .subcellular_normalize_locations(final)
  list(prediction = pred, scores = scores)
}

#' Predict plant protein localization with WoLF PSORT
#'
#' @param sequence Protein sequence.
#' @param id Protein identifier.
#' @param timeout Request timeout in seconds.
#' @return Provider result used by the consensus localization module.
#' @export
predict_wolfpsort <- function(sequence, id = "query_protein", timeout = 90) {
  sequence <- .plant_mploc_clean_sequence(sequence)
  if (nchar(sequence) < 30L) {
    return(.subcellular_provider_result(
      "WoLF PSORT", "unavailable",
      error = "WoLF PSORT requires at least 30 amino acids."
    ))
  }
  fasta <- .plant_mploc_fasta(sequence, id)
  tryCatch({
    ses <- .subcellular_submit_generic(
      .subcellular_provider_urls[["WoLF PSORT"]],
      fasta, "WoLF PSORT", timeout = timeout
    )
    txt <- .subcellular_session_text(ses)
    parsed <- .subcellular_parse_wolfpsort(txt)
    if (!length(parsed$prediction)) {
      return(.subcellular_provider_result(
        "WoLF PSORT", "unavailable", raw_result = txt, result_url = ses$url,
        error = "WoLF PSORT did not return a scored result before the timeout. Check the raw response and result URL."
      ))
    }
    .subcellular_provider_result(
      "WoLF PSORT", "success", parsed$prediction,
      raw_result = txt, result_url = ses$url, score_table = parsed$scores
    )
  }, error = function(e) {
    .subcellular_provider_result("WoLF PSORT", "unavailable", error = conditionMessage(e))
  })
}

.subcellular_parse_busca <- function(text) {
  if (is.null(text) || !nzchar(text)) return(character())
  lines <- trimws(unlist(strsplit(text, "[\r\n]+")))
  terms <- lines[grepl("^C:[[:space:]]*|^cellular component:", lines, ignore.case = TRUE)]
  .subcellular_normalize_locations(terms)
}

.subcellular_parse_busca_json <- function(value) {
  leaves <- unlist(value, recursive = TRUE, use.names = FALSE)
  terms <- as.character(leaves[!is.na(leaves)])
  terms <- terms[grepl("^C:[[:space:]]*", terms, ignore.case = TRUE)]
  .subcellular_normalize_locations(terms)
}

.subcellular_busca_json <- function(result_url, timeout = 30) {
  if (!grepl("^https://busca\\.biocomp\\.unibo\\.it/[A-Za-z0-9-]+/showresult/?$", result_url)) {
    stop("BUSCA did not return a recognized job URL.", call. = FALSE)
  }
  url <- sub("showresult/?$", "getjson/", result_url)
  response <- httr::GET(url, httr::timeout(min(as.numeric(timeout), 30)))
  httr::stop_for_status(response)
  raw <- httr::content(response, as = "text", encoding = "UTF-8")
  list(prediction = .subcellular_parse_busca_json(
    jsonlite::fromJSON(raw, simplifyVector = FALSE)
  ), raw = raw)
}

.subcellular_parse_busca_page <- function(page) {
  headers <- tryCatch(
    rvest::html_text2(rvest::html_elements(page, "#resultdata thead th")),
    error = function(e) character()
  )
  term_column <- which(tolower(trimws(headers)) == "go-term")
  rows <- tryCatch(rvest::html_elements(page, "#resultdata tbody tr"), error = function(e) list())
  if (!length(term_column) || !length(rows)) return(character())
  terms <- vapply(rows, function(row) {
    cells <- rvest::html_elements(row, "td")
    if (length(cells) < term_column[[1L]]) return("")
    rvest::html_text2(cells[[term_column[[1L]]]])
  }, character(1))
  .subcellular_normalize_locations(terms)
}

#' Predict plant protein localization with BUSCA
#'
#' @param sequence Protein sequence.
#' @param id Protein identifier.
#' @param timeout Request timeout in seconds.
#' @return Provider result used by the consensus localization module.
#' @export
predict_busca <- function(sequence, id = "query_protein", timeout = 120) {
  sequence <- .plant_mploc_clean_sequence(sequence)
  fasta <- .plant_mploc_fasta(sequence, id)
  tryCatch({
    ses <- .subcellular_submit_generic(
      .subcellular_provider_urls[["BUSCA"]],
      fasta, "BUSCA", timeout = timeout
    )
    txt <- .subcellular_session_text(ses)
    complete <- grepl("/showresult/", ses$url, fixed = TRUE) &&
      !grepl("queued|pending|running|refresh every minute", txt, ignore.case = TRUE)
    pred <- if (complete) .subcellular_parse_busca_page(ses) else character()
    raw <- txt
    if (complete && !length(pred)) {
      job <- tryCatch(.subcellular_busca_json(ses$url, timeout), error = function(e) NULL)
      if (!is.null(job)) {
        pred <- job$prediction
        raw <- job$raw
      }
    }
    if (!length(pred)) {
      queued <- grepl("queued|pending|running|refresh every minute", txt, ignore.case = TRUE)
      return(.subcellular_provider_result(
        "BUSCA", if (queued) "queued" else "unavailable",
        raw_result = txt, result_url = ses$url,
        error = if (queued) "BUSCA is still processing this job; use its result URL to check later."
                else "BUSCA did not return a recognizable completed result."
      ))
    }
    .subcellular_provider_result(
      "BUSCA", "success", pred, raw_result = raw, result_url = ses$url
    )
  }, error = function(e) {
    .subcellular_provider_result("BUSCA", "unavailable", error = conditionMessage(e))
  })
}

.subcellular_run_plant_mploc <- function(sequence, id, timeout, ignore_ssl) {
  tryCatch({
    x <- predict_plant_mploc(
      sequence = sequence, id = id, timeout = timeout,
      ignore_ssl = ignore_ssl, verbose = FALSE
    )
    .subcellular_provider_result(
      "Plant-mPLoc", "success", x$prediction,
      raw_result = x$raw_result, result_url = x$result_url
    )
  }, error = function(e) {
    .subcellular_provider_result("Plant-mPLoc", "unavailable", error = conditionMessage(e))
  })
}

.subcellular_consensus <- function(provider_results) {
  ok <- provider_results[vapply(provider_results, function(x) identical(x$status, "success"), logical(1))]
  if (!length(ok)) {
    return(list(
      prediction = character(), confidence = "Unavailable",
      agreement = 0L, successful = 0L, requested = length(provider_results),
      votes = data.frame()
    ))
  }
  locs <- unlist(lapply(ok, function(x) unique(x$prediction)), use.names = FALSE)
  tab <- sort(table(locs), decreasing = TRUE)
  max_vote <- max(tab)
  prediction <- names(tab)[tab == max_vote]
  successful <- length(ok)
  confidence <- if (successful >= 3L && max_vote >= 3L) {
    "High"
  } else if (max_vote >= 2L) {
    "Moderate"
  } else {
    "Low"
  }
  votes <- data.frame(
    location = names(tab),
    votes = as.integer(tab),
    agreement = sprintf("%d/%d", as.integer(tab), successful),
    stringsAsFactors = FALSE
  )
  list(
    prediction = prediction,
    confidence = confidence,
    agreement = as.integer(max_vote),
    successful = successful,
    requested = length(provider_results),
    votes = votes
  )
}

#' Run consensus plant subcellular-localization prediction
#'
#' Providers are independent: a failed web service is retained as unavailable
#' evidence and does not invalidate successful results from other providers.
#'
#' @param sequence Protein sequence.
#' @param id Protein identifier.
#' @param providers Character vector of providers.
#' @param timeout Timeout per provider.
#' @param ignore_ssl Allow the legacy Plant-mPLoc certificate.
#' @return A consensus localization result.
#' @export
predict_subcellular_localization <- function(
    sequence,
    id = "query_protein",
    providers = c("Plant-mPLoc", "WoLF PSORT", "BUSCA"),
    timeout = 120,
    ignore_ssl = TRUE) {
  sequence <- .plant_mploc_clean_sequence(sequence)
  id <- trimws(as.character(id %||% "query_protein")[[1L]])
  if (!nzchar(id)) id <- "query_protein"
  providers <- intersect(unique(providers), names(.subcellular_provider_urls))
  if (!length(providers)) stop("Select at least one prediction source.", call. = FALSE)

  results <- lapply(providers, function(provider) {
    switch(
      provider,
      "Plant-mPLoc" = .subcellular_run_plant_mploc(sequence, id, timeout, ignore_ssl),
      "WoLF PSORT" = predict_wolfpsort(sequence, id, timeout),
      "BUSCA" = predict_busca(sequence, id, timeout)
    )
  })
  names(results) <- providers
  consensus <- .subcellular_consensus(results)
  structure(
    list(
      protein_id = id,
      sequence = sequence,
      length = nchar(sequence),
      providers = results,
      consensus = consensus,
      submitted_at = Sys.time()
    ),
    class = "ProtVis_subcellular_localization"
  )
}

.subcellular_deeploc_locations <- c(
  "Cytoplasm", "Nucleus", "Extracellular", "Cell membrane",
  "Mitochondrion", "Plastid", "Endoplasmic reticulum",
  "Lysosome/Vacuole", "Golgi apparatus", "Peroxisome"
)
.subcellular_deeploc_membranes <- c("Peripheral", "Transmembrane", "Lipid anchor", "Soluble")

.subcellular_deeploc_labels <- function(value, choices, field, required = FALSE) {
  value <- paste(value %||% "", collapse = ";")
  parts <- trimws(unlist(strsplit(gsub("\n", ";", value, fixed = TRUE), "[;,]+")))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    if (required) stop(paste("Enter DeepLoc", field, "labels."), call. = FALSE)
    return(character())
  }
  matched <- choices[match(tolower(parts), tolower(choices))]
  if (anyNA(matched)) {
    stop(paste("Unknown DeepLoc", field, "label:", paste(parts[is.na(matched)], collapse = "; ")), call. = FALSE)
  }
  unique(matched)
}

.subcellular_deeploc_probabilities <- function(value) {
  value <- trimws(value %||% "")
  if (!nzchar(value)) return(data.frame(type = character(), label = character(), probability = numeric()))
  lines <- strsplit(value, "\n", fixed = TRUE)[[1L]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  rows <- lapply(lines, function(line) {
    parts <- trimws(strsplit(line, "[,\t;]+")[[1L]])
    if (length(parts) != 2L) stop("Probability rows must be label,probability.", call. = FALSE)
    labels <- c(.subcellular_deeploc_locations, .subcellular_deeploc_membranes)
    requested <- if (tolower(parts[[1L]]) == "chloroplast") "Plastid" else parts[[1L]]
    idx <- match(tolower(requested), tolower(labels))
    if (is.na(idx)) stop(paste("Unknown probability label:", parts[[1L]]), call. = FALSE)
    score <- suppressWarnings(as.numeric(parts[[2L]]))
    if (!is.finite(score) || score < 0 || score > 1) {
      stop(paste("Probability must be from 0 to 1 for", parts[[1L]]), call. = FALSE)
    }
    data.frame(type = if (idx <= length(.subcellular_deeploc_locations)) "Localization" else "Membrane association",
               label = labels[[idx]], probability = score)
  })
  out <- do.call(rbind, rows)
  if (anyDuplicated(out$label)) stop("Each DeepLoc probability label must appear once.", call. = FALSE)
  expected <- c(.subcellular_deeploc_locations, .subcellular_deeploc_membranes)
  if (!setequal(out$label, expected)) {
    stop("DeepLoc probability table needs all 10 localization and 4 membrane scores.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}

.subcellular_deeploc_importance <- function(file) {
  if (is.null(file)) return(NULL)
  tab <- utils::read.csv(file$datapath, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab)) stop("Sorting importance CSV is empty.", call. = FALSE)
  keys <- tolower(gsub("[^a-z0-9]", "", names(tab)))
  pos <- match(TRUE, keys %in% c("position", "pos", "residueindex", "index", "seqpos"))
  val <- match(TRUE, keys %in% c("importance", "score", "attention", "sortingsignalimportance"))
  if (is.na(pos) || is.na(val)) {
    stop("Sorting importance CSV needs position and importance columns.", call. = FALSE)
  }
  positions <- suppressWarnings(as.integer(tab[[pos]]))
  scores <- suppressWarnings(as.numeric(tab[[val]]))
  if (anyNA(positions) || any(positions < 1L) || anyDuplicated(positions) ||
      any(!is.finite(scores)) || any(scores < 0 | scores > 1)) {
    stop("Sorting importance must have unique positive positions and scores from 0 to 1.", call. = FALSE)
  }
  if (any(positions > 100000L)) stop("Sorting importance position exceeds the supported range.", call. = FALSE)
  data.frame(position = positions,
             residue = if ("residue" %in% keys) as.character(tab[[match("residue", keys)]]) else NA_character_,
             importance = scores, stringsAsFactors = FALSE)
}

.subcellular_import_deeploc <- function(location, membrane = "", signals = "",
                                        probabilities = "", importance = NULL,
                                        sequence = NULL) {
  locations <- .subcellular_deeploc_labels(
    location, c(.subcellular_deeploc_locations, "Chloroplast"), "localization", TRUE
  )
  membranes <- .subcellular_deeploc_labels(membrane, .subcellular_deeploc_membranes, "membrane")
  signals <- trimws(unlist(strsplit(gsub("\n", ";", signals %||% "", fixed = TRUE), "[;,]+")))
  signals <- unique(signals[nzchar(signals)])
  scores <- .subcellular_deeploc_probabilities(probabilities)
  if (nrow(scores) && !all(ifelse(locations == "Chloroplast", "Plastid", locations) %in% scores$label)) {
    stop("Provide a probability for every predicted localization, or leave probabilities empty.", call. = FALSE)
  }
  if (!is.null(importance) && !is.null(sequence)) {
    if (any(importance$position > nchar(sequence))) {
      stop("Sorting importance contains a position beyond the protein sequence.", call. = FALSE)
    }
    valid_residues <- !is.na(importance$residue) & nzchar(importance$residue)
    if (any(valid_residues) &&
        any(toupper(importance$residue[valid_residues]) !=
            substring(sequence, importance$position[valid_residues],
                      importance$position[valid_residues]))) {
      stop("Sorting importance residues do not match the current protein sequence.", call. = FALSE)
    }
  }
  provider <- .subcellular_provider_result(
    "DeepLoc 2.1", "success", locations,
    details = "Entered from the DeepLoc 2.1 web result by the user",
    result_url = "https://services.healthtech.dtu.dk/services/DeepLoc-2.1/",
    score_table = scores
  )
  provider$membrane_types <- membranes
  provider$signals <- signals
  provider$sorting_importance <- importance
  provider
}

.subcellular_deeploc_summary <- function(file, protein_id) {
  if (is.null(file)) return(NULL)
  tab <- utils::read.csv(file$datapath, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab)) stop("DeepLoc summary CSV is empty.", call. = FALSE)
  keys <- tolower(gsub("[^a-z0-9]", "", names(tab)))
  id_col <- match(TRUE, keys %in% c("proteinid", "protein", "id", "name", "sequenceid", "entry"))
  if (nrow(tab) > 1L) {
    if (is.na(id_col)) stop("Summary has multiple proteins but no protein ID column.", call. = FALSE)
    hit <- which(as.character(tab[[id_col]]) == protein_id)
    if (length(hit) != 1L) stop("Select a single DeepLoc summary row matching the current protein ID.", call. = FALSE)
    tab <- tab[hit, , drop = FALSE]
  } else if (!is.na(id_col) && nzchar(as.character(tab[[id_col]][[1L]])) &&
             !identical(as.character(tab[[id_col]][[1L]]), protein_id)) {
    stop("DeepLoc summary protein ID does not match the current protein ID.", call. = FALSE)
  }
  col_value <- function(aliases) {
    idx <- match(TRUE, keys %in% aliases)
    if (is.na(idx)) "" else as.character(tab[[idx]][[1L]])
  }
  locations <- col_value(c("predictedlocalizations", "localizations", "localization",
                           "predictedlocation", "prediction"))
  membranes <- col_value(c("predictedmembraneassociation", "predictedmembranetypes",
                           "membraneassociation", "membranetypes", "membranelabels"))
  signals <- col_value(c("predictedsignals", "sortingsignals", "signals", "signal"))
  labels <- c(.subcellular_deeploc_locations, .subcellular_deeploc_membranes)
  scores <- character()
  for (label in labels) {
    target <- tolower(gsub("[^a-z0-9]", "", label))
    aliases <- if (identical(label, "Plastid")) c("plastid", "chloroplast") else target
    idx <- match(TRUE, keys %in% unlist(lapply(aliases, function(key) {
      c(key, paste0(key, "probability"), paste0("probability", key), paste0(key, "score"))
    })))
    if (!is.na(idx) && !is.na(tab[[idx]][[1L]]) && nzchar(as.character(tab[[idx]][[1L]]))) {
      scores <- c(scores, paste(label, tab[[idx]][[1L]], sep = ","))
    }
  }
  if (!nzchar(locations)) {
    stop("No predicted localizations column found in DeepLoc summary. Enter labels manually.", call. = FALSE)
  }
  list(location = locations, membrane = membranes, signals = signals,
       probabilities = paste(scores, collapse = "\n"))
}

.subcellular_deeploc_example <- function() {
  paste(c(
    "Cytoplasm,0.2085", "Nucleus,0.1770", "Extracellular,0.5127",
    "Cell membrane,0.1617", "Mitochondrion,0.0995", "Plastid,0.0019",
    "Endoplasmic reticulum,0.7058", "Lysosome/Vacuole,0.4825",
    "Golgi apparatus,0.4283", "Peroxisome,0.0013",
    "Peripheral,0.4020", "Transmembrane,0.0870",
    "Lipid anchor,0.1280", "Soluble,0.8140"
  ), collapse = "\n")
}

.subcellular_deeploc_fasta <- function(sequence, id) {
  sequence <- .plant_mploc_clean_sequence(sequence)
  if (nchar(sequence) < 10L) {
    stop("DeepLoc 2.1 requires at least 10 amino acids.", call. = FALSE)
  }
  id <- gsub("[^A-Za-z0-9_.-]+", "_", trimws(as.character(id %||% "")))
  if (!nzchar(id)) id <- "query_protein"
  .plant_mploc_fasta(sequence, id)
}

.subcellular_evidence_table <- function(value) {
  rows <- lapply(value$providers, function(x) {
    data.frame(
      source = x$source,
      status = switch(x$status, success = "Available", queued = "Queued", "Unavailable"),
      prediction = if (length(x$prediction)) paste(x$prediction, collapse = "; ") else "—",
      detail = if (!is.null(x$error) && nzchar(x$error)) x$error else x$details %||% "",
      result_url = x$result_url %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rbind(
    out,
    data.frame(
      source = "Consensus",
      status = if (length(value$consensus$prediction)) value$consensus$confidence else "Unavailable",
      prediction = if (length(value$consensus$prediction)) {
        paste(value$consensus$prediction, collapse = "; ")
      } else "—",
      detail = sprintf(
        "%d of %d selected providers returned usable results",
        value$consensus$successful, value$consensus$requested
      ),
      result_url = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}

#' ProtVis subcellular localization UI
#'
#' @param id Module namespace.
#' @return Shiny UI.
#' @export
subcellular_localization_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360, open = "open", gap = "12px",
      shiny::div(
        shiny::h4("Subcellular localization", style = "font-weight:700; margin-bottom:6px;"),
        shiny::p(
          "Consensus plant-protein localization from independent prediction services.",
          style = "color:#667085; font-size:13px;"
        )
      ),
      bslib::accordion(
        open = c("Protein input", "Prediction sources", "Run"),
        bslib::accordion_panel(
          "Protein input",
          shiny::textInput(ns("protein_id"), "Protein ID", value = "ZmProtein_demo"),
          shiny::textAreaInput(
            ns("sequence"), "Protein sequence",
            value = .plant_mploc_demo_sequence, rows = 12,
            placeholder = "Paste a complete protein sequence using one-letter amino-acid codes"
          ),
          shiny::actionButton(
            ns("load_demo"), "LOAD DEMO SEQUENCE",
            icon = bsicons::bs_icon("stars"),
            class = "btn-outline-primary w-100 mb-2"
          ),
          shiny::uiOutput(ns("sequence_status"))
        ),
        bslib::accordion_panel(
          "Prediction sources",
          shiny::checkboxGroupInput(
            ns("providers"), NULL,
            choices = c("Plant-mPLoc", "WoLF PSORT", "BUSCA"),
            selected = "WoLF PSORT"
          ),
          shiny::div(
            style = "font-size:12px;color:#667085;",
            "WoLF PSORT is the default. Plant-mPLoc and BUSCA are optional web services; BUSCA jobs can remain queued beyond this session."
          )
        ),
        bslib::accordion_panel(
          "Advanced",
          shiny::numericInput(ns("timeout"), "Timeout per source (seconds)", 120, min = 10, max = 600),
          shiny::checkboxInput(ns("ignore_ssl"), "Allow legacy Plant-mPLoc certificate", TRUE)
        ),
        bslib::accordion_panel(
          "Run",
          shiny::div(
            class = "d-grid gap-2",
            shiny::actionButton(
              ns("run"), "PREDICT LOCALIZATION",
              icon = bsicons::bs_icon("play-fill"), class = "btn-primary fw-bold"
            ),
            shiny::actionButton(
              ns("clear"), "CLEAR",
              icon = bsicons::bs_icon("x-circle"), class = "btn-outline-secondary"
            )
          ),
          shiny::hr(),
          shiny::div(style = "font-size:12px;font-weight:700;color:#667085;margin-bottom:6px;", "EXTERNAL RESOURCES"),
          shiny::tags$a(
            href = .subcellular_provider_urls[["Plant-mPLoc"]], target = "_blank",
            class = "btn btn-sm btn-outline-info w-100 mb-1",
            bsicons::bs_icon("box-arrow-up-right"), " Plant-mPLoc"
          ),
          shiny::tags$a(
            href = .subcellular_provider_urls[["WoLF PSORT"]], target = "_blank",
            class = "btn btn-sm btn-outline-info w-100 mb-1",
            bsicons::bs_icon("box-arrow-up-right"), " WoLF PSORT"
          ),
          shiny::tags$a(
            href = .subcellular_provider_urls[["BUSCA"]], target = "_blank",
            class = "btn btn-sm btn-outline-info w-100",
            bsicons::bs_icon("box-arrow-up-right"), " BUSCA"
          )
        )
      )
    ),
    shiny::div(
      class = "p-3",
      bslib::layout_column_wrap(
        width = 1 / 3,
        bslib::card(
          bslib::card_body(
            shiny::div(style = "font-size:12px;color:#667085;font-weight:700;", "STATUS"),
            shiny::uiOutput(ns("run_status"))
          )
        ),
        bslib::card(
          bslib::card_body(
            shiny::div(style = "font-size:12px;color:#667085;font-weight:700;", "PROTEIN"),
            shiny::uiOutput(ns("protein_summary"))
          )
        ),
        bslib::card(
          bslib::card_body(
            shiny::div(style = "font-size:12px;color:#667085;font-weight:700;", "PREDICTED LOCATION"),
            shiny::uiOutput(ns("location_summary"))
          )
        )
      ),
      bslib::card(
        bslib::card_header("DeepLoc 2.1 website result"),
        bslib::card_body(
          shiny::p(
            "If automatic web sources are unavailable, download this protein as FASTA, run DeepLoc 2.1 on its website, and enter the predicted localization labels below.",
            style = "font-size:13px;color:#667085;"
          ),
          shiny::div(
            style = "display:flex;gap:8px;flex-wrap:wrap;",
            shiny::downloadButton(ns("deeploc_fasta"), "DOWNLOAD FASTA", class = "btn-outline-primary"),
            shiny::tags$a(
              href = "https://services.healthtech.dtu.dk/services/DeepLoc-2.1/",
              target = "_blank", rel = "noopener noreferrer",
              class = "btn btn-outline-info",
              bsicons::bs_icon("box-arrow-up-right"), " Open DeepLoc 2.1"
            )
          ),
          shiny::fileInput(ns("deeploc_summary"), "DeepLoc CSV Summary (optional)",
                           accept = c(".csv", "text/csv")),
          shiny::p("A selected Summary CSV supplies labels and probabilities in place of the text fields. For multiple proteins, set the current protein ID to the matching CSV row.",
                   style = "font-size:12px;color:#667085;"),
          shiny::textAreaInput(
            ns("deeploc_location"), "Predicted localizations (semicolon separated)",
            rows = 2, placeholder = "For example: Endoplasmic reticulum"
          ),
          shiny::textInput(ns("deeploc_membrane"), "Predicted membrane association",
                           placeholder = "For example: Soluble"),
          shiny::textInput(ns("deeploc_signals"), "Predicted sorting signals",
                           placeholder = "For example: Signal peptide"),
          shiny::textAreaInput(
            ns("deeploc_probabilities"), "All 10 localization and 4 membrane probabilities (label,probability per line)",
            rows = 6, placeholder = "Endoplasmic reticulum,0.7058\nSoluble,0.8140"
          ),
          shiny::fileInput(ns("deeploc_importance"), "Sorting signal importance CSV (optional)",
                           accept = c(".csv", "text/csv")),
          shiny::downloadButton(ns("deeploc_example"), "DOWNLOAD EXAMPLE CSV",
                                class = "btn-outline-secondary"),
          shiny::tags$details(
            shiny::tags$summary("Example data from the screenshot"),
            shiny::p("Rabbit CASQ1 example: localization Endoplasmic reticulum; membrane Soluble; signal Signal peptide. These values belong to that protein only."),
            shiny::tags$pre(.subcellular_deeploc_example()),
            shiny::p("Sorting importance CSV columns: position,residue,importance. Use DeepLoc's CSV export when available; the screenshot alone does not provide numeric residue scores.")
          ),
          shiny::actionButton(
            ns("import_deeploc"), "ADD DEEPLOC RESULT",
            icon = bsicons::bs_icon("check-circle"), class = "btn-primary"
          )
        )
      ),
      bslib::card(
        full_screen = TRUE, min_height = 650,
        bslib::card_header("Subcellular Localization Results"),
        bslib::card_body(
          bslib::navset_card_tab(
            bslib::nav_panel("Prediction", shiny::uiOutput(ns("prediction_panel"))),
            bslib::nav_panel("Evidence", DT::DTOutput(ns("evidence_table"))),
            bslib::nav_panel("Result table", DT::DTOutput(ns("result_table"))),
            bslib::nav_panel(
              "DeepLoc details",
              shiny::uiOutput(ns("deeploc_details")),
              DT::DTOutput(ns("deeploc_scores")),
              shiny::plotOutput(ns("deeploc_importance_plot"), height = "240px"),
              DT::DTOutput(ns("deeploc_importance_table"))
            ),
            bslib::nav_panel(
              "Raw responses",
              shiny::div(
                style = "max-height:520px; overflow:auto; white-space:pre-wrap;",
                shiny::verbatimTextOutput(ns("raw_result"))
              )
            ),
            bslib::nav_panel(
              "Method",
              shiny::p(
                "ProtVis queries selected services independently and harmonizes their labels. Consensus is based on provider agreement; provider failures are retained as evidence rather than treated as a whole-analysis failure."
              ),
              shiny::tags$ul(
                shiny::tags$li("Plant-mPLoc: plant-specific single- and multi-location prediction."),
                shiny::tags$li("WoLF PSORT: sequence-based eukaryotic localization prediction with a Plant model."),
                shiny::tags$li("BUSCA: integrative localization and targeting-feature prediction for plant proteins.")
              )
            )
          )
        ),
        bslib::card_footer(
          shiny::downloadButton(ns("download"), "DOWNLOAD EVIDENCE", class = "btn-outline-primary"),
          shiny::downloadButton(ns("download_deeploc_scores"), "DOWNLOAD DEEPLOC PROBABILITIES",
                                class = "btn-outline-primary"),
          shiny::downloadButton(ns("download_deeploc_importance"), "DOWNLOAD SORTING IMPORTANCE",
                                class = "btn-outline-primary")
        )
      )
    )
  )
}

#' ProtVis subcellular localization server
#'
#' @param id Module namespace.
#' @param shared_state Optional ProtVis shared state.
#' @return Reactive localization result.
#' @export
subcellular_localization_server <- function(id, shared_state = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    result <- shiny::reactiveVal(NULL)
    status <- shiny::reactiveVal("Ready")
    error_message <- shiny::reactiveVal(NULL)
    running <- shiny::reactiveVal(FALSE)

    output$deeploc_example <- shiny::downloadHandler(
      filename = function() "deeploc_2_1_CASQ1_example.csv",
      content = function(file) {
        probabilities <- .subcellular_deeploc_probabilities(.subcellular_deeploc_example())
        row <- as.list(setNames(rep("", length(c(.subcellular_deeploc_locations,
                                                 .subcellular_deeploc_membranes))),
                                c(.subcellular_deeploc_locations, .subcellular_deeploc_membranes)))
        for (i in seq_len(nrow(probabilities))) {
          row[[probabilities$label[[i]]]] <- probabilities$probability[[i]]
        }
        utils::write.csv(data.frame(
          "Protein ID" = "sp_P07221_CASQ1_RABIT",
          "Predicted localizations" = "Endoplasmic reticulum",
          "Predicted membrane association" = "Soluble",
          "Predicted signals" = "Signal peptide",
          row, check.names = FALSE
        ), file, row.names = FALSE)
      }
    )

    output$deeploc_fasta <- shiny::downloadHandler(
      filename = function() {
        id <- gsub("[^A-Za-z0-9_.-]+", "_", input$protein_id %||% "")
        if (!nzchar(id)) id <- "query_protein"
        paste0(id, "_deeploc.fasta")
      },
      content = function(file) {
        writeLines(.subcellular_deeploc_fasta(input$sequence, input$protein_id), file)
      },
      contentType = "text/plain"
    )

    shiny::observeEvent(input$sequence, {
      result(NULL)
      status("Ready")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$protein_id, {
      result(NULL)
      status("Ready")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$import_deeploc, {
      tryCatch({
        sequence <- .plant_mploc_clean_sequence(input$sequence)
        .subcellular_deeploc_fasta(sequence, input$protein_id)
        summary <- .subcellular_deeploc_summary(input$deeploc_summary, input$protein_id)
        imported <- .subcellular_import_deeploc(
          if (is.null(summary)) input$deeploc_location else summary$location,
          if (is.null(summary)) input$deeploc_membrane else summary$membrane,
          if (is.null(summary)) input$deeploc_signals else summary$signals,
          if (is.null(summary)) input$deeploc_probabilities else summary$probabilities,
          .subcellular_deeploc_importance(input$deeploc_importance),
          sequence = sequence
        )
        value <- result()
        if (is.null(value) || !identical(value$sequence, sequence) ||
            !identical(value$protein_id, input$protein_id)) {
          value <- list(
            protein_id = input$protein_id, sequence = sequence,
            length = nchar(sequence), providers = list(),
            submitted_at = Sys.time()
          )
        }
        value$providers[["DeepLoc 2.1"]] <- imported
        value$consensus <- .subcellular_consensus(value$providers)
        class(value) <- "ProtVis_subcellular_localization"
        evidence <- .subcellular_evidence_table(value)
        .protvis_record_shared_run(
          shared_state,
          module = "subcellular_localization",
          method = "DeepLoc 2.1 web result",
          category = "subcellular_localization",
          parameters = list(
            protein_id = value$protein_id,
            website = imported$result_url,
            entry = "user-entered result"
          ),
          tables = list(evidence = evidence, consensus_votes = value$consensus$votes,
                        deeploc_probabilities = imported$score_table,
                        deeploc_sorting_importance = imported$sorting_importance),
          statistics = list(sequence = sequence)
        )
        result(value)
        error_message(NULL)
        status("Imported")
        shiny::showNotification("DeepLoc 2.1 result added.", type = "message")
      }, error = function(e) {
        shiny::showNotification(
          paste("DeepLoc 2.1:", conditionMessage(e)),
          type = "error", duration = 8
        )
      })
    }, ignoreInit = TRUE)

    output$sequence_status <- shiny::renderUI({
      seq <- toupper(gsub("\\s+", "", input$sequence %||% ""))
      if (!nzchar(seq)) {
        return(shiny::div(style = "font-size:12px;color:#667085;", "No sequence entered"))
      }
      valid <- grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", seq)
      shiny::div(
        style = paste0("font-size:12px;color:", if (valid) "#157347" else "#b42318", ";"),
        if (valid) paste(nchar(seq), "aa; valid alphabet") else "Invalid characters detected"
      )
    })

    shiny::observeEvent(input$clear, {
      shiny::updateTextAreaInput(session, "sequence", value = "")
      result(NULL); error_message(NULL); status("Ready")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$load_demo, {
      shiny::updateTextInput(session, "protein_id", value = "ZmProtein_demo")
      shiny::updateTextAreaInput(session, "sequence", value = .plant_mploc_demo_sequence)
      result(NULL); error_message(NULL); status("Ready")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$run, {
      if (isTRUE(running())) return(invisible(NULL))
      if (!length(input$providers %||% character())) {
        shiny::showNotification("Select at least one prediction source.", type = "warning")
        return(invisible(NULL))
      }
      running(TRUE)
      shinyjs::disable("run")
      on.exit({
        running(FALSE)
        shinyjs::enable("run")
      }, add = TRUE)

      status("Running")
      error_message(NULL)
      previous <- result()
      result(NULL)

      tryCatch({
        seq <- .plant_mploc_clean_sequence(input$sequence)
        manual <- if (!is.null(previous) &&
                      identical(previous$sequence, seq) &&
                      identical(previous$protein_id, input$protein_id)) {
          previous$providers[["DeepLoc 2.1"]]
        } else NULL
        value <- shiny::withProgress(
          message = "Predicting subcellular localization...", value = 0.1,
          {
            ans <- predict_subcellular_localization(
              sequence = seq,
              id = input$protein_id,
              providers = input$providers,
              timeout = input$timeout,
              ignore_ssl = isTRUE(input$ignore_ssl)
            )
            shiny::incProgress(0.9, detail = "Consensus assembled")
            ans
          }
        )
        if (!is.null(manual)) {
          value$providers[["DeepLoc 2.1"]] <- manual
          value$consensus <- .subcellular_consensus(value$providers)
        }
        result(value)
        status(if (value$consensus$successful == value$consensus$requested) {
          "Completed"
        } else if (value$consensus$successful > 0L) {
          "Completed with warnings"
        } else {
          "Unavailable"
        })

        evidence <- .subcellular_evidence_table(value)
        .protvis_record_shared_run(
          shared_state,
          module = "subcellular_localization",
          method = "Consensus localization",
          category = "subcellular_localization",
          parameters = list(
            protein_id = input$protein_id,
            providers = input$providers,
            timeout = input$timeout
          ),
          tables = list(evidence = evidence, consensus_votes = value$consensus$votes,
                        deeploc_probabilities = manual$score_table,
                        deeploc_sorting_importance = manual$sorting_importance)
        )
      }, error = function(e) {
        error_message(conditionMessage(e))
        status("Failed")
        shiny::showNotification(conditionMessage(e), type = "error", duration = 8)
      })
    }, ignoreInit = TRUE)

    output$run_status <- shiny::renderUI({
      s <- status()
      col <- switch(
        s,
        "Completed" = "#157347",
        "Imported" = "#157347",
        "Completed with warnings" = "#b7791f",
        "Unavailable" = "#b7791f",
        "Failed" = "#b42318",
        "#0b6fa4"
      )
      shiny::div(
        style = paste0("font-size:20px;font-weight:700;color:", col, ";margin-top:14px;"),
        s
      )
    })

    output$protein_summary <- shiny::renderUI({
      value <- result()
      if (is.null(value)) {
        return(shiny::div(style = "font-size:20px;color:#0b6fa4;margin-top:14px;", "—"))
      }
      shiny::div(
        style = "margin-top:12px;",
        shiny::div(style = "font-weight:700;", value$protein_id),
        shiny::div(style = "font-size:12px;color:#667085;", paste(value$length, "aa"))
      )
    })

    output$location_summary <- shiny::renderUI({
      value <- result()
      if (is.null(value) || !length(value$consensus$prediction)) {
        return(shiny::div(style = "font-size:20px;color:#0b6fa4;margin-top:14px;", "—"))
      }
      shiny::div(
        style = "margin-top:10px;",
        shiny::div(
          style = "font-size:18px;font-weight:700;color:#0b6fa4;",
          paste(value$consensus$prediction, collapse = "; ")
        ),
        shiny::div(
          style = "font-size:12px;color:#667085;margin-top:3px;",
          paste(value$consensus$confidence, "confidence ·",
                value$consensus$successful, "/", value$consensus$requested, "sources available")
        )
      )
    })

    output$prediction_panel <- shiny::renderUI({
      value <- result()
      if (is.null(value)) {
        msg <- error_message()
        if (!is.null(msg)) {
          return(shiny::div(
            class = "alert alert-danger",
            shiny::strong("Analysis failed"), shiny::br(), msg
          ))
        }
        return(shiny::div(
          style = "padding:28px;color:#667085;",
          "Run the analysis to assemble localization evidence."
        ))
      }

      if (!length(value$consensus$prediction)) {
        failures <- vapply(value$providers, function(x) {
          paste0(x$source, ": ", x$error %||% "No usable result")
        }, character(1))
        return(shiny::div(
          class = "alert alert-warning",
          shiny::strong("No automatic source returned a usable prediction."),
          shiny::tags$ul(lapply(failures, shiny::tags$li)),
          "Use the DeepLoc 2.1 website result panel above to add a verified prediction."
        ))
      }

      provider_cards <- lapply(value$providers, function(x) {
        ok <- identical(x$status, "success")
        shiny::div(
          style = paste0(
            "border:1px solid ", if (ok) "#cfe8f6" else "#ead7a5",
            ";border-radius:10px;padding:12px;margin-top:8px;background:#fff;"
          ),
          shiny::div(style = "font-weight:700;", x$source),
          shiny::div(
            style = paste0("font-size:13px;color:", if (ok) "#157347" else "#9a6700", ";"),
            if (ok) paste(x$prediction, collapse = "; ") else paste("Unavailable —", x$error %||% "")
          )
        )
      })

      shiny::tagList(
        shiny::div(
          style = paste0(
            "padding:18px;border-radius:12px;background:#f7fbfe;",
            "border:1px solid #9bd4f5;color:#075985;"
          ),
          shiny::div(style = "font-size:12px;font-weight:700;", "CONSENSUS PREDICTION"),
          shiny::div(
            style = "font-size:28px;font-weight:750;margin-top:5px;",
            paste(value$consensus$prediction, collapse = "; ")
          ),
          shiny::div(
            style = "font-size:13px;margin-top:5px;",
            paste0(
              value$consensus$confidence, " confidence · ",
              value$consensus$agreement, "/", value$consensus$successful,
              " available providers agree on the leading location"
            )
          )
        ),
        provider_cards
      )
    })

    output$evidence_table <- DT::renderDT({
      value <- result()
      shiny::req(value)
      DT::datatable(
        .subcellular_evidence_table(value),
        rownames = FALSE,
        options = list(dom = "t", pageLength = 10)
      )
    })

    output$result_table <- DT::renderDT({
      value <- result()
      shiny::req(value)
      x <- value$consensus$votes
      if (!nrow(x)) {
        x <- data.frame(location = character(), votes = integer(), agreement = character())
      }
      DT::datatable(x, rownames = FALSE, options = list(dom = "t", pageLength = 12))
    })

    deeploc_result <- shiny::reactive({
      value <- result()
      shiny::req(value)
      value$providers[["DeepLoc 2.1"]]
    })

    output$deeploc_details <- shiny::renderUI({
      x <- deeploc_result()
      if (is.null(x)) return(shiny::p("Add a DeepLoc 2.1 result to view its full output."))
      shiny::tagList(
        shiny::p(shiny::strong("Predicted localizations: "), paste(x$prediction, collapse = "; ")),
        shiny::p(shiny::strong("Membrane association: "),
                 if (length(x$membrane_types)) paste(x$membrane_types, collapse = "; ") else "Not supplied"),
        shiny::p(shiny::strong("Sorting signals: "),
                 if (length(x$signals)) paste(x$signals, collapse = "; ") else "Not supplied"),
        shiny::p("Probabilities are the model's per-label scores. Consensus votes use predicted labels, not probability values.")
      )
    })

    output$deeploc_scores <- DT::renderDT({
      x <- deeploc_result()
      shiny::req(x)
      tab <- x$score_table
      if (is.null(tab)) tab <- .subcellular_deeploc_probabilities("")
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 14))
    })

    output$deeploc_importance_plot <- shiny::renderPlot({
      x <- deeploc_result()
      shiny::req(x, x$sorting_importance)
      d <- x$sorting_importance[order(x$sorting_importance$position), , drop = FALSE]
      graphics::plot(d$position, d$importance, type = "l", ylim = c(0, 1),
                     xlab = "Residue position", ylab = "Sorting signal importance",
                     col = "#168fd0", lwd = 1.5)
    })

    output$deeploc_importance_table <- DT::renderDT({
      x <- deeploc_result()
      shiny::req(x)
      d <- x$sorting_importance
      if (is.null(d)) d <- data.frame(position = integer(), residue = character(),
                                      importance = numeric())
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })

    output$raw_result <- shiny::renderText({
      value <- result()
      shiny::req(value)
      paste(vapply(value$providers, function(x) {
        paste0(
          "===== ", x$source, " =====\n",
          "Status: ", x$status, "\n",
          if (!is.null(x$error) && nzchar(x$error)) paste0("Error: ", x$error, "\n") else "",
          if (nzchar(x$raw_result %||% "")) x$raw_result else "(no raw response retained)"
        )
      }, character(1)), collapse = "\n\n")
    })

    output$download <- shiny::downloadHandler(
      filename = function() {
        paste0(
          "subcellular_localization_",
          gsub("[^A-Za-z0-9_.-]", "_", input$protein_id %||% "protein"),
          ".csv"
        )
      },
      content = function(file) {
        value <- result()
        shiny::req(value)
        utils::write.csv(.subcellular_evidence_table(value), file, row.names = FALSE)
      }
    )

    output$download_deeploc_scores <- shiny::downloadHandler(
      filename = function() "deeploc_probabilities.csv",
      content = function(file) {
        x <- deeploc_result()
        shiny::req(x, x$score_table)
        utils::write.csv(x$score_table, file, row.names = FALSE)
      }
    )
    output$download_deeploc_importance <- shiny::downloadHandler(
      filename = function() "deeploc_sorting_importance.csv",
      content = function(file) {
        x <- deeploc_result()
        shiny::req(x, x$sorting_importance)
        utils::write.csv(x$sorting_importance, file, row.names = FALSE)
      }
    )

    result
  })
}
