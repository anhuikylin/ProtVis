# Experimental plant PTM benchmarks. Source identifications are not re-scored.
.protvis_arabidopsis_glyco_target <- function() {
  root <- system.file("extdata", "ptm", "arabidopsis_glyco", package = "ProtVis")
  if (!nzchar(root)) {
    root <- file.path("inst", "extdata", "ptm", "arabidopsis_glyco")
  }
  meta <- utils::read.delim(file.path(root, "metadata.tsv"),
                          stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  values <- stats::setNames(as.list(meta$value), meta$field)
  list(benchmark_id = values$benchmark_id, ptm_type = "glycosylation",
       organism = values$organism, project = values$project,
       protein = values$protein, sequence = values$sequence,
       modified_sequence = "N[HexNAc]VTHAPRPGGFSSSVVSGLSQGSGEYFTR",
       fragmentation = "ETD", position = as.integer(values$modified_position),
       mass = as.numeric(values$mass_shift_da),
       precursor_mz = as.numeric(values$precursor_mz),
       precursor_charge = as.integer(values$precursor_charge),
       spectrum_title = values$spectrum_title, provenance = values,
       mgf = file.path(root, "arabidopsis_N117_HexNAc.mgf"))
}

.protvis_arabidopsis_glyco_bundle <- function() {
  .protvis_plant_ptm_bundle(.protvis_arabidopsis_glyco_target())
}

.protvis_arabidopsis_methyl_bundle <- function() {
  .protvis_plant_ptm_bundle(.protvis_arabidopsis_methyl_target())
}

.protvis_plant_ptm_bundle <- function(target) {
  .protvis_ptm_require_spectrum_packages()
  psm <- data.frame(sequence = target$sequence, spectrumID = "index=0",
                    chargeState = target$precursor_charge, passThreshold = TRUE,
                    experimentalMassToCharge = target$precursor_mz,
                    calculatedMassToCharge = target$precursor_mz,
                    spectrum.title = target$spectrum_title,
                    stringsAsFactors = FALSE, check.names = FALSE)
  psm$DatabaseAccess <- I(list(target$protein))
  psm$modLocation <- I(list(target$position))
  psm$modMass <- I(list(target$mass))
  psm$modName <- I(list(target$modification %||% "HexNAc"))
  if (identical(target$ptm_type, "methylation")) {
    psm$Mascot.score <- as.numeric(target$provenance$score)
  }
  spectra <- Spectra::Spectra(target$mgf, source = MsBackendMgf::MsBackendMgf())
  if (length(spectra) != 1L) stop("Expected one experimental plant PTM spectrum.")
  list(psm = psm, catalog = .protvis_ptm_psm_catalog(psm), spectra = spectra,
       metadata = as.data.frame(Spectra::spectraData(spectra), optional = TRUE),
       source = target$source %||% "Author MSViewer experimental ETD scan 3863; source-assigned HexNAc at Asn117",
       benchmark = target)
}

# c = b + NH3; z radical = y - NH2. Charge conversion uses proton mass.
# ETD cleavage N-terminal to proline is excluded. No glycan-loss inference.
.protvis_ptm_etd_theoretical <- function(sequence, modifications = NULL) {
  theory <- .protvis_ptm_theoretical(sequence, modifications)
  proton <- 1.007276466621
  c_mz <- theory$b_mz + 17.026549101
  z_mz <- theory$y_mz - 16.018724068
  n <- nchar(sequence)
  residues <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  valid <- residues[2:n] != "P"
  c_mz[!valid] <- NA_real_
  z_mz[!rev(valid)] <- NA_real_
  candidates <- list()
  for (charge in 1:3) {
    suffix <- if (charge == 1L) "" else paste(rep("+", charge), collapse = "")
    for (series in c("c", "z")) {
      masses <- if (series == "c") c_mz else z_mz
      candidates[[length(candidates) + 1L]] <- data.frame(
        label = paste0(series, seq_along(masses), suffix),
        mz = (masses + (charge - 1L) * proton) / charge,
        series = series, charge = charge, neutral = FALSE,
        priority = charge, stringsAsFactors = FALSE)
    }
  }
  ions <- do.call(rbind, candidates)
  ions <- ions[is.finite(ions$mz) & ions$mz > 0, , drop = FALSE]
  theory$candidates <- ions[order(ions$mz), , drop = FALSE]
  theory$key_ions <- ions[ions$charge == 1L, c("label", "mz", "series"), drop = FALSE]
  c_full <- c(c_mz, NA_real_)
  z_full <- c(NA_real_, rev(z_mz))
  theory$fragment_table <- data.frame(
    B = seq_len(n), `C Ions` = c_full, `C+2H` = (c_full + proton) / 2,
    `C+3H` = (c_full + 2 * proton) / 3, AA = theory$fragment_table$AA,
    `Z Ions` = z_full, `Z+2H` = (z_full + proton) / 2,
    `Z+3H` = (z_full + 2 * proton) / 3, Y = n:1L, check.names = FALSE)
  numeric <- vapply(theory$fragment_table, is.numeric, logical(1L))
  theory$fragment_table[numeric] <- lapply(theory$fragment_table[numeric], round, 1L)
  theory$fragmentation <- "ETD"
  theory
}

.protvis_arabidopsis_methyl_target <- function() {
  root <- system.file("extdata", "ptm", "arabidopsis_methyl", package = "ProtVis")
  if (!nzchar(root)) root <- file.path("inst", "extdata", "ptm", "arabidopsis_methyl")
  meta <- utils::read.delim(file.path(root, "metadata.tsv"), stringsAsFactors = FALSE, quote = "")
  values <- stats::setNames(as.list(meta$value), meta$field)
  list(benchmark_id = values$benchmark_id, ptm_type = "methylation",
       organism = values$organism, project = values$project,
       protein = values$protein, sequence = values$sequence,
       modified_sequence = "GGR[Dimethyl]GYGQPPQQQQQYGGPQEYQGR",
       modification = "Dimethyl", fragmentation = "HCD",
       position = as.integer(values$modified_position), mass = as.numeric(values$mass_shift_da),
       precursor_mz = as.numeric(values$precursor_mz),
       precursor_charge = as.integer(values$precursor_charge),
       spectrum_title = values$spectrum_title, provenance = values,
       source = "PRIDE PXD043460 experimental HCD scan 18622; rank 1, passThreshold=true",
       mgf = file.path(root, "arabidopsis_AGO1_R62_Dimethyl.mgf"))
}
