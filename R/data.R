#' Read BPC files from a folder
#'
#' Finds the needed BPC files in a folder by name and reads them.
#'
#' @param dir Folder holding the BPC csv files for one cohort.
#' @param patterns Named file-name patterns for each table.
#' @return A named list of data frames: `ca_dx`, `regimens`, `imaging`,
#'   `med_onc`. Missing files come back as `NULL` with a warning.
#' @examples
#' \dontrun{
#' raw <- read_bpc("~/genie/NSCLC_2.0-public")
#' x <- lot_data(raw$ca_dx, raw$regimens, raw$imaging, raw$med_onc)
#' }
#' @export
read_bpc <- function(dir,
                     patterns = c(ca_dx = "cancer_level_dataset_index",
                                  regimens = "regimen_cancer_level",
                                  imaging = "imaging_level",
                                  med_onc = "med_onc_note_level")) {
  files <- list.files(dir, full.names = TRUE)
  lapply(patterns, function(p) {
    f <- files[grepl(p, basename(files), ignore.case = TRUE)]
    if (!length(f)) {
      warning("No file matching '", p, "' in ", dir, call. = FALSE)
      return(NULL)
    }
    utils::read.csv(f[1], stringsAsFactors = FALSE, check.names = FALSE)
  })
}

#' Prepare BPC tables for LOT calculation
#'
#' Turns the raw BPC tables into three tidy tables: one row per cancer case,
#' one row per drug taken, and one row per progression event. All times are
#' days since diagnosis of the index cancer.
#'
#' @param ca_dx Cancer diagnosis table (`cancer_level_dataset_index`).
#' @param regimens Regimen table (`regimen_cancer_level_dataset`).
#' @param imaging Optional imaging table (`imaging_level_dataset`).
#' @param med_onc Optional oncologist note table (`med_onc_note_level_dataset`).
#' @param cols Column mapping from [bpc_columns()].
#' @return An object of class `lot_data`.
#' @examples
#' sim <- simulate_bpc(50)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' x
#' @export
lot_data <- function(ca_dx, regimens, imaging = NULL, med_onc = NULL,
                     cols = bpc_columns()) {
  # cases ---------------------------------------------------------------
  rid <- as.character(get_col(ca_dx, cols$id, "diagnosis"))
  seq <- if (cols$ca_seq %in% names(ca_dx)) as.character(ca_dx[[cols$ca_seq]]) else rep("0", length(rid))
  cases <- data.frame(
    case_id = paste(rid, seq, sep = "_"),
    record_id = rid,
    ca_seq = seq,
    cohort = as.character(get_col(ca_dx, cols$cohort, "diagnosis", FALSE)),
    stage = as.character(get_col(ca_dx, cols$stage, "diagnosis")),
    dmets_day = to_days(get_col(ca_dx, cols$dmets_time, "diagnosis", FALSE), cols$dmets_unit),
    followup_day = to_days(get_col(ca_dx, cols$os_time, "diagnosis", FALSE), cols$os_unit),
    stringsAsFactors = FALSE
  )
  cases$early_stage <- is_early_stage(cases$stage)
  cases <- cases[!duplicated(cases$case_id), , drop = FALSE]

  # drugs (one row per drug per regimen) -----------------------------------
  rr <- as.character(get_col(regimens, cols$id, "regimen"))
  rs <- if (cols$ca_seq %in% names(regimens)) as.character(regimens[[cols$ca_seq]]) else rep("0", length(rr))
  rkey <- paste(rr, rs, sep = "_")
  rnum <- get_col(regimens, cols$reg_number, "regimen", FALSE)
  rdrugs <- as.character(get_col(regimens, cols$reg_drugs, "regimen"))
  rstart <- to_days(get_col(regimens, cols$reg_start, "regimen"), cols$reg_unit)
  rend <- to_days(get_col(regimens, cols$reg_end, "regimen", FALSE), cols$reg_unit)
  parts <- lapply(rdrugs, function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(character(0))
    p <- trimws(strsplit(s, cols$reg_sep, fixed = TRUE)[[1]])
    unique(p[nzchar(p)])
  })
  n <- lengths(parts)
  drugs <- data.frame(
    case_id = rep(rkey, n),
    regimen_number = rep(rnum, n),
    drug = unlist(parts, use.names = FALSE),
    start = rep(rstart, n),
    end = rep(rend, n),
    stringsAsFactors = FALSE
  )
  drugs <- drugs[drugs$case_id %in% cases$case_id & !is.na(drugs$start), , drop = FALSE]
  drugs$end[is.na(drugs$end) | drugs$end < drugs$start] <-
    drugs$start[is.na(drugs$end) | drugs$end < drugs$start]
  drugs <- drugs[order(drugs$case_id, drugs$start, drugs$drug), , drop = FALSE]
  rownames(drugs) <- NULL

  # progression events --------------------------------------------------------
  prog <- rbind(
    prog_events(imaging, cols, cols$img_time, cols$img_status, cols$img_unit, "imaging"),
    prog_events(med_onc, cols, cols$md_time, cols$md_status, cols$md_unit, "med_onc")
  )
  if (is.null(prog)) {
    prog <- data.frame(case_id = character(), day = numeric(), source = character())
  } else {
    m <- merge(prog, cases[, c("case_id", "record_id", "ca_seq")], by = "record_id")
    m <- m[is.na(m$ev_seq) | m$ev_seq == m$ca_seq, , drop = FALSE]
    prog <- unique(m[order(m$case_id, m$day), c("case_id", "day", "source")])
    rownames(prog) <- NULL
  }

  structure(list(cases = cases, drugs = drugs, prog = prog, cols = cols),
            class = "lot_data")
}

prog_events <- function(df, cols, tcol, scol, unit, source) {
  if (is.null(df) || !nrow(df)) return(NULL)
  rid <- as.character(get_col(df, cols$id, source))
  t <- to_days(get_col(df, tcol, source), unit)
  s <- as.character(get_col(df, scol, source))
  sq <- if (cols$ca_seq %in% names(df)) as.character(df[[cols$ca_seq]]) else rep(NA_character_, length(rid))
  hit <- !is.na(s) & !is.na(t) & grepl(cols$progression_pattern, s, ignore.case = TRUE)
  data.frame(record_id = rid[hit], ev_seq = sq[hit], day = t[hit],
             source = rep(source, sum(hit)), stringsAsFactors = FALSE)
}

#' @export
print.lot_data <- function(x, ...) {
  cs <- x$cases
  cat("<lot_data>\n")
  cat(sprintf("  Cancer cases:        %d (early stage I-III: %d, stage IV: %d)\n",
              nrow(cs), sum(cs$early_stage %in% TRUE), sum(cs$early_stage %in% FALSE)))
  cat(sprintf("  Drug records:        %d across %d regimens\n",
              nrow(x$drugs), nrow(unique(x$drugs[, c("case_id", "start", "end")]))))
  cat(sprintf("  Progression events:  %d (imaging %d, oncologist notes %d)\n",
              nrow(x$prog), sum(x$prog$source == "imaging"), sum(x$prog$source == "med_onc")))
  invisible(x)
}
