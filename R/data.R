#' Read BPC files from a folder
#'
#' Finds the needed BPC files in a folder by name and reads them.
#'
#' @param dir Folder holding the BPC csv files for one cohort.
#' @param patterns Named file-name patterns for each table.
#' @return A named list of data frames: `ca_dx`, `regimens`, `imaging`,
#'   `med_onc`, `pathology`. Missing files come back as `NULL` with a warning.
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
                                  med_onc = "med_onc_note_level",
                                  pathology = "pathology_report_level")) {
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
#' @param therapy Optional surgery or radiation events, from
#'   [therapy_events()]: one row per event with `record_id`, `day` (days since
#'   diagnosis) and `type`. Used by the therapy index events.
#' @param cols Column mapping from [bpc_columns()].
#' @return An object of class `lot_data`.
#' @examples
#' sim <- simulate_bpc(50)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' x
#' @export
lot_data <- function(ca_dx, regimens, imaging = NULL, med_onc = NULL,
                     therapy = NULL, cols = bpc_columns()) {
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
  end_source <- ifelse(is.na(rend), NA_character_, "recorded")
  if (!is.null(cols$reg_end_lastdose) && !is.na(cols$reg_end_lastdose)) {
    ld_cols <- grep(cols$reg_end_lastdose, names(regimens), value = TRUE)
    if (length(ld_cols)) {
      ld <- sapply(ld_cols, function(cl) to_days(regimens[[cl]], cols$reg_unit))
      ld <- matrix(ld, nrow = nrow(regimens))
      last_dose <- suppressWarnings(apply(ld, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)))
      fill <- is.na(rend) & !is.na(last_dose)
      rend[fill] <- last_dose[fill]; end_source[fill] <- "last dose"
    }
  }
  if (!is.null(cols$reg_end_fallback) && !is.na(cols$reg_end_fallback) &&
      cols$reg_end_fallback %in% names(regimens)) {
    fb <- to_days(regimens[[cols$reg_end_fallback]], cols$reg_unit)
    fill <- is.na(rend) & !is.na(fb)
    rend[fill] <- fb[fill]; end_source[fill] <- "first drug end"
  }
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
    end_source = rep(end_source, n),
    stringsAsFactors = FALSE
  )
  drugs <- drugs[drugs$case_id %in% cases$case_id & !is.na(drugs$start), , drop = FALSE]
  bad_end <- is.na(drugs$end) | drugs$end < drugs$start
  drugs$end_source[bad_end] <- "missing (start used)"
  drugs$end[bad_end] <- drugs$start[bad_end]
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

  # surgery / radiation events -----------------------------------------------
  th <- data.frame(case_id = character(), day = numeric(), type = character(),
                   stringsAsFactors = FALSE)
  if (!is.null(therapy) && nrow(therapy)) {
    tt <- therapy
    tt$ev_seq <- if ("ca_seq" %in% names(tt)) as.character(tt$ca_seq) else NA_character_
    tt$record_id <- as.character(tt$record_id)
    m <- merge(tt[, c("record_id", "ev_seq", "day", "type")],
               cases[, c("case_id", "record_id", "ca_seq")], by = "record_id")
    m <- m[(is.na(m$ev_seq) | m$ev_seq == m$ca_seq) & !is.na(m$day), , drop = FALSE]
    th <- unique(m[order(m$case_id, m$day), c("case_id", "day", "type")])
    rownames(th) <- NULL
  }

  structure(list(cases = cases, drugs = drugs, prog = prog, therapy = th, cols = cols),
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
  reg <- unique(x$drugs[, c("case_id", "start", "end", "end_source")])
  tab <- table(factor(reg$end_source, levels = c("recorded", "last dose", "first drug end", "missing (start used)")))
  cat(sprintf("  Regimen end dates:   %s\n",
              paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = ", ")))
  if (!is.null(x$therapy) && nrow(x$therapy)) {
    tt <- table(x$therapy$type)
    cat(sprintf("  Other therapy:       %s\n", paste(sprintf("%s %d", names(tt), as.integer(tt)), collapse = ", ")))
  }
  cat(sprintf("  Progression events:  %d (imaging %d, oncologist notes %d)\n",
              nrow(x$prog), sum(x$prog$source == "imaging"), sum(x$prog$source == "med_onc")))
  invisible(x)
}

#' Build surgery or radiation events from any BPC table
#'
#' Turns rows of a table (for example the pathology report table, or a
#' radiation table) into therapy events for [lot_data()].
#'
#' @param df The raw table.
#' @param day_col Column with days since diagnosis of the event.
#' @param type Label for these events, e.g. "surgery" or "radiation".
#' @param filter_col Optional column to filter on (e.g. the procedure type).
#' @param pattern Regular expression; rows whose `filter_col` matches are kept
#'   (case-insensitive), e.g. "resection|excision|ectomy".
#' @param id_col Patient id column.
#' @param unit Unit of `day_col`: "days", "months" or "years".
#' @return A data frame with `record_id`, `ca_seq` (if present), `day`, `type`.
#' @examples
#' path <- data.frame(record_id = c("A", "A", "B"),
#'                    path_proc = c("Biopsy", "Surgical resection", "Lobectomy"),
#'                    dx_path_proc_days = c(5, 40, 30))
#' therapy_events(path, "dx_path_proc_days", "surgery",
#'                filter_col = "path_proc", pattern = "resection|ectomy")
#' @export
therapy_events <- function(df, day_col, type, filter_col = NULL, pattern = NULL,
                           id_col = "record_id", unit = "days") {
  keep <- rep(TRUE, nrow(df))
  if (!is.null(filter_col)) {
    v <- as.character(get_col(df, filter_col, type))
    keep <- !is.na(v) & grepl(pattern, v, ignore.case = TRUE)
  }
  out <- data.frame(record_id = as.character(get_col(df, id_col, type))[keep],
                    day = to_days(get_col(df, day_col, type), unit)[keep],
                    type = rep(type, sum(keep)), stringsAsFactors = FALSE)
  if ("ca_seq" %in% names(df)) out$ca_seq <- as.character(df$ca_seq[keep])
  out[!is.na(out$day), , drop = FALSE]
}

#' Surgery events from the BPC pathology report table
#'
#' Uses pathology reports of type "Surgical pathology" as a proxy for surgery
#' dates. Some surgical pathology reports come from diagnostic procedures,
#' so treat this as an approximation.
#'
#' @param path The pathology report table (`pathology_report_level_dataset`).
#' @param type_col Procedure type column.
#' @param day_col Days-since-diagnosis column.
#' @param pattern Which procedure types count as surgery.
#' @return Events for the `therapy` argument of [lot_data()].
#' @examples
#' path <- data.frame(record_id = c("A", "A"), path_proc_type = c("Cytology", "Surgical pathology"),
#'                    dx_path_proc_days = c(3, 25))
#' bpc_surgery_events(path)
#' @export
bpc_surgery_events <- function(path, type_col = "path_proc_type",
                               day_col = "dx_path_proc_days", pattern = "^surgical") {
  therapy_events(path, day_col, "surgery", filter_col = type_col, pattern = pattern)
}
