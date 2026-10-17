#' Compare your LOT definition with naive counting from diagnosis
#'
#' Naive counting starts at diagnosis, so early-stage patients' adjuvant or
#' neoadjuvant treatment is wrongly counted as line 1. This shows how often
#' that happens.
#'
#' @param x A [lot_data()] object.
#' @param index Index dates from [set_index()].
#' @param rules Rules from [lot_rules()].
#' @return A list with `summary` (headline numbers) and `cases` (one row per
#'   indexed patient: their correct line 1, what naive counting called it, and
#'   how many pre-index lines were excluded).
#' @examples
#' sim <- simulate_bpc(100)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' compare_lot(x)$summary
#' @export
compare_lot <- function(x, index = set_index(x), rules = lot_rules()) {
  correct <- calculate_lot(x, index, rules)
  naive <- calculate_lot(x, set_index(x, "diagnosis"), rules)
  cs <- x$cases
  idx <- as.data.frame(index)
  keep <- idx$case_id[!is.na(idx$index_day)]
  rows <- lapply(keep, function(id) {
    c1 <- correct[correct$case_id == id & correct$line == 1, ]
    nv <- naive[naive$case_id == id, ]
    iday <- idx$index_day[idx$case_id == id]
    naive_l1 <- nv[nv$line == 1, ]
    # which naive line number lands on the true line 1 start?
    naive_label <- if (nrow(c1)) {
      hit <- nv$line[nv$start <= c1$start & nv$end >= c1$start]
      if (!length(hit)) hit <- nv$line[which.min(abs(nv$start - c1$start))]
      max(hit)
    } else NA_integer_
    data.frame(case_id = id,
               early_stage = cs$early_stage[cs$case_id == id],
               index_event = idx$index_event[idx$case_id == id],
               index_day = iday,
               correct_line1 = if (nrow(c1)) c1$regimen else NA_character_,
               correct_line1_start = if (nrow(c1)) c1$start else NA_real_,
               naive_line1 = if (nrow(naive_l1)) naive_l1$regimen else NA_character_,
               naive_line1_start = if (nrow(naive_l1)) naive_l1$start else NA_real_,
               naive_label_of_true_line1 = naive_label,
               pre_index_lines_excluded = sum(nv$start < iday - rules$lookback),
               stringsAsFactors = FALSE)
  })
  tab <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(tab)) {
    tab$mislabeled <- !is.na(tab$correct_line1) &
      (is.na(tab$naive_line1_start) | tab$naive_line1_start != tab$correct_line1_start)
  }
  treated <- tab[!is.na(tab$correct_line1), , drop = FALSE]
  early <- treated[treated$early_stage %in% TRUE &
                     !grepl("^diagnosis", treated$index_event), , drop = FALSE]
  pct <- function(a, b) if (b) round(100 * a / b, 1) else NA_real_
  summary <- data.frame(
    group = c("All indexed, treated patients", "Early stage (I-III) recurrences"),
    patients = c(nrow(treated), nrow(early)),
    mislabeled = c(sum(treated$mislabeled), sum(early$mislabeled)),
    pct_mislabeled = c(pct(sum(treated$mislabeled), nrow(treated)),
                       pct(sum(early$mislabeled), nrow(early))),
    pre_index_lines_excluded = c(sum(treated$pre_index_lines_excluded),
                                 sum(early$pre_index_lines_excluded)),
    stringsAsFactors = FALSE
  )
  structure(list(summary = summary, cases = tab, correct = correct, naive = naive),
            class = "lot_comparison")
}

#' @export
print.lot_comparison <- function(x, ...) {
  cat("<lot_comparison> your index vs naive counting from diagnosis\n\n")
  print(x$summary, row.names = FALSE)
  e <- x$summary[2, ]
  if (!is.na(e$pct_mislabeled)) {
    cat(sprintf("\nNaive counting gives the wrong line 1 for %s%% of early-stage patients who recurred.\n",
                e$pct_mislabeled))
  }
  invisible(x)
}

#' Check your index against a recurrence flag in the data
#'
#' BPC includes `dmets_stage_i_iii`, a flag for stage I-III patients who later
#' developed distant metastasis. This compares it with the index events you
#' chose, for early-stage patients.
#'
#' @param x A [lot_data()] object.
#' @param index Index dates from [set_index()].
#' @param ca_dx The raw diagnosis table that holds the flag.
#' @param flag_col Name of the flag column (1 = recurred with distant spread).
#' @param cols Column mapping from [bpc_columns()].
#' @return A list with `table` (counts) and `summary` (percentages).
#' @examples
#' \dontrun{
#' validate_index(x, set_index(x, "recurrence"), raw$ca_dx, cols = cols)
#' }
#' @export
validate_index <- function(x, index, ca_dx, flag_col = "dmets_stage_i_iii",
                           cols = x$cols) {
  key <- paste(ca_dx[[cols$id]],
               if (cols$ca_seq %in% names(ca_dx)) ca_dx[[cols$ca_seq]] else "0",
               sep = "_")
  if (!flag_col %in% names(ca_dx)) stop("Column '", flag_col, "' not found.", call. = FALSE)
  idx <- as.data.frame(index)
  early <- x$cases$case_id[x$cases$early_stage %in% TRUE]
  v <- idx[idx$case_id %in% early, ]
  v$flag <- ca_dx[[flag_col]][match(v$case_id, key)]
  v$ours <- ifelse(is.na(v$index_day), "no index",
                   ifelse(grepl("upstaged", v$index_event), "upstaged (index at diagnosis)",
                          "recurrence found"))
  tab <- table(ours = v$ours, flag = v$flag, useNA = "ifany")
  f1 <- v$flag %in% c(1, "1", "Yes")
  f0 <- v$flag %in% c(0, "0", "No")
  handled <- sum(f1 & v$ours != "no index")
  agree <- sum(f1 & v$ours != "no index") + sum(f0 & v$ours == "no index")
  summary <- data.frame(
    measure = c("Flagged patients handled (recurrence or upstaged)",
                "Overall agreement with flag",
                "Extra recurrences found (local or progression only)"),
    value = c(sprintf("%d of %d (%.1f%%)", handled, sum(f1), 100 * handled / max(1, sum(f1))),
              sprintf("%d of %d (%.1f%%)", agree, sum(f1 | f0), 100 * agree / max(1, sum(f1 | f0))),
              sprintf("%d", sum(f0 & v$ours == "recurrence found"))),
    stringsAsFactors = FALSE)
  list(table = tab, summary = summary)
}

#' One-row summary of a cohort, for comparing cancers side by side
#'
#' @param x A [lot_data()] object.
#' @param index Index dates from [set_index()].
#' @param rules Rules from [lot_rules()] or [lot_preset()].
#' @param ca_dx Optional raw diagnosis table with `dmets_stage_i_iii`, for the
#'   validation column.
#' @param label Cohort name. Default: the cohort column in the data.
#' @return A one-row data frame. Bind several with `rbind()`.
#' @examples
#' sim <- simulate_bpc(80)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' cohort_summary(x, set_index(x), lot_preset("NSCLC"), sim$ca_dx)
#' @export
cohort_summary <- function(x, index = set_index(x), rules = lot_rules(),
                           ca_dx = NULL, label = NULL) {
  cmp <- compare_lot(x, index, rules)
  sm <- summarize_lot(cmp$correct)
  val <- if (!is.null(ca_dx) && "dmets_stage_i_iii" %in% names(ca_dx)) {
    v <- validate_index(x, index, ca_dx)$summary$value[1]
    regmatches(v, regexpr("[0-9.]+%", v))
  } else NA_character_
  if (is.null(label)) label <- paste(unique(stats::na.omit(x$cases$cohort)), collapse = "/")
  data.frame(
    cohort = label,
    cancer_cases = nrow(x$cases),
    early_stage = sum(x$cases$early_stage %in% TRUE),
    stage_iv = sum(x$cases$early_stage %in% FALSE),
    patients_with_lines = cmp$summary$patients[1],
    early_recurrences_treated = cmp$summary$patients[2],
    naive_line1_wrong = cmp$summary$mislabeled[2],
    pct_naive_line1_wrong = cmp$summary$pct_mislabeled[2],
    pct_reach_line2 = if (nrow(sm) >= 2) sm$pct_of_indexed[2] else 0,
    median_days_line1 = if (nrow(sm)) sm$median_days_on_line[1] else NA,
    flagged_recurrences_handled = val,
    stringsAsFactors = FALSE
  )
}
