#' Sensitivity analysis across LOT definitions
#'
#' Runs the same data under several reasonable definitions and puts the
#' headline results side by side. Standard practice in real-world evidence:
#' it shows whether a conclusion depends on the definition chosen.
#'
#' @param x A [lot_data()] object.
#' @param scenarios Named list; each element a list with `index` (from
#'   [set_index()]) and `rules` (from [lot_rules()]). Default: a standard set.
#' @param k Mask counts below this (see [suppress_small()]); `0` turns it off.
#' @return A data frame, one row per scenario.
#' @examples
#' sim <- simulate_bpc(80)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' lot_sensitivity(x)
#' @export
lot_sensitivity <- function(x, scenarios = default_scenarios(x), k = 10) {
  rows <- lapply(names(scenarios), function(nm) {
    sc <- scenarios[[nm]]
    cmp <- compare_lot(x, sc$index, sc$rules)
    sm <- summarize_lot(cmp$correct)
    pct_line <- function(l) if (nrow(sm) >= l) sm$pct_of_indexed[l] else 0
    data.frame(
      scenario = nm,
      patients_with_lines = cmp$summary$patients[1],
      early_recurrences = cmp$summary$patients[2],
      pct_naive_line1_wrong = cmp$summary$pct_mislabeled[2],
      pct_reach_line2 = pct_line(2),
      pct_reach_line3 = pct_line(3),
      median_days_line1 = if (nrow(sm)) sm$median_days_on_line[1] else NA,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (k > 0) out <- suppress_small(out, cols = c("patients_with_lines", "early_recurrences"), k = k)
  out
}

#' The standard set of definitions used by [lot_sensitivity()]
#'
#' @param x A [lot_data()] object.
#' @param rules Base rules. Default [lot_rules()].
#' @return A named list of scenarios.
#' @export
default_scenarios <- function(x, rules = lot_rules()) {
  r <- function(...) { a <- unclass(rules); b <- list(...); a[names(b)] <- b; do.call(lot_rules, a[names(formals(lot_rules))]) }
  list(
    "Base: recurrence, 28-day window, 120-day gap" =
      list(index = set_index(x, "recurrence"), rules = rules),
    "Challenge definition: progression after initial therapy" =
      list(index = set_index(x, "progression_after_therapy"), rules = rules),
    "Strict: distant metastasis only" =
      list(index = set_index(x, "metastasis", stage_iv = "diagnosis"), rules = rules),
    "Progression confirmed by scan AND oncologist note" =
      list(index = set_index(x, "recurrence", prog_source = "both"), rules = rules),
    "Gap 60 days" = list(index = set_index(x, "recurrence"), rules = r(gap = 60)),
    "Gap 180 days" = list(index = set_index(x, "recurrence"), rules = r(gap = 180)),
    "Window 42 days" = list(index = set_index(x, "recurrence"), rules = r(window = 42)),
    "Trial drugs excluded" = list(index = set_index(x, "recurrence"),
                                  rules = r(ignore_drugs = "Investigational Drug"))
  )
}

#' Cohort attrition table
#'
#' Counts patients at each step from all cancer cases to those with a
#' line of therapy after the index, as reviewers expect in real-world
#' evidence studies.
#'
#' @param x A [lot_data()] object.
#' @param index Index dates from [set_index()].
#' @param rules Rules from [lot_rules()].
#' @param k Mask counts below this; `0` turns it off.
#' @return A data frame of steps and counts.
#' @examples
#' sim <- simulate_bpc(80)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' lot_attrition(x)
#' @export
lot_attrition <- function(x, index = set_index(x), rules = lot_rules(), k = 10) {
  cs <- x$cases
  idx <- as.data.frame(index)
  lots <- calculate_lot(x, index, rules)
  early <- cs$case_id[cs$early_stage %in% TRUE]
  iv <- cs$case_id[cs$early_stage %in% FALSE]
  has_idx <- idx$case_id[!is.na(idx$index_day)]
  treated <- unique(lots$case_id)
  up <- idx$case_id[grepl("upstaged", idx$index_event)]
  steps <- data.frame(
    step = c("All cancer cases",
             "  Stage IV at diagnosis",
             "  Early stage (I-III)",
             "    Upstaged (spread found within the first days)",
             "    Recurred (index found)",
             "    No recurrence found (no lines counted)",
             "  Unknown stage",
             "Cases with an index date",
             "Cases with at least one line after the index",
             "  Early-stage recurrences with a line"),
    patients = c(nrow(cs), length(iv), length(early),
                 length(intersect(early, up)),
                 length(setdiff(intersect(early, has_idx), up)),
                 length(setdiff(early, has_idx)),
                 sum(is.na(cs$early_stage)),
                 length(has_idx), length(treated),
                 length(intersect(setdiff(intersect(early, has_idx), up), treated))),
    stringsAsFactors = FALSE)
  if (k > 0) steps <- suppress_small(steps, cols = "patients", k = k)
  steps
}
