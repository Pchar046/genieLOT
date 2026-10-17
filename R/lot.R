#' Line-of-therapy rules
#'
#' LOT definitions differ by disease and by user. Each rule is a setting here.
#'
#' How a new line starts (checked in order for each drug after the index date):
#' 1. Drugs that start within `window` days of the line's first drug belong
#'    to that line (they make up the regimen).
#' 2. After the window, a drug **not already in the line** starts a new line,
#'    unless it is in `add_without_new_line` or swaps with a line drug
#'    (`swaps`).
#' 3. A break longer than `gap` days after all the line's drugs stopped starts
#'    a new line, even if the same drugs restart.
#' 4. Optional: if `progression_new_line = TRUE`, any drug started after a
#'    progression event starts a new line.
#'
#' Dropping a drug (e.g. carboplatin stopped, pemetrexed maintenance continues)
#' never starts a new line.
#'
#' @param window Days after a line starts during which drugs join it.
#' @param gap Days without treatment after which a restart is a new line.
#'   Use `Inf` to switch off.
#' @param swaps List of drug groups that may replace each other in the same
#'   line, e.g. `list(c("Carboplatin", "Cisplatin"))`.
#' @param add_without_new_line Drugs that can be added later without a new
#'   line (e.g. "Bevacizumab" for some definitions).
#' @param ignore_drugs Drugs to drop entirely (e.g. "Investigational Drug").
#' @param progression_new_line Start a new line after each progression event.
#' @param lookback Count drugs started up to this many days before the index
#'   date (treatment started just before the index event was recorded).
#' @param max_lines Keep at most this many lines.
#' @return A list of class `lot_rules`.
#' @examples
#' lot_rules()
#' lot_rules(window = 42, gap = 90, add_without_new_line = "Bevacizumab")
#' @export
lot_rules <- function(window = 28,
                      gap = 120,
                      swaps = list(c("Carboplatin", "Cisplatin")),
                      add_without_new_line = character(),
                      ignore_drugs = character(),
                      progression_new_line = FALSE,
                      lookback = 0,
                      max_lines = Inf) {
  stopifnot(is.numeric(window), window >= 0, is.numeric(gap), gap >= 0,
            is.list(swaps), is.numeric(lookback), lookback >= 0)
  structure(list(window = window, gap = gap, swaps = swaps,
                 add_without_new_line = add_without_new_line,
                 ignore_drugs = ignore_drugs,
                 progression_new_line = progression_new_line,
                 lookback = lookback, max_lines = max_lines),
            class = "lot_rules")
}

#' @export
print.lot_rules <- function(x, ...) {
  cat("<lot_rules>\n")
  cat(sprintf("  Regimen window:        %s days\n", x$window))
  cat(sprintf("  Gap for new line:      %s days\n", x$gap))
  cat(sprintf("  Swaps (same line):     %s\n",
              if (length(x$swaps)) paste(vapply(x$swaps, paste, "", collapse = "/"), collapse = "; ") else "none"))
  cat(sprintf("  Add-ons (same line):   %s\n",
              if (length(x$add_without_new_line)) paste(x$add_without_new_line, collapse = ", ") else "none"))
  cat(sprintf("  Ignored drugs:         %s\n",
              if (length(x$ignore_drugs)) paste(x$ignore_drugs, collapse = ", ") else "none"))
  cat(sprintf("  New line on progression: %s\n", x$progression_new_line))
  cat(sprintf("  Lookback before index: %s days\n", x$lookback))
  invisible(x)
}

# TRUE where a drug name starts with any of `names` (case-insensitive), so
# "Goserelin" matches "Goserelin Acetate".
drug_in <- function(drug, names) {
  names <- tolower(trimws(names))
  names <- names[nzchar(names)]
  d <- tolower(trimws(drug))
  if (!length(names) || !length(d)) return(rep(FALSE, length(d)))
  Reduce(`|`, lapply(names, function(n) startsWith(d, n)))
}

standardize_drugs <- function(drug, swaps) {
  out <- tolower(trimws(drug))
  for (g in swaps) {
    out[drug_in(drug, g)] <- paste(sort(tolower(trimws(g))), collapse = "/")
  }
  out
}

#' Calculate lines of therapy
#'
#' @param x A [lot_data()] object.
#' @param index Index dates from [set_index()].
#' @param rules Rules from [lot_rules()].
#' @return A data frame of class `lot_result`, one row per patient line:
#'   `case_id`, `line`, `regimen`, `start`, `end`, `days_from_index`,
#'   `duration`, `time_to_next_line`. The drug-level assignment is kept in
#'   `attr(result, "drugs")`.
#' @examples
#' sim <- simulate_bpc(50)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' lots <- calculate_lot(x, set_index(x, "recurrence"), lot_rules())
#' head(lots)
#' @export
calculate_lot <- function(x, index = set_index(x), rules = lot_rules()) {
  stopifnot(inherits(x, "lot_data"), inherits(rules, "lot_rules"))
  idx <- as.data.frame(index)[, c("case_id", "index_day", "index_event")]
  d <- merge(x$drugs, idx, by = "case_id")
  d <- d[!is.na(d$index_day) & d$start >= d$index_day - rules$lookback, , drop = FALSE]
  d <- d[!drug_in(d$drug, rules$ignore_drugs), , drop = FALSE]
  d$drug_std <- standardize_drugs(d$drug, rules$swaps)
  d$is_addon <- drug_in(d$drug, rules$add_without_new_line)

  empty_lines <- data.frame(case_id = character(), line = integer(), regimen = character(),
                            start = numeric(), end = numeric(), index_day = numeric(),
                            index_event = character(), days_from_index = numeric(),
                            duration = numeric(), time_to_next_line = numeric(),
                            stringsAsFactors = FALSE)
  if (!nrow(d)) {
    return(structure(empty_lines, class = c("lot_result", "data.frame"),
                     drugs = d, index = index, rules = rules))
  }

  progs <- split(x$prog$day, x$prog$case_id)
  pieces <- lapply(split(d, d$case_id), function(dd) {
    assign_lines(dd, rules, progs[[dd$case_id[1]]])
  })
  d <- do.call(rbind, pieces)
  d <- d[d$line <= rules$max_lines, , drop = FALSE]
  rownames(d) <- NULL

  lines <- do.call(rbind, lapply(split(d, list(d$case_id, d$line), drop = TRUE), function(z) {
    data.frame(case_id = z$case_id[1], line = z$line[1],
               regimen = paste(sort(unique(z$drug)), collapse = " + "),
               start = min(z$start), end = max(z$end),
               index_day = z$index_day[1], index_event = z$index_event[1],
               stringsAsFactors = FALSE)
  }))
  lines <- lines[order(lines$case_id, lines$line), , drop = FALSE]
  lines$days_from_index <- lines$start - lines$index_day
  lines$duration <- lines$end - lines$start
  nxt <- stats::ave(lines$start, lines$case_id, FUN = function(s) c(s[-1], NA))
  lines$time_to_next_line <- nxt - lines$start
  rownames(lines) <- NULL

  structure(lines, class = c("lot_result", "data.frame"),
            drugs = d, index = index, rules = rules)
}

assign_lines <- function(dd, rules, prog_days) {
  dd <- dd[order(dd$start, dd$drug_std), , drop = FALSE]
  if (is.null(prog_days)) prog_days <- numeric(0)
  line <- integer(nrow(dd))
  cur <- 0L
  line_start <- NA_real_
  win_end <- -Inf
  line_end <- -Inf
  line_drugs <- character(0)
  for (i in seq_len(nrow(dd))) {
    s <- dd$start[i]
    g <- dd$drug_std[i]
    if (cur == 0L) {
      new <- TRUE
    } else if (s <= win_end) {
      new <- FALSE
    } else {
      new_drug <- !(g %in% line_drugs) && !isTRUE(dd$is_addon[i])
      gap_break <- (s - line_end) > rules$gap
      prog_break <- isTRUE(rules$progression_new_line) &&
        any(prog_days > win_end & prog_days <= s)
      new <- new_drug || gap_break || prog_break
    }
    if (new) {
      cur <- cur + 1L
      line_start <- s
      win_end <- s + rules$window
      line_drugs <- g
      line_end <- dd$end[i]
    } else {
      line_drugs <- union(line_drugs, g)
      line_end <- max(line_end, dd$end[i])
    }
    line[i] <- cur
  }
  dd$line <- line
  dd
}

#' @export
print.lot_result <- function(x, n = 10, ...) {
  n_cases <- length(unique(x$case_id))
  cat(sprintf("<lot_result> %d lines for %d patients\n", nrow(x), n_cases))
  if (nrow(x)) {
    tab <- table(factor(x$line, levels = sort(unique(x$line))))
    cat("  Patients reaching each line: ",
        paste(sprintf("L%s=%d", names(tab), as.integer(tab)), collapse = "  "), "\n\n")
    print(utils::head(as.data.frame(x)[, c("case_id", "line", "regimen", "start",
                                          "end", "days_from_index")], n), row.names = FALSE)
  }
  invisible(x)
}

#' Summarize lines of therapy
#'
#' @param lots Result of [calculate_lot()].
#' @param top_n Number of most common regimens to list per line.
#' @return A data frame with one row per line number.
#' @examples
#' sim <- simulate_bpc(80)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' summarize_lot(calculate_lot(x))
#' @export
summarize_lot <- function(lots, top_n = 3) {
  if (!nrow(lots)) return(data.frame())
  idx <- attr(lots, "index")
  n_index <- if (!is.null(idx)) sum(!is.na(idx$index_day)) else length(unique(lots$case_id))
  out <- do.call(rbind, lapply(split(lots, lots$line), function(z) {
    tr <- sort(table(z$regimen), decreasing = TRUE)
    tr <- utils::head(tr, top_n)
    data.frame(line = z$line[1],
               patients = nrow(z),
               pct_of_indexed = round(100 * nrow(z) / n_index, 1),
               median_days_on_line = stats::median(z$duration),
               median_days_to_next_line = suppressWarnings(stats::median(z$time_to_next_line, na.rm = TRUE)),
               top_regimens = paste(sprintf("%s (%d)", names(tr), as.integer(tr)), collapse = "; "),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
