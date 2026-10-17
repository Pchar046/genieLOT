#' Choose the index date for each patient
#'
#' The index date is where line counting starts. Treatments before it are
#' ignored. This is how early-stage patients' adjuvant or neoadjuvant
#' treatment is kept out of the line count.
#'
#' Events:
#' * `"recurrence"` (default): stage I-III patients start at the first sign the
#'   cancer came back (distant metastasis or a progression event on or after
#'   `min_days`). Stage IV patients start at diagnosis (see `stage_iv`).
#'   Early-stage patients who never recur get no index and no lines.
#' * `"diagnosis"`: day 0 for everyone (the naive approach).
#' * `"metastasis"`: date of first distant metastasis.
#' * `"progression"`: first progression event on or after `min_days`.
#' * `"custom"`: your own dates, given in `custom`.
#'
#' @param x A [lot_data()] object.
#' @param event Index event, see Details.
#' @param prog_source Which records count as progression: `"either"` source,
#'   `"imaging"` only, `"med_onc"` only, or `"both"` (imaging and an oncologist
#'   note agreeing within `confirm_days`).
#' @param confirm_days Window for `prog_source = "both"`.
#' @param min_days Ignore progression events earlier than this many days after
#'   diagnosis. Default 90 for `"recurrence"` (skips the initial work-up and
#'   treatment period), 0 otherwise.
#' @param stage_iv What to do with stage IV patients: `"diagnosis"` (index at
#'   day 0) or `"same"` (use the same event as everyone else). Default
#'   `"diagnosis"` for `"recurrence"`, `"same"` otherwise.
#' @param upstaged Treat early-stage patients whose distant metastasis was
#'   found within `min_days` of diagnosis as stage IV (index at diagnosis).
#'   Their cancer had already spread when found; it is not a recurrence.
#'   Applies when `stage_iv = "diagnosis"`. Default `TRUE`.
#' @param custom For `event = "custom"`: a data frame with `record_id` (or
#'   `case_id`) and `index_day` (days since diagnosis).
#' @return A data frame of class `lot_index` with `case_id`, `index_day`,
#'   `index_event`.
#' @examples
#' sim <- simulate_bpc(50)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' head(set_index(x, "recurrence"))
#' @export
set_index <- function(x,
                      event = c("recurrence", "diagnosis", "metastasis",
                                "progression", "custom"),
                      prog_source = c("either", "imaging", "med_onc", "both"),
                      confirm_days = 30,
                      min_days = NULL,
                      stage_iv = NULL,
                      upstaged = TRUE,
                      custom = NULL) {
  stopifnot(inherits(x, "lot_data"))
  event <- match.arg(event)
  prog_source <- match.arg(prog_source)
  if (is.null(min_days)) min_days <- if (event == "recurrence") 90 else 0
  if (is.null(stage_iv)) stage_iv <- if (event == "recurrence") "diagnosis" else "same"
  stage_iv <- match.arg(stage_iv, c("diagnosis", "same"))
  cs <- x$cases

  first_prog <- function() {
    p <- progression_days(x$prog, prog_source, confirm_days)
    p <- p[p$day >= min_days, , drop = FALSE]
    fp <- if (nrow(p)) tapply(p$day, p$case_id, min) else numeric(0)
    unname(fp[cs$case_id])
  }
  dmets <- ifelse(!is.na(cs$dmets_day) & cs$dmets_day >= 0, cs$dmets_day, NA)

  day <- switch(event,
    diagnosis = rep(0, nrow(cs)),
    metastasis = dmets,
    progression = first_prog(),
    recurrence = {
      fp <- first_prog()
      dm <- ifelse(!is.na(dmets) & dmets >= min_days, dmets, NA)
      suppressWarnings(pmin(fp, dm, na.rm = TRUE))
    },
    custom = {
      if (is.null(custom) || !"index_day" %in% names(custom)) {
        stop("For event = 'custom', give `custom` with columns record_id (or case_id) and index_day.",
             call. = FALSE)
      }
      key <- if ("case_id" %in% names(custom)) cs$case_id else cs$record_id
      kc <- if ("case_id" %in% names(custom)) custom$case_id else custom$record_id
      as.numeric(custom$index_day)[match(key, as.character(kc))]
    }
  )
  day[is.infinite(day)] <- NA
  label <- rep(event, nrow(cs))
  if (stage_iv == "diagnosis" && event != "diagnosis") {
    iv <- cs$early_stage %in% FALSE
    day[iv] <- 0
    label[iv] <- "diagnosis (stage IV)"
    if (isTRUE(upstaged)) {
      up <- cs$early_stage %in% TRUE & !is.na(dmets) & dmets < min_days
      day[up] <- 0
      label[up] <- "diagnosis (upstaged)"
    }
  }
  structure(
    data.frame(case_id = cs$case_id, index_day = day, index_event = label,
               stringsAsFactors = FALSE),
    class = c("lot_index", "data.frame"),
    settings = list(event = event, prog_source = prog_source,
                    confirm_days = confirm_days, min_days = min_days,
                    stage_iv = stage_iv, upstaged = upstaged)
  )
}

progression_days <- function(prog, source, confirm_days) {
  if (source %in% c("imaging", "med_onc")) {
    return(prog[prog$source == source, c("case_id", "day")])
  }
  if (source == "either") return(prog[, c("case_id", "day")])
  img <- prog[prog$source == "imaging", ]
  md <- prog[prog$source == "med_onc", ]
  if (!nrow(img) || !nrow(md)) return(data.frame(case_id = character(), day = numeric()))
  m <- merge(img, md, by = "case_id", suffixes = c("_img", "_md"))
  m <- m[abs(m$day_img - m$day_md) <= confirm_days, ]
  data.frame(case_id = m$case_id, day = pmin(m$day_img, m$day_md))
}
