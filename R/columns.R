#' Column mapping for GENIE BPC files
#'
#' BPC variable names can differ slightly between cohorts and releases.
#' Every column the package reads is listed here, so a mismatch is fixed in
#' one place: `bpc_columns(stage = "stage_dx")`.
#'
#' Time columns are converted to days since diagnosis of the index cancer.
#' Set the matching `*_unit` to "days", "months" or "years".
#'
#' @param ... Named overrides of the defaults.
#' @return A list of class `bpc_columns`.
#' @examples
#' bpc_columns()
#' bpc_columns(stage = "stage_dx", dmets_time = "dx_to_dmets_mos",
#'             dmets_unit = "months")
#' @export
bpc_columns <- function(...) {
  cols <- list(
    id = "record_id",
    ca_seq = "ca_seq",
    cohort = "cohort",
    stage = "stage_dx_iv",
    dmets_time = "dx_to_dmets_yrs", dmets_unit = "years",
    os_time = "tt_os_dx_yrs", os_unit = "years",
    reg_number = "regimen_number",
    reg_drugs = "regimen_drugs", reg_sep = ",",
    reg_start = "dx_reg_start_int", reg_end = "dx_reg_end_all_int",
    reg_unit = "days",
    img_time = "dx_scan_days", img_status = "image_ca", img_unit = "days",
    md_time = "dx_md_visit_days", md_status = "md_ca_status", md_unit = "days",
    progression_pattern = "progress"
  )
  over <- list(...)
  bad <- setdiff(names(over), names(cols))
  if (length(bad)) {
    stop("Unknown column setting(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  cols[names(over)] <- over
  structure(cols, class = "bpc_columns")
}

#' @export
print.bpc_columns <- function(x, ...) {
  cat("<bpc_columns>\n")
  for (n in names(x)) cat(sprintf("  %-20s %s\n", n, x[[n]]))
  invisible(x)
}

to_days <- function(x, unit) {
  x <- suppressWarnings(as.numeric(x))
  switch(unit,
         days = x,
         months = x * 30.4375,
         years = x * 365.25,
         stop("unit must be 'days', 'months' or 'years'", call. = FALSE))
}

get_col <- function(df, col, table, required = TRUE) {
  if (is.null(col) || !col %in% names(df)) {
    if (required) {
      stop(sprintf(paste0("Column '%s' not found in the %s table. ",
                          "Set the right name with bpc_columns()."),
                   col, table), call. = FALSE)
    }
    return(rep(NA, nrow(df)))
  }
  df[[col]]
}

#' Is a stage value early stage (I-III)?
#'
#' @param stage Character vector of stage labels, e.g. "Stage II".
#' @return Logical vector. `NA` when the stage is unknown.
#' @examples
#' is_early_stage(c("Stage I", "Stage IIIB", "Stage IV", NA))
#' @export
is_early_stage <- function(stage) {
  s <- toupper(trimws(as.character(stage)))
  out <- ifelse(grepl("IV|\\b4\\b", s), FALSE,
                ifelse(grepl("STAGE\\s*(I|II|III|1|2|3)|^(I|II|III)[ABC]?$|^[123][ABC]?$", s),
                       TRUE, NA))
  out[is.na(stage)] <- NA
  out
}
