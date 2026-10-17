#' Hide small counts in a summary table
#'
#' Replaces counts below `k` with "<k" so small groups can't point to
#' individual patients. Also masks "(n)" counts inside text such as
#' "Osimertinib (4)". Use before sharing or showing results.
#'
#' @param df A data frame (e.g. from [summarize_lot()] or [compare_lot()]).
#' @param cols Count columns to check. Default: integer-valued numeric columns
#'   named like patients, mislabeled, n or count.
#' @param k Smallest count that may be shown. Default 10.
#' @return The data frame with small counts masked (count columns become text).
#' @examples
#' suppress_small(data.frame(line = 1:3, patients = c(120, 40, 6)))
#' @export
suppress_small <- function(df, cols = NULL, k = 10) {
  if (!nrow(df)) return(df)
  if (is.null(cols)) {
    cols <- grep("^(patients|mislabeled|n|count|pre_index_lines_excluded)$", names(df), value = TRUE)
  }
  for (cl in cols) {
    v <- df[[cl]]
    df[[cl]] <- ifelse(!is.na(v) & v > 0 & v < k, paste0("<", k), as.character(v))
  }
  for (cl in names(df)[vapply(df, is.character, logical(1))]) {
    df[[cl]] <- mask_counts_in_text(df[[cl]], k)
  }
  df
}

mask_counts_in_text <- function(x, k) {
  vapply(x, function(s) {
    if (is.na(s)) return(s)
    m <- gregexpr("\\((\\d+)\\)", s)
    hits <- regmatches(s, m)[[1]]
    if (!length(hits)) return(s)
    nums <- as.integer(gsub("\\D", "", hits))
    repl <- ifelse(nums < k, paste0("(<", k, ")"), hits)
    regmatches(s, m) <- list(repl)
    s
  }, character(1), USE.NAMES = FALSE)
}
