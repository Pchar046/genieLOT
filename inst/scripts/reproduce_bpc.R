# Reproduce the genieLOT results on AACR Project GENIE BPC data
# ----------------------------------------------------------------
# Needs approved GENIE BPC access (public releases are self-serve on Synapse).
# Put these 4 files from one cohort (e.g. NSCLC 2.0-public) in one folder:
#   cancer_level_dataset_index.csv, regimen_cancer_level_dataset.csv,
#   imaging_level_dataset.csv, med_onc_note_level_dataset.csv
#
# Run one cohort:
#   source(system.file("scripts", "reproduce_bpc.R", package = "genieLOT"))
#   reproduce_bpc("C:/genie_data")
# Compare several cancers (one folder per cohort):
#   reproduce_cohorts(c(NSCLC = "C:/genie_data", CRC = "C:/genie_crc", BrCa = "C:/genie_brca"))
#
# Output: summary tables only (counts under 10 masked), saved to ./genieLOT_results.
# No patient-level rows are written.

library(genieLOT)

reproduce_bpc <- function(data_dir, out_dir = "genieLOT_results",
                          preset = "NSCLC", min_count = 10) {
  raw <- read_bpc(data_dir)
  cols <- bpc_columns(img_status = "image_overall",
                      dmets_time = "dx_to_dmets_days", dmets_unit = "days",
                      os_time = "tt_os_dx_days", os_unit = "days")
  x <- lot_data(raw$ca_dx, raw$regimens, raw$imaging, raw$med_onc, cols = cols)
  idx <- set_index(x, "recurrence")
  rules <- lot_preset(preset)

  cmp <- compare_lot(x, idx, rules)
  val <- validate_index(x, idx, raw$ca_dx, cols = cols)
  lines_all <- summarize_lot(cmp$correct)
  lines_no_trial <- summarize_lot(calculate_lot(
    x, idx, lot_rules(window = rules$window, gap = rules$gap, swaps = rules$swaps,
                      add_without_new_line = rules$add_without_new_line,
                      ignore_drugs = "Investigational Drug")))

  dir.create(out_dir, showWarnings = FALSE)
  save <- function(df, name) {
    df <- suppress_small(df, k = min_count)
    utils::write.csv(df, file.path(out_dir, name), row.names = FALSE)
    df
  }
  out <- list(
    cohort = data.frame(cancer_cases = nrow(x$cases),
                        early_stage = sum(x$cases$early_stage %in% TRUE),
                        stage_iv = sum(x$cases$early_stage %in% FALSE)),
    naive_vs_correct = save(cmp$summary, "naive_vs_correct.csv"),
    validation = save(val$summary, "validation.csv"),
    lines = save(lines_all, "lines_of_therapy.csv"),
    lines_without_trial_drugs = save(lines_no_trial, "lines_without_trial_drugs.csv")
  )
  utils::write.csv(out$cohort, file.path(out_dir, "cohort.csv"), row.names = FALSE)
  print(out)
  message("Summary tables saved to ", normalizePath(out_dir))
  invisible(out)
}

bpc_cols_default <- function() {
  bpc_columns(img_status = "image_overall",
              dmets_time = "dx_to_dmets_days", dmets_unit = "days",
              os_time = "tt_os_dx_days", os_unit = "days")
}

reproduce_cohorts <- function(dirs, out_dir = "genieLOT_results", min_count = 10) {
  rows <- lapply(names(dirs), function(nm) {
    message("Running ", nm, " ...")
    raw <- read_bpc(dirs[[nm]])
    cols <- bpc_cols_default()
    x <- lot_data(raw$ca_dx, raw$regimens, raw$imaging, raw$med_onc, cols = cols)
    preset <- if (nm %in% c("NSCLC", "CRC", "BrCa")) nm else "Custom"
    cohort_summary(x, set_index(x, "recurrence"), lot_preset(preset), raw$ca_dx, label = nm)
  })
  tab <- do.call(rbind, rows)
  tab <- suppress_small(tab, cols = c("patients_with_lines", "early_recurrences_treated",
                                      "naive_line1_wrong"), k = min_count)
  dir.create(out_dir, showWarnings = FALSE)
  utils::write.csv(tab, file.path(out_dir, "cohort_comparison.csv"), row.names = FALSE)
  print(tab)
  invisible(tab)
}
