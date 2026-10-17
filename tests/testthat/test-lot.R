make_data <- function(regs, stage = "Stage II", dmets_yrs = NA, img = NULL) {
  ca <- data.frame(record_id = "P1", ca_seq = 0, cohort = "NSCLC", stage_dx_iv = stage,
                   dx_to_dmets_yrs = dmets_yrs, tt_os_dx_yrs = 5)
  regimens <- data.frame(record_id = "P1", ca_seq = 0,
                         regimen_number = seq_len(nrow(regs)),
                         regimen_drugs = regs$drugs,
                         dx_reg_start_int = regs$start,
                         dx_reg_end_all_int = regs$end)
  lot_data(ca, regimens, imaging = img)
}

test_that("drugs starting within the window form one line", {
  x <- make_data(data.frame(drugs = c("Carboplatin, Pemetrexed Disodium", "Pembrolizumab"),
                            start = c(100, 110), end = c(184, 300)))
  l <- calculate_lot(x, set_index(x, "diagnosis"))
  expect_equal(nrow(l), 1)
  expect_equal(l$regimen, "Carboplatin + Pembrolizumab + Pemetrexed Disodium")
})

test_that("dropping a drug (maintenance) does not start a new line", {
  x <- make_data(data.frame(drugs = c("Carboplatin, Pemetrexed Disodium", "Pemetrexed Disodium"),
                            start = c(100, 185), end = c(184, 400)))
  expect_equal(max(calculate_lot(x, set_index(x, "diagnosis"))$line), 1)
})

test_that("adding a new drug after the window starts a new line", {
  x <- make_data(data.frame(drugs = c("Osimertinib", "Docetaxel"),
                            start = c(100, 400), end = c(390, 500)))
  expect_equal(calculate_lot(x, set_index(x, "diagnosis"))$line, c(1, 2))
})

test_that("add-on drugs can be exempted", {
  x <- make_data(data.frame(drugs = c("Carboplatin, Paclitaxel", "Bevacizumab"),
                            start = c(100, 160), end = c(190, 300)))
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"))), 2)
  r <- lot_rules(add_without_new_line = "Bevacizumab")
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"), r)), 1)
})

test_that("swapping cisplatin for carboplatin stays in the same line", {
  x <- make_data(data.frame(drugs = c("Cisplatin, Pemetrexed Disodium", "Carboplatin"),
                            start = c(100, 160), end = c(150, 200)))
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"))), 1)
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"), lot_rules(swaps = list()))), 2)
})

test_that("a long gap restarts as a new line, a short one does not", {
  regs <- data.frame(drugs = c("Osimertinib", "Osimertinib"),
                     start = c(100, 0), end = c(200, 0))
  regs$start[2] <- 260; regs$end[2] <- 400
  x <- make_data(regs)
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"), lot_rules(gap = 120))), 1)
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"), lot_rules(gap = 30))), 2)
})

test_that("early-stage adjuvant therapy is excluded with a recurrence index", {
  img <- data.frame(record_id = "P1", dx_scan_days = c(60, 150, 900),
                    image_ca = c("Stable/No change", "Stable/No change",
                                 "Progressing/Worsening/Enlarging"))
  x <- make_data(data.frame(drugs = c("Cisplatin, Vinorelbine Tartrate",
                                      "Carboplatin, Pemetrexed Disodium, Pembrolizumab",
                                      "Docetaxel"),
                            start = c(45, 920, 1230), end = c(129, 1200, 1350)),
                 img = img)
  idx <- set_index(x, "recurrence")
  expect_equal(idx$index_day, 900)
  l <- calculate_lot(x, idx)
  expect_equal(l$regimen[1], "Carboplatin + Pembrolizumab + Pemetrexed Disodium")
  expect_equal(nrow(l), 2)
  naive <- calculate_lot(x, set_index(x, "diagnosis"))
  expect_equal(naive$regimen[1], "Cisplatin + Vinorelbine Tartrate")
  cmp <- compare_lot(x)
  expect_true(cmp$cases$mislabeled)
  expect_equal(cmp$cases$naive_label_of_true_line1, 2)
})

test_that("early-stage patients who never recur get no lines", {
  x <- make_data(data.frame(drugs = "Cisplatin, Vinorelbine Tartrate", start = 45, end = 129))
  expect_true(is.na(set_index(x, "recurrence")$index_day))
  expect_equal(nrow(calculate_lot(x, set_index(x, "recurrence"))), 0)
})

test_that("stage IV patients index at diagnosis", {
  x <- make_data(data.frame(drugs = "Osimertinib", start = 20, end = 400),
                 stage = "Stage IV", dmets_yrs = 0)
  expect_equal(set_index(x, "recurrence")$index_day, 0)
})

test_that("distant metastasis date can trigger recurrence", {
  x <- make_data(data.frame(drugs = c("Cisplatin", "Osimertinib"),
                            start = c(45, 760), end = c(129, 900)),
                 dmets_yrs = 2)
  expect_equal(round(set_index(x, "recurrence")$index_day), 730)
})

test_that("custom index dates work", {
  x <- make_data(data.frame(drugs = c("Osimertinib", "Docetaxel"),
                            start = c(100, 400), end = c(390, 500)))
  idx <- set_index(x, "custom", custom = data.frame(record_id = "P1", index_day = 395))
  l <- calculate_lot(x, idx)
  expect_equal(l$regimen, "Docetaxel")
})

test_that("column names can be remapped", {
  ca <- data.frame(record_id = "P1", stage_dx = "Stage IV")
  reg <- data.frame(record_id = "P1", regimen_drugs = "Osimertinib",
                    dx_reg_start_int = 10, dx_reg_end_all_int = 50)
  expect_error(lot_data(ca, reg), "stage_dx_iv")
  x <- lot_data(ca, reg, cols = bpc_columns(stage = "stage_dx"))
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"))), 1)
})

test_that("simulated data runs end to end", {
  sim <- simulate_bpc(80)
  x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
  cmp <- compare_lot(x)
  expect_gt(cmp$summary$pct_mislabeled[2], 0)
  expect_s3_class(plot_patient(x, cmp$correct, "SIM-0001", naive = cmp$naive), "ggplot")
})

test_that("early-stage patients with spread found soon after diagnosis are upstaged", {
  x <- make_data(data.frame(drugs = "Osimertinib", start = 30, end = 300),
                 dmets_yrs = 20 / 365.25)
  idx <- set_index(x, "recurrence")
  expect_equal(idx$index_day, 0)
  expect_equal(idx$index_event, "diagnosis (upstaged)")
  expect_true(is.na(set_index(x, "recurrence", upstaged = FALSE)$index_day))
  expect_equal(nrow(calculate_lot(x, idx)), 1)
})

test_that("validate_index counts agreement with a flag", {
  sim <- simulate_bpc(80)
  sim$ca_dx$dmets_stage_i_iii <- as.integer(!is.na(sim$ca_dx$dx_to_dmets_yrs) &
                                              sim$ca_dx$stage_dx_iv != "Stage IV")
  x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
  val <- validate_index(x, set_index(x), sim$ca_dx)
  expect_match(val$summary$value[1], "100.0%")
})

test_that("drug names match by prefix, ignoring case", {
  expect_equal(genieLOT:::drug_in(c("Goserelin Acetate", "Letrozole", "goserelin"), "Goserelin"),
               c(TRUE, FALSE, TRUE))
})

test_that("CRC preset keeps bevacizumab add-on and 5-FU/capecitabine swap in one line", {
  x <- make_data(data.frame(drugs = c("Fluorouracil, Oxaliplatin, Leucovorin Calcium",
                                      "Bevacizumab", "Capecitabine"),
                            start = c(100, 160, 200), end = c(190, 300, 350)),
                 stage = "Stage IV", dmets_yrs = 0)
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"), lot_preset("CRC"))), 1)
  expect_equal(nrow(calculate_lot(x, set_index(x, "diagnosis"), lot_rules(swaps = list()))), 3)
})

test_that("small counts are masked", {
  df <- suppress_small(data.frame(line = 1:3, patients = c(120, 40, 6),
                                  top = c("A (50)", "B (4)", "C (12); D (3)")))
  expect_equal(df$patients, c("120", "40", "<10"))
  expect_equal(df$top, c("A (50)", "B (<10)", "C (12); D (<10)"))
})

test_that("the app server computes every output", {
  skip_if_not_installed("shiny")
  app <- run_lot_app()
  shiny::testServer(app, {
    session$setInputs(preset = "NSCLC", event = "recurrence", src = "either", min_days = 90,
                      window = 28, gap = 120, swaps = "Carboplatin/Cisplatin", addons = "",
                      trial = TRUE, prog_new = FALSE, privacy = TRUE, pid = "SIM-0001_0")
    expect_true(cmp()$summary$pct_mislabeled[2] > 0)
    expect_match(output$cards$html, "AACR-flagged")
    expect_match(val()$summary$value[1], "100.0%")
    session$setInputs(trial = FALSE, privacy = FALSE)
    expect_true(nrow(smry()) > 0)
    expect_match(output$attrition, "All cancer cases")
    session$setInputs(run_sens = 1)
    expect_equal(nrow(sens()), 8)
  })
})

test_that("missing end dates are filled from the last dose", {
  ca <- data.frame(record_id = "P1", stage_dx_iv = "Stage IV")
  reg <- data.frame(record_id = "P1", regimen_drugs = "Carboplatin, Pemetrexed Disodium",
                    dx_reg_start_int = 10, dx_reg_end_all_int = NA,
                    dx_drug_end_or_lastadm_int_1 = 90, dx_drug_end_or_lastadm_int_2 = 150,
                    dx_reg_end_any_int = 90)
  x <- lot_data(ca, reg)
  expect_equal(unique(x$drugs$end), 150)
  expect_equal(unique(x$drugs$end_source), "last dose")
  reg$dx_drug_end_or_lastadm_int_1 <- NA; reg$dx_drug_end_or_lastadm_int_2 <- NA
  expect_equal(unique(lot_data(ca, reg)$drugs$end_source), "first drug end")
})

test_that("progression after initial therapy follows the challenge definition", {
  img <- data.frame(record_id = "P1", dx_scan_days = c(20, 400, 900),
                    image_ca = "Progressing/Worsening/Enlarging")
  x <- make_data(data.frame(drugs = c("Cisplatin, Vinorelbine Tartrate", "Osimertinib"),
                            start = c(300, 920), end = c(380, 1100)), img = img)
  # surgery at day 30: progression on day 20 is before any therapy, so skipped
  surg <- data.frame(record_id = "P1", day = 30, type = "surgery")
  x2 <- lot_data(data.frame(record_id = "P1", ca_seq = 0, stage_dx_iv = "Stage II",
                            dx_to_dmets_yrs = NA, tt_os_dx_yrs = 5),
                 data.frame(record_id = "P1", ca_seq = 0,
                            regimen_drugs = c("Cisplatin, Vinorelbine Tartrate", "Osimertinib"),
                            dx_reg_start_int = c(300, 920), dx_reg_end_all_int = c(380, 1100)),
                 imaging = img, therapy = surg)
  expect_equal(set_index(x2, "initial_therapy")$index_day, 30)
  expect_equal(set_index(x2, "progression_after_therapy")$index_day, 400)
  expect_equal(set_index(x2, "progression_after_therapy", therapy_types = "drug")$index_day, 400)
  # without the surgery record, first therapy is the day-300 drug: same answer
  expect_equal(set_index(x, "progression_after_therapy")$index_day, 400)
  expect_equal(calculate_lot(x2, set_index(x2, "progression_after_therapy"))$regimen, "Osimertinib")
})

test_that("therapy_events filters rows by pattern", {
  path <- data.frame(record_id = c("A", "A", "B"),
                     path_proc = c("Biopsy", "Surgical resection", "Lobectomy"),
                     dx_path_proc_days = c(5, 40, 30))
  ev <- therapy_events(path, "dx_path_proc_days", "surgery", "path_proc", "resection|ectomy")
  expect_equal(ev$day, c(40, 30))
})

test_that("cohort_summary gives one comparable row", {
  sim <- simulate_bpc(80)
  x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
  row <- cohort_summary(x, set_index(x), lot_preset("NSCLC"), sim$ca_dx)
  expect_equal(nrow(row), 1)
  expect_equal(row$flagged_recurrences_handled, "100.0%")
})

test_that("sensitivity and attrition tables run", {
  sim <- simulate_bpc(80)
  x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
  s <- lot_sensitivity(x)
  expect_equal(nrow(s), 8)
  a <- lot_attrition(x, k = 0)
  expect_equal(a$patients[1], 80)
  expect_equal(a$patients[2] + a$patients[3] + a$patients[7], 80)
})

test_that("surgery events come from surgical pathology reports", {
  path <- data.frame(record_id = c("A", "A", "A"),
                     path_proc_type = c("Cytology", "Surgical pathology", "Other"),
                     dx_path_proc_days = c(3, 25, 40))
  ev <- bpc_surgery_events(path)
  expect_equal(ev$day, 25); expect_equal(ev$type, "surgery")
})
