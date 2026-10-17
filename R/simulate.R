#' Simulate fake BPC-style data
#'
#' Makes fake patients with the same table layout and column names as GENIE
#' BPC, so the package can be tried without data access. No real patient data.
#'
#' Patient `SIM-0001` is a worked example: stage II, adjuvant chemo after
#' surgery, recurrence about 2.5 years later, then two lines of therapy.
#'
#' @param n Number of patients.
#' @param cohort Cohort label.
#' @param seed Random seed.
#' @return A list of data frames: `ca_dx`, `regimens`, `imaging`, `med_onc`.
#' @examples
#' sim <- simulate_bpc(20)
#' str(sim, max.level = 1)
#' @export
simulate_bpc <- function(n = 200, cohort = "NSCLC", seed = 42) {
  if (!is.null(seed)) set.seed(seed)
  ca <- vector("list", n)
  reg <- list()
  img <- list()
  md <- list()
  add_reg <- function(id, drugs, s, e) {
    reg[[length(reg) + 1]] <<- data.frame(
      record_id = id, ca_seq = 0, regimen_drugs = paste(drugs, collapse = ", "),
      dx_reg_start_int = round(s), dx_reg_end_all_int = round(e),
      stringsAsFactors = FALSE)
  }
  first_lines <- list(
    induction = c("Carboplatin", "Pemetrexed Disodium", "Pembrolizumab"),
    egfr = c("Osimertinib"),
    doublet = c("Carboplatin", "Paclitaxel")
  )
  later_lines <- list(c("Docetaxel", "Ramucirumab"), c("Docetaxel"),
                      c("Gemcitabine Hydrochloride"), c("Nivolumab"))

  for (i in seq_len(n)) {
    id <- sprintf("SIM-%04d", i)
    demo <- i == 1
    early <- if (demo) TRUE else stats::runif(1) < 0.6
    stage <- if (demo) "Stage II" else if (early) sample(c("Stage I", "Stage II", "Stage III"), 1) else "Stage IV"
    fu <- if (demo) 1600 else stats::runif(1, 700, 2600)
    progs <- numeric(0)
    dmets <- NA_real_
    recur <- FALSE

    if (early) {
      if (demo || stats::runif(1) < 0.7) {
        s <- if (demo) 45 else stats::runif(1, 35, 75)
        add_reg(id, c(if (demo || stats::runif(1) < 0.7) "Cisplatin" else "Carboplatin", "Vinorelbine Tartrate"), s, s + 84)
        if (stage == "Stage III" && stats::runif(1) < 0.5) add_reg(id, "Durvalumab", s + 110, s + 110 + 330)
      }
      recur <- demo || stats::runif(1) < 0.55
      if (recur) {
        t0 <- if (demo) 900 else stats::runif(1, 300, min(1400, fu - 200))
        dmets <- t0
        progs <- t0
      }
    } else {
      t0 <- 0
      dmets <- 0
    }

    if (!early || recur) {
      nlines <- if (demo) 2 else sample(1:3, 1, prob = c(0.4, 0.35, 0.25))
      cur <- t0 + (if (demo) 20 else stats::runif(1, 10, 35))
      l1 <- if (demo) "induction" else sample(names(first_lines), 1)
      if (l1 == "induction") {
        add_reg(id, first_lines$induction, cur, cur + 84)
        dur <- if (demo) 280 else stats::runif(1, 150, 420)
        add_reg(id, c("Pemetrexed Disodium", "Pembrolizumab"), cur + 85, cur + dur)
      } else {
        dur <- if (demo) 280 else stats::runif(1, 120, 500)
        add_reg(id, first_lines[[l1]], cur, cur + dur)
      }
      end <- cur + dur
      for (k in seq_len(nlines - 1)) {
        p <- end + (if (demo) 10 else stats::runif(1, 0, 20))
        if (p > fu) break
        progs <- c(progs, p)
        cur <- p + (if (demo) 20 else stats::runif(1, 10, 30))
        dur <- if (demo) 120 else stats::runif(1, 60, 240)
        drugs <- if (demo) "Docetaxel" else later_lines[[sample(length(later_lines), 1)]]
        add_reg(id, drugs, cur, cur + dur)
        end <- cur + dur
      }
      # occasional treatment break and restart of the same drug (tests gap rule)
      if (!demo && stats::runif(1) < 0.1) add_reg(id, utils::tail(strsplit(reg[[length(reg)]]$regimen_drugs, ", ")[[1]], 1), end + 200, end + 300)
    }

    scans <- seq(30, fu, by = 90)
    st <- "Stable/No change"
    img[[i]] <- rbind(
      data.frame(record_id = id, dx_scan_days = round(scans), image_ca = st, stringsAsFactors = FALSE),
      if (length(progs)) data.frame(record_id = id, dx_scan_days = round(progs),
                                    image_ca = "Progressing/Worsening/Enlarging", stringsAsFactors = FALSE)
    )
    visits <- seq(20, fu, by = 60)
    md[[i]] <- rbind(
      data.frame(record_id = id, dx_md_visit_days = round(visits), md_ca_status = "Stable/No change",
                 stringsAsFactors = FALSE),
      if (length(progs)) data.frame(record_id = id, dx_md_visit_days = round(progs + 7),
                                    md_ca_status = "Progressing/Worsening/Enlarging", stringsAsFactors = FALSE)
    )
    ca[[i]] <- data.frame(record_id = id, ca_seq = 0, cohort = cohort, stage_dx_iv = stage,
                          dx_to_dmets_yrs = round(dmets / 365.25, 3),
                          dmets_stage_i_iii = as.integer(early && recur),
                          tt_os_dx_yrs = round(fu / 365.25, 3), stringsAsFactors = FALSE)
  }
  regimens <- do.call(rbind, reg)
  regimens <- regimens[order(regimens$record_id, regimens$dx_reg_start_int), ]
  regimens$regimen_number <- stats::ave(seq_len(nrow(regimens)), regimens$record_id, FUN = seq_along)
  imaging <- do.call(rbind, img)
  med_onc <- do.call(rbind, md)
  rownames(regimens) <- rownames(imaging) <- rownames(med_onc) <- NULL
  list(ca_dx = do.call(rbind, ca), regimens = regimens,
       imaging = imaging[order(imaging$record_id, imaging$dx_scan_days), ],
       med_onc = med_onc[order(med_onc$record_id, med_onc$dx_md_visit_days), ])
}
