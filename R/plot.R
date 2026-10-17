#' Plot one patient's treatment timeline
#'
#' Shows every drug as a bar, coloured by line of therapy. Treatments before
#' the index date are grey ("excluded"). Progression events are marked.
#'
#' @param x A [lot_data()] object.
#' @param lots Result of [calculate_lot()].
#' @param case_id Patient case to plot (`record_id` also works).
#' @param naive Optional naive result (from `calculate_lot(x, set_index(x,
#'   "diagnosis"))`) to show a side-by-side "naive vs correct" view.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_bpc(30)
#' x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
#' lots <- calculate_lot(x)
#' plot_patient(x, lots, "SIM-0001")
#' @export
plot_patient <- function(x, lots, case_id, naive = NULL) {
  id <- case_id
  if (!id %in% x$cases$case_id) id <- x$cases$case_id[x$cases$record_id == case_id][1]
  if (is.na(id)) stop("Patient not found: ", case_id, call. = FALSE)

  build <- function(res, label) {
    dd <- x$drugs[x$drugs$case_id == id, ]
    asg <- attr(res, "drugs")
    asg <- asg[asg$case_id == id, c("drug", "start", "line")]
    dd <- merge(dd, asg, by = c("drug", "start"), all.x = TRUE)
    idx <- as.data.frame(attr(res, "index"))
    dd$line_label <- ifelse(is.na(dd$line), "Before index (excluded)", paste("Line", dd$line))
    dd$panel <- label
    dd$index_day <- idx$index_day[idx$case_id == id]
    dd
  }
  dd <- build(lots, "Your definition")
  if (!is.null(naive)) dd <- rbind(build(naive, "Naive (from diagnosis)"), dd)
  dd$panel <- factor(dd$panel, levels = unique(dd$panel))
  dd$drug <- factor(dd$drug, levels = rev(unique(dd$drug[order(dd$start)])))
  lv <- c("Before index (excluded)", paste("Line", sort(unique(stats::na.omit(dd$line)))))
  dd$line_label <- factor(dd$line_label, levels = lv)
  pal <- c("grey75", grDevices::hcl.colors(max(1, length(lv) - 1), "Dark 3"))
  names(pal) <- lv

  pr <- x$prog[x$prog$case_id == id, ]
  pr <- unique(pr[, c("day", "source")])
  ix <- unique(dd[, c("panel", "index_day")])
  cs <- x$cases[x$cases$case_id == id, ]

  p <- ggplot2::ggplot(dd) +
    ggplot2::geom_segment(ggplot2::aes(x = start, xend = pmax(end, start + 7),
                                       y = drug, yend = drug, colour = line_label),
                          linewidth = 6, lineend = "butt") +
    ggplot2::geom_vline(data = ix[!is.na(ix$index_day), ],
                        ggplot2::aes(xintercept = index_day),
                        linetype = "dashed", colour = "black") +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
    ggplot2::labs(x = "Days since diagnosis", y = NULL,
                  title = sprintf("%s  |  %s", cs$record_id, cs$stage),
                  subtitle = "Dashed line = index date. Triangles = progression events.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  if (nrow(pr)) {
    p <- p + ggplot2::geom_point(data = pr, ggplot2::aes(x = day, y = 0.4), shape = 17,
                                 size = 3, colour = "#B22222", inherit.aes = FALSE)
  }
  if (!is.null(naive)) p <- p + ggplot2::facet_wrap(~panel, ncol = 1)
  p
}

#' Plot how many patients reach each line
#'
#' @param lots Result of [calculate_lot()].
#' @return A ggplot object.
#' @export
plot_lines <- function(lots) {
  s <- summarize_lot(lots)
  s$label <- paste0("Line ", s$line)
  s$label <- factor(s$label, levels = s$label)
  ggplot2::ggplot(s, ggplot2::aes(x = label, y = patients)) +
    ggplot2::geom_col(fill = "#2E6F9E", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(patients, " (", pct_of_indexed, "%)")),
                       vjust = -0.4, size = 4) +
    ggplot2::labs(x = NULL, y = "Patients", title = "Patients reaching each line of therapy") +
    ggplot2::theme_minimal(base_size = 12)
}

utils::globalVariables(c("start", "end", "drug", "line_label", "index_day", "day",
                         "label", "patients", "pct_of_indexed", "panel"))
