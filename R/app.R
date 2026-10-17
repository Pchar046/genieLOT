#' Interactive LOT explorer
#'
#' A Shiny app: choose a disease preset, the index event and the line rules,
#' and see headline numbers, lines reached, the naive-vs-correct comparison,
#' validation against the BPC recurrence flag, and patient timelines update
#' live.
#'
#' Privacy mode (on by default) hides patient timelines and masks any count
#' under 10, so the app can be shown on a projector or shared screen.
#'
#' @param x A [lot_data()] object. Defaults to simulated data.
#' @param ca_dx Optional raw diagnosis table holding `dmets_stage_i_iii`, to
#'   show the validation tab.
#' @param privacy Start with privacy mode on.
#' @param cols Column mapping used to build `x` (for validation).
#' @examples
#' \dontrun{
#' run_lot_app()                                   # fake data
#' shiny::runApp(run_lot_app(x, ca_dx = raw$ca_dx, cols = cols))
#' }
#' @export
run_lot_app <- function(x = NULL, ca_dx = NULL, privacy = TRUE, cols = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Install shiny first: install.packages('shiny')", call. = FALSE)
  }
  if (is.null(x)) {
    sim <- simulate_bpc(200)
    x <- lot_data(sim$ca_dx, sim$regimens, sim$imaging, sim$med_onc)
    ca_dx <- sim$ca_dx
  }
  if (is.null(cols)) cols <- x$cols
  has_flag <- !is.null(ca_dx) && "dmets_stage_i_iii" %in% names(ca_dx)
  ids <- unique(x$drugs$case_id)
  cohort <- unique(stats::na.omit(x$cases$cohort))
  start_preset <- if (length(cohort) && cohort[1] %in% c("NSCLC", "CRC", "BrCa")) cohort[1] else "NSCLC"
  fmt_swaps <- function(sw) paste(vapply(sw, paste, "", collapse = "/"), collapse = "; ")
  parse_list <- function(s, sep) { p <- trimws(strsplit(s, sep)[[1]]); p[nzchar(p)] }

  css <- "
    body { background: #f6f7f9; }
    .cards { display: flex; gap: 14px; flex-wrap: wrap; margin: 6px 0 18px; }
    .card { flex: 1 1 180px; background: #fff; border-radius: 10px; padding: 14px 18px;
            box-shadow: 0 1px 3px rgba(0,0,0,.08); border-left: 5px solid #2E6F9E; }
    .card.alert { border-left-color: #C0392B; }
    .card.ok { border-left-color: #1E8449; }
    .card .v { font-size: 30px; font-weight: 700; line-height: 1.1; color: #1b2631; }
    .card .l { font-size: 13px; color: #5d6d7e; margin-top: 4px; }
    .privacy { background: #eaf2f8; border-radius: 8px; padding: 8px 12px; font-size: 13px;
               color: #1f4e79; margin-bottom: 10px; }
    .well { background: #fff; }
    h2 { font-weight: 700; }
  "
  card <- function(value, label, cls = "") {
    shiny::div(class = paste("card", cls), shiny::div(class = "v", value),
               shiny::div(class = "l", label))
  }

  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(css))),
    shiny::titlePanel("genieLOT: line of therapy explorer"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 3,
        shiny::selectInput("preset", "Disease preset",
                           c("Lung (NSCLC)" = "NSCLC", "Colorectal (CRC)" = "CRC",
                             "Breast (BrCa)" = "BrCa", "Custom" = "Custom"),
                           selected = start_preset),
        shiny::selectInput("event", "Index event (where counting starts)",
                           c("Recurrence (early stage) / diagnosis (stage IV)" = "recurrence",
                             "Diagnosis (naive)" = "diagnosis",
                             "Distant metastasis" = "metastasis",
                             "First progression" = "progression",
                             "Progression after initial therapy (challenge definition)" = "progression_after_therapy",
                             "Initial therapy" = "initial_therapy")),
        shiny::selectInput("src", "Progression source",
                           c("Imaging or oncologist note" = "either", "Imaging only" = "imaging",
                             "Oncologist note only" = "med_onc", "Both agree" = "both")),
        shiny::sliderInput("min_days", "Ignore progression in first N days", 0, 365, 90, step = 15),
        shiny::hr(),
        shiny::sliderInput("window", "Regimen window (days)", 0, 90, 28),
        shiny::sliderInput("gap", "Treatment break that starts a new line (days)", 30, 365, 120, step = 15),
        shiny::textInput("swaps", "Swaps that stay in the same line (A/B; C/D)", ""),
        shiny::textInput("addons", "Add-on drugs that don't start a new line (comma separated)", ""),
        shiny::checkboxInput("trial", "Count investigational (trial) drugs as lines", TRUE),
        shiny::checkboxInput("prog_new", "Start a new line after every progression", FALSE),
        shiny::hr(),
        shiny::checkboxInput("privacy", "Privacy mode (hide patients, mask counts under 10)", privacy),
        shiny::conditionalPanel("!input.privacy",
          shiny::selectInput("pid", "Patient timeline", ids, selected = ids[1]))
      ),
      shiny::mainPanel(width = 9,
        shiny::uiOutput("privacy_note"),
        shiny::uiOutput("cards"),
        shiny::tabsetPanel(
          shiny::tabPanel("Overview",
            shiny::br(),
            shiny::plotOutput("lines_plot", height = "300px"),
            shiny::h4("Naive counting from diagnosis vs your definition"),
            shiny::tableOutput("cmp")),
          shiny::tabPanel("Lines",
            shiny::br(),
            shiny::tableOutput("summary"),
            shiny::downloadButton("dl", "Download summary (csv)")),
          shiny::tabPanel("Cohort",
            shiny::br(),
            shiny::p("Who is in the analysis, step by step, under the current index and rules."),
            shiny::tableOutput("attrition")),
          shiny::tabPanel("Sensitivity",
            shiny::br(),
            shiny::p("The same data under 8 standard definitions, side by side. Takes a moment."),
            shiny::actionButton("run_sens", "Run sensitivity analysis"),
            shiny::br(), shiny::br(),
            shiny::tableOutput("sens")),
          shiny::tabPanel("Validation",
            shiny::br(),
            shiny::uiOutput("val_ui")),
          shiny::tabPanel("Patient timeline",
            shiny::br(),
            shiny::uiOutput("patient_ui")),
          shiny::tabPanel("Rules",
            shiny::br(),
            shiny::verbatimTextOutput("rules_txt"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    shiny::observeEvent(input$preset, {
      r <- lot_preset(input$preset)
      shiny::updateSliderInput(session, "window", value = r$window)
      shiny::updateSliderInput(session, "gap", value = r$gap)
      shiny::updateTextInput(session, "swaps", value = fmt_swaps(r$swaps))
      shiny::updateTextInput(session, "addons", value = paste(r$add_without_new_line, collapse = ", "))
    })
    rules <- shiny::reactive({
      sw <- lapply(parse_list(input$swaps, ";"), function(g) parse_list(g, "/"))
      lot_rules(window = input$window, gap = input$gap,
                swaps = sw[lengths(sw) > 1],
                add_without_new_line = parse_list(input$addons, ","),
                ignore_drugs = if (isTRUE(input$trial)) character() else "Investigational Drug",
                progression_new_line = input$prog_new)
    })
    idx <- shiny::reactive(set_index(x, input$event, prog_source = input$src,
                                     min_days = input$min_days))
    cmp <- shiny::reactive(compare_lot(x, idx(), rules()))
    smry <- shiny::reactive(summarize_lot(cmp()$correct))
    val <- shiny::reactive(if (has_flag) validate_index(x, idx(), ca_dx, cols = cols) else NULL)
    mask <- function(df) if (isTRUE(input$privacy)) suppress_small(df) else df
    show_n <- function(n) if (isTRUE(input$privacy) && n > 0 && n < 10) "<10" else format(n, big.mark = ",")

    output$privacy_note <- shiny::renderUI({
      if (isTRUE(input$privacy)) {
        shiny::div(class = "privacy", "Privacy mode is on: patient timelines are hidden and counts under 10 are masked.")
      }
    })

    output$cards <- shiny::renderUI({
      s <- cmp()$summary
      sm <- smry()
      l2 <- if (nrow(sm) >= 2) sm$pct_of_indexed[2] else 0
      v <- val()
      vpct <- if (!is.null(v)) regmatches(v$summary$value[1], regexpr("[0-9.]+%", v$summary$value[1])) else "n/a"
      pct <- s$pct_mislabeled[2]
      shiny::div(class = "cards",
        card(show_n(s$patients[1]), "Patients with a line after the index"),
        card(if (is.na(pct)) "n/a" else paste0(pct, "%"),
             sprintf("Early-stage recurrences where naive counting gets Line 1 wrong (%s patients)",
                     show_n(s$patients[2])), "alert"),
        card(paste0(l2, "%"), "Reach Line 2"),
        card(vpct, "AACR-flagged recurrences we handle", "ok"))
    })

    output$lines_plot <- shiny::renderPlot({
      sm <- smry()
      if (!nrow(sm)) return(NULL)
      sm <- utils::head(sm, 8)
      sm$label <- factor(paste("Line", sm$line), levels = paste("Line", sm$line))
      sm$txt <- ifelse(isTRUE(input$privacy) & sm$patients < 10, "<10",
                       paste0(sm$patients, " (", sm$pct_of_indexed, "%)"))
      ggplot2::ggplot(sm, ggplot2::aes(x = label, y = patients)) +
        ggplot2::geom_col(fill = "#2E6F9E", width = 0.65) +
        ggplot2::geom_text(ggplot2::aes(label = txt), vjust = -0.4, size = 4.5) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
        ggplot2::labs(x = NULL, y = "Patients", title = "Patients reaching each line of therapy") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    })

    pretty <- function(df, labels) { names(df) <- labels[seq_along(df)]; df }
    output$cmp <- shiny::renderTable(pretty(mask(cmp()$summary),
      c("Group", "Patients", "Line 1 wrong (naive)", "% wrong", "Early treatments excluded")))
    output$summary <- shiny::renderTable(pretty(mask(smry()),
      c("Line", "Patients", "% of indexed", "Median days on line", "Median days to next line",
        "Most common regimens")))
    output$dl <- shiny::downloadHandler(
      filename = function() "lot_summary.csv",
      content = function(f) utils::write.csv(mask(smry()), f, row.names = FALSE))

    output$attrition <- shiny::renderTable({
      a <- lot_attrition(x, idx(), rules(), k = if (isTRUE(input$privacy)) 10 else 0)
      names(a) <- c("Step", "Cases"); a
    })
    sens <- shiny::eventReactive(input$run_sens, {
      shiny::withProgress(message = "Running 8 definitions", value = 0.3, {
        lot_sensitivity(x, default_scenarios(x, rules()), k = if (isTRUE(input$privacy)) 10 else 0)
      })
    })
    output$sens <- shiny::renderTable({
      pretty(sens(), c("Definition", "Patients with lines", "Early-stage relapses",
                       "% naive Line 1 wrong", "% reach Line 2", "% reach Line 3",
                       "Median days on Line 1"))
    })

    output$val_ui <- shiny::renderUI({
      if (!has_flag) return(shiny::p("Pass the raw diagnosis table (ca_dx) with dmets_stage_i_iii to see validation."))
      shiny::tagList(
        shiny::p("Early-stage patients: our index event vs BPC's own flag for later distant metastasis (dmets_stage_i_iii)."),
        shiny::tableOutput("val_sum"),
        shiny::tableOutput("val_tab"))
    })
    output$val_sum <- shiny::renderTable(val()$summary)
    output$val_tab <- shiny::renderTable({
      t <- as.data.frame.matrix(val()$table)
      t <- cbind(ours = rownames(t), t)
      names(t)[-1] <- paste("flag =", names(t)[-1])
      if (isTRUE(input$privacy)) suppress_small(t, cols = names(t)[-1]) else t
    })

    output$patient_ui <- shiny::renderUI({
      if (isTRUE(input$privacy)) {
        return(shiny::p("Hidden in privacy mode. Turn privacy mode off to view a patient's timeline."))
      }
      shiny::plotOutput("patient", height = "480px")
    })
    output$patient <- shiny::renderPlot(
      plot_patient(x, cmp()$correct, input$pid, naive = cmp()$naive))

    output$rules_txt <- shiny::renderPrint({
      print(rules())
      cat("\nIndex settings:\n")
      utils::str(attr(idx(), "settings"))
    })
  }
  shiny::shinyApp(ui, server)
}

utils::globalVariables(c("txt"))
