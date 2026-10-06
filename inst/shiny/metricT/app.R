# Metric-T graphical interface.  Launch with metricT::run_metricT_app()
library(shiny)
library(metricT)
options(shiny.maxRequestSize = 4 * 1024^3)

fmt_table <- function(t) {
  data.frame(
    Comparison = t$comparison, Band = t$band,
    `DC(wPLI) %` = sprintf("%.1f", t$DC_wPLI), `DC(MSC) %` = sprintf("%.1f", t$DC_MSC),
    `T (pp)` = sprintf("%+.1f", t$T), Reversal = ifelse(t$reversal, "YES", ""),
    `p (raw)` = sprintf("%.3f", t$p_raw), `min attainable p` = sprintf("%.3f", t$p_min), `p (Bonferroni)` = sprintf("%.3f", t$p_bonferroni),
    `p (Holm)` = sprintf("%.3f", t$p_holm), `q (BH)` = sprintf("%.3f", t$q_BH),
    `p (Westfall-Young)` = ifelse(is.na(t$p_WY), "", sprintf("%.3f", t$p_WY)),
    check.names = FALSE, stringsAsFactors = FALSE)
}

ui <- fluidPage(
  titlePanel("Metric-T: is the direction of your EEG connectivity result metric-dependent?"),
  sidebarLayout(
    sidebarPanel(width = 4,
      h4("1. Data"),
      radioButtons("src", NULL, c("Demo (simulated EEG)" = "demo",
                                  "Feature file (.npz from Python pipeline, or .csv)" = "feat",
                                  "Preprocessed EEG (one CSV per subject)" = "eeg")),
      conditionalPanel("input.src == 'demo'",
        sliderInput("n_per", "Subjects per group", 8, 60, 20),
        sliderInput("zlA", "Zero-lag fraction, group A", 0, 1, 0.58, 0.01),
        sliderInput("zlB", "Zero-lag fraction, group B", 0, 1, 0.42, 0.01),
        sliderInput("stB", "Coupling strength of B relative to A", 0.2, 1.5, 0.9, 0.05),
        sliderInput("gsd", "Channel-specific variability (0 = all pairs vary in unison)", 0, 1.5, 0.7, 0.1),
        sliderInput("n_ch", "Channels", 3, 19, 8)),
      conditionalPanel("input.src == 'feat'",
        fileInput("featfile", "features.npz or features CSV", accept = c(".npz", ".csv"))),
      conditionalPanel("input.src == 'eeg'",
        fileInput("eegfiles", "EEG files (samples x channels, header = channel names)",
                  multiple = TRUE, accept = ".csv"),
        fileInput("partfile", "Participants table (.tsv/.csv: participant_id, group)",
                  accept = c(".tsv", ".csv", ".txt")),
        numericInput("sfreq", "Sampling frequency (Hz)", 500, min = 1),
        fluidRow(column(6, numericInput("eplen", "Epoch (s)", 2, min = 0.1, step = 0.5)),
                 column(6, numericInput("overlap", "Overlap", 0.5, min = 0, max = 0.95, step = 0.05)))),
      actionButton("load", "Load data"),
      hr(),
      h4("2. Analysis"),
      uiOutput("grp_ui"),
      fluidRow(column(6, numericInput("nperm", "Permutations", 10000, min = 100, step = 1000)),
               column(6, numericInput("seed", "Seed", 42))),
      checkboxInput("exact", "Exact: enumerate every label assignment (two groups, small samples)", FALSE),
      actionButton("run", "Run Metric-T", class = "btn-primary")
    ),
    mainPanel(width = 8,
      tabsetPanel(id = "tabs",
        tabPanel("Data", br(), verbatimTextOutput("info")),
        tabPanel("Results", br(), tableOutput("tab"), textOutput("revtxt"), br(),
                 downloadButton("dl_csv", "Download table (CSV)")),
        tabPanel("DC plot", plotOutput("dcplot", height = "640px"),
                 downloadButton("dl_png", "Download figure (PNG)")),
        tabPanel("Null distribution", br(), uiOutput("cond_ui"), plotOutput("nullplot", height = "420px")),
        tabPanel("R code", br(), p("Script that reproduces this analysis without the GUI:"),
                 verbatimTextOutput("code")),
        tabPanel("About", br(),
          p("Metric-T compares two connectivity metrics computed from the same cross-spectrum.",
            "Direction consistency (DC) is the percentage of channel-pair features for which the first group",
            "exceeds the second; T = DC(wPLI) - DC(MSC). A directional reversal is flagged when the two DC",
            "values fall on opposite sides of 50%."),
          p("Significance: two-sided label permutation shared by all conditions; Bonferroni, Holm,",
            "Benjamini-Hochberg and Westfall-Young max-T adjustments are reported."),
          p(strong("Reading the p-value."), "DC keeps only the sign of each mean difference. When channel pairs",
            "vary in unison across subjects, relabelled data also give DC near 0% or 100%, and even a complete",
            "reversal (|T| = 100 pp) cannot become significant. The column 'min attainable p' shows the smallest",
            "p-value possible in each condition; if it exceeds 0.05, a large p is not evidence that the direction is stable.",
            "Try the demo with channel-specific variability set to 0."),
          p("Reference: Kimura A (2026). Neuroscience Informatics 6:100286. doi:10.1016/j.neuri.2026.100286"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(dat = NULL, res = NULL, code = NULL)

  observeEvent(input$load, {
    rv$res <- NULL
    dat <- tryCatch({
      if (input$src == "demo") {
        withProgress(message = "Simulating EEG and computing connectivity", {
          s <- mt_simulate(n = c(A = input$n_per, B = input$n_per), zero_lag = c(input$zlA, input$zlB),
                           strength = c(1, input$stB), gain_sd = input$gsd, n_ch = input$n_ch, seed = 1)
        })
        s$label <- "simulated demo data"
        s$code <- sprintf("dat <- mt_simulate(n = c(A = %d, B = %d), zero_lag = c(%.2f, %.2f), strength = c(1, %.2f), gain_sd = %.1f, n_ch = %d, seed = 1)",
                          input$n_per, input$n_per, input$zlA, input$zlB, input$stB, input$gsd, input$n_ch)
        s
      } else if (input$src == "feat") {
        req(input$featfile)
        f <- input$featfile
        s <- if (grepl("\\.npz$", f$name, ignore.case = TRUE)) mt_read_npz(f$datapath) else mt_read_features_csv(f$datapath)
        s$label <- f$name
        s$code <- sprintf("dat <- %s(\"%s\")", if (grepl("\\.npz$", f$name, ignore.case = TRUE)) "mt_read_npz" else "mt_read_features_csv", f$name)
        s
      } else {
        req(input$eegfiles, input$partfile)
        pf <- input$partfile
        part <- if (grepl("\\.csv$", pf$name, ignore.case = TRUE)) read.csv(pf$datapath, stringsAsFactors = FALSE)
                else read.delim(pf$datapath, stringsAsFactors = FALSE)
        names(part) <- tolower(names(part))
        idc <- intersect(c("participant_id", "subject", "id"), names(part))[1]
        if (is.na(idc) || !("group" %in% names(part))) stop("Participants table needs 'participant_id' and 'group' columns.")
        ids <- sub("\\.[^.]*$", "", input$eegfiles$name)
        m <- match(ids, part[[idc]])
        if (anyNA(m)) stop("No group found for: ", paste(ids[is.na(m)], collapse = ", "))
        withProgress(message = "Computing wPLI and MSC", value = 0, {
          feats <- mt_features_from_eeg_csv(input$eegfiles$datapath, ids, input$sfreq, input$eplen, input$overlap,
                                            progress = function(s) setProgress(s / length(ids)))
        })
        list(features = feats, group = as.character(part$group[m]),
             label = sprintf("%d EEG files", length(ids)),
             code = sprintf(paste0("files <- list.files(\"<folder>\", pattern = \"\\\\.csv$\", full.names = TRUE)\n",
                                   "part  <- read.delim(\"%s\")\n",
                                   "feats <- mt_features_from_eeg_csv(files, sfreq = %g, epoch_len = %g, overlap = %g)\n",
                                   "dat   <- list(features = feats, group = part$group[match(rownames(feats[[1]]$wpli), part$participant_id)])"),
                            pf$name, input$sfreq, input$eplen, input$overlap))
      }
    }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 10); NULL })
    rv$dat <- dat
    updateTabsetPanel(session, "tabs", selected = "Data")
  })

  output$grp_ui <- renderUI({
    req(rv$dat)
    lev <- unique(rv$dat$group)
    prs <- combn(lev, 2, simplify = FALSE)
    lab <- c(vapply(prs, paste, "", collapse = " vs "), vapply(prs, function(p) paste(rev(p), collapse = " vs "), ""))
    tagList(
      selectInput("control", "Reference (control) group", c("(none)", lev),
                  selected = if ("CN" %in% lev) "CN" else "(none)"),
      selectizeInput("cmps", "Comparisons (first vs second)", lab, selected = NULL, multiple = TRUE,
                     options = list(placeholder = "default: each group vs control, then the rest")))
  })

  output$info <- renderText({
    if (is.null(rv$dat)) return("Choose a data source and press 'Load data'.")
    d <- rv$dat; tb <- table(d$group)
    paste0("Source: ", d$label, "\nSubjects: ", length(d$group), "\nGroups: ",
           paste0(names(tb), " (n=", tb, ")", collapse = ", "), "\nBands: ",
           paste(names(d$features), collapse = ", "), "\nChannel-pair features per band: ",
           paste(vapply(d$features, function(f) ncol(f$wpli), 1L), collapse = ", "),
           "\n\nwPLI range: ", paste(signif(range(unlist(lapply(d$features, function(f) range(f$wpli)))), 3), collapse = " to "),
           "\nMSC range:  ", paste(signif(range(unlist(lapply(d$features, function(f) range(f$msc)))), 3), collapse = " to "))
  })

  observeEvent(input$run, {
    req(rv$dat)
    ctrl <- if (is.null(input$control) || input$control == "(none)") NULL else input$control
    cmps <- if (length(input$cmps)) strsplit(input$cmps, " vs ", fixed = TRUE) else NULL
    res <- tryCatch(
      withProgress(message = "Permuting group labels", value = 0, {
        metric_t(rv$dat$features, rv$dat$group, comparisons = cmps, control = ctrl,
                 n_perm = input$nperm, seed = input$seed, exact = isTRUE(input$exact),
                 progress = function(f) setProgress(f))
      }), error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 10); NULL })
    rv$res <- res
    cm <- if (is.null(cmps)) "" else paste0(", comparisons = list(",
            paste(vapply(cmps, function(p) sprintf("c(\"%s\", \"%s\")", p[1], p[2]), ""), collapse = ", "), ")")
    rv$code <- paste0("library(metricT)\n", rv$dat$code, "\n",
                      "res <- metric_t(dat$features, dat$group", cm,
                      if (!is.null(ctrl)) sprintf(", control = \"%s\"", ctrl) else "",
                      if (isTRUE(input$exact)) ", exact = TRUE)\n" else
                        sprintf(", n_perm = %d, seed = %d)\n", as.integer(input$nperm), as.integer(input$seed)),
                      "print(res)\nplot(res)\nwrite.csv(as.data.frame(res), \"metricT_results.csv\", row.names = FALSE)")
    if (!is.null(res)) updateTabsetPanel(session, "tabs", selected = "Results")
  })

  output$tab <- renderTable({ req(rv$res); fmt_table(rv$res$table) }, striped = TRUE, spacing = "s")
  output$revtxt <- renderText({
    req(rv$res); t <- rv$res$table
    sprintf("Directional reversals: %d of %d conditions. Smallest raw p = %.3f; %d condition(s) with Holm-adjusted p < 0.05.",
            sum(t$reversal), nrow(t), min(t$p_raw), sum(t$p_holm < 0.05))
  })
  output$dcplot <- renderPlot({ req(rv$res); plot(rv$res) })
  output$cond_ui <- renderUI({
    req(rv$res); t <- rv$res$table
    selectInput("cond", "Condition", stats::setNames(seq_len(nrow(t)), paste(t$comparison, t$band, sep = " | ")))
  })
  output$nullplot <- renderPlot({ req(rv$res, input$cond); mt_null_plot(rv$res, as.integer(input$cond)) })
  output$code <- renderText({ if (is.null(rv$code)) "Run an analysis first." else rv$code })
  output$dl_csv <- downloadHandler("metricT_results.csv",
    function(file) write.csv(rv$res$table, file, row.names = FALSE))
  output$dl_png <- downloadHandler("metricT_DC.png", function(file) {
    grDevices::png(file, width = 1600, height = 520 * length(unique(rv$res$table$comparison)), res = 200)
    plot(rv$res); grDevices::dev.off()
  })
}

shinyApp(ui, server)
