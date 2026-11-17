## =====================================================================
##  CUBpanel -- shiny front end
##  ---------------------------------------------------------------------
##  Run locally with:   shiny::runApp(".")
##  Required: shiny.  Optional: readxl (for .xlsx upload).
##
##  The front end owns no statistics.  Every number it shows comes from
##  cubpanel_analyse() in R/cubpanel-engine.R and is formatted by
##  R/cubpanel-report.R, so the shiny path and the command-line path of
##  scripts/reproduce-paper.R produce identical output.
## =====================================================================

library(shiny)

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE)))
  source(f, local = FALSE)

options(shiny.maxRequestSize = 64 * 1024^2, digits = 6, width = 160)

HAS_READXL <- requireNamespace("readxl", quietly = TRUE)

## ---------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------

cub_read_upload <- function(path, name, sep = "auto", sheet = 1L) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) {
    if (!HAS_READXL)
      stop("reading .xlsx requires the readxl package: ",
           'install.packages("readxl")')
    as.data.frame(readxl::read_excel(path, sheet = as.integer(sheet)),
                  stringsAsFactors = FALSE)
  } else if (ext == "rds") {
    as.data.frame(readRDS(path))
  } else {
    s <- if (identical(sep, "auto")) {
      l1 <- readLines(path, n = 1L, warn = FALSE)
      cnt <- c(`,` = lengths(regmatches(l1, gregexpr(",", l1))),
               `;` = lengths(regmatches(l1, gregexpr(";", l1))),
               `\t` = lengths(regmatches(l1, gregexpr("\t", l1))))
      names(cnt)[which.max(cnt)]
    } else sep
    utils::read.table(path, header = TRUE, sep = s, quote = "\"",
                      stringsAsFactors = FALSE, check.names = TRUE,
                      comment.char = "", fill = TRUE)
  }
}

cub_guess_role <- function(nms, patterns, fallback = NA_character_) {
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  if (is.na(fallback)) nms[1] else fallback
}

cub_mk_formula <- function(vars) {
  if (!length(vars)) return(~ 1)
  stats::reformulate(vars)
}

.pre <- function(txt) tags$pre(paste(txt, collapse = "\n"))

## memoised panel generators: the shiny session should not pay twice
.PANEL_CACHE <- new.env(parent = emptyenv())
cub_demo_panel <- function(which) {
  key <- paste0("panel", which)
  if (!is.null(.PANEL_CACHE[[key]])) return(.PANEL_CACHE[[key]])
  d <- if (which == "A") cub_sim_panelA() else cub_sim_panelB()
  .PANEL_CACHE[[key]] <- d
  d
}

FIG_SIZE <- list(
  waves       = c(2200, 1100), within      = c(1700, 1200),
  dispersion  = c(1700, 1250), weights     = c(1700, 1150),
  ladder      = c(1700, 1250), profiles    = c(2200, 1100),
  diagnostics = c(2200, 1700), memory      = c(2200, 1100))

## ---------------------------------------------------------------------
## user interface
## ---------------------------------------------------------------------

ui <- fluidPage(
  title = "CUBpanel",
  tags$head(tags$style(HTML("
    body { font-size: 13px; }
    pre  { font-size: 11.5px; line-height: 1.32; background: #fbfbfb;
           border: 1px solid #e3e3e3; padding: 9px; }
    h4   { margin-top: 4px; font-weight: 600; }
    .wellish { background:#f7f7f7; border:1px solid #e3e3e3;
               border-radius:3px; padding:8px 10px; margin-bottom:8px; }
    .card { display:inline-block; border:1px solid #ddd; border-radius:3px;
            padding:6px 11px; margin:0 7px 7px 0; background:#fff;
            min-width:118px; }
    .card .lab { font-size:10.5px; color:#666; text-transform:uppercase;
                 letter-spacing:.4px; }
    .card .val { font-size:17px; font-weight:600; }
    .flag-ok   { color:#1a7f37; } .flag-bad { color:#a40e26; }
    .nav-tabs > li > a { padding: 5px 9px; font-size: 12px; }
    .shiny-input-container { margin-bottom: 7px; }
  "))),

  titlePanel(div(
    span("CUBpanel", style = "font-weight:700;"),
    span(" \u2014 panel CUB models for repeated ordinal ratings, with a free shelter effect",
         style = "font-weight:400; color:#555; font-size:15px;"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      div(class = "wellish",
          h4("1. Data"),
          radioButtons("src", NULL, inline = FALSE,
                       choices = c("Simulated panel B (all three components)" = "B",
                                   "Simulated panel A (collapsing mixture)" = "A",
                                   "Upload a file" = "up"),
                       selected = "B"),
          conditionalPanel(
            "input.src == 'up'",
            fileInput("file", "CSV, TSV, XLSX or RDS",
                      accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls", ".rds")),
            fluidRow(
              column(6, selectInput("sep", "Separator",
                                    c("auto" = "auto", "comma" = ",",
                                      "semicolon" = ";", "tab" = "\t"))),
              column(6, numericInput("sheet", "Sheet", 1, min = 1, step = 1))),
            if (!HAS_READXL)
              helpText(HTML("<span class='flag-bad'>readxl is not installed; .xlsx upload is unavailable.</span>"))),
          actionButton("load", "Load data", class = "btn-primary btn-sm",
                       width = "100%")),

      div(class = "wellish",
          h4("2. Column roles"),
          uiOutput("ui_roles")),

      div(class = "wellish",
          h4("3. Scale and shelter"),
          uiOutput("ui_scale")),

      div(class = "wellish",
          h4("4. The three equations"),
          uiOutput("ui_eqs")),

      div(class = "wellish",
          h4("5. Numerical controls"),
          fluidRow(
            column(6, numericInput("Q", "Q (estimation)", 15, 5, 41, 2)),
            column(6, numericInput("Qpred", "Q (prediction)", 31, 5, 61, 2))),
          fluidRow(
            column(6, numericInput("nstart", "Pooled starts", 3, 1, 10, 1)),
            column(6, numericInput("chibarB", "chi-bar draws", 5000, 200,
                                   50000, 100))),
          checkboxInput("selfcheck", "Run self-checks before fitting", TRUE),
          checkboxInput("rhoaudit", "Correlation audit (profile)", TRUE),
          checkboxInput("prof", "Profile interval for sigma^2_1", TRUE),
          checkboxInput("profgrid", "Profile GRID for sigma^2_1 (slower)",
                        FALSE),
          checkboxInput("oos", "Out-of-sample on the held-out final wave",
                        TRUE),
          numericInput("budget", "Seconds per profile budget", 150, 0, 3600,
                       10)),

      actionButton("run", "Run the analysis", class = "btn-success",
                   width = "100%", icon = icon("play")),
      br(), br(),
      downloadButton("dl_report", "Full report (.txt)", class = "btn-sm"),
      downloadButton("dl_rds", "Analysis object (.rds)", class = "btn-sm"),
      br(), br(),
      helpText(HTML(paste0(
        "<small>The interface owns no statistics: every number comes from ",
        "<code>cubpanel_analyse()</code>. Reproduce the manuscript with ",
        "<code>Rscript scripts/reproduce-paper.R</code>.</small>")))
    ),

    mainPanel(
      width = 9,
      uiOutput("cards"),
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("Data",
                 br(), plotOutput("p_waves", height = "330px"),
                 fluidRow(
                   column(6, downloadButton("dl_waves", "PNG", class = "btn-xs")),
                   column(6, downloadButton("dl_within", "PNG (within)",
                                            class = "btn-xs"))),
                 br(), plotOutput("p_within", height = "330px"),
                 br(), uiOutput("t_data")),

        tabPanel("Dispersion",
                 br(), plotOutput("p_disp", height = "420px"),
                 downloadButton("dl_disp", "PNG", class = "btn-xs"),
                 br(), br(), uiOutput("t_disp"),
                 br(), plotOutput("p_weights", height = "380px"),
                 downloadButton("dl_weights", "PNG", class = "btn-xs"),
                 br(), br(), uiOutput("t_verify")),

        tabPanel("Self-checks", br(), uiOutput("t_checks")),

        tabPanel("Ladder",
                 br(), plotOutput("p_ladder", height = "400px"),
                 downloadButton("dl_ladder", "PNG", class = "btn-xs"),
                 br(), br(), uiOutput("t_ladder"),
                 br(), uiOutput("t_chains")),

        tabPanel("Correlation",
                 br(), plotOutput("p_prof", height = "380px"),
                 downloadButton("dl_prof", "PNG", class = "btn-xs"),
                 br(), br(), uiOutput("t_rho"),
                 br(), uiOutput("t_profile")),

        tabPanel("Coefficients",
                 br(), uiOutput("ui_coefpick"), uiOutput("t_coef")),

        tabPanel("Tests", br(), uiOutput("t_tests")),

        tabPanel("Fit",
                 br(), plotOutput("p_diag", height = "640px"),
                 downloadButton("dl_diag", "PNG", class = "btn-xs"),
                 br(), br(), uiOutput("t_fit")),

        tabPanel("Quadrature", br(), uiOutput("t_quad")),

        tabPanel("Subject effects", br(), uiOutput("t_eb")),

        tabPanel("Out-of-sample", br(), uiOutput("t_oos")),

        tabPanel("Memory",
                 br(), plotOutput("p_mem", height = "380px"),
                 downloadButton("dl_mem", "PNG", class = "btn-xs"),
                 br(), br(), uiOutput("t_mem")),

        tabPanel("Export",
                 br(),
                 fluidRow(
                   column(4, textInput("prefix", "Key prefix", "X")),
                   column(8, br(),
                          downloadButton("dl_tex", "LaTeX bodies (.tex)",
                                         class = "btn-sm"),
                          downloadButton("dl_csv", "Ladder (.csv)",
                                         class = "btn-sm"),
                          downloadButton("dl_coefcsv", "Coefficients (.csv)",
                                         class = "btn-sm"))),
                 uiOutput("t_export")),

        tabPanel("About", br(), uiOutput("t_about"))
      )
    )
  )
)

## ---------------------------------------------------------------------
## server
## ---------------------------------------------------------------------

server <- function(input, output, session) {

  RAW <- reactiveVal(NULL)
  AN  <- reactiveVal(NULL)
  MSG <- reactiveVal(NULL)

  ## ---- load ---------------------------------------------------------
  observeEvent(input$load, {
    AN(NULL); MSG(NULL)
    d <- tryCatch({
      if (input$src == "up") {
        req(input$file)
        cub_read_upload(input$file$datapath, input$file$name,
                        input$sep, input$sheet)
      } else {
        withProgress(message = "Generating the simulated panel", value = 0.5, {
          cub_demo_panel(input$src)
        })
      }
    }, error = function(e) { MSG(conditionMessage(e)); NULL })
    if (is.null(d)) return(NULL)
    for (j in names(d)) if (is.character(d[[j]])) d[[j]] <- factor(d[[j]])
    attr(d, "source") <- if (input$src == "up") input$file$name
    else paste("simulated panel", input$src)
    RAW(d)
  }, ignoreInit = TRUE)

  ## load panel B once at startup so the interface is never empty
  observeEvent(TRUE, {
    d <- tryCatch(cub_demo_panel("B"), error = function(e) NULL)
    if (!is.null(d)) {
      attr(d, "source") <- "simulated panel B"
      RAW(d)
    }
  }, once = TRUE)

  ## ---- role / scale / equation inputs -------------------------------
  output$ui_roles <- renderUI({
    d <- RAW(); if (is.null(d)) return(helpText("Load data first."))
    nms <- names(d)
    tagList(
      selectInput("id", "Respondent identifier", nms,
                  cub_guess_role(nms, c("^respond", "^id$", "subject", "unit"))),
      selectInput("time", "Wave identifier", nms,
                  cub_guess_role(nms, c("^period", "^wave", "^time", "^year"))),
      selectInput("resp", "Ordinal response", nms,
                  cub_guess_role(nms, c("^response", "^rating", "^score",
                                        "^y$"))))
  })

  output$ui_scale <- renderUI({
    d <- RAW(); req(d, input$resp)
    y <- suppressWarnings(as.integer(d[[input$resp]]))
    mmax <- if (all(is.na(y))) 7L else max(y, na.rm = TRUE)
    tagList(
      checkboxInput("auto_m", "Resolve m from max(response)", TRUE),
      conditionalPanel("!input.auto_m",
                       numericInput("m", "Scale categories m", mmax, 4, 40, 1)),
      selectInput("cs", "Declared shelter category c",
                  choices = seq_len(mmax),
                  selected = if (mmax >= 7) 4 else ceiling(mmax / 2)),
      helpText(HTML(paste0(
        "<small>c is <b>declared</b>, never chosen by maximising the ",
        "likelihood: such a search leaves no admissible reference ",
        "distribution for the tests that follow. Observed support: {",
        paste(sort(unique(y[!is.na(y)])), collapse = ","), "}.</small>"))))
  })

  output$ui_eqs <- renderUI({
    d <- RAW(); req(d, input$id, input$time, input$resp)
    cand <- setdiff(names(d), c(input$id, input$time, input$resp))
    sel <- function(wanted) intersect(wanted, cand)
    tagList(
      selectizeInput("f_feel", "Feeling (phi)", cand, multiple = TRUE,
                     selected = sel(c("Gender", "Experience", "Inflation"))),
      selectizeInput("f_unc", "Uncertainty (upsilon)", cand, multiple = TRUE,
                     selected = sel(c("Education", "Delivery_Speed",
                                      "Customer_Service"))),
      selectizeInput("f_shel", "Shelter (delta)", cand, multiple = TRUE,
                     selected = sel("Inflation")),
      helpText(HTML(paste0(
        "<small>Each component owns a predictor of its own, so a covariate ",
        "may raise indecision while leaving evaluation untouched. An empty ",
        "selection leaves that equation with an intercept only.</small>"))))
  })

  ## ---- run ----------------------------------------------------------
  observeEvent(input$run, {
    d <- RAW(); req(d)
    MSG(NULL)
    res <- tryCatch(withProgress(message = "CUBpanel", value = 0, {
      D <- cub_prepare(
        d, id = input$id, time = input$time, response = input$resp,
        m = if (isTRUE(input$auto_m)) NULL else input$m,
        feeling = cub_mk_formula(input$f_feel),
        uncertainty = cub_mk_formula(input$f_unc),
        shelter = cub_mk_formula(input$f_shel),
        shelter_cat = as.integer(input$cs))

      ctl <- cubpanel_control(
        Q = as.integer(input$Q), Qpred = as.integer(input$Qpred),
        nstart_pool = as.integer(input$nstart),
        selfcheck = isTRUE(input$selfcheck),
        rho_audit = isTRUE(input$rhoaudit),
        rho_budget = as.numeric(input$budget),
        profile = isTRUE(input$prof),
        profile_budget = as.numeric(input$budget),
        profile_grid = isTRUE(input$profgrid),
        profile_grid_budget = as.numeric(input$budget),
        chibar_B = as.integer(input$chibarB),
        oos = isTRUE(input$oos),
        progress = function(frac, msg)
          setProgress(value = frac, detail = msg))

      truth <- NULL; re_true <- NULL
      if (identical(attr(d, "panel"), "B")) {
        tr <- attr(d, "truth")
        truth <- cub_truth_vector(D, tr)
        attr(truth, "truth_list") <- tr
        re_true <- attr(d, "re")[D$subj, , drop = FALSE]
      }
      cubpanel_analyse(D, ctl, truth = truth, re_true = re_true)
    }), error = function(e) { MSG(conditionMessage(e)); NULL })
    if (!is.null(res)) {
      AN(res)
      updateTabsetPanel(session, "tabs", selected = "Ladder")
    }
  }, ignoreInit = TRUE)

  ## ---- status cards -------------------------------------------------
  output$cards <- renderUI({
    if (!is.null(MSG()))
      return(div(class = "wellish",
                 HTML(sprintf("<span class='flag-bad'><b>Error:</b> %s</span>",
                              htmltools::htmlEscape(MSG())))))
    an <- AN()
    if (is.null(an)) {
      d <- RAW()
      return(div(class = "wellish", HTML(sprintf(
        "Loaded: <b>%s</b>%s. Set the column roles, declare the shelter category, then press <b>Run the analysis</b>.",
        if (is.null(d)) "nothing yet" else attr(d, "source"),
        if (is.null(d)) "" else sprintf(" (%d rows, %d columns)",
                                        nrow(d), ncol(d))))))
    }
    card <- function(lab, val, cls = "") div(
      class = "card", div(class = "lab", lab),
      div(class = paste("val", cls), HTML(val)))
    R <- an$floor[["R"]]
    tagList(
      card("Respondents", an$D$N),
      card("Est. rows", an$D$n),
      card("Scale m", an$D$m),
      card("Shelter c", an$D$cs),
      card("Dispersion R", sprintf("%.3f", R),
           if (R < 1) "flag-bad" else "flag-ok"),
      card("BIC_N selects", an$fits[[an$selected$bicN]]$label),
      card("AIC selects", an$fits[[an$selected$aic]]$label),
      card("Quadrature",
           if (an$quad$adequate) "certified" else
             if (isTRUE(an$quad$refitted)) "refitted" else "insufficient",
           if (an$quad$adequate || isTRUE(an$quad$refitted))
             "flag-ok" else "flag-bad"),
      card("Collapse", if (isTRUE(an$verif$collapsed)) "YES" else "no",
           if (isTRUE(an$verif$collapsed)) "flag-bad" else "flag-ok"),
      card("Elapsed", sprintf("%.0f s", an$secs)))
  })

  ## ---- text blocks --------------------------------------------------
  blk <- function(fun, ...) renderUI({
    an <- AN(); if (is.null(an)) return(helpText("Run the analysis first."))
    .pre(fun(an, ...))
  })
  output$t_data    <- blk(cub_report_data)
  output$t_disp    <- blk(cub_report_dispersion)
  output$t_verify  <- blk(cub_report_verify)
  output$t_checks  <- blk(cub_report_checks)
  output$t_ladder  <- blk(cub_report_ladder)
  output$t_chains  <- blk(cub_report_chains)
  output$t_rho     <- blk(cub_report_rho)
  output$t_profile <- blk(cub_report_profile)
  output$t_tests   <- blk(cub_report_tests)
  output$t_fit     <- blk(cub_report_fit)
  output$t_quad    <- blk(cub_report_quad)
  output$t_eb      <- blk(cub_report_eb)
  output$t_oos     <- blk(cub_report_oos)
  output$t_mem     <- blk(cub_report_memory)

  output$ui_coefpick <- renderUI({
    an <- AN(); req(an)
    radioButtons("coefpick", "Rung", inline = TRUE,
                 choices = setNames(names(an$coef),
                                    vapply(an$coef, `[[`, "", "label")),
                 selected = an$primary_coef)
  })
  output$t_coef <- renderUI({
    an <- AN(); if (is.null(an)) return(helpText("Run the analysis first."))
    .pre(cub_report_coef(an, if (is.null(input$coefpick)) NULL
                         else input$coefpick))
  })
  output$t_export <- renderUI({
    an <- AN(); if (is.null(an)) return(helpText("Run the analysis first."))
    pfx <- if (nzchar(input$prefix %||% "")) input$prefix else "X"
    .pre(c(cub_report_export(an, pfx), "", cub_latex_ladder(an), "",
           cub_latex_coef(an), "", cub_latex_tests(an)))
  })
  `%||%` <- function(a, b) if (is.null(a)) b else a

  output$t_about <- renderUI(.pre(c(
    .rule("About", "#"),
    "CUBpanel estimates a three-component mixture for repeated ordinal",
    "ratings.  A shifted binomial carries evaluation, a discrete uniform",
    "carries indecision, and a point mass at an analyst-declared category",
    "carries retreat to a refuge.  Each component owns a logistic predictor",
    "of its own, so a covariate may raise indecision while leaving",
    "evaluation untouched, and the software reports that as such.",
    "",
    "Six commitments are visible in every output panel.",
    "",
    "  1. Dimensions are read, not assumed.  The declared m is printed",
    "     beside the observed support, because the two can differ and the",
    "     difference changes the model.",
    "  2. Quadrature is fixed, Cholesky-scaled and CERTIFIED rather than",
    "     adaptive, with the tolerance stated on the scale of the inference",
    "     the quadrature supports: N * err(Q), not err(Q).",
    "  3. The shelter category is declared, never selected by maximising",
    "     the likelihood over the m candidates.",
    "  4. Boundary honesty.  No variance component is forced away from zero",
    "     and no Wald statistic is printed for a parameter whose null value",
    "     lies on the boundary.  When a frontier flag is raised, the",
    "     chi-bar-squared reference is WITHHELD rather than printed with a",
    "     disclaimer: such a p-value has no known relation to the",
    "     hypothesis, and printing it beside valid ones is the most",
    "     effective way to mislead a reader.",
    "  5. The dispersion floor is a PRE-FIT gate.  Both extra components",
    "     can only add dispersion, so one weighted least-squares solve",
    "     decides whether the family can represent the data at all.",
    "  6. Two marginal effects, never one.  A covariate's effect on the",
    "     shape of a response distribution and its effect on the location",
    "     of that distribution are different statements and routinely",
    "     carry opposite signs.",
    "",
    "Sigma is parameterised through its Cholesky factor, which makes",
    "rho_12 = 0 an INTERIOR point and the correlation testable by an",
    "ordinary chi^2_1 likelihood ratio rather than by a cone argument.",
    "",
    "Dependencies: base, stats, parallel, shiny.  Optional: readxl.",
    sprintf("Session: R %s on %s.", getRversion(), R.version$platform))))

  ## ---- plots --------------------------------------------------------
  pl <- function(which) renderPlot({
    an <- AN()
    if (is.null(an)) {
      d <- RAW()
      if (which == "waves" && !is.null(d)) {
        op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op))
        plot.new()
        text(0.5, 0.5, "Press 'Run the analysis'", cex = 1.2, col = grey(0.4))
        return(invisible())
      }
      plot.new()
      text(0.5, 0.5, "Press 'Run the analysis'", cex = 1.2, col = grey(0.4))
      return(invisible())
    }
    cub_draw(an, which)
  })
  output$p_waves   <- pl("waves")
  output$p_within  <- pl("within")
  output$p_disp    <- pl("dispersion")
  output$p_weights <- pl("weights")
  output$p_ladder  <- pl("ladder")
  output$p_prof    <- pl("profiles")
  output$p_diag    <- pl("diagnostics")
  output$p_mem     <- pl("memory")

  dlfig <- function(which) downloadHandler(
    filename = function() sprintf("CUBpanel-%s.png", which),
    content = function(file) {
      an <- AN(); req(an)
      sz <- FIG_SIZE[[which]]
      cub_png(cub_draw(an, which), file, sz[1], sz[2])
    })
  output$dl_waves   <- dlfig("waves")
  output$dl_within  <- dlfig("within")
  output$dl_disp    <- dlfig("dispersion")
  output$dl_weights <- dlfig("weights")
  output$dl_ladder  <- dlfig("ladder")
  output$dl_prof    <- dlfig("profiles")
  output$dl_diag    <- dlfig("diagnostics")
  output$dl_mem     <- dlfig("memory")

  ## ---- downloads ----------------------------------------------------
  output$dl_report <- downloadHandler(
    filename = function()
      sprintf("CUBpanel-report-%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      an <- AN(); req(an)
      writeLines(cubpanel_report(an, prefix = input$prefix %||% "X"), file)
    })
  output$dl_rds <- downloadHandler(
    filename = function() "CUBpanel-analysis.rds",
    content = function(file) { an <- AN(); req(an); saveRDS(an, file) })
  output$dl_tex <- downloadHandler(
    filename = function() "CUBpanel-tables.tex",
    content = function(file) {
      an <- AN(); req(an)
      writeLines(c(cub_latex_ladder(an), "", cub_latex_coef(an), "",
                   cub_latex_tests(an)), file)
    })
  output$dl_csv <- downloadHandler(
    filename = function() "CUBpanel-ladder.csv",
    content = function(file) {
      an <- AN(); req(an)
      utils::write.csv(an$ladder$table, file, row.names = TRUE)
    })
  output$dl_coefcsv <- downloadHandler(
    filename = function() "CUBpanel-coefficients.csv",
    content = function(file) {
      an <- AN(); req(an)
      nm <- if (is.null(input$coefpick)) an$primary_coef else input$coefpick
      utils::write.csv(an$coef[[nm]]$tab, file, row.names = TRUE)
    })
}

shinyApp(ui, server)