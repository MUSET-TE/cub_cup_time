## =====================================================================
##  Reproduce every table and figure of the manuscript.
##
##    Rscript scripts/reproduce-paper.R            demo, no replication
##    Rscript scripts/reproduce-paper.R mc 500     demo + 500 replications
##
##  Writes:  output/CUBpanel-fig{1..6}-*.png
##           output/CUBpanel-report-panel{A,B}.txt
##           output/CUBpanel-tables.tex
##           output/CUBpanel-pending.txt
## =====================================================================

rm(list = ls()); invisible(gc(FALSE))
options(digits = 6, width = 135, warn = 1, stringsAsFactors = FALSE)
T0 <- Sys.time()
tic <- function(...) cat(sprintf("[%7.1fs] %s\n",
                                 as.numeric(difftime(Sys.time(), T0,
                                                     units = "secs")),
                                 paste0(...)))

root <- if (file.exists("R/cubpanel-engine.R")) "." else ".."
for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$",
                          full.names = TRUE))) source(f)
OUT <- file.path(root, "output")
dir.create(OUT, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
RUN_MC <- length(args) && args[1] == "mc"
R_MC   <- if (RUN_MC && length(args) > 1) as.integer(args[2]) else 500L
NCORES <- max(1L, min(parallel::detectCores(), 4L))

cat(sprintf("CUBpanel reproduction | R %s | %d cores | MC: %s\n",
            getRversion(), NCORES, if (RUN_MC) paste("R =", R_MC) else "off"))

CTL <- cubpanel_control(
  Q = 15L, Qpred = 31L, Qcheck = c(9L, 15L, 21L, 27L),
  nstart_pool = 3L, selfcheck = TRUE,
  rho_audit = TRUE, rho_budget = 600,
  profile = TRUE, profile_budget = 600,
  profile_grid = TRUE, profile_grid_budget = 600, profile_grid_n = 7L,
  chibar_B = 5000L, oos = TRUE, cold_start_audit = TRUE,
  progress = function(frac, msg) tic(sprintf("  [%3.0f%%] %s", 100 * frac, msg)))

## ---- panel A ---------------------------------------------------------
tic("panel A: generation (the calibration VERIFIES the collapse)")
datA <- cub_sim_panelA()
cat(sprintf("  calibration: %s\n",
            paste(sprintf("%s = %.3g", names(attr(datA, "calibration")),
                          attr(datA, "calibration")), collapse = " | ")))
DA <- cub_prepare(datA, shelter_cat = attr(datA, "shelter_cat"))
tic("panel A: analysis")
anA <- cubpanel_analyse(DA, CTL)
writeLines(cubpanel_report(anA, prefix = "A"),
           file.path(OUT, "CUBpanel-report-panelA.txt"))

## ---- panel B ---------------------------------------------------------
tic("panel B: generation")
datB <- cub_sim_panelB()
DB <- cub_prepare(datB, shelter_cat = attr(datB, "shelter_cat"))
truthB <- cub_truth_vector(DB, attr(datB, "truth"))
attr(truthB, "truth_list") <- attr(datB, "truth")
tic("panel B: analysis")
anB <- cubpanel_analyse(DB, CTL, truth = truthB,
                        re_true = attr(datB, "re")[DB$subj, , drop = FALSE])
writeLines(cubpanel_report(anB, prefix = "B"),
           file.path(OUT, "CUBpanel-report-panelB.txt"))

## ---- figures ---------------------------------------------------------
tic("figures")
RA <- anA$floor[["R"]]; RB <- anB$floor[["R"]]

cub_png(cub_plot_waves_pair(DA, RA, "Panel A", DB, RB, "Panel B"),
        file.path(OUT, "CUBpanel-fig1-data.png"), 2200, 1100)

cub_png(cub_plot_dispersion(list(
  list(label = "panel B", D = DB, mu = anB$stats_selected$mu, R = RB,
       pch_obs = 16, pch_floor = 1, col = grey(0.10)),
  list(label = "panel A", D = DA, mu = anA$stats_selected$mu, R = RA,
       pch_obs = 17, pch_floor = 2, col = grey(0.42)))),
  file.path(OUT, "CUBpanel-fig2-dispersion.png"), 1700, 1250)

cub_png(cub_plot_weights(DA, anA$verif$post0,
                         "Panel A: CUBpanel(0) mixing weights at the numerical floor"),
        file.path(OUT, "CUBpanel-fig3-collapse.png"), 1700, 1150)

cub_png(cub_plot_profiles(anB$rho, anB$profile_s1_grid, anB$profile_s1),
        file.path(OUT, "CUBpanel-fig4-profiles.png"), 2200, 1100)

cub_png(cub_plot_diagnostics(DB, anB$stats_selected, anB$fit_selected$label),
        file.path(OUT, "CUBpanel-fig5-diagnostics.png"), 2200, 1700)

fdB <- anB$fits[[anB$spec_dyn]]
cub_png(cub_plot_memory(anB$memory, DB$m,
                        unname(fdB$th["xi:lag_f"]),
                        unname(fdB$th["ups:lag_u"])),
        file.path(OUT, "CUBpanel-fig6-memory.png"), 2200, 1100)

## two extra figures the shiny panels show and the paper may use
cub_png(cub_plot_ladder(anB$ladder),
        file.path(OUT, "CUBpanel-figA1-ladder.png"), 1700, 1250)
cub_png(cub_plot_within(DB),
        file.path(OUT, "CUBpanel-figA2-within.png"), 1700, 1200)

## ---- LaTeX bodies and pending values ---------------------------------
tic("LaTeX bodies and pending values")
writeLines(c("%% ===== panel A =====", cub_latex_ladder(anA), "",
             cub_latex_coef(anA), "", cub_latex_tests(anA), "",
             "%% ===== panel B =====", cub_latex_ladder(anB), "",
             cub_latex_coef(anB), "",
             cub_latex_coef(anB, if ("CRE-LAG" %in% names(anB$coef))
               "CRE-LAG" else anB$primary_coef), "",
             cub_latex_tests(anB)),
           file.path(OUT, "CUBpanel-tables.tex"))
PA <- cubpanel_export(anA, "A"); PB <- cubpanel_export(anB, "B")
writeLines(c("# values for \\pending{...}",
             sprintf("%-34s = %s", names(PA), unlist(PA)), "",
             sprintf("%-34s = %s", names(PB), unlist(PB))),
           file.path(OUT, "CUBpanel-pending.txt"))

## ---- console summary -------------------------------------------------
cat("\n"); cat(paste(cub_report_ladder(anA), collapse = "\n"), "\n")
cat("\n"); cat(paste(cub_report_verify(anA), collapse = "\n"), "\n")
cat("\n"); cat(paste(cub_report_ladder(anB), collapse = "\n"), "\n")
cat("\n"); cat(paste(cub_report_rho(anB), collapse = "\n"), "\n")
cat("\n"); cat(paste(cub_report_oos(anB), collapse = "\n"), "\n")

cat(sprintf(paste0("\nSUMMARY\n",
                   "  panel A: R = %.3f | BIC_N -> %s | AIC -> %s | collapse: %s\n",
                   "  panel B: R = %.3f | BIC_N -> %s | AIC -> %s | rho_12: %s\n"),
            anA$floor[["R"]], anA$fits[[anA$selected$bicN]]$label,
            anA$fits[[anA$selected$aic]]$label,
            if (isTRUE(anA$verif$collapsed)) "EXACT" else "none",
            anB$floor[["R"]], anB$fits[[anB$selected$bicN]]$label,
            anB$fits[[anB$selected$aic]]$label,
            if (is.null(anB$rho)) "not audited"
            else if (anB$rho$informative) "informative"
            else "uninformative corner -> Sigma diagonal"))

## ---- replication study -----------------------------------------------
if (RUN_MC) {
  tic(sprintf("replication study: R = %d, %d cores", R_MC, NCORES))
  mc <- cubpanel_mc(R = R_MC, ncores = NCORES,
                    cache = file.path(OUT, "mc-cache.rds"))
  txt <- cub_report_mc(mc)
  writeLines(txt, file.path(OUT, "CUBpanel-montecarlo.txt"))
  cat("\n"); cat(paste(txt, collapse = "\n"), "\n")
} else {
  tic("replication cost calibration (two replications)")
  mc <- cubpanel_mc(R = NCORES, ncores = NCORES)
  cat("\n"); cat(paste(cub_report_mc(mc), collapse = "\n"), "\n")
}

tic(sprintf("DONE.  Output in %s", normalizePath(OUT)))