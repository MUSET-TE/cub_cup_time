## =====================================================================
##  CUBpanel -- report layer.  Every function returns a character vector,
##  so the same text appears in the console, in the shiny panels and in
##  the downloadable log.  No computation happens here.
## =====================================================================

.rule <- function(title, ch = "=") {
  paste0(ch, ch, " ", title, " ", strrep(ch, max(3, 66 - nchar(title))))
}
.pr <- function(x, digits = 4, row.names = TRUE) {
  utils::capture.output(print(format(x, digits = digits),
                              row.names = row.names))
}
.kv <- function(k, v, w = 44) sprintf("%-*s %s", w, paste0(k, ":"), v)

## ---------------------------------------------------------------------
## block 1: data as detected
## ---------------------------------------------------------------------
cub_report_data <- function(an) {
  D <- an$D; i <- D$info; dd <- an$decomp
  c(.rule("Data as detected"),
    .kv("Respondent / wave / response columns",
        paste(D$roles, collapse = ", ")),
    .kv("Raw records", i$n_raw),
    .kv("Removed by listwise deletion", i$n_del),
    .kv("Interior gaps in the wave sequence",
        sprintf("%d%s", i$n_gap, if (i$n_gap) " (dropped)" else "")),
    .kv("Estimation rows (first wave reserved)", i$rows_est),
    .kv("Respondents N", D$N),
    .kv("Waves per respondent",
        sprintf("%d..%d  (%s)", i$minT, i$maxT,
                if (i$balanced) "balanced" else "unbalanced")),
    .kv("Scale m (declared)", D$m),
    .kv("Scale m resolved from max(response)", D$m_resolved),
    .kv("Observed support",
        sprintf("{%s}%s", paste(i$support, collapse = ","),
                if (length(i$unused))
                  sprintf("   categories never chosen: {%s}",
                          paste(i$unused, collapse = ",")) else "")),
    .kv("Declared shelter category c",
        sprintf("%d  (%s of the scale)", D$cs,
                if (D$cs == (D$m + 1) / 2) "midpoint"
                else if (D$cs %in% c(1L, D$m)) "endpoint" else "interior")),
    .kv("Lagged specifications available",
        if (i$lag_ok) "yes (max T_i >= 3)" else "NO: max T_i < 3"),
    "",
    "Additive between-within decomposition on the estimation rows:",
    sprintf("  total %.3f = between %.3f + within %.3f   (between share %.3f)",
            dd[["total"]], dd[["between"]], dd[["within"]],
            dd[["between_share"]]),
    sprintf("  moment-based intraclass correlation %.3f  (descriptive only)",
            dd[["icc"]]),
    "",
    "Formulas:",
    sprintf("  feeling     %s", deparse(D$formulas$feeling)),
    sprintf("  uncertainty %s", deparse(D$formulas$uncertainty)),
    sprintf("  shelter     %s", deparse(D$formulas$shelter)),
    "",
    "Lag centring constants (the estimator targets the CENTRED intercepts):",
    sprintf("  mean(lag_f) = %.4f   mean(lag_u) = %.4f",
            D$ctr[["lag_f"]], D$ctr[["lag_u"]]),
    "",
    sprintf("Correlated-random-effects block (%s):", D$cre$mode),
    sprintf("  feeling     -> {%s}", paste(colnames(D$Xm), collapse = ", ")),
    sprintf("  uncertainty -> {%s}", paste(colnames(D$Wm), collapse = ", ")),
    "  (a covariate with within-share below 0.01 contributes no mean,",
    "   because that mean would be collinear with the covariate itself)",
    "",
    "Within-respondent share of variation by covariate:",
    .pr(data.frame(within_share = round(D$wsh, 4))),
    "",
    "Relative frequencies by wave and category:",
    .pr(round(D$tab_wave, 4)),
    "",
    "Wave-to-wave transition counts (pooled):",
    utils::capture.output(print(D$tab_trans)))
}

## ---------------------------------------------------------------------
## block 2: dispersion diagnostic
## ---------------------------------------------------------------------
cub_report_dispersion <- function(an) {
  f <- an$floor; D <- an$D
  verdict <- if (f[["R"]] < 1)
    "R < 1: NO member of the family with positive mixing weights can be optimal."
  else "R >= 1: the family has room to work."
  c(.rule("Pre-fit dispersion diagnostic"),
    "Both extra components can only ADD dispersion at a fixed conditional",
    "mean, so the kernel variance (mu-1)(m-mu)/(m-1) is a FLOOR on what the",
    "family can represent.  One weighted least-squares solve decides it.",
    "",
    .kv("Mean dispersion floor  E[V_SB(mu_hat)]", sprintf("%.4f", f[["floor"]])),
    .kv("Conditional residual variance", sprintf("%.4f", f[["resid_var"]])),
    .kv("Dispersion ratio R", sprintf("%.4f", f[["R"]])),
    .kv("E|c - mu_hat|",
        sprintf("%.4f   (c = %d)", f[["dist_cs"]], D$cs)),
    "",
    paste0("VERDICT: ", verdict),
    if (f[["dist_cs"]] > (D$m - 1) / 3 && f[["R"]] < 1)
      paste0("A shelter declared far from the fitted mean compounds this: it ",
             "contributes (c - mu)^2 to the mixture variance.") else "")
}

## ---------------------------------------------------------------------
## block 3: self-checks
## ---------------------------------------------------------------------
cub_report_checks <- function(an) {
  if (is.null(an$checks)) return(c(.rule("Self-checks"), "not run"))
  k <- an$checks
  tb <- k$table
  tb$value <- sprintf("%10.3e", tb$value)
  tb$tolerance <- sprintf("%7.1e", tb$tolerance)
  tb$verdict <- ifelse(tb$pass, "OK",
                       ifelse(tb$gate, "** FAIL **", "see verdict below"))
  tb$gate <- ifelse(tb$gate, "gate", "diagnostic")
  tb$pass <- NULL
  c(.rule("Self-checks"),
    "Each row compares two computations of the same quantity that agree only",
    "if the implementation is correct.  C1-C6 and C8 are gates; C7 is a",
    "reported diagnostic whose verdict is anchored to the scale of the",
    "inference the quadrature supports.",
    "",
    .pr(tb, row.names = FALSE),
    "",
    sprintf("C7: per-respondent error of GH(Q) against a deterministic %s product",
            "301 x 301"),
    .pr(data.frame(Q = k$quad$Q,
                   err_per_subject = signif(k$quad$err_per_subject, 3),
                   implied_aggregate = signif(k$quad$implied_aggregate, 3)),
        row.names = FALSE),
    sprintf("  monotone in Q: %s", if (k$monotone) "yes" else "NO - INVESTIGATE"),
    sprintf("  N * err(Q = %d) = %.2e against a tolerance of %.0e",
            an$control$Q, k$implied_aggregate, k$Qtol),
    sprintf("  VERDICT: %s%s",
            if (k$Q_adequate) "Q is adequate"
            else sprintf("Q is NOT adequate; Q = %d is the smallest certified order",
                         k$Q_recommended),
            if (!k$Q_adequate)
              "   (the pessimistic bound; see the Quadrature block for the realised effect)"
            else ""),
    sprintf("  Gates passed: %s", if (k$ok) "all" else "NOT ALL"))
}

## ---------------------------------------------------------------------
## blocks 4-5: ladder and chains
## ---------------------------------------------------------------------
cub_report_ladder <- function(an) {
  L <- an$ladder
  out <- c(.rule("Specification ladder"),
           sprintf("All rungs fitted to an identical estimation sample (%d rows, N = %d): %s",
                   an$D$n, an$D$N,
                   if (L$same_sample) "verified by hash" else "** MISMATCH **"),
           "",
           .pr(L$table, digits = 6),
           "",
           .kv("BIC_N selects", an$fits[[L$best_bicN]]$label),
           .kv("AIC selects", an$fits[[L$best_aic]]$label))
  fr <- L$table$frontier != "-"
  if (any(fr))
    out <- c(out, "",
             "Frontier flags were raised.  Where they are, the simulated",
             "chi-bar-squared reference is withheld and the delta-method",
             "standard errors for sigma_12 and rho_12 are suppressed.")
  gains <- c()
  if (all(c("0", "RE") %in% names(an$fits)))
    gains <- c(gains, sprintf("  heterogeneity over the pooled rung: %+.1f log-likelihood units",
                              an$fits[["RE"]]$ll - an$fits[["0"]]$ll))
  if (all(c("RE", "RE-LAG") %in% names(an$fits)))
    gains <- c(gains, sprintf("  dynamics over the RE rung:          %+.1f",
                              an$fits[["RE-LAG"]]$ll - an$fits[["RE"]]$ll))
  if (all(c("RE-LAG", "CRE-LAG") %in% names(an$fits)))
    gains <- c(gains, sprintf("  CRE block over RE-LAG:              %+.1f",
                              an$fits[["CRE-LAG"]]$ll - an$fits[["RE-LAG"]]$ll))
  if (length(gains)) out <- c(out, "", "Increments:", gains)
  out
}

cub_report_chains <- function(an) {
  one <- function(ch, ttl) {
    c(ttl, .pr(ch$table, digits = 6),
      sprintf("  smallest step %+.2e | audit at 1e-6: %s | audit at %.1e: %s",
              ch$min_step,
              if (ch$audit_strict) "pass" else "FLAGGED",
              ch$nest_tol,
              if (ch$audit_lr) "pass" else "FLAGGED"), "")
  }
  out <- c(.rule("Nested chains and the nesting audit"),
           "A nested model attaining a higher log-likelihood than the model",
           "nesting it is an optimisation failure.  The audit is reported at a",
           "strict tolerance of 1e-6 and at one anchored to the scale of the",
           "inference, max(1e-2, N * err(Q)).  On a ladder that is flat to",
           "within the quadrature error the strict audit fires on differences",
           "that are not failures, which is why both are printed.", "",
           one(an$chain1, "Chain I: heterogeneity before dynamics"))
  if (!is.null(an$chain2))
    out <- c(out, one(an$chain2, "Chain II: dynamics before heterogeneity"))
  if (!is.null(an$chain2)) {
    d1 <- diff(an$chain1$table$loglik); d2 <- diff(an$chain2$table$loglik)
    out <- c(out,
             "The two chains decompose the SAME total differently.  Observed",
             "persistence and permanent heterogeneity are confounded at these",
             "panel lengths; the branch does not resolve the confounding, it",
             "makes the size of it visible:",
             sprintf("  chain I  increments: %s",
                     paste(sprintf("%+.1f", d1), collapse = "  ")),
             sprintf("  chain II increments: %s",
                     paste(sprintf("%+.1f", d2), collapse = "  ")))
  }
  out
}

## ---------------------------------------------------------------------
## block 6: correlation audit
## ---------------------------------------------------------------------
cub_report_rho <- function(an) {
  if (is.null(an$rho)) return(c(.rule("Correlation audit"), "not run"))
  r <- an$rho
  out <- c(.rule("Correlation audit: H0 rho_12 = 0"),
           "Because Sigma is parameterised through its Cholesky factor,",
           "rho_12 = 0 is equivalent to L21 = 0, an INTERIOR point.  The",
           "restriction is therefore an ordinary chi^2_1 likelihood ratio.",
           "",
           .kv("rho_12 estimate", sprintf("%+.4f", r$rho)),
           .kv("LR against the diagonal restriction",
               sprintf("%.4f on 1 df, p = %s", r$lr, cub_pfmt(r$p))))
  if (!is.null(r$prof)) {
    out <- c(out, "", "Profile (every other parameter re-optimised):",
             .pr(r$prof, digits = 5, row.names = FALSE),
             sprintf("  95%% profile set on the evaluated grid: [%s, %s]",
                     cub_num(r$ci[1], 2), cub_num(r$ci[2], 2)))
    if (isTRUE(r$truncated))
      out <- c(out,
               sprintf("  the grid was truncated by the time budget: %d of %d nodes evaluated",
                       nrow(r$prof), length(r$grid)))
  }
  c(out, "",
    paste0("VERDICT: ", if (r$informative)
      "rho_12 is informative; the unrestricted Sigma is reported."
      else paste0("UNINFORMATIVE CORNER.  The diagonal specification is ",
                  "implied, not because the correlation is zero but because ",
                  "nothing has been learned about it.")))
}

## ---------------------------------------------------------------------
## block 7: coefficients
## ---------------------------------------------------------------------
cub_report_coef <- function(an, which = NULL) {
  if (!length(an$coef)) return(c(.rule("Coefficients"), "not computed"))
  if (is.null(which)) which <- names(an$coef)
  out <- c(.rule("Coefficients, respondent-clustered standard errors"),
           "A positive coefficient means more of the quantity named in its own",
           "equation.  AME is the average derivative of that quantity on the",
           "(0,1) scale; AME_ER is the average derivative of E(R_it) on the",
           "category scale.  BOTH are reported because they routinely disagree",
           "in sign: raising upsilon pulls the conditional mean towards the",
           "scale midpoint, so the sign of the mean-scale effect depends on",
           "which side of the midpoint the fitted mean lies.", "")
  for (nm in which) {
    cc <- an$coef[[nm]]
    out <- c(out, sprintf("-- %s --", cc$label),
             .pr(cc$tab, digits = 4),
             sprintf("  cond(H) = %.4g | H positive definite: %s | max|sum_i score_i| = %.2e | KKT informative: %s",
                     cc$kappa, cc$pd, cc$score_sum, cc$kkt_informative))
    if (!is.null(cc$vcn) && isTRUE(cc$vcn$boundary))
      out <- c(out,
               "  NOTE: a variance component is at the frontier, so sigma_12",
               "  and rho_12 are not identified and their standard errors are",
               "  withheld rather than printed with a disclaimer.")
    if (is.finite(cc$cover))
      out <- c(out, sprintf(
        "  %d of %d comparable parameters have nominal 95%% intervals containing the generating value.",
        cc$cover, cc$ncomparable),
        "  This is a COUNT ON ONE REALISATION, not coverage.  An interval wide",
        "  enough covers by construction; see the replication driver.")
    if (nm == "CRE-LAG")
      out <- c(out,
               "  NOTE: in this rung the variance components are CONDITIONAL on",
               "  (y_i1, xbar), so they are not comparable with marginal",
               "  generating values and the True column is left blank.")
    out <- c(out, "")
  }
  out
}

## ---------------------------------------------------------------------
## block 8: profile interval
## ---------------------------------------------------------------------
cub_report_profile <- function(an) {
  if (is.null(an$profile_s1))
    return(c(.rule("Profile interval"), "not computed"))
  p <- an$profile_s1
  out <- c(.rule("Profile interval for sigma^2_1"),
           sprintf("  %s   (%.0f s)", cub_fmt_profile(p), p$secs))
  if (!p$hi_ok || !p$lo_ok)
    out <- c(out,
             "  A truncated profile is reported as a one-sided limit rather",
             "  than closed by extrapolation.")
  out <- c(out,
           "  sigma^2_2 and rho_12: delta-method intervals from the",
           "  coefficient table (an INTERIOR optimum in the Cholesky",
           "  parameterisation); the rho_12 profile is in the audit above.")
  if (!is.null(an$profile_s1_grid))
    out <- c(out, "", "Profile grid:",
             .pr(an$profile_s1_grid, digits = 5, row.names = FALSE))
  out
}

## ---------------------------------------------------------------------
## block 9: tests
## ---------------------------------------------------------------------
cub_report_tests <- function(an) {
  tb <- an$tests
  pr <- data.frame(
    hypothesis = tb$hypothesis, reference = tb$reference,
    statistic = ifelse(is.finite(tb$statistic),
                       sprintf("%.3f", tb$statistic), "---"),
    df = ifelse(is.finite(tb$df), as.character(tb$df), "---"),
    p = vapply(seq_len(nrow(tb)), function(i)
      if (grepl("WITHHELD", tb$reference[i])) "withheld"
      else cub_pfmt(tb$p[i],
                    if (grepl("B = ", tb$reference[i]))
                      an$control$chibar_B else NA_integer_), ""),
    note = tb$note, stringsAsFactors = FALSE)
  c(.rule("Test battery"),
    .pr(pr, row.names = FALSE),
    "",
    "Reading notes.",
    "  Wald statistics use the respondent-clustered covariance; no Wald",
    "  statistic is ever applied to a variance component.",
    "  H0: Sigma = 0 lies on the boundary of the cone of positive",
    "  semi-definite matrices with a third parameter unidentified under it.",
    "  The simulated chi-bar-squared reference is obtained by exact",
    "  projection onto that cone and is WITHHELD when the plug-in",
    "  information comes from a boundary solution, because such a p-value",
    "  has no known relation to the hypothesis.",
    "  H0: rho_f = rho_u is not scale-free in magnitude, since the two",
    "  features carry different units; a rejection in which the estimates",
    "  carry OPPOSITE SIGNS is scale-free, and that is what to look for.",
    "  The CRE block is tested jointly AND in halves, because the Mundlak",
    "  means test random-effects exogeneity while the initial-condition",
    "  terms test the conditioning on the first wave.  A joint rejection",
    "  does not say which failed.")
}

## ---------------------------------------------------------------------
## blocks 10-11: fit and quadrature
## ---------------------------------------------------------------------
cub_report_fit <- function(an) {
  g <- an$stats_selected; D <- an$D
  tab <- rbind(observed = g$obs, fitted = g$fit, difference = g$obs - g$fit)
  colnames(tab) <- seq_len(D$m)
  c(.rule("Fit of the selected rung"),
    .kv("Rung", an$fit_selected$label),
    .kv("Predictive quadrature order", g$Qpred),
    .kv("max |sum_r p_itr - 1|", sprintf("%.1e", g$rowsum)),
    .kv("Pooled normalised dissimilarity", sprintf("%.4f", g$diss)),
    .kv("Mean log-score", sprintf("%.4f", g$logscore)),
    "",
    .pr(round(tab, 4)),
    "",
    "Dissimilarity by wave:",
    paste0("  ", paste(sprintf("wave %s: %.4f", names(g$diss_wave),
                               g$diss_wave), collapse = "   ")),
    "",
    sprintf("Variance account: observed %.3f vs implied %.3f = E[Var] %.3f + Var[E] %.3f (ratio %.3f)",
            g$var_obs, g$EV + g$VE, g$EV, g$VE, (g$EV + g$VE) / g$var_obs),
    "",
    "Binned residuals over bins of the fitted conditional mean:",
    .pr(g$binres, digits = 3, row.names = FALSE),
    sprintf("  %d of %d bins outside +-2/sqrt(n_b)", g$out_band, g$nbins))
}

cub_report_quad <- function(an) {
  q <- an$quad
  out <- c(.rule("Quadrature order"),
           "The integral is evaluated on a FIXED Cholesky-scaled",
           "Gauss-Hermite product grid.  Adaptation is deliberately not used:",
           "the conditional log-integrand flattens as the uncertainty weight",
           "approaches its boundary, and an adaptive rule then rescales the",
           "grid by a near-singular curvature, which is the one error that",
           "destroys the comparability of log-likelihoods along a ladder.",
           "",
           "Log-likelihood at the retained parameter vector, over Q:",
           .pr(q$table, digits = 7, row.names = FALSE),
           sprintf("  |ll(Q = %d) - ll(Q = %d)| = %.2e against a tolerance of %.0e",
                   an$control$Q, max(q$table$Q), q$drift, q$tol),
           sprintf("  VERDICT: %s", if (q$adequate) "Q is sufficient"
                   else "Q is insufficient"))
  if (isTRUE(q$refitted))
    out <- c(out,
             sprintf("  The selected rung was REFITTED at Q = %d: log-likelihood moved from %.4f to %.4f (%+.4f).",
                     q$Q_used, q$ll_before, q$ll_after, q$ll_change),
             "  The pessimistic N*err bound is not realised at the optimum;",
             "  both numbers are printed so that the gap is visible.")
  c(out, "",
    "Recommended practice, which the software supports: certify first, fit",
    "once.  Certifying after fitting leaves a discrepancy of the size above.")
}

## ---------------------------------------------------------------------
## block 12: subject effects
## ---------------------------------------------------------------------
cub_report_eb <- function(an) {
  e <- an$eb
  rows <- rbind(`theta estimated` = e$estimated)
  if (!is.null(e$oracle)) rows <- rbind(rows, `theta ORACLE` = e$oracle)
  out <- c(.rule("Subject effects (empirical Bayes)"),
           .pr(round(rows, 3)),
           "",
           sprintf("sigma_1 estimated %.3f%s | sigma_2 estimated %.3f%s",
                   e$sigma[["sigma_1"]],
                   if (!is.null(e$sigma_true))
                     sprintf(" vs true %.3f", e$sigma_true[["sigma_1"]]) else "",
                   e$sigma[["sigma_2"]],
                   if (!is.null(e$sigma_true))
                     sprintf(" vs true %.3f", e$sigma_true[["sigma_2"]]) else ""),
           "",
           "Reading notes.",
           "  sd(posterior means) < sigma is SHRINKAGE, not error.",
           "  The posterior means are an ORDERING device, not estimates of",
           "  individual effects, and no p-value is attached to any of them.")
  if (!is.null(e$ll_oracle))
    out <- c(out, sprintf(
      "  ll(oracle) = %.2f vs ll(ML) = %.2f: the ML value must be the larger, and is (%s).",
      e$ll_oracle, e$ll_ml, e$ll_ml >= e$ll_oracle),
      "  Compare the two rows above before reading any shrinkage as a defect:",
      "  if the oracle row shows the same pattern, it is a property of the",
      "  design and not of the estimator.")
  out
}

## ---------------------------------------------------------------------
## block 13: out-of-sample
## ---------------------------------------------------------------------
cub_report_oos <- function(an) {
  if (is.null(an$oos))
    return(c(.rule("Out-of-sample"),
             paste0("not reported: ",
                    if (!is.null(an$oos_note)) an$oos_note else "disabled")))
  o <- an$oos
  c(.rule("Out-of-sample: held-out final wave"),
    "Each respondent's final wave is held out and every model is refitted on",
    "the remainder.  Two precautions: the lagged features and the Mundlak",
    "means are recentred and recomputed on the TRAINING rows only, and both",
    "benchmarks receive the union of the regressors of all three equations.",
    "",
    sprintf("Test rows: %d.  Common information set (%d regressors):",
            o$n_test, length(o$regressors)),
    paste0("  ", paste(o$regressors, collapse = ", ")),
    "",
    .pr(round(o$scores, 4)),
    "",
    sprintf("Training log-likelihoods: %s",
            paste(sprintf("%s %.2f (p = %d)", names(o$ll), o$ll, o$p),
                  collapse = " | ")),
    sprintf("Metric winners against the proportional-odds benchmark: %s",
            paste(sprintf("%s = %s", names(o$winner), o$winner),
                  collapse = " | ")),
    sprintf("VERDICT: %s (%d of 4 metrics; the designated arbiter, RPS, favours %s)",
            o$verdict, o$n_win, o$arbiter),
    "",
    "The case for this family is not a predictive case.  A proportional-odds",
    "model with the same information set predicts the next rating about as",
    "well; the extra structure is spent on DECOMPOSING the conditional",
    "distribution rather than on sharpening it.")
}

## ---------------------------------------------------------------------
## block 14: memory profile and collapse audit
## ---------------------------------------------------------------------
cub_report_memory <- function(an) {
  if (is.null(an$memory))
    return(c(.rule("Memory profile"),
             "unavailable: the fitted rung carries no lagged features"))
  fd <- an$fits[[an$spec_dyn]]
  rf <- if ("xi:lag_f" %in% names(fd$th)) unname(fd$th["xi:lag_f"]) else NA_real_
  ru <- if ("ups:lag_u" %in% names(fd$th)) unname(fd$th["ups:lag_u"]) else NA_real_
  c(.rule("Memory profile"),
    "The feeling feature is a signed deviation from the scale midpoint and is",
    "monotone in the previous rating.  The uncertainty feature is",
    "ECCENTRICITY: distance from the midpoint in half-range units.  It is",
    "V-shaped and takes equal values at r and m+1-r, so a single coefficient",
    "cannot summarise what the uncertainty dynamics do, and rho_u < 0 means",
    "that a CENTRAL previous answer predicts higher uncertainty now.",
    "",
    .pr(an$memory, digits = 4, row.names = FALSE),
    "",
    sprintf("Logit shift from the lowest to the highest previous rating, feeling:      %+.3f",
            2 * rf),
    sprintf("Logit shift from a midpoint to an endpoint previous rating, uncertainty: %+.3f",
            ru),
    "Report the coefficient, build the discussion on the profile.")
}

cub_report_verify <- function(an) {
  v <- an$verif
  out <- c(.rule("Collapse audit"),
           "These identities are available whichever regime the data fall in.",
           "V1 is computed by stats::glm.fit and uses nothing from the mixture",
           "code, so it is an independent confirmation.",
           "",
           .pr(v$table, row.names = FALSE),
           "",
           sprintf("From CUBpanel(0): mean upsilon %.3e (max %.3e) | mean delta %.3e (max %.3e)",
                   v$mean_ups, v$max_ups, v$mean_del, v$max_del),
           sprintf("Model-implied total variance %.3f against an observed %.3f (ratio %.3f, sd ratio %.3f)",
                   v$implied_total, v$implied_total / v$V4, v$V4, v$sd_ratio))
  if (isTRUE(v$collapsed))
    out <- c(out, "",
             "VERDICT: the mixture has COLLAPSED onto the boundary of the",
             "family.  The model that fitted these data is a shifted-binomial",
             "regression with the mixing parameters wasted.  The weights are at",
             "the numerical floor, not merely small, and the collapse is exact",
             "rather than approximate.  Do not report the uncertainty or",
             "shelter coefficients of a collapsed rung: they are slopes on",
             "components of negligible weight, their Wald statistics test nulls",
             "no analyst holds, and their marginal effects are zero.")
  else
    out <- c(out, "",
             "VERDICT: no collapse.  All three components carry weight.")
  if (!is.null(v$cold))
    out <- c(out, "",
             sprintf("Cold-start audit: %s attained %.4f from moment-based starts against %.4f warm-started.",
                     v$cold$label, v$cold$ll,
                     an$fits[[if ("CRE-LAG" %in% names(an$fits)) "CRE-LAG"
                              else an$spec_dyn]]$ll),
             "A warm start that happens to sit in a basin the cold start",
             "escapes is exactly what an audit should detect, and comparing",
             "warm against warm cannot detect it.")
  out
}

## ---------------------------------------------------------------------
## block 15: export
## ---------------------------------------------------------------------

#' Every reported scalar as a key-value list, for a manuscript
cubpanel_export <- function(an, prefix = "X") {
  D <- an$D; dd <- an$decomp; f <- an$floor; v <- an$verif
  P <- list()
  put <- function(k, val) P[[paste0(prefix, ":", k)]] <<- as.character(val)
  put("N", D$N); put("T", D$info$maxT); put("rows", D$n)
  put("m", D$m); put("cs", D$cs)
  put("support", sprintf("$\\{%s\\}$", paste(range(D$info$support),
                                             collapse = ",\\dots,")))
  put("ndel", D$info$n_del); put("ngap", D$info$n_gap)
  put("vtot", cub_num(dd[["total"]])); put("vb", cub_num(dd[["between"]]))
  put("vw", cub_num(dd[["within"]]))
  put("bshare", cub_num(dd[["between_share"]]))
  put("icc", cub_num(dd[["icc"]]))
  put("R", cub_num(f[["R"]])); put("floor", cub_num(f[["floor"]]))
  put("residvar", cub_num(f[["resid_var"]]))
  put("distcs", cub_num(f[["dist_cs"]]))
  put("ctrf", cub_num(D$ctr[["lag_f"]], 4))
  put("ctru", cub_num(D$ctr[["lag_u"]], 4))
  for (nm in names(an$fits)) {
    sfx <- tolower(gsub("[^A-Za-z0-9]", "", nm))
    put(paste0("ll", sfx), cub_num(an$fits[[nm]]$ll, 2))
    put(paste0("aic", sfx), cub_num(an$fits[[nm]]$aic, 2))
    put(paste0("bic", sfx), cub_num(an$fits[[nm]]$bicN, 2))
    put(paste0("diss", sfx), cub_num(an$diss[[nm]], 4))
    put(paste0("p", sfx), an$fits[[nm]]$p)
  }
  put("selbicN", an$fits[[an$selected$bicN]]$label)
  put("selaic", an$fits[[an$selected$aic]]$label)
  if (!is.null(an$rho)) {
    put("rhohat", cub_num(an$rho$rho, 3)); put("rholr", cub_num(an$rho$lr, 3))
    put("rhop", cub_ptex(an$rho$p))
    put("rhoverdict", if (an$rho$informative) "informative"
        else "uninformative corner")
  }
  if (length(an$coef)) {
    cc <- an$coef[[an$primary_coef]]
    put("cond", sprintf("$%.3g$", cc$kappa))
    tb <- cc$tab
    for (j in seq_len(nrow(tb))) {
      key <- gsub("[^A-Za-z0-9]", "", rownames(tb)[j])
      if ("True" %in% names(tb))
        put(paste0("t", key), if (is.finite(tb$True[j]))
          sprintf("$%.3f$", tb$True[j]) else "---")
      put(paste0("e", key), sprintf("$%.3f$", tb$est[j]))
      put(paste0("s", key), if (is.finite(tb$se[j]))
        sprintf("$%.3f$", tb$se[j]) else "$\\dagger$")
      if ("AME" %in% names(tb)) {
        put(paste0("a", key), if (is.finite(tb$AME[j]))
          sprintf("$%.4f$", tb$AME[j]) else "---")
        put(paste0("q", key), if (is.finite(tb$p_AME[j]))
          cub_ptex(tb$p_AME[j]) else "---")
        put(paste0("y", key), if (is.finite(tb$AME_ER[j]))
          sprintf("$%.4f$", tb$AME_ER[j]) else "---")
      }
    }
  }
  tt <- an$tests
  for (j in seq_len(nrow(tt))) {
    key <- gsub("[^A-Za-z0-9]", "", tt$hypothesis[j])
    put(paste0("w", key), cub_num(tt$statistic[j], 3))
    put(paste0("qq", key), cub_ptex(tt$p[j]))
  }
  g <- an$stats_selected
  put("dissfit", cub_num(g$diss, 4)); put("nbins", g$nbins)
  put("outband", g$out_band); put("logscore", cub_num(g$logscore, 4))
  put("impliedtot", cub_num(g$EV + g$VE))
  put("varratio", cub_num(v$V4)); put("sdratio", cub_num(v$sd_ratio))
  put("V1", sprintf("%.1e", v$V1)); put("V2", sprintf("%.0e", v$V2))
  put("V3", sprintf("%.0e", v$V3))
  put("V5", if (isTRUE(v$V5)) "yes" else if (isFALSE(v$V5)) "no" else "n/a")
  if (!is.null(an$oos)) {
    put("rps", cub_num(an$oos$scores[1, "RPS"], 4))
    put("rpspo", cub_num(an$oos$scores[2, "RPS"], 4))
    put("rpssb", cub_num(an$oos$scores[3, "RPS"], 4))
    put("oosverdict", an$oos$verdict)
  }
  if (!is.null(an$profile_s1)) {
    put("s1est", cub_num(an$profile_s1$est, 4))
    put("s1lo", cub_num(an$profile_s1$lo, 4))
    put("s1hi", sprintf("%s%s", if (an$profile_s1$hi_ok) "" else ">",
                        cub_num(an$profile_s1$hi, 4)))
  }
  P
}

cub_report_export <- function(an, prefix = "X") {
  P <- cubpanel_export(an, prefix)
  c(.rule("Export: every reported scalar as a key-value list"),
    sprintf("  %-28s = %s", names(P), unlist(P)))
}

## LaTeX table bodies, ready to paste -----------------------------------
cub_latex_ladder <- function(an) {
  tb <- an$ladder$table
  c("% ladder body",
    vapply(seq_len(nrow(tb)), function(i) sprintf(
      "%-34s & %2d & $%.2f$ & $%.2f$ & $%.2f$ & $%.4f$ \\\\",
      tb$spec[i], tb$p[i], tb$loglik[i], tb$AIC[i], tb$BIC_N[i], tb$Diss[i]),
      ""))
}

cub_latex_coef <- function(an, which = NULL) {
  if (is.null(which)) which <- an$primary_coef
  tb <- an$coef[[which]]$tab
  has_true <- "True" %in% names(tb)
  c(sprintf("%% coefficient body for %s", an$coef[[which]]$label),
    vapply(seq_len(nrow(tb)), function(i) sprintf(
      "%-26s &%s $%.3f$ & %s & %s & %s & %s \\\\",
      rownames(tb)[i],
      if (has_true) paste0(" ", if (is.finite(tb$True[i]))
        sprintf("$%.3f$", tb$True[i]) else "---", " &") else "",
      tb$est[i],
      if (is.finite(tb$se[i])) sprintf("$%.3f$", tb$se[i]) else "$\\dagger$",
      if (!is.null(tb$AME) && is.finite(tb$AME[i]))
        sprintf("$%.4f$", tb$AME[i]) else "---",
      if (!is.null(tb$p_AME) && is.finite(tb$p_AME[i]))
        cub_ptex(tb$p_AME[i]) else "---",
      if (!is.null(tb$AME_ER) && is.finite(tb$AME_ER[i]))
        sprintf("$%.4f$", tb$AME_ER[i]) else "---"), ""))
}

cub_latex_tests <- function(an) {
  tb <- an$tests
  c("% test battery body",
    vapply(seq_len(nrow(tb)), function(i) sprintf(
      "%-34s & %s & $%.3f$ & %s & %s \\\\",
      tb$hypothesis[i], tb$reference[i], tb$statistic[i],
      if (is.finite(tb$df[i])) sprintf("$%d$", as.integer(tb$df[i])) else "---",
      cub_ptex(tb$p[i])), ""))
}

## ---------------------------------------------------------------------
## the whole report
## ---------------------------------------------------------------------
cubpanel_report <- function(an, blocks = NULL, prefix = "X") {
  all_blocks <- c(
    data = "cub_report_data", dispersion = "cub_report_dispersion",
    checks = "cub_report_checks", ladder = "cub_report_ladder",
    chains = "cub_report_chains", rho = "cub_report_rho",
    coef = "cub_report_coef", profile = "cub_report_profile",
    tests = "cub_report_tests", fit = "cub_report_fit",
    quad = "cub_report_quad", eb = "cub_report_eb",
    oos = "cub_report_oos", memory = "cub_report_memory",
    verify = "cub_report_verify")
  if (is.null(blocks)) blocks <- names(all_blocks)
  out <- c(
    .rule("CUBpanel analysis report", "#"),
    sprintf("  %s | R %s | %s", format(an$stamp$when, "%Y-%m-%d %H:%M:%S"),
            an$stamp$R, an$stamp$platform),
    sprintf("  elapsed %.1f s | estimation quadrature Q = %d | predictive Q = %d",
            an$secs, an$control$Q, an$control$Qpred),
    "")
  for (b in blocks)
    out <- c(out, do.call(all_blocks[[b]], list(an)), "")
  out <- c(out, cub_report_export(an, prefix), "",
           .rule("LaTeX bodies", "#"), "",
           cub_latex_ladder(an), "",
           cub_latex_coef(an), "",
           cub_latex_tests(an))
  out
}

cub_report_mc <- function(mcraw) {
  if (is.null(mcraw))
    return(c(.rule("Replication study"), "no usable replications"))
  s <- cub_mc_summary(mcraw)
  ach <- max(s$param$R_ok)
  secs <- attr(mcraw, "secs"); nc <- attr(mcraw, "ncores")
  per <- secs / max(ach, 1) * nc
  Rmin <- attr(mcraw, "R_min")
  if (ach < Rmin)
    return(c(.rule("Replication study: NOT REPORTED"),
             sprintf("%d usable replications against a hard floor of %d.",
                     ach, Rmin),
             "",
             "A table of bias, root mean squared error and coverage computed on",
             "a handful of replications is not a weak version of a simulation",
             "study, it is a different and misleading object.  The floor exists",
             "to remove the temptation to print it anyway.",
             "",
             sprintf("Observed cost: %.0f s per replication per core, %d cores.",
                     per, nc),
             sprintf("  R = %d: %.1f core-hours, %.1f h wall-clock.",
                     Rmin, per * Rmin / 3600, per * Rmin / nc / 3600),
             sprintf("  R = 500: %.1f core-hours, %.1f h wall-clock.",
                     per * 500 / 3600, per * 500 / nc / 3600)))
  pa <- s$param
  for (cc in c("bias", "rmse", "sd_est", "mean_se", "mean_est"))
    pa[[cc]] <- round(pa[[cc]], 4)
  pa$cover95 <- round(pa$cover95, 3)
  c(.rule("Replication study"),
    sprintf("%d usable replications in %.0f s (%.0f s per replication per core).",
            ach, secs, per),
    "", .pr(pa, row.names = FALSE), "",
    "Diagnostics:", .pr(s$diag, row.names = FALSE), "",
    "Reading notes.",
    "  cover95 near 0.95 validates the clustered standard errors.",
    "  frontier_rate is how often Sigma reaches the boundary.",
    "  reject_rho12_5pct is the POWER of the chi^2_1 test, L21 = 0 being",
    "  interior.  The sigma^2 and rho_12 cells for CRE-LAG are NA because",
    "  there the variances are conditional on (y_i1, xbar).",
    sprintf("  Design: covariates were %s across replications.",
            if (isTRUE(attr(mcraw, "fixed_design"))) "held FIXED"
            else "re-generated"))
}