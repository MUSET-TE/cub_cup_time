## =====================================================================
##  CUBpanel -- figures.  Base graphics only; every function takes
##  explicit arguments, so the same code backs the shiny panels, the
##  manuscript figures and any user script.
##  Greyscale throughout: the figures must survive a print journal.
## =====================================================================

CUB_GREY <- grey(c(0.15, 0.35, 0.50, 0.62, 0.72, 0.80, 0.86))

cub_par <- function(...) {
  graphics::par(family = "sans", mar = c(4.1, 4.3, 2.4, 1.1),
                mgp = c(2.5, 0.7, 0), tcl = -0.3, cex.axis = 0.88,
                cex.lab = 0.95, cex.main = 0.98, bty = "l", las = 1, ...)
}

#' Write a figure to PNG at publication resolution
cub_png <- function(expr, file, width = 2200, height = 1100, res = 180) {
  grDevices::png(file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  cub_par()
  force(expr)
  invisible(file)
}

.cub_wave_cols <- function(nw) {
  if (nw <= length(CUB_GREY)) CUB_GREY[seq_len(nw)]
  else grey(seq(0.12, 0.86, length.out = nw))
}

## ---------------------------------------------------------------------
## F1  wave profiles, before any model is fitted
## ---------------------------------------------------------------------
cub_plot_waves <- function(D, R = NA_real_, main = NULL) {
  cub_par()
  tw <- D$tab_wave
  waves <- rownames(tw); nw <- length(waves)
  cols <- .cub_wave_cols(nw)
  if (is.null(main))
    main <- sprintf("Observed frequencies by wave:  N = %d,  T = %d%s",
                    D$N, D$info$maxT,
                    if (is.finite(R)) sprintf(",  R = %.3f", R) else "")
  graphics::matplot(seq_len(D$m), t(tw), type = "o", pch = 16, cex = 0.7,
                    lty = 1, lwd = 1.5, col = cols, xaxt = "n",
                    ylim = c(0, max(tw) * 1.14),
                    xlab = "Response category", ylab = "Relative frequency",
                    main = main)
  graphics::axis(1, at = seq_len(D$m))
  graphics::abline(v = D$cs, lty = 3, col = "black", lwd = 1.4)
  graphics::text(D$cs, max(tw) * 1.11, sprintf("declared c = %d", D$cs),
                 pos = if (D$cs > D$m / 2) 2 else 4, cex = 0.78)
  if (length(D$info$unused)) {
    graphics::points(D$info$unused, rep(0, length(D$info$unused)),
                     pch = 4, cex = 1.3, lwd = 2)
    graphics::text(D$info$unused, 0, "never chosen", pos = 4, cex = 0.75,
                   offset = 0.6)
  }
  graphics::legend("topright", legend = paste("wave", waves), col = cols,
                   lty = 1, pch = 16, bty = "n", cex = 0.72, pt.cex = 0.7,
                   ncol = if (nw > 4) 2 else 1)
  invisible(NULL)
}

#' Two panels side by side (manuscript Figure 1)
cub_plot_waves_pair <- function(D1, R1, lab1, D2, R2, lab2) {
  op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
  cub_plot_waves(D1, R1, sprintf("%s:  N = %d,  T = %d,  R = %.3f",
                                 lab1, D1$N, D1$info$maxT, R1))
  cub_plot_waves(D2, R2, sprintf("%s:  N = %d,  T = %d,  R = %.3f",
                                 lab2, D2$N, D2$info$maxT, R2))
}

## ---------------------------------------------------------------------
## F2  within-respondent share of variation in each covariate
## ---------------------------------------------------------------------
cub_plot_within <- function(D) {
  cub_par(mar = c(4.1, 10.5, 2.4, 1.1))
  w <- sort(D$wsh)
  bp <- graphics::barplot(w, horiz = TRUE, las = 1, col = grey(0.55),
                          border = NA, xlim = c(0, max(c(w, 0.6)) * 1.18),
                          xlab = "within-respondent share of variation",
                          main = "Covariates: within versus between identification")
  graphics::text(w, bp, sprintf(" %.3f", w), pos = 4, cex = 0.8, xpd = NA)
  graphics::abline(v = 0.01, lty = 2)
  graphics::text(0.01, max(bp), "0.01: Mundlak-mean threshold", pos = 4,
                 cex = 0.74, srt = 0)
  invisible(NULL)
}

## ---------------------------------------------------------------------
## F3  dispersion floor of the family, with panels overlaid
## ---------------------------------------------------------------------

#' @param series list of list(label=, D=, mu=, R=, pch_obs=, pch_floor=, col=)
cub_plot_dispersion <- function(series, m = NULL, k = 8L) {
  cub_par()
  if (!length(series)) { graphics::plot.new(); return(invisible(NULL)) }
  if (is.null(m)) m <- series[[1]]$D$m
  oct <- function(y, mu, mm) {
    br <- unique(stats::quantile(mu, seq(0, 1, length.out = k + 1L)))
    g <- cut(mu, br, include.lowest = TRUE)
    data.frame(mu = tapply(mu, g, mean),
               obs = tapply((y - mu)^2, g, mean),
               flo = tapply((mu - 1) * (mm - mu) / (mm - 1), g, mean))
  }
  O <- lapply(series, function(s) oct(s$D$y, s$mu, s$D$m))
  mu <- seq(1.05, m - 0.05, length.out = 400)
  flo <- (mu - 1) * (m - mu) / (m - 1)
  ymax <- max(c(flo, unlist(lapply(O, function(o) o$obs))), na.rm = TRUE) * 1.08
  graphics::plot(mu, flo, type = "l", lwd = 2.4, col = "black",
                 ylim = c(0, ymax),
                 xlab = expression("fitted conditional mean  " * hat(mu)[it]),
                 ylab = "conditional dispersion",
                 main = bquote("Dispersion floor " * V[SB](mu) *
                                 " of the family at m = " * .(m)))
  graphics::lines(mu, rep((m^2 - 1) / 12, length(mu)), lty = 2,
                  col = grey(0.55), lwd = 1.3)
  graphics::text(mu[20], (m^2 - 1) / 12, "uniform component variance",
                 pos = 3, cex = 0.74, col = grey(0.35))
  leg_txt <- expression(V[SB](mu) * ", the floor")
  leg_pch <- NA; leg_lty <- 1; leg_lwd <- 2.4; leg_col <- "black"
  for (j in seq_along(series)) {
    s <- series[[j]]; o <- O[[j]]
    col <- if (!is.null(s$col)) s$col else grey(0.10 + 0.32 * (j - 1))
    po <- if (!is.null(s$pch_obs)) s$pch_obs else c(16, 17, 15)[j]
    pf <- if (!is.null(s$pch_floor)) s$pch_floor else c(1, 2, 0)[j]
    graphics::segments(o$mu, o$flo, o$mu, o$obs, col = col, lwd = 1)
    graphics::points(o$mu, o$obs, pch = po, cex = 1.25, col = col)
    graphics::points(o$mu, o$flo, pch = pf, cex = 1.25, col = col)
    leg_txt <- c(leg_txt,
                 sprintf("%s observed  (R = %.3f)", s$label, s$R),
                 sprintf("%s floor, by octile", s$label))
    leg_pch <- c(leg_pch, po, pf); leg_lty <- c(leg_lty, NA, NA)
    leg_lwd <- c(leg_lwd, NA, NA); leg_col <- c(leg_col, col, col)
  }
  graphics::legend("topleft", bty = "n", cex = 0.80, legend = leg_txt,
                   pch = leg_pch, lty = leg_lty, lwd = leg_lwd, col = leg_col)
  invisible(NULL)
}

## ---------------------------------------------------------------------
## F4  fitted mixing weights on a log10 axis (the collapse figure)
## ---------------------------------------------------------------------
cub_plot_weights <- function(D, post, main = NULL) {
  cub_par()
  u <- pmax(post$ups, 1e-300); d <- pmax(post$del, 1e-300)
  o <- order(D$y, u)
  if (is.null(main))
    main <- "Fitted mixing weights, CUBpanel(0), every estimation row"
  ylim <- c(min(log10(c(u, d))) - 0.6, 0)
  graphics::plot(seq_along(u), log10(u[o]), pch = 16, cex = 0.8,
                 col = grey(0.12), ylim = ylim,
                 xlab = "estimation row (ordered by response, then by weight)",
                 ylab = expression(log[10] * " fitted mixing weight"),
                 main = main)
  graphics::points(seq_along(d), log10(d[o]), pch = 1, cex = 0.8,
                   col = grey(0.45))
  graphics::abline(h = -2, lty = 2, lwd = 1.5)
  graphics::text(length(u) * 0.015, -2, "one per cent", pos = 3, cex = 0.78)
  graphics::abline(h = log10(1 / D$m), lty = 3, lwd = 1.3, col = grey(0.4))
  graphics::text(length(u) * 0.015, log10(1 / D$m), expression(1/m),
                 pos = 3, cex = 0.78, col = grey(0.3))
  graphics::legend("bottomright", bty = "n", cex = 0.82,
                   legend = c(sprintf("uncertainty   max = %.1e", max(post$ups)),
                              sprintf("shelter       max = %.1e", max(post$del))),
                   pch = c(16, 1), col = c(grey(0.12), grey(0.45)))
  invisible(NULL)
}

## ---------------------------------------------------------------------
## F5  profile deviances: rho_12 (interior) and sigma^2_1 (boundary)
## ---------------------------------------------------------------------
cub_plot_profiles <- function(rho = NULL, s1_grid = NULL, s1 = NULL) {
  op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
  cut1 <- stats::qchisq(0.95, 1)

  cub_par()
  if (!is.null(rho) && !is.null(rho$prof) && nrow(rho$prof) >= 2L) {
    pr <- rho$prof
    graphics::plot(pr$rho, pr$deviance, type = "o", pch = 16, lwd = 1.8,
                   cex = 0.95, col = grey(0.12),
                   ylim = c(0, max(c(pr$deviance, cut1)) * 1.14),
                   xlab = expression(rho[12]), ylab = "profile deviance",
                   main = expression("(a)  " * rho[12] *
                                       ":  interior restriction, " * chi[1]^2))
    graphics::abline(h = cut1, lty = 2, lwd = 1.5)
    graphics::abline(v = 0, lty = 3, col = grey(0.5))
    if (is.finite(rho$rho))
      graphics::abline(v = rho$rho, lty = 1, lwd = 1, col = grey(0.55))
    ok <- pr$deviance <= cut1
    if (any(ok)) graphics::rug(range(pr$rho[ok]), side = 1, lwd = 3,
                               ticksize = 0.045)
    graphics::legend("topright", bty = "n", cex = 0.78, legend = c(
      sprintf("LR = %.3f,  p = %.3f", rho$lr, rho$p),
      sprintf("95%% set on grid: [%.2f, %.2f]", rho$ci[1], rho$ci[2]),
      if (isTRUE(rho$truncated)) "grid truncated by time budget" else ""))
  } else {
    graphics::plot.new()
    graphics::title(main = "(a)  correlation profile not computed")
  }

  cub_par()
  if (!is.null(s1_grid) && nrow(s1_grid) >= 2L) {
    graphics::plot(s1_grid$sigma2_1, s1_grid$deviance, type = "o", pch = 16,
                   lwd = 1.8, cex = 0.95, col = grey(0.12),
                   ylim = c(0, max(c(s1_grid$deviance, cut1)) * 1.14),
                   xlab = expression(sigma[1]^2), ylab = "profile deviance",
                   main = expression("(b)  " * sigma[1]^2 *
                                       ":  boundary at zero"))
    graphics::abline(h = cut1, lty = 2, lwd = 1.5)
    if (!is.null(s1)) {
      graphics::abline(v = s1$est, lty = 1, lwd = 1, col = grey(0.55))
      graphics::legend("topright", bty = "n", cex = 0.78, legend = c(
        sprintf("point estimate %.3f", s1$est),
        sprintf("95%%: [%.3f, %s]", s1$lo,
                if (s1$hi_ok) sprintf("%.3f", s1$hi)
                else sprintf("> %.2f", s1$hi)),
        if (!s1$hi_ok) "upper limit truncated" else ""))
    }
  } else {
    graphics::plot.new()
    graphics::title(main = expression("(b)  " * sigma[1]^2 *
                                        " profile grid not computed"))
  }
  invisible(NULL)
}

## ---------------------------------------------------------------------
## F6  four-panel fit diagnostics
## ---------------------------------------------------------------------
cub_plot_diagnostics <- function(D, g, label = "") {
  op <- graphics::par(mfrow = c(2, 2)); on.exit(graphics::par(op), add = TRUE)

  ## (a) dissimilarity by wave
  cub_par()
  dw <- as.numeric(g$diss_wave); wv <- seq_along(dw)
  graphics::plot(wv, dw, type = "o", pch = 16, lwd = 1.8, col = grey(0.12),
                 xaxt = "n", ylim = c(0, max(c(dw, 0.05)) * 1.3),
                 xlab = "wave", ylab = "normalised dissimilarity",
                 main = "(a)  Dissimilarity index by wave")
  graphics::axis(1, at = wv, labels = names(g$diss_wave))
  graphics::abline(h = 0.05, lty = 2)
  graphics::abline(h = g$diss, lty = 3, lwd = 1.6)
  graphics::text(max(wv), 0.05, "0.05 reference", pos = 3, cex = 0.75, adj = 1)
  graphics::text(min(wv), g$diss, sprintf("pooled %.4f", g$diss), pos = 3,
                 cex = 0.75)

  ## (b) observed against fitted
  cub_par()
  mx <- max(c(g$obs, g$fit))
  bp <- graphics::barplot(rbind(g$obs, g$fit), beside = TRUE,
                          names.arg = seq_len(D$m),
                          col = c(grey(0.25), grey(0.72)), border = NA,
                          ylim = c(0, mx * 1.25),
                          xlab = "response category",
                          ylab = "relative frequency",
                          main = "(b)  Observed against fitted distribution")
  graphics::legend("topright", c("observed", "fitted"),
                   fill = c(grey(0.25), grey(0.72)), border = NA,
                   bty = "n", cex = 0.82)
  graphics::text(mean(bp[, D$cs]), mx * 1.16, sprintf("c = %d", D$cs),
                 cex = 0.78)

  ## (c) binned residuals
  cub_par()
  br <- g$binres
  graphics::plot(br$mu, br$resid, pch = 16, cex = 1.1, col = grey(0.12),
                 ylim = range(c(br$resid, br$band, -br$band)) * 1.25,
                 xlab = expression("bin mean of " * hat(mu)[it]),
                 ylab = "mean residual",
                 main = sprintf("(c)  Binned residuals (%d of %d outside)",
                                g$out_band, g$nbins))
  graphics::abline(h = 0, lty = 1, col = grey(0.6))
  graphics::lines(br$mu, br$band, lty = 2)
  graphics::lines(br$mu, -br$band, lty = 2)

  ## (d) empirical Bayes subject effects
  cub_par()
  Up <- g$post$Up; Qu <- g$pred$Uq
  if (is.null(Qu) || all(Qu == 0)) {
    graphics::plot.new()
    graphics::title(main = "(d)  no subject effects in this rung")
  } else {
    o <- order(Up)
    graphics::plot(seq_along(Up), Up[o], pch = 16, cex = 0.55,
                   col = grey(0.12), ylim = range(Qu),
                   xlab = "respondent, ordered by posterior mean",
                   ylab = expression("posterior mean of " * u[i]),
                   main = "(d)  Feeling subject effects, 95% intervals")
    graphics::segments(seq_along(Up), Qu[o, 1], seq_along(Up), Qu[o, 2],
                       col = grey(0.65), lwd = 0.8)
    graphics::points(seq_along(Up), Up[o], pch = 16, cex = 0.55,
                     col = grey(0.12))
    graphics::abline(h = 0, lty = 2)
    graphics::legend("topleft", bty = "n", cex = 0.76, legend = c(
      sprintf("sd(posterior means) = %.3f", stats::sd(Up)),
      "shrinkage, not error"))
  }
  invisible(NULL)
}

## ---------------------------------------------------------------------
## F7  memory profile: eccentricity is V-shaped, so rho_u needs a picture
## ---------------------------------------------------------------------
cub_plot_memory <- function(MP, m, rho_f = NA_real_, rho_u = NA_real_) {
  if (is.null(MP)) {
    graphics::plot.new()
    graphics::title(main = "memory profile unavailable: no lagged features")
    return(invisible(NULL))
  }
  op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)

  cub_par()
  graphics::plot(MP$previous, MP$phi, type = "o", pch = 16, lwd = 2,
                 cex = 1.05, col = grey(0.12), ylim = c(0, 1), xaxt = "n",
                 xlab = "previous observed rating", ylab = "fitted weight",
                 main = expression("(a)  " * varphi *
                                     " is monotone, " * upsilon * " is not"))
  graphics::axis(1, at = seq_len(m))
  graphics::lines(MP$previous, MP$upsilon, type = "o", pch = 1, lwd = 2,
                  cex = 1.05, col = grey(0.42))
  graphics::lines(MP$previous, MP$delta, lty = 3, lwd = 1.6, col = grey(0.6))
  graphics::abline(v = (m + 1) / 2, lty = 3, col = grey(0.6))
  graphics::legend("topleft", bty = "n", cex = 0.80,
                   legend = c(expression(varphi * "  (feeling)"),
                              expression(upsilon * "  (uncertainty)"),
                              expression(delta * "  (shelter, no lag)")),
                   pch = c(16, 1, NA), lty = c(1, 1, 3),
                   col = c(grey(0.12), grey(0.42), grey(0.6)))

  cub_par()
  graphics::plot(MP$previous, MP$ER, type = "o", pch = 16, lwd = 2, cex = 1.05,
                 col = grey(0.12), xaxt = "n",
                 ylim = range(MP$ER) + c(-0.15, 0.15),
                 xlab = "previous observed rating",
                 ylab = expression(E(R[it])),
                 main = "(b)  implied mean of the next rating")
  graphics::axis(1, at = seq_len(m))
  graphics::abline(h = (m + 1) / 2, lty = 2, col = grey(0.5))
  graphics::text(1, (m + 1) / 2, "scale midpoint", pos = 3, cex = 0.75,
                 col = grey(0.35))
  if (is.finite(rho_f) && is.finite(rho_u))
    graphics::legend("bottomright", bty = "n", cex = 0.76, legend = c(
      sprintf("feeling shift, lowest to highest: %+.3f", 2 * rho_f),
      sprintf("uncertainty shift, mid to end:   %+.3f", rho_u)))
  invisible(NULL)
}

## ---------------------------------------------------------------------
## F8  the ladder, as a picture
## ---------------------------------------------------------------------
cub_plot_ladder <- function(ladder) {
  cub_par(mar = c(8.5, 4.3, 2.4, 1.1))
  tb <- ladder$table
  k <- seq_len(nrow(tb))
  rng <- range(c(tb$AIC, tb$BIC_N))
  graphics::plot(k, tb$AIC, type = "o", pch = 16, lwd = 1.8, col = grey(0.12),
                 xaxt = "n", ylim = rng + c(-1, 1) * diff(rng) * 0.08,
                 xlab = "", ylab = "information criterion",
                 main = "Specification ladder")
  graphics::lines(k, tb$BIC_N, type = "o", pch = 1, lwd = 1.8,
                  col = grey(0.45))
  graphics::axis(1, at = k, labels = tb$spec, las = 2, cex.axis = 0.78)
  graphics::points(which.min(tb$AIC), min(tb$AIC), pch = 16, cex = 2.1,
                   col = grey(0.12))
  graphics::points(which.min(tb$BIC_N), min(tb$BIC_N), pch = 1, cex = 2.1,
                   lwd = 2, col = grey(0.45))
  graphics::legend("topright", bty = "n", cex = 0.82,
                   legend = c("AIC", expression(BIC[N])),
                   pch = c(16, 1), lty = 1,
                   col = c(grey(0.12), grey(0.45)))
  invisible(NULL)
}

## ---------------------------------------------------------------------
## convenience: every figure an analysis object can produce
## ---------------------------------------------------------------------
cub_figure_list <- function() c(
  waves       = "Wave profiles (observed)",
  within      = "Within-respondent share by covariate",
  dispersion  = "Dispersion floor of the family",
  weights     = "Fitted mixing weights, log10 scale",
  ladder      = "Specification ladder",
  profiles    = "Profile deviances",
  diagnostics = "Fit diagnostics (four panels)",
  memory      = "Memory profile")

cub_draw <- function(an, which = "waves") {
  D <- an$D
  switch(which,
         waves = cub_plot_waves(D, an$floor[["R"]]),
         within = cub_plot_within(D),
         dispersion = cub_plot_dispersion(list(list(
           label = "this panel", D = D,
           mu = an$stats_selected$mu, R = an$floor[["R"]]))),
         weights = cub_plot_weights(D, an$verif$post0),
         ladder = cub_plot_ladder(an$ladder),
         profiles = cub_plot_profiles(an$rho, an$profile_s1_grid,
                                      an$profile_s1),
         diagnostics = cub_plot_diagnostics(D, an$stats_selected,
                                            an$fit_selected$label),
         memory = {
           fd <- an$fits[[an$spec_dyn]]
           cub_plot_memory(an$memory, D$m,
                           if ("xi:lag_f" %in% names(fd$th))
                             unname(fd$th["xi:lag_f"]) else NA_real_,
                           if ("ups:lag_u" %in% names(fd$th))
                             unname(fd$th["ups:lag_u"]) else NA_real_)
         },
         stop("unknown figure: ", which))
}