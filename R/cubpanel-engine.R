## =====================================================================
##  CUBpanel -- estimation engine for panel CUB models with shelter
##  ---------------------------------------------------------------------
##  Dependencies: base, stats, parallel.  No I/O, no external packages.
##
##  Sigma is parameterised through its Cholesky factor:
##      u_i = L11 z1 ,  v_i = L21 z1 + L22 z2
##      L11 = exp(a) , L22 = exp(b) , L21 free
##      sigma1^2 = L11^2 , sigma2^2 = L21^2 + L22^2 , sigma12 = L11 L21
##      rho12    = L21 / sqrt(L21^2 + L22^2)
##  Consequently rho12 = 0 <=> L21 = 0, an INTERIOR point, so the
##  correlation is testable by an ordinary chi^2_1 likelihood ratio.
##
##  Specification rungs:
##      "SB"       degenerate shifted binomial (lower boundary)
##      "0"        pooled three-component mixture, Sigma = 0
##      "LAG"      pooled, observation-driven dynamics
##      "RE"       correlated subject effects, no dynamics
##      "RE-LAG"   subject effects and dynamics
##      "CRE-LAG"  adds the correlated-random-effects block
## =====================================================================

## ---------------------------------------------------------------------
## 0.  numerical utilities
## ---------------------------------------------------------------------

.el <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

.pl <- function(x) {
  d <- dim(x); r <- stats::plogis(x); if (!is.null(d)) dim(r) <- d; r
}
.pll <- function(x) {
  d <- dim(x); r <- stats::plogis(x, log.p = TRUE); if (!is.null(d)) dim(r) <- d; r
}
.lse  <- function(v) { mx <- max(v); mx + log(sum(exp(v - mx))) }
.rsum <- function(x, g) rowsum(as.matrix(x), g, reorder = TRUE)

#' Format a p-value, flagging simulation-limited tails
cub_pfmt <- function(p, B = NA_integer_) {
  if (length(p) != 1L || !is.finite(p)) return("NA")
  if (is.finite(B) && B > 0 && p <= 1 / (B + 1) + 1e-12)
    return(sprintf("< %.1e", 1 / (B + 1)))
  if (p < 1e-4) sprintf("%.3e", p) else sprintf("%.4g", p)
}
cub_ptex <- function(p) {
  if (!is.finite(p)) "---"
  else if (p < 1e-4) "$<10^{-4}$" else sprintf("$%.3f$", p)
}
cub_num <- function(x, d = 3) {
  if (!is.finite(x)) "---" else formatC(x, format = "f", digits = d)
}

## Golub-Welsch construction of the Gauss-Hermite rule from the
## symmetric Jacobi matrix of the Hermite recursion.  No tables.
cub_gh1 <- function(Q) {
  Q <- as.integer(Q)
  k <- seq_len(Q - 1L)
  J <- matrix(0, Q, Q)
  J[cbind(k, k + 1L)] <- sqrt(k / 2)
  J[cbind(k + 1L, k)] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(z = sqrt(2) * e$values[o],
       w = (e$vectors[1L, o]^2) / sum(e$vectors[1L, ]^2))
}

.GH_CACHE <- new.env(parent = emptyenv())

cub_gh2 <- function(Q) {
  key <- paste0("Q", Q)
  if (!is.null(.GH_CACHE[[key]])) return(.GH_CACHE[[key]])
  g <- cub_gh1(Q)
  w <- g$w / sum(g$w)
  o <- list(z1 = rep(g$z, each = Q), z2 = rep(g$z, times = Q),
            w  = rep(w, each = Q) * rep(w, times = Q), Q = Q)
  o$lw <- log(o$w)
  .GH_CACHE[[key]] <- o
  o
}

## ---------------------------------------------------------------------
## 1.  additive variance decomposition (between + within = total)
## ---------------------------------------------------------------------

cub_decomp <- function(y, gi) {
  gi <- as.integer(factor(gi)); n <- length(y); N <- max(gi)
  mb <- stats::ave(y, gi); gm <- mean(y)
  SSb <- sum((mb - gm)^2); SSw <- sum((y - mb)^2); SSt <- SSb + SSw
  MSb <- SSb / max(N - 1L, 1L); MSw <- SSw / max(n - N, 1L)
  s2u <- max((MSb - MSw) / (n / N), 0)
  c(total = SSt / (n - 1), between = SSb / (n - 1), within = SSw / (n - 1),
    MS_between = MSb, MS_within = MSw, sigma2_u = s2u, sigma2_e = MSw,
    icc = if (s2u + MSw > 0) s2u / (s2u + MSw) else 0,
    between_share = SSb / SSt, within_share = SSw / SSt)
}

## Fraction of a covariate's total sum of squares that survives
## subtraction of its respondent mean.  Exactly zero for time-invariant
## covariates, which is the quantity the Mundlak block is selected on.
cub_wshare <- function(x, gi) {
  if (is.factor(x) || is.character(x)) x <- as.numeric(as.integer(as.factor(x)))
  SSt <- sum((x - mean(x))^2)
  if (SSt <= 0) 0 else sum((x - stats::ave(x, gi))^2) / SSt
}

## ---------------------------------------------------------------------
## 2.  data preparation
## ---------------------------------------------------------------------

#' Prepare a panel of ordinal ratings for CUBpanel estimation
#'
#' Sorts by respondent and wave, applies listwise deletion, reserves each
#' respondent's first observation as the initial condition, builds the two
#' lagged features from the previously OBSERVED rating (not from calendar
#' wave t-1), centres them and stores the centring constants, assembles the
#' three design matrices and the correlated-random-effects block.
#'
#' @param m number of scale categories.  Resolved from max(response) when
#'   NULL, but the observed support is reported beside the resolved value
#'   because the two can differ and the difference changes the model.
#' @param shelter_cat the refuge category, declared by the analyst.  There
#'   is no default: selecting it by maximising the likelihood over the m
#'   candidates leaves no admissible reference distribution for the tests
#'   that follow.
cub_prepare <- function(data, id = "Respondent", time = "Period",
                        response = "Response", m = NULL,
                        feeling     = ~ Gender + Experience + Inflation,
                        uncertainty = ~ Education + Delivery_Speed + Customer_Service,
                        shelter     = ~ Inflation,
                        shelter_cat = NULL,
                        cre = c("wooldridge", "mundlak"),
                        drop_gaps = TRUE) {

  cre <- match.arg(cre)
  cov_vars <- unique(c(all.vars(feeling), all.vars(uncertainty), all.vars(shelter)))
  vars <- unique(c(id, time, response, cov_vars))
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop("columns absent from the data: ", paste(miss, collapse = ", "))

  d  <- data[, vars, drop = FALSE]
  n0 <- nrow(d)
  d  <- d[stats::complete.cases(d), , drop = FALSE]
  n_del <- n0 - nrow(d)
  if (!nrow(d)) stop("no complete cases remain after listwise deletion")

  d <- d[order(d[[id]], d[[time]]), , drop = FALSE]
  y <- as.integer(d[[response]])
  if (any(is.na(y))) stop("the response must be coercible to integer")
  m_res <- max(y)
  if (is.null(m)) m <- m_res else m <- as.integer(m)
  if (min(y) < 1L || max(y) > m)
    stop(sprintf("the response must lie in 1..m; observed range is %d..%d with m = %d",
                 min(y), max(y), m))
  if (m < 4L)
    stop("identification: m >= 4 is required for the two-component mixture")
  if (is.null(shelter_cat))
    stop("a shelter component requires shelter_cat to be declared explicitly")
  shelter_cat <- as.integer(shelter_cat)
  if (m < 5L)
    stop("identification: m >= 5 is required when a shelter component is requested")
  if (!(shelter_cat %in% seq_len(m)))
    stop("shelter_cat lies outside the scale")

  gi0 <- as.integer(factor(d[[id]]))
  wv  <- as.integer(factor(d[[time]], levels = sort(unique(d[[time]]))))
  wv_lab <- as.character(d[[time]])
  Ti0 <- as.integer(table(gi0)); maxT <- max(Ti0)
  first <- !duplicated(gi0)

  prev <- c(NA_integer_, y[-length(y)]);  prev[first] <- NA_integer_
  pw   <- c(NA_integer_, wv[-length(wv)]); pw[first]  <- NA_integer_
  gap  <- !first & (wv - pw) != 1L
  n_gap <- sum(gap)
  keep <- !first & !(drop_gaps & gap)
  if (sum(keep) < 2L)
    stop("too few estimation rows remain after reserving the initial condition")

  y_first <- y[first]
  subj <- sort(unique(gi0[keep]))
  gi   <- match(gi0[keep], subj)
  y1s  <- y_first[subj]
  N    <- length(subj)
  if (N < 2L) stop("fewer than two respondents have T_i >= 2")

  ctr <- (m + 1) / 2; hw <- (m - 1) / 2
  lf_raw <- (prev[keep] - ctr) / hw                 # signed deviation
  lu_raw <- abs(prev[keep] - ctr) / hw              # eccentricity
  c_lf <- mean(lf_raw); c_lu <- mean(lu_raw)
  lf <- lf_raw - c_lf; lu <- lu_raw - c_lu
  if1 <- ((y1s - ctr) / hw)[gi];    if1 <- if1 - mean(if1)
  iu1 <- (abs(y1s - ctr) / hw)[gi]; iu1 <- iu1 - mean(iu1)

  dk <- d[keep, , drop = FALSE]
  for (v in cov_vars) if (is.character(dk[[v]])) dk[[v]] <- factor(dk[[v]])

  mm <- function(f) {
    M <- stats::model.matrix(f, data = dk)
    attr(M, "assign") <- NULL; attr(M, "contrasts") <- NULL
    M
  }
  Xb <- mm(feeling); Wb <- mm(uncertainty); Z <- mm(shelter)
  Xl <- cbind(lag_f = lf); Wl <- cbind(lag_u = lu)

  ## Mundlak means only for covariates with within-respondent variation:
  ## a time-invariant covariate's mean is collinear with the covariate.
  mk_cre <- function(M, initcol, initname) {
    cand <- setdiff(colnames(M), "(Intercept)")
    sel  <- cand[vapply(cand, function(cc) cub_wshare(M[, cc], gi) > 0.01, TRUE)]
    O <- NULL
    if (length(sel)) {
      O <- vapply(sel, function(cc) stats::ave(M[, cc], gi), numeric(nrow(M)))
      O <- as.matrix(O); colnames(O) <- paste0(sel, "_bar")
      O <- sweep(O, 2, colMeans(O))
    }
    E <- matrix(initcol - mean(initcol), ncol = 1,
                dimnames = list(NULL, initname))
    cbind(O, E)
  }
  if (identical(cre, "wooldridge")) {
    Xm <- mk_cre(Xb, if1, "y1_f"); Wm <- mk_cre(Wb, iu1, "y1_u")
  } else {
    Xm <- mk_cre(Xb, stats::ave(lf, gi), "lag_f_bar")
    Wm <- mk_cre(Wb, stats::ave(lu, gi), "lag_u_bar")
  }
  bar <- function(M) grep("_bar$", colnames(M), value = TRUE)

  wave_est <- wv_lab[keep]
  tab_wave <- prop.table(table(Wave = wave_est,
                               Category = factor(y[keep], levels = seq_len(m))), 1)
  prev_est <- prev[keep]
  tab_trans <- table(From = factor(prev_est, levels = seq_len(m)),
                     To   = factor(y[keep],  levels = seq_len(m)))

  out <- list(
    y = as.integer(y[keep]), gi = gi, wave = wave_est,
    wave_idx = wv[keep], n = sum(keep), N = N, m = m, m_resolved = m_res,
    cs = shelter_cat, subj = subj,
    Xb = Xb, Xl = Xl, Wb = Wb, Wl = Wl, Z = Z, Xm = Xm, Wm = Wm,
    ctr = c(lag_f = c_lf, lag_u = c_lu),
    cre = list(mode = cre,
               xbar = bar(Xm), xinit = setdiff(colnames(Xm), bar(Xm)),
               wbar = bar(Wm), winit = setdiff(colnames(Wm), bar(Wm))),
    formulas = list(feeling = feeling, uncertainty = uncertainty,
                    shelter = shelter),
    roles = c(id = id, time = time, response = response),
    info = list(n_raw = n0, n_del = n_del, Ti = Ti0,
                balanced = length(unique(Ti0)) == 1L,
                minT = min(Ti0), maxT = maxT,
                lag_ok = maxT >= 3L, oos_ok = min(Ti0) >= 3L,
                n_gap = n_gap, rows_est = sum(keep),
                support = sort(unique(y[keep])),
                unused = setdiff(seq_len(m), sort(unique(y[keep])))),
    wsh = vapply(cov_vars, function(v) cub_wshare(dk[[v]], gi), 0),
    tab_wave = tab_wave, tab_trans = tab_trans,
    hash = paste(sum(y[keep]), sum(keep), N, sep = "|"))
  out
}

cub_subset <- function(D, keep) {
  for (nm in c("Xb", "Xl", "Wb", "Wl", "Z", "Xm", "Wm"))
    if (!is.null(D[[nm]])) D[[nm]] <- D[[nm]][keep, , drop = FALSE]
  g <- D$gi[keep]; subj <- sort(unique(g))
  D$gi <- match(g, subj); D$subj <- D$subj[subj]
  D$y <- D$y[keep]; D$wave <- D$wave[keep]; D$wave_idx <- D$wave_idx[keep]
  D$n <- length(D$y); D$N <- length(subj)
  D$hash <- paste(sum(D$y), D$n, D$N, sep = "|")
  D
}

#' Pre-fit dispersion diagnostic
#'
#' Both components that distinguish the mixture from its shifted-binomial
#' kernel can only ADD dispersion at a fixed conditional mean, so the kernel
#' variance (mu-1)(m-mu)/(m-1) is a floor on what the family can represent.
#' This function compares the conditional residual variance with that floor
#' at the cost of one weighted least-squares solve.  R < 1 predicts that no
#' member of the family with positive mixing weights can be optimal.
cub_floor_diag <- function(D) {
  M <- D$Xb
  fit <- suppressWarnings(stats::glm.fit(
    M, (D$y - 1) / (D$m - 1), weights = rep(D$m - 1, D$n),
    family = stats::binomial()))
  pr <- stats::plogis(as.vector(M %*% fit$coefficients))
  mu <- 1 + (D$m - 1) * pr
  fl <- mean((mu - 1) * (D$m - mu) / (D$m - 1))
  rv <- sum((D$y - mu)^2) / max(D$n - ncol(M), 1L)
  dd <- cub_decomp(D$y, D$gi)
  c(floor = fl, resid_var = rv, R = rv / fl,
    MS_within = dd[["MS_within"]], R_within = dd[["MS_within"]] / fl,
    dist_cs = mean(abs(D$cs - mu)),
    ll_sb_glm = sum(stats::dbinom(D$y - 1L, D$m - 1L, pr, log = TRUE)))
}

## ---------------------------------------------------------------------
## 3.  core: marginal likelihood, analytic gradient, respondent scores
## ---------------------------------------------------------------------

cub_core <- function(D, spec, rho_free = TRUE) {
  stopifnot(spec %in% c("SB", "0", "LAG", "RE", "RE-LAG", "CRE-LAG"))
  lag <- spec %in% c("LAG", "RE-LAG", "CRE-LAG")
  mk  <- spec == "CRE-LAG"
  re  <- spec %in% c("RE", "RE-LAG", "CRE-LAG")
  use_u <- use_s <- spec != "SB"
  if (lag && !D$info$lag_ok)
    stop("lagged specifications are unavailable: max T_i < 3")
  if (!re) rho_free <- FALSE

  X <- D$Xb; if (lag) X <- cbind(X, D$Xl); if (mk) X <- cbind(X, D$Xm)
  W <- NULL; Z <- NULL
  if (use_u) { W <- D$Wb; if (lag) W <- cbind(W, D$Wl); if (mk) W <- cbind(W, D$Wm) }
  if (use_s) Z <- D$Z

  n <- D$n; N <- D$N; m <- D$m; y <- D$y; gi <- D$gi; cs <- D$cs
  Ind <- as.numeric(y == cs)
  lck <- lchoose(m - 1L, y - 1L); y1 <- y - 1L; my <- m - y

  kx <- ncol(X)
  kw <- if (use_u) ncol(W) else 0L
  kz <- if (use_s) ncol(Z) else 0L
  nv <- if (re) (if (rho_free) 3L else 2L) else 0L
  ib <- seq_len(kx)
  ig <- if (kw) kx + seq_len(kw) else integer(0)
  ia <- if (kz) kx + kw + seq_len(kz) else integer(0)
  iv <- if (re) kx + kw + kz + seq_len(nv) else integer(0)
  p  <- kx + kw + kz + nv
  nms <- c(paste0("xi:", colnames(X)),
           if (kw) paste0("ups:", colnames(W)),
           if (kz) paste0("del:", colnames(Z)),
           if (re) c("log.L11", "log.L22", if (rho_free) "L21"))

  chol_of <- function(th) {
    list(L11 = exp(th[iv[1]]), L22 = exp(th[iv[2]]),
         L21 = if (rho_free) th[iv[3]] else 0)
  }

  ev <- function(th, Q = 15L, grad = TRUE, scores = FALSE, post = FALSE,
                 pred = FALSE, quant = FALSE) {
    if (quant) pred <- TRUE
    if (pred)  post <- TRUE
    th <- as.numeric(th)
    stopifnot(length(th) == p, all(is.finite(th)))

    lam0 <- as.vector(X %*% th[ib])
    eta0 <- if (use_u) as.vector(W %*% th[ig]) else NULL
    del  <- if (use_s) .pl(as.vector(Z %*% th[ia])) else NULL

    if (re) {
      G <- cub_gh2(Q); L <- chol_of(th)
      uq <- L$L11 * G$z1
      vq <- L$L21 * G$z1 + L$L22 * G$z2
      lw <- G$lw
    } else { G <- NULL; L <- NULL; uq <- 0; vq <- 0; lw <- 0 }
    nQ <- length(uq)

    LAM <- outer(lam0, uq, "+")
    lp <- .pll(LAM); lq <- .pll(-LAM); P <- exp(lp)
    B  <- exp(lck + y1 * lp + my * lq)
    UP <- if (use_u) .pl(outer(eta0, vq, "+")) else NULL
    Gm <- if (use_u) (1 - UP) * B + UP / m else B
    F_ <- if (use_s) del * Ind + (1 - del) * Gm else Gm
    F_ <- pmax(F_, 1e-300); dim(F_) <- c(n, nQ)
    LF <- log(F_); dim(LF) <- c(n, nQ)

    A  <- .rsum(LF, gi) + rep(lw, each = N)
    mx <- A[cbind(seq_len(N), max.col(A, ties.method = "first"))]
    Lc <- mx + log(rowSums(exp(A - mx)))
    ll <- sum(Lc)
    o <- list(ll = ll, nll = if (is.finite(ll)) -ll else 1e10)
    if (!(grad || scores || post)) return(o)
    Wq <- exp(A - Lc)

    if (grad || scores) {
      Wr <- Wq[gi, , drop = FALSE]; iF <- 1 / F_
      od <- if (use_s) (1 - del) else 1
      ou <- if (use_u) (1 - UP) else 1
      dlam <- iF * (od * ou * B * (y1 - (m - 1) * P))
      glam <- rowSums(Wr * dlam)
      g <- numeric(p)
      g[ib] <- as.vector(crossprod(X, glam))
      deta <- NULL; geta <- NULL; gzet <- NULL
      if (use_u) {
        deta <- iF * (od * UP * (1 - UP) * (1 / m - B))
        geta <- rowSums(Wr * deta)
        g[ig] <- as.vector(crossprod(W, geta))
      }
      if (use_s) {
        dzet <- iF * (del * (1 - del) * (Ind - Gm))
        gzet <- rowSums(Wr * dzet)
        g[ia] <- as.vector(crossprod(Z, gzet))
      }
      if (re) {
        cl <- colSums(Wr * dlam); ce <- colSums(Wr * deta)
        g[iv[1]] <- sum(cl * uq)                  # d u / d log L11 = u
        g[iv[2]] <- sum(ce * (L$L22 * G$z2))      # d v / d log L22 = L22 z2
        if (rho_free) g[iv[3]] <- sum(ce * G$z1)  # d v / d L21     = z1
      }
      o$grad <- -g
      if (scores) {
        Sc <- matrix(0, N, p)
        Sc[, ib] <- .rsum(X * glam, gi)
        if (use_u) Sc[, ig] <- .rsum(W * geta, gi)
        if (use_s) Sc[, ia] <- .rsum(Z * gzet, gi)
        if (re) {
          Dl <- .rsum(Wr * dlam, gi); De <- .rsum(Wr * deta, gi)
          Sc[, iv[1]] <- as.vector(Dl %*% uq)
          Sc[, iv[2]] <- as.vector(De %*% (L$L22 * G$z2))
          if (rho_free) Sc[, iv[3]] <- as.vector(De %*% G$z1)
        }
        o$scores <- Sc
      }
    }

    if (post) {
      Up <- if (re) as.vector(Wq %*% uq) else rep(0, N)
      Vp <- if (re) as.vector(Wq %*% vq) else rep(0, N)
      o$post <- list(
        Wq = Wq, uq = uq, vq = vq, wpri = if (re) G$w else 1,
        lam0 = lam0, eta0 = if (use_u) eta0 else rep(NA_real_, n),
        del = if (use_s) del else rep(0, n), Up = Up, Vp = Vp,
        xi  = .pl(lam0 + Up[gi]),
        ups = if (use_u) .pl(eta0 + Vp[gi]) else rep(0, n),
        sdU = if (re) sqrt(pmax(as.vector(Wq %*% uq^2) - Up^2, 0)) else rep(0, N),
        sdV = if (re) sqrt(pmax(as.vector(Wq %*% vq^2) - Vp^2, 0)) else rep(0, N))
    }

    if (pred) {
      wpri <- if (re) G$w else 1
      Pm <- matrix(0, n, m)
      rv <- 0:(m - 1L); lchv <- lchoose(m - 1L, rv)
      for (q in seq_len(nQ)) {
        lpr <- .pll(lam0 + uq[q]); lqr <- .pll(-(lam0 + uq[q]))
        BM  <- exp(outer(lpr, rv) + outer(lqr, rev(rv)) + rep(lchv, each = n))
        upq <- if (use_u) .pl(eta0 + vq[q]) else 0
        Pm  <- Pm + wpri[q] * ((1 - upq) * BM + upq / m)
      }
      if (use_s) { Pm <- (1 - del) * Pm; Pm[, cs] <- Pm[, cs] + del }
      o$pred <- list(Pm = Pm, Qpred = Q, Uq = NULL, Vq = NULL)
      if (quant) {
        qt <- function(nd, probs = c(0.025, 0.975)) {
          if (!re) return(cbind(rep(0, N), rep(0, N)))
          o2 <- order(nd); nds <- nd[o2]; Wo <- Wq[, o2, drop = FALSE]
          t(apply(Wo, 1, function(w) {
            sw <- sum(w)
            if (!is.finite(sw) || sw <= 0) return(range(nds))
            cw <- cumsum(w) / sw
            cwm <- cw - w / (2 * sw)                    # Hazen positions
            ok <- c(TRUE, diff(cwm) > 0)
            stats::approx(cwm[ok], nds[ok], xout = probs,
                          rule = 2, ties = "ordered")$y
          }))
        }
        o$pred$Uq <- qt(uq); o$pred$Vq <- qt(vq)
      }
    }
    o
  }

  ## Transparent looped reference, used only by self-check C5.
  ev_ref <- function(th, Q) {
    stopifnot(re, use_u, use_s)
    lam0 <- as.vector(X %*% th[ib]); eta0 <- as.vector(W %*% th[ig])
    del  <- stats::plogis(as.vector(Z %*% th[ia]))
    G <- cub_gh2(Q); L <- chol_of(th)
    uq <- L$L11 * G$z1; vq <- L$L21 * G$z1 + L$L22 * G$z2
    tot <- 0
    for (i in seq_len(N)) {
      k <- which(gi == i); acc <- numeric(length(uq))
      for (q in seq_along(uq)) {
        pr <- stats::plogis(lam0[k] + uq[q])
        up <- stats::plogis(eta0[k] + vq[q])
        bb <- stats::dbinom(y[k] - 1L, m - 1L, pr)
        acc[q] <- sum(log(del[k] * (y[k] == cs) +
                            (1 - del[k]) * ((1 - up) * bb + up / m))) + G$lw[q]
      }
      tot <- tot + .lse(acc)
    }
    tot
  }

  list(ev = ev, ev_ref = ev_ref, chol_of = chol_of, p = p, nms = nms,
       re = re, rho_free = rho_free, spec = spec, use_u = use_u,
       use_s = use_s, lag = lag, mk = mk, kx = kx, kw = kw, kz = kz,
       ib = ib, ig = ig, ia = ia, iv = iv, X = X, W = W, Z = Z, D = D)
}

## ---------------------------------------------------------------------
## 4.  stage 1: EM on the pooled likelihood (three weighted GLMs)
## ---------------------------------------------------------------------

cub_em_pooled <- function(co, iter = 300L, tol = 1e-11) {
  D <- co$D; y <- D$y; m <- D$m; cs <- D$cs; n <- D$n
  X <- co$X; W <- co$W; Z <- co$Z
  Ind <- as.numeric(y == cs); lck <- lchoose(m - 1L, y - 1L)
  sfit <- function(M, r, w) {
    v <- tryCatch(suppressWarnings(stats::glm.fit(
      M, r, weights = w, family = stats::binomial())$coefficients),
      error = function(e) rep(0, ncol(M)))
    v[!is.finite(v)] <- 0; v
  }
  b <- sfit(X, pmin(pmax((y - 1) / (m - 1), .02), .98), rep(m - 1, n))
  g <- if (co$use_u) { v <- numeric(ncol(W)); v[1] <- stats::qlogis(0.25); v } else NULL
  a <- if (co$use_s) {
    v <- numeric(ncol(Z))
    v[1] <- stats::qlogis(min(max(mean(Ind) - 1 / m, 0.02), 0.5)); v
  } else NULL
  llo <- -Inf
  for (it in seq_len(iter)) {
    lp <- .pll(as.vector(X %*% b)); lq <- .pll(-as.vector(X %*% b))
    B  <- exp(lck + (y - 1) * lp + (m - y) * lq)
    up <- if (co$use_u) .pl(as.vector(W %*% g)) else rep(0, n)
    dl <- if (co$use_s) .pl(as.vector(Z %*% a)) else rep(0, n)
    f  <- pmax(dl * Ind + (1 - dl) * ((1 - up) * B + up / m), 1e-300)
    ll <- sum(log(f)); if (abs(ll - llo) < tol) break; llo <- ll
    ts <- dl * Ind / f
    tb <- (1 - dl) * (1 - up) * B / f
    tu <- (1 - dl) * up / (m * f)
    b <- sfit(X, (y - 1) / (m - 1), pmax(tb, 1e-8) * (m - 1))
    if (co$use_u) { wt <- pmax(tb + tu, 1e-8); g <- sfit(W, tu / wt, wt) }
    if (co$use_s) a <- sfit(Z, pmin(pmax(ts, 1e-8), 1 - 1e-8), rep(1, n))
  }
  c(b, g, a)
}

## ---------------------------------------------------------------------
## 5.  stage 2: marginal likelihood maximisation; restriction maps
## ---------------------------------------------------------------------

cub_proj_gmax <- function(g, x, lo, up, eps = 1e-6) {
  gg <- g; al <- x <= lo + eps; au <- x >= up - eps
  if (any(al)) gg[al] <- pmin(gg[al], 0)
  if (any(au)) gg[au] <- pmax(gg[au], 0)
  max(abs(gg))
}

cub_bounds <- function(co) {
  k <- co$p - length(co$iv)
  list(lo = c(rep(-20, k), if (co$re) c(-6, -6, if (co$rho_free) -6)),
       up = c(rep( 20, k), if (co$re) c(1.5, 1.5, if (co$rho_free) 6)))
}

## Frontier flags are read off the NATURAL parameters, not off the
## Cholesky coordinates in which the optimiser works.  When a flag is
## raised, the plug-in chi-bar-squared reference and the delta-method
## standard errors for sigma_12 and rho_12 are withheld.
cub_nat_bnd <- function(th, co, lo) {
  if (!co$re) return(character(0))
  i <- co$iv; L <- co$chol_of(th)
  s2 <- sqrt(L$L21^2 + L$L22^2)
  rho <- if (s2 > 0) L$L21 / s2 else 0
  out <- character(0)
  if (th[i[1]] <= lo[i[1]] + 1e-3) out <- c(out, "sigma_1 -> 0")
  if (s2 <= 2e-3) out <- c(out, "sigma_2 -> 0")
  else if (co$rho_free && abs(rho) >= 0.995)
    out <- c(out, sprintf("|rho_12| -> 1 (rho = %+.4f)", rho))
  out
}

cub_fit <- function(D, spec, rho_free = TRUE, Q = 15L, nstart = 1L,
                    start = NULL, tmap = NULL, label = NULL, seed = 11L,
                    gtol = 1e-3, maxit = 600L) {
  co <- cub_core(D, spec, rho_free); p <- co$p
  if (is.null(label)) label <- spec
  bb <- cub_bounds(co); lo <- bb$lo; up <- bb$up

  th0 <- if (!is.null(start)) as.numeric(start) else
    c(cub_em_pooled(co),
      if (co$re) c(log(0.45), log(0.45), if (co$rho_free) 0))
  stopifnot(length(th0) == p)
  th0 <- pmin(pmax(th0, lo + 1e-6), up - 1e-6)

  if (is.null(tmap)) {
    fo <- function(x) x; jacf <- NULL; nd <- p
    x0 <- th0; flo <- lo; fup <- up
  } else {
    fo <- tmap$fo; jacf <- tmap$jac; nd <- tmap$nd
    x0 <- tmap$x0; flo <- tmap$lo; fup <- tmap$up
    x0 <- pmin(pmax(x0, flo + 1e-6), fup - 1e-6)
  }

  fn <- function(x) {
    v <- co$ev(fo(x), Q, grad = FALSE)$nll
    if (is.finite(v)) v else 1e10
  }
  gr <- function(x) {
    g <- co$ev(fo(x), Q)$grad
    g <- if (is.null(jacf)) g else
      as.vector(crossprod(if (is.function(jacf)) jacf(x) else jacf, g))
    g[!is.finite(g)] <- 0; g
  }
  run <- function(x, factr)
    stats::optim(x, fn, gr, method = "L-BFGS-B", lower = flo, upper = fup,
                 control = list(maxit = maxit, factr = factr,
                                pgtol = gtol / 20))

  set.seed(seed); best <- NULL
  for (s in seq_len(max(1L, nstart))) {
    x <- if (s == 1L) x0 else
      pmin(pmax(x0 + stats::rnorm(nd, 0, 0.30), flo + 1e-6), fup - 1e-6)
    r <- try(run(x, 1e5), silent = TRUE)
    if (inherits(r, "try-error")) next
    if (is.null(best) || r$value < best$value) best <- r
  }
  if (is.null(best))
    stop(sprintf("cub_fit(%s): every start failed", label))

  ## tightening restarts until the projected gradient certifies KKT
  for (k in 1:3) {
    if (cub_proj_gmax(gr(best$par), best$par, flo, fup) < gtol) break
    r2 <- try(run(best$par, 1e2), silent = TRUE)
    if (inherits(r2, "try-error") || r2$value >= best$value - 1e-12) break
    best <- r2
  }

  th <- fo(best$par); names(th) <- co$nms
  gfull <- co$ev(th, Q, grad = TRUE)$grad
  gx <- if (is.null(jacf)) gfull else
    as.vector(crossprod(if (is.function(jacf)) jacf(best$par) else jacf, gfull))
  gpr <- cub_proj_gmax(gx, best$par, flo, fup)
  bnd <- if (is.null(tmap))
    c(paste0(co$nms[th <= lo + 1e-3], "[lo]"),
      paste0(co$nms[th >= up - 1e-3], "[up]")) else character(0)
  bnat <- if (is.null(tmap)) cub_nat_bnd(th, co, lo) else character(0)

  out <- list(label = label, spec = spec, rho_free = rho_free, th = th,
              ll = -best$value, p = nd, p_full = p, Q = Q, co = co, D = D,
              n = D$n, N = D$N, tmap = tmap, lo = lo, up = up,
              conv = best$convergence, gmax = gpr, boundary = bnd,
              boundary_nat = bnat,
              kkt_informative = (gpr < gtol) && !length(bnat),
              hash = D$hash, ok = (gpr < gtol),
              nest_ok = NA, nest_ok_lr = NA, nest_gap = NA_real_,
              aic  = 2 * best$value + 2 * nd,
              bic  = 2 * best$value + log(D$n) * nd,
              bicN = 2 * best$value + log(D$N) * nd)
  out
}

cub_fit_line <- function(f) sprintf(
  "  %-30s p=%2d  ll=%13.4f  AIC=%10.2f  BIC_N=%10.2f  |g|=%.1e %-7s%s%s",
  f$label, f$p, f$ll, f$aic, f$bicN, f$gmax,
  if (f$ok) "KKT-ok" else "NOCONV",
  if (length(f$boundary_nat))
    paste0(" | FRONTIER: ", paste(f$boundary_nat, collapse = ", ")) else "",
  if (isFALSE(f$nest_ok))
    sprintf(" | NESTING AUDIT: %s (%.1e)",
            if (isFALSE(f$nest_ok_lr)) "FAIL" else "strict-only",
            f$nest_gap) else "")

cub_pad_start <- function(f_prev, D, spec, rho_free = TRUE, sigma0 = 0.45) {
  co <- cub_core(D, spec, rho_free)
  v <- numeric(co$p); names(v) <- co$nms
  if (co$re) { v["log.L11"] <- log(sigma0); v["log.L22"] <- log(sigma0) }
  k <- intersect(names(f_prev$th), co$nms)
  v[k] <- f_prev$th[k]
  as.numeric(v)
}

#' Fit one rung of the ladder, warm-started from every rung it nests
#'
#' A nested model attaining a higher log-likelihood than the model nesting
#' it is an optimisation failure.  The audit is reported at two tolerances:
#' a strict 1e-6, and one anchored to the scale of the inference the
#' quadrature supports, max(1e-2, N * err(Q)), supplied via `nest_tol`.
cub_fit_rung <- function(D, spec, label, fits, parents = character(0),
                         rho_free = TRUE, Q = 15L, nstart = 1L,
                         nest_tol = 1e-2, gtol = 1e-3) {
  llref <- if (length(parents))
    max(vapply(fits[parents], `[[`, 0, "ll")) else -Inf
  starts <- if (!length(parents)) list(NULL) else c(
    lapply(parents, function(nm) cub_pad_start(fits[[nm]], D, spec, rho_free)),
    lapply(c(0.20, 0.75), function(s)
      cub_pad_start(fits[[parents[1]]], D, spec, rho_free, s)))
  best <- NULL; used <- 0L
  for (st in starts) {
    used <- used + 1L
    f <- try(cub_fit(D, spec, rho_free, Q = Q, nstart = nstart, start = st,
                     label = label, seed = 11L + 7L * used, gtol = gtol),
             silent = TRUE)
    if (inherits(f, "try-error")) next
    if (is.null(best) || f$ll > best$ll) best <- f
    if (best$ll >= llref - 1e-6 && best$ok) break
  }
  if (is.null(best)) stop(sprintf("cub_fit_rung(%s): total failure", label))
  gap <- if (is.finite(llref)) best$ll - llref else NA_real_
  best$nest_gap   <- gap
  best$nest_ok    <- if (is.finite(gap)) gap >= -1e-6 else NA
  best$nest_ok_lr <- if (is.finite(gap)) gap >= -nest_tol else NA
  best$n_starts   <- used
  best
}

## restriction maps -----------------------------------------------------

cub_map_fix <- function(f, k, val) {                  # fix coordinate k
  p <- f$co$p; idx <- setdiff(seq_len(p), k)
  M <- matrix(0, p, p - 1L); M[cbind(idx, seq_len(p - 1L))] <- 1
  e <- numeric(p); e[k] <- val
  list(fo = function(x) as.vector(M %*% x) + e, jac = M, nd = p - 1L,
       x0 = f$th[idx], lo = f$lo[idx], up = f$up[idx])
}
cub_map_equal <- function(f, k1, k2) {                # theta[k2] = theta[k1]
  p <- f$co$p; idx <- setdiff(seq_len(p), k2)
  M <- matrix(0, p, p - 1L); M[cbind(idx, seq_len(p - 1L))] <- 1
  M[k2, match(k1, idx)] <- 1
  list(fo = function(x) as.vector(M %*% x), jac = M, nd = p - 1L,
       x0 = f$th[idx], lo = f$lo[idx], up = f$up[idx])
}
cub_map_rho <- function(f, r) {            # fix rho_12 = r, free log sigma_2
  p <- f$co$p; iv <- f$co$iv; i2 <- iv[2]; i3 <- iv[3]
  idx <- setdiff(seq_len(p), c(i2, i3)); nf <- length(idx)
  L <- f$co$chol_of(f$th); s2 <- sqrt(L$L21^2 + L$L22^2)
  c0 <- log(max(s2, 1e-3)); k <- 0.5 * log(max(1 - r^2, 1e-12))
  list(
    fo = function(x) {
      t <- numeric(p); t[idx] <- x[1:nf]
      t[i2] <- x[nf + 1L] + k; t[i3] <- r * exp(x[nf + 1L]); t
    },
    jac = function(x) {
      J <- matrix(0, p, nf + 1L)
      J[cbind(idx, seq_len(nf))] <- 1
      J[i2, nf + 1L] <- 1; J[i3, nf + 1L] <- r * exp(x[nf + 1L]); J
    },
    nd = nf + 1L, x0 = c(f$th[idx], c0),
    lo = c(f$lo[idx], -6), up = c(f$up[idx], 1.5))
}

#' Profile interval for sigma^2_1, truncation reported as an inequality
cub_profile_s1 <- function(f, lev = 0.95, nexp = 5L, maxit = 7L,
                           budget = 150) {
  t0 <- Sys.time(); k <- f$co$iv[1]; cutv <- stats::qchisq(lev, 1) / 2
  pfx <- function(val) {
    r <- try(cub_fit(f$D, f$spec, f$rho_free, Q = f$Q, nstart = 1L,
                     tmap = cub_map_fix(f, k, val), label = "", maxit = 300L),
             silent = TRUE)
    if (inherits(r, "try-error")) -Inf else r$ll
  }
  side <- function(dir) {
    lo_ <- f$th[k]; step <- 0.4; hi_ <- NA_real_
    for (j in seq_len(nexp)) {
      if (.el(t0) > budget) return(c(f$th[k] + dir * step, 0))
      cand <- f$th[k] + dir * step
      if (f$ll - pfx(cand) >= cutv) { hi_ <- cand; break }
      lo_ <- cand; step <- step * 1.8
    }
    if (is.na(hi_)) return(c(f$th[k] + dir * step, 0))
    for (j in seq_len(maxit)) {
      if (.el(t0) > budget) break
      mid <- (lo_ + hi_) / 2
      if (f$ll - pfx(mid) < cutv) lo_ <- mid else hi_ <- mid
    }
    c((lo_ + hi_) / 2, 1)
  }
  Lo <- side(-1); Up <- side(+1)
  list(par = "sigma2_1", est = exp(2 * f$th[k]),
       lo = exp(2 * Lo[1]), hi = exp(2 * Up[1]),
       lo_ok = Lo[2] == 1, hi_ok = Up[2] == 1, secs = .el(t0))
}

cub_profile_s1_grid <- function(f, lo = 0.28, hi = 1.25, ngrid = 7L,
                                budget = 180) {
  k <- f$co$iv[1]; t0 <- Sys.time()
  g <- seq(log(lo), log(hi), length.out = ngrid)
  rows <- list()
  for (v in g) {
    if (.el(t0) > budget) break
    r <- try(cub_fit(f$D, f$spec, f$rho_free, Q = f$Q, nstart = 1L,
                     tmap = cub_map_fix(f, k, v), label = "", maxit = 300L),
             silent = TRUE)
    if (inherits(r, "try-error")) next
    rows[[length(rows) + 1L]] <- data.frame(
      sigma2_1 = exp(2 * v), loglik = r$ll, deviance = 2 * (f$ll - r$ll))
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

cub_fmt_profile <- function(p) sprintf(
  "%s = %.4f   95%% profile [%s, %s]", p$par, p$est,
  if (p$lo_ok) sprintf("%.4f", p$lo) else sprintf("< %.3g", p$lo),
  if (p$hi_ok) sprintf("%.4f", p$hi) else sprintf("> %.3g", p$hi))

## ---------------------------------------------------------------------
## 6.  observed information, clustered sandwich, natural components
## ---------------------------------------------------------------------

## Central differences of the ANALYTIC score: one order of accuracy lost
## instead of two.
cub_hess <- function(f, rel = 1e-4) {
  th <- f$th; p <- length(th); h <- rel * pmax(abs(th), 1)
  H <- matrix(0, p, p)
  for (i in seq_len(p)) {
    e <- numeric(p); e[i] <- h[i]
    H[, i] <- (f$co$ev(th + e, f$Q)$grad - f$co$ev(th - e, f$Q)$grad) / (2 * h[i])
  }
  (H + t(H)) / 2
}

cub_vcov <- function(f, cluster = TRUE, rel = 1e-4) {
  if (!is.null(f$tmap))
    stop("cub_vcov: restricted fit; standard errors are not reported")
  H <- cub_hess(f, rel)
  e <- eigen(H, symmetric = TRUE); pdef <- min(e$values) > 0
  Hi <- if (pdef) chol2inv(chol(H)) else {
    d <- e$values; d[abs(d) < max(abs(d)) * 1e-10] <- Inf
    e$vectors %*% diag(1 / d) %*% t(e$vectors)
  }
  V <- Hi; ssum <- NA_real_
  if (cluster) {
    Sc <- f$co$ev(f$th, f$Q, grad = FALSE, scores = TRUE)$scores
    ssum <- max(abs(colSums(Sc))); Nc <- nrow(Sc)
    V <- Hi %*% (crossprod(Sc) * Nc / max(Nc - 1, 1)) %*% Hi
  }
  dimnames(V) <- list(f$co$nms, f$co$nms)
  attr(V, "kappa") <- max(abs(e$values)) / min(abs(e$values))
  attr(V, "pd") <- pdef; attr(V, "score_sum") <- ssum
  V
}

cub_vc_natural <- function(f, V = NULL) {
  co <- f$co; if (!co$re) return(NULL)
  i <- co$iv; L <- co$chol_of(f$th)
  L11 <- L$L11; L22 <- L$L22; L21 <- L$L21
  s21 <- L11^2; s22 <- L21^2 + L22^2; s12 <- L11 * L21
  s2 <- sqrt(s22); rho <- if (s2 > 1e-12) L21 / s2 else NA_real_
  est <- c(sigma2_1 = s21, sigma2_2 = s22, sigma_12 = s12, rho_12 = rho)
  J <- if (co$rho_free)
    rbind(c(2 * s21, 0, 0),
          c(0, 2 * L22^2, 2 * L21),
          c(s12, 0, L11),
          c(0,
            if (s2 > 1e-12) -L21 * L22^2 / s2^3 else 0,
            if (s2 > 1e-12)  L22^2 / s2^3 else 0))
  else rbind(c(2 * s21, 0), c(0, 2 * L22^2), c(0, 0), c(0, 0))
  Vn <- NULL; se <- rep(NA_real_, 4L)
  if (!is.null(V)) {
    Vn <- J %*% V[i, i, drop = FALSE] %*% t(J)
    se <- sqrt(pmax(diag(Vn), 0))
  }
  bad <- length(cub_nat_bnd(f$th, co, f$lo)) > 0L
  if (bad) se[3:4] <- NA_real_
  if (!co$rho_free) se[3:4] <- NA_real_
  list(est = est, se = se, J = J, Vnat = Vn,
       V12 = if (!is.null(Vn)) Vn[1:3, 1:3, drop = FALSE] else NULL,
       boundary = bad)
}

## ---------------------------------------------------------------------
## 7.  average marginal effects, on the own scale AND on E(R)
## ---------------------------------------------------------------------

#' Two average marginal effects per coefficient
#'
#' AME is the average derivative of the quantity named in the coefficient's
#' own equation, on the (0,1) scale.  AME_ER is the average derivative of
#' E(R_it) on the category scale.  They routinely disagree in sign: raising
#' the uncertainty weight pulls the conditional mean towards (m+1)/2, so the
#' sign depends on which side of the midpoint the fitted mean lies.
cub_ame <- function(f, V = NULL, Q = f$Q) {
  co <- f$co; m <- co$D$m; cs <- co$D$cs; n <- co$D$n; p <- co$p
  amef <- function(th) {
    lam0 <- as.vector(co$X %*% th[co$ib])
    eta0 <- if (co$use_u) as.vector(co$W %*% th[co$ig]) else NULL
    del  <- if (co$use_s) .pl(as.vector(co$Z %*% th[co$ia])) else rep(0, n)
    if (co$re) {
      G <- cub_gh2(Q); L <- co$chol_of(th)
      uq <- L$L11 * G$z1; vq <- L$L21 * G$z1 + L$L22 * G$z2; w <- G$w
    } else { uq <- 0; vq <- 0; w <- 1 }
    XI  <- .pl(outer(lam0, uq, "+"))
    UPS <- if (co$use_u) .pl(outer(eta0, vq, "+")) else matrix(0, n, length(uq))
    mub <- 1 + (m - 1) * XI
    muG <- (1 - UPS) * mub + UPS * (m + 1) / 2
    d_xi_o <- mean(as.vector((XI * (1 - XI)) %*% w))
    d_up_o <- if (co$use_u) mean(as.vector((UPS * (1 - UPS)) %*% w)) else 0
    d_dl_o <- mean(del * (1 - del))
    d_xi_y <- mean(as.vector(((1 - del) * ((1 - UPS) * (m - 1) * XI * (1 - XI))) %*% w))
    d_up_y <- if (co$use_u)
      mean(as.vector(((1 - del) * (((m + 1) / 2 - mub) * UPS * (1 - UPS))) %*% w)) else 0
    d_dl_y <- mean(del * (1 - del) * (cs - as.vector(muG %*% w)))
    o <- numeric(2 * p)
    o[co$ib] <- d_xi_o * th[co$ib]; o[p + co$ib] <- d_xi_y * th[co$ib]
    if (co$use_u) { o[co$ig] <- d_up_o * th[co$ig]; o[p + co$ig] <- d_up_y * th[co$ig] }
    if (co$use_s) { o[co$ia] <- d_dl_o * th[co$ia]; o[p + co$ia] <- d_dl_y * th[co$ia] }
    o
  }
  a <- amef(f$th)
  out <- data.frame(AME = a[seq_len(p)], se_AME = NA_real_, p_AME = NA_real_,
                    AME_ER = a[p + seq_len(p)], row.names = co$nms)
  if (is.null(V)) return(out)
  J <- vapply(seq_len(p), function(i) {
    h <- 1e-5 * max(abs(f$th[i]), 1); e <- numeric(p); e[i] <- h
    (amef(f$th + e) - amef(f$th - e)) / (2 * h)
  }, numeric(2 * p))
  se <- sqrt(pmax(diag(J %*% V %*% t(J)), 0))[seq_len(p)]
  out$se_AME <- se
  out$p_AME <- 2 * stats::pnorm(-abs(out$AME / se))
  out
}

## ---------------------------------------------------------------------
## 8.  coefficient table
## ---------------------------------------------------------------------

cub_coeftab <- function(f, V = NULL, A = NULL, truth = NULL) {
  co <- f$co
  if (is.null(V)) V <- cub_vcov(f)
  se <- sqrt(pmax(diag(V), 0))
  bd <- co$nms %in% sub("\\[(lo|up)\\]$", "", f$boundary)
  z  <- ifelse(bd, NA_real_, f$th / se)
  tb <- data.frame(est = f$th, se = ifelse(bd, NA_real_, se), z = z,
                   p = 2 * stats::pnorm(-abs(z)), row.names = co$nms)
  if (!is.null(A)) {
    tb$AME <- A$AME; tb$p_AME <- A$p_AME; tb$AME_ER <- A$AME_ER
    ic <- grepl("\\(Intercept\\)$", co$nms)
    tb$AME[ic] <- NA; tb$p_AME[ic] <- NA; tb$AME_ER[ic] <- NA
    if (length(co$iv)) {
      tb$AME[co$iv] <- NA; tb$p_AME[co$iv] <- NA; tb$AME_ER[co$iv] <- NA
    }
  }
  vn <- cub_vc_natural(f, V)
  if (!is.null(vn)) {
    nat <- data.frame(est = vn$est, se = vn$se, z = NA_real_, p = NA_real_,
                      row.names = names(vn$est))
    if (!is.null(A)) {
      nat$AME <- NA_real_; nat$p_AME <- NA_real_; nat$AME_ER <- NA_real_
    }
    tb <- rbind(tb, nat)
  }
  cover <- NA_integer_; ncmp <- NA_integer_
  if (!is.null(truth)) {
    tb <- cbind(True = truth[match(rownames(tb), names(truth))], tb)
    cmp <- is.finite(tb$True) & is.finite(tb$se)
    ncmp <- sum(cmp)
    cover <- sum(abs(tb$True[cmp] - tb$est[cmp]) / tb$se[cmp] < 1.96)
  }
  list(tab = tb, kappa = attr(V, "kappa"), pd = attr(V, "pd"),
       score_sum = attr(V, "score_sum"), V = V, vcn = vn,
       cover = cover, ncomparable = ncmp, label = f$label,
       kkt_informative = f$kkt_informative)
}

## ---------------------------------------------------------------------
## 9.  fit assessment, ladder, nested chains, tests
## ---------------------------------------------------------------------

cub_fit_stats <- function(f, Qpred = f$Q, quant = FALSE) {
  E <- f$co$ev(f$th, Qpred, grad = FALSE, pred = TRUE, quant = quant)
  Pm <- E$pred$Pm; D <- f$D; m <- D$m
  obs  <- as.numeric(table(factor(D$y, levels = seq_len(m)))) / D$n
  fitp <- colMeans(Pm)
  mu <- as.vector(Pm %*% seq_len(m))
  v2 <- as.vector(Pm %*% (seq_len(m)^2)) - mu^2
  dw <- tapply(seq_len(D$n), D$wave, function(k)
    0.5 * sum(abs(as.numeric(table(factor(D$y[k], levels = seq_len(m)))) / length(k) -
                    colMeans(Pm[k, , drop = FALSE]))))
  nb_target <- max(3L, min(10L, floor(D$n / 15)))
  brk <- unique(stats::quantile(mu, seq(0, 1, length.out = nb_target + 1L)))
  q <- cut(mu, brk, include.lowest = TRUE)
  nb <- as.integer(table(q))
  br <- data.frame(bin = levels(q), n = nb,
                   mu = as.numeric(tapply(mu, q, mean)),
                   resid = as.numeric(tapply(D$y - mu, q, mean)),
                   band = 2 * as.numeric(tapply(sqrt(v2), q, mean)) / sqrt(pmax(nb, 1)))
  br <- br[!is.na(br$mu), ]
  list(diss = 0.5 * sum(abs(obs - fitp)), obs = obs, fit = fitp,
       diss_wave = dw, rowsum = max(abs(rowSums(Pm) - 1)),
       binres = br, mu = mu, nbins = nrow(br),
       out_band = sum(abs(br$resid) > br$band),
       logscore = -mean(log(pmax(Pm[cbind(seq_len(D$n), D$y)], 1e-300))),
       EV = mean(v2), VE = mean((mu - mean(mu))^2),
       var_obs = stats::var(D$y),
       floor_var = mean((mu - 1) * (m - mu) / (m - 1)),
       post = E$post, pred = E$pred, Qpred = Qpred)
}

cub_compare <- function(fits, diss = NULL) {
  tb <- data.frame(
    spec = vapply(fits, `[[`, "", "label"),
    p = vapply(fits, `[[`, 0, "p"),
    loglik = vapply(fits, `[[`, 0, "ll"),
    AIC = vapply(fits, `[[`, 0, "aic"),
    BIC_N = vapply(fits, `[[`, 0, "bicN"),
    Diss = if (is.null(diss)) NA_real_ else as.numeric(diss),
    KKT = vapply(fits, function(f) if (f$ok) "ok" else "NO", ""),
    frontier = vapply(fits, function(f)
      if (length(f$boundary_nat)) paste(f$boundary_nat, collapse = "; ") else "-", ""),
    stringsAsFactors = FALSE)
  tb$best_AIC   <- ifelse(seq_len(nrow(tb)) == which.min(tb$AIC), "<<", "")
  tb$best_BIC_N <- ifelse(seq_len(nrow(tb)) == which.min(tb$BIC_N), "<<", "")
  rownames(tb) <- names(fits)
  list(table = tb,
       same_sample = length(unique(vapply(fits, `[[`, "", "hash"))) == 1L,
       best_bicN = names(fits)[which.min(tb$BIC_N)],
       best_aic  = names(fits)[which.min(tb$AIC)])
}

cub_chain <- function(fits, diss = NULL, nest_tol = 1e-2) {
  ll <- vapply(fits, `[[`, 0, "ll"); p <- vapply(fits, `[[`, 0, "p")
  tb <- data.frame(spec = vapply(fits, `[[`, "", "label"), p = p, loglik = ll,
                   Diss = if (is.null(diss)) NA_real_ else as.numeric(diss),
                   LR = c(NA, pmax(2 * diff(ll), 0)), df = c(NA, diff(p)),
                   stringsAsFactors = FALSE)
  tb$p_LR <- ifelse(is.na(tb$LR), NA,
                    stats::pchisq(tb$LR, pmax(tb$df, 1), lower.tail = FALSE))
  rownames(tb) <- names(fits)
  list(table = tb, min_step = min(diff(ll)),
       audit_strict = all(diff(ll) > -1e-6),
       audit_lr = all(diff(ll) > -nest_tol),
       nest_tol = nest_tol,
       same_sample = length(unique(vapply(fits, `[[`, "", "hash"))) == 1L)
}

cub_wald <- function(f, V, nms, what = "") {
  miss <- setdiff(nms, f$co$nms)
  if (length(miss) || !length(nms)) return(NULL)
  i <- match(nms, f$co$nms)
  st <- tryCatch(as.numeric(t(f$th[i]) %*% solve(V[i, i], f$th[i])),
                 error = function(e) NA_real_)
  c(W = st, df = length(i),
    p = stats::pchisq(st, length(i), lower.tail = FALSE))
}

#' chi-bar-squared p-value for H0: Sigma = 0
#'
#' The reference law is obtained by exact projection of draws from the
#' estimated sampling distribution of (sigma^2_1, sigma^2_2, sigma_12) onto
#' the cone of positive semi-definite 2x2 matrices, in the metric of the
#' estimated information.  The textbook (1/4, 1/2, 1/4) weights do not apply
#' here: they are for two independent variances, not for a matrix with a
#' correlation.
cub_chibar <- function(V12, lr, B = 5000L, seed = 5L, ngrid = 721L) {
  bad <- c(p = NA_real_, mean_stat = NA_real_, B = B)
  if (is.null(V12) || any(!is.finite(V12))) return(bad)
  V12 <- (V12 + t(V12)) / 2
  e <- eigen(V12, symmetric = TRUE)
  if (min(e$values) <= 0) return(bad)
  A <- e$vectors %*% diag(1 / e$values) %*% t(e$vectors)
  L <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
  an <- seq(0, pi, length.out = ngrid + 1L)[-(ngrid + 1L)]
  Sv <- rbind(cos(an)^2, sin(an)^2, cos(an) * sin(an))
  AS <- A %*% Sv; den <- colSums(Sv * AS)
  set.seed(seed)
  stat <- vapply(seq_len(B), function(b) {
    z <- as.vector(L %*% stats::rnorm(3))
    q0 <- as.numeric(t(z) %*% A %*% z)
    if (z[1] >= 0 && z[2] >= 0 && z[1] * z[2] - z[3]^2 >= 0) return(q0)
    nmr <- as.vector(crossprod(z, AS)); t_ <- pmax(nmr / den, 0)
    max(q0 - min(q0 - 2 * t_ * nmr + t_^2 * den), 0)
  }, 0)
  c(p = (1 + sum(stat >= lr - 1e-9)) / (1 + B), mean_stat = mean(stat), B = B)
}

## ---------------------------------------------------------------------
## 10. correlation audit: chi^2_1 test plus a profile
## ---------------------------------------------------------------------

cub_rho_audit <- function(f_free, f_diag,
                          grid = c(-0.3, 0, 0.3, 0.5, 0.7, 0.9),
                          budget = 160) {
  stopifnot(f_free$co$rho_free, !f_diag$co$rho_free)
  lr <- max(2 * (f_free$ll - f_diag$ll), 0)
  pv <- stats::pchisq(lr, 1, lower.tail = FALSE)
  L <- f_free$co$chol_of(f_free$th); s2 <- sqrt(L$L21^2 + L$L22^2)
  rh <- if (s2 > 0) L$L21 / s2 else NA_real_
  t0 <- Sys.time(); prof <- NULL; truncated <- FALSE
  if (is.finite(budget) && budget > 0) {
    rows <- list()
    for (r in grid) {
      if (.el(t0) > budget) { truncated <- TRUE; break }
      fr <- try(cub_fit(f_free$D, f_free$spec, TRUE, Q = f_free$Q, nstart = 1L,
                        tmap = cub_map_rho(f_free, r), label = "",
                        maxit = 300L), silent = TRUE)
      if (inherits(fr, "try-error")) next
      rows[[length(rows) + 1L]] <- data.frame(
        rho = r, loglik = fr$ll, deviance = 2 * (f_free$ll - fr$ll))
    }
    if (length(rows)) prof <- do.call(rbind, rows)
  }
  ci <- c(NA_real_, NA_real_)
  if (!is.null(prof) && nrow(prof) >= 3L) {
    ok <- prof$deviance <= stats::qchisq(0.95, 1)
    if (any(ok)) ci <- range(prof$rho[ok])
  }
  list(lr = lr, p = pv, rho = rh, prof = prof, ci = ci,
       truncated = truncated, grid = grid,
       informative = lr >= stats::qchisq(0.95, 1), secs = .el(t0))
}

## ---------------------------------------------------------------------
## 11. self-checks.  C1-C6 and C8 are gates; C7 is a reported diagnostic
##     whose verdict is anchored to the scale of the inference.
## ---------------------------------------------------------------------

cub_self_checks <- function(D, spec = NULL, Q = 15L, Qtol = 1e-2,
                            Qcheck = c(9L, 15L, 21L, 31L)) {
  if (is.null(spec)) spec <- if (D$info$lag_ok) "RE-LAG" else "RE"
  rows <- list()
  add <- function(id, what, val, tol, gate = TRUE) {
    rows[[length(rows) + 1L]] <<- data.frame(
      id = id, identity = what, value = val, tolerance = tol,
      gate = gate, pass = is.finite(val) && val < tol,
      stringsAsFactors = FALSE)
  }

  G <- cub_gh2(15L)
  add("C1", "bivariate GH moments to order 4, and orthogonality",
      max(abs(c(sum(G$w) - 1, sum(G$w * G$z1), sum(G$w * G$z1^2) - 1,
                sum(G$w * G$z1^4) - 3, sum(G$w * G$z2^2) - 1,
                sum(G$w * G$z1 * G$z2)))), 1e-10)

  co <- cub_core(D, spec, TRUE)
  set.seed(4)
  th <- c(cub_em_pooled(co), log(.65), log(.55), 0.25) +
    stats::rnorm(co$p, 0, .12)

  add("C2", "sum_r p_itr = 1 at every row",
      max(abs(rowSums(co$ev(th, 15L, grad = FALSE, pred = TRUE)$pred$Pm) - 1)),
      1e-11)

  lim <- function(sp_re, sp_pool) {
    cr <- cub_core(D, sp_re, TRUE); cp <- cub_core(D, sp_pool)
    stopifnot(identical(cp$nms, cr$nms[seq_len(cp$p)]))
    set.seed(5); t0 <- cub_em_pooled(cp) + stats::rnorm(cp$p, 0, .12)
    abs(cr$ev(c(t0, -10, -10, 0), 15L, grad = FALSE)$ll -
          cp$ev(t0, 1L, grad = FALSE)$ll)
  }
  add("C3a", "Sigma -> 0: RE rung equals the pooled rung", lim("RE", "0"), 1e-6)
  if (D$info$lag_ok)
    add("C3b", "Sigma -> 0: RE-LAG rung equals the LAG rung",
        lim("RE-LAG", "LAG"), 1e-6)

  g_an <- co$ev(th, 11L)$grad
  g_nu <- vapply(seq_along(th), function(i) {
    h <- 1e-5 * max(abs(th[i]), 1); e <- numeric(length(th)); e[i] <- h
    (co$ev(th + e, 11L, grad = FALSE)$nll -
       co$ev(th - e, 11L, grad = FALSE)$nll) / (2 * h)
  }, 0)
  add("C4", "analytic gradient against central differences",
      max(abs(g_an - g_nu) / pmax(abs(g_nu), 1)), 1e-6)
  add("C5", "vectorised likelihood against a looped reference",
      abs(co$ev(th, 9L, grad = FALSE)$ll - co$ev_ref(th, 9L)), 1e-9)
  add("C6", "sum_i score_i equals the total gradient",
      max(abs(-colSums(co$ev(th, 11L, scores = TRUE)$scores) - g_an)), 1e-9)

  L <- co$chol_of(th); G8 <- cub_gh2(21L)
  u8 <- L$L11 * G8$z1; v8 <- L$L21 * G8$z1 + L$L22 * G8$z2
  add("C8", "covariance implied by the scaled nodes equals Sigma",
      max(abs(c(sum(G8$w * u8^2) - L$L11^2,
                sum(G8$w * v8^2) - (L$L21^2 + L$L22^2),
                sum(G8$w * u8 * v8) - L$L11 * L$L21))), 1e-10)

  ## C7: GH(Q) against a deterministic 301 x 301 product on a few subjects
  lam0 <- as.vector(co$X %*% th[co$ib])
  eta0 <- as.vector(co$W %*% th[co$ig])
  del  <- .pl(as.vector(co$Z %*% th[co$ia])); m <- D$m
  set.seed(9); ss <- sample(D$N, min(5L, D$N))
  Gq <- 301L; a <- seq(-6, 6, length.out = Gq); hh <- a[2] - a[1]
  Vg <- L$L21 * matrix(a, Gq, Gq) + L$L22 * matrix(a, Gq, Gq, byrow = TRUE)
  lpri <- outer(stats::dnorm(a, log = TRUE), stats::dnorm(a, log = TRUE), "+") +
    2 * log(hh)
  exact <- vapply(ss, function(i) {
    k <- which(D$gi == i); S <- lpri
    for (tt in k) {
      bb <- stats::dbinom(D$y[tt] - 1L, m - 1L,
                          stats::plogis(lam0[tt] + L$L11 * a))
      UPq <- stats::plogis(eta0[tt] + Vg)
      S <- S + log(pmax(del[tt] * (D$y[tt] == D$cs) +
                          (1 - del[tt]) * ((1 - UPq) * bb + UPq / m), 1e-300))
    }
    .lse(as.vector(S))
  }, 0)
  gh_ll <- function(i, Qx) {
    G <- cub_gh2(Qx); k <- which(D$gi == i)
    uq <- L$L11 * G$z1; vq <- L$L21 * G$z1 + L$L22 * G$z2
    U <- outer(lam0[k], uq, "+"); Vm <- outer(eta0[k], vq, "+")
    Bk <- exp(lchoose(m - 1L, D$y[k] - 1L) + (D$y[k] - 1) * .pll(U) +
                (m - D$y[k]) * .pll(-U))
    Fk <- del[k] * (D$y[k] == D$cs) +
      (1 - del[k]) * ((1 - .pl(Vm)) * Bk + .pl(Vm) / m)
    .lse(colSums(log(pmax(Fk, 1e-300))) + G$lw)
  }
  Qs <- sort(unique(as.integer(c(Qcheck, Q))))
  tabQ <- data.frame(
    Q = Qs,
    err_per_subject = vapply(Qs, function(Qx)
      max(abs(vapply(ss, function(i) gh_ll(i, Qx), 0) - exact)), 0))
  tabQ$implied_aggregate <- D$N * tabQ$err_per_subject
  eQ <- tabQ$err_per_subject[tabQ$Q == Q]
  monotone <- all(diff(tabQ$err_per_subject) <= 0)
  add("C7", sprintf("quadrature accuracy at Q = %d (aggregate N*err)", Q),
      D$N * eQ, Qtol, gate = FALSE)

  tab <- do.call(rbind, rows)
  list(table = tab, quad = tabQ, monotone = monotone,
       err_per_subject = eQ, implied_aggregate = D$N * eQ,
       Qtol = Qtol, spec = spec,
       ok = all(tab$pass[tab$gate]) && monotone,
       Q_adequate = D$N * eQ < Qtol,
       Q_recommended = {
         ok <- tabQ$Q[tabQ$implied_aggregate < Qtol]
         if (length(ok)) min(ok) else max(tabQ$Q)
       })
}

## ---------------------------------------------------------------------
## 12. out-of-sample: leakage-free split, benchmark on the SAME information
## ---------------------------------------------------------------------

cub_oos_prep <- function(D) {
  last <- stats::ave(D$wave_idx, D$gi, FUN = max) == D$wave_idx
  tr <- !last
  ## recentre every constructed feature on the training rows only
  D$Xl[, 1] <- D$Xl[, 1] - mean(D$Xl[tr, 1])
  D$Wl[, 1] <- D$Wl[, 1] - mean(D$Wl[tr, 1])
  srcv <- function(j) {
    b <- sub("_bar$", "", j)
    if (b == "lag_f") return(D$Xl[, 1])
    if (b == "lag_u") return(D$Wl[, 1])
    if (b %in% colnames(D$Xb)) return(D$Xb[, b])
    if (b %in% colnames(D$Wb)) return(D$Wb[, b])
    NULL
  }
  for (nm in c("Xm", "Wm")) {
    M <- D[[nm]]; if (is.null(M) || !ncol(M)) next
    for (j in colnames(M)) {
      v <- if (grepl("_bar$", j)) srcv(j) else NULL
      if (!is.null(v)) {
        mu <- tapply(v[tr], D$gi[tr], mean)
        w <- as.numeric(mu[as.character(D$gi)])
        M[, j] <- w - mean(w[tr])
      } else M[, j] <- M[, j] - mean(M[tr, j])
    }
    D[[nm]] <- M
  }
  list(train = cub_subset(D, tr), test = cub_subset(D, last))
}

## exactly the regressors the mixture sees, intercept excluded
cub_info_set <- function(D, spec, rho_free = TRUE) {
  co <- cub_core(D, spec, rho_free)
  M <- cbind(co$X, if (co$use_u) co$W, if (co$use_s) co$Z)
  M <- M[, colnames(M) != "(Intercept)", drop = FALSE]
  M[, !duplicated(colnames(M)), drop = FALSE]
}

cub_pred_cond <- function(f_tr, D_te, Qp) {
  co_te <- cub_core(D_te, f_tr$spec, f_tr$rho_free)
  E  <- f_tr$co$ev(f_tr$th, Qp, grad = FALSE, post = TRUE)$post
  th <- f_tr$th; m <- D_te$m; n <- D_te$n
  lam0 <- as.vector(co_te$X %*% th[co_te$ib])
  eta0 <- if (co_te$use_u) as.vector(co_te$W %*% th[co_te$ig]) else NULL
  del  <- if (co_te$use_s) .pl(as.vector(co_te$Z %*% th[co_te$ia])) else rep(0, n)
  P <- matrix(0, n, m); rv <- 0:(m - 1L); lchv <- lchoose(m - 1L, rv)
  for (q in seq_along(E$uq)) {
    lp <- .pll(lam0 + E$uq[q]); lq <- .pll(-(lam0 + E$uq[q]))
    BM <- exp(outer(lp, rv) + outer(lq, rev(rv)) + rep(lchv, each = n))
    up <- if (co_te$use_u) .pl(eta0 + E$vq[q]) else 0
    P  <- P + E$Wq[D_te$gi, q] * ((1 - up) * BM + up / m)
  }
  P <- (1 - del) * P
  if (co_te$use_s) P[, D_te$cs] <- P[, D_te$cs] + del
  P / rowSums(P)
}

## random-effects proportional-odds benchmark, reimplemented so that it
## receives exactly the mixture's information set
cub_po_fit <- function(y, X, gi, m, Q = 20L) {
  g <- cub_gh1(Q); k <- ncol(X); N <- max(gi)
  nll <- function(th) {
    cut <- c(-Inf, cumsum(c(th[1], exp(th[2:(m - 1)]))), Inf)
    xb <- as.vector(X %*% th[m - 1 + seq_len(k)])
    s <- exp(th[m - 1 + k + 1])
    A <- matrix(0, N, Q)
    for (q in seq_len(Q)) {
      e <- xb + s * g$z[q]
      pr <- stats::plogis(cut[y + 1L] - e) - stats::plogis(cut[y] - e)
      A[, q] <- .rsum(log(pmax(pr, 1e-300)), gi)[, 1] + log(g$w[q])
    }
    mx <- A[cbind(seq_len(N), max.col(A, ties.method = "first"))]
    -sum(mx + log(rowSums(exp(A - mx))))
  }
  th0 <- c(stats::qlogis(min(max(mean(y <= 1) + .02, .01), .5)),
           rep(log(0.8), m - 2), numeric(k), log(0.5))
  r <- stats::optim(th0, nll, method = "BFGS",
                    control = list(maxit = 800L, reltol = 1e-11))
  r2 <- stats::optim(r$par, nll, method = "BFGS",
                     control = list(maxit = 800L, reltol = 1e-12))
  if (r2$value < r$value) r <- r2
  list(th = r$par, ll = -r$value, Q = Q, m = m, k = k, g = g)
}

cub_po_pred <- function(fit, ytr, Xtr, gtr, Xte, gte) {
  th <- fit$th; m <- fit$m; k <- fit$k; g <- fit$g; Q <- fit$Q; N <- max(gtr)
  cut <- c(-Inf, cumsum(c(th[1], exp(th[2:(m - 1)]))), Inf)
  s <- exp(th[m - 1 + k + 1]); bt <- th[m - 1 + seq_len(k)]
  xtr <- as.vector(Xtr %*% bt); A <- matrix(0, N, Q)
  for (q in seq_len(Q)) {
    e <- xtr + s * g$z[q]
    pr <- stats::plogis(cut[ytr + 1L] - e) - stats::plogis(cut[ytr] - e)
    A[, q] <- .rsum(log(pmax(pr, 1e-300)), gtr)[, 1] + log(g$w[q])
  }
  mx <- A[cbind(seq_len(N), max.col(A, ties.method = "first"))]
  W <- exp(A - mx); W <- W / rowSums(W)
  xte <- as.vector(Xte %*% bt); P <- matrix(0, length(xte), m)
  for (q in seq_len(Q)) {
    e <- xte + s * g$z[q]
    for (r in seq_len(m))
      P[, r] <- P[, r] + W[gte, q] *
        (stats::plogis(cut[r + 1L] - e) - stats::plogis(cut[r] - e))
  }
  P / rowSums(P)
}

cub_score_pred <- function(P, y) {
  m <- ncol(P); CP <- t(apply(P, 1, cumsum))
  I <- outer(y, seq_len(m), "<=") * 1
  c(logscore = -mean(log(pmax(P[cbind(seq_along(y), y)], 1e-300))),
    RPS = mean(rowSums((CP - I)^2)),
    hit = mean(apply(P, 1, which.max) == y),
    diss = 0.5 * sum(abs(
      as.numeric(table(factor(y, levels = seq_len(m)))) / length(y) -
        colMeans(P))))
}

## ---------------------------------------------------------------------
## 13. memory profile: the eccentricity feature is V-shaped, so the
##     coefficient rho_u does not summarise what the dynamics do
## ---------------------------------------------------------------------

cub_memory_profile <- function(f, Q = f$Q) {
  co <- f$co; D <- co$D; m <- D$m; th <- f$th
  if (!co$lag) return(NULL)
  ctr <- (m + 1) / 2; hw <- (m - 1) / 2
  xb <- colMeans(D$Xb); wb <- colMeans(D$Wb); zb <- colMeans(D$Z)
  nxb <- length(xb); nwb <- length(wb)
  bb <- th[co$ib]; gg <- if (co$use_u) th[co$ig] else NULL
  aa <- if (co$use_s) th[co$ia] else NULL
  G <- cub_gh2(Q); L <- co$chol_of(th)
  uq <- if (co$re) L$L11 * G$z1 else 0
  vq <- if (co$re) L$L21 * G$z1 + L$L22 * G$z2 else 0
  wq <- if (co$re) G$w else 1
  del <- if (co$use_s) stats::plogis(sum(zb * aa)) else 0
  out <- t(vapply(seq_len(m), function(r) {
    lf <- (r - ctr) / hw - D$ctr[["lag_f"]]
    lu <- abs(r - ctr) / hw - D$ctr[["lag_u"]]
    lam0 <- sum(xb * bb[seq_len(nxb)]) + lf * bb[nxb + 1L]
    phi <- .pl(lam0 + uq)
    ups <- if (co$use_u)
      .pl(sum(wb * gg[seq_len(nwb)]) + lu * gg[nwb + 1L] + vq) else 0
    muG <- (1 - ups) * (1 + (m - 1) * phi) + ups * (m + 1) / 2
    c(previous = r, eccentricity = abs(r - ctr) / hw,
      phi = sum(wq * phi), upsilon = sum(wq * ups), delta = del,
      ER = (1 - del) * sum(wq * muG) + del * D$cs,
      share_evaluation = (1 - del) * (1 - sum(wq * ups)),
      share_indifference = (1 - del) * sum(wq * ups),
      share_shelter = del)
  }, numeric(9)))
  as.data.frame(out)
}

## ---------------------------------------------------------------------
## 14. panel generation
## ---------------------------------------------------------------------

cub_sim_covariates <- function(N, Tn, seed) {
  set.seed(seed)
  gen <- factor(sample(c("Female", "Male"), N, TRUE),
                levels = c("Female", "Male"))
  edu <- factor(sample(c("Low", "Medium", "High"), N, TRUE),
                levels = c("Low", "Medium", "High"))
  cv <- function() rep(stats::rnorm(N), each = Tn) +
    stats::rnorm(N * Tn, 0, 1.05)
  data.frame(Respondent = rep(seq_len(N), each = Tn),
             Period = rep(seq_len(Tn), times = N),
             Gender = rep(gen, each = Tn), Education = rep(edu, each = Tn),
             Experience = cv(), Inflation = cv(),
             Delivery_Speed = cv(), Customer_Service = cv())
}

## Panel B: all three components carry non-negligible weight; the shelter
## sits at the scale midpoint, which is the harder configuration.
CUB_TRUTH_B <- list(
  beta = c("(Intercept)" = 0.00, "GenderMale" = 0.30,
           "Experience" = 0.25, "Inflation" = -0.20),
  rhof = 0.35,
  gam  = c("(Intercept)" = -0.50, "EducationMedium" = -0.30,
           "EducationHigh" = -0.60, "Delivery_Speed" = -0.35,
           "Customer_Service" = -0.30),
  rhou = -0.90,
  alp  = c("(Intercept)" = -1.50, "Inflation" = 0.40),
  s1 = 0.49, s2 = 0.36, rho12 = 0.50)

cub_sim_panelB <- function(N = 300L, Tn = 7L, m = 7L, cs = 4L,
                           seed = 20240701L, seed_cov = 20240101L,
                           truth = CUB_TRUTH_B) {
  tr <- truth
  d <- cub_sim_covariates(N, Tn, seed_cov)
  set.seed(seed)
  S <- matrix(c(tr$s1, tr$rho12 * sqrt(tr$s1 * tr$s2),
                tr$rho12 * sqrt(tr$s1 * tr$s2), tr$s2), 2, 2)
  UV <- matrix(stats::rnorm(2 * N), N, 2) %*% chol(S)
  Xb <- cbind(1, as.integer(d$Gender) - 1L, d$Experience, d$Inflation)
  Wb <- cbind(1, as.integer(d$Education == "Medium"),
              as.integer(d$Education == "High"),
              d$Delivery_Speed, d$Customer_Service)
  Zb <- cbind(1, d$Inflation)
  lam0 <- as.vector(Xb %*% tr$beta); eta0 <- as.vector(Wb %*% tr$gam)
  del  <- stats::plogis(as.vector(Zb %*% tr$alp))
  ctr <- (m + 1) / 2; hw <- (m - 1) / 2
  Y <- matrix(NA_integer_, N, Tn); idx <- matrix(seq_len(N * Tn), Tn, N)
  for (t in seq_len(Tn)) {
    k <- idx[t, ]
    if (t == 1L) { lf <- 0; lu <- 0 } else {
      lf <- (Y[, t - 1L] - ctr) / hw; lu <- abs(Y[, t - 1L] - ctr) / hw
    }
    lam <- lam0[k] + tr$rhof * lf + UV[, 1]
    eta <- eta0[k] + tr$rhou * lu + UV[, 2]
    up <- stats::plogis(eta)
    u1 <- stats::runif(N); u2 <- stats::runif(N)
    yb <- 1L + stats::rbinom(N, m - 1L, stats::plogis(lam))
    Y[, t] <- as.integer(ifelse(u1 < del[k], cs,
                                ifelse(u2 < up,
                                       1L + floor(stats::runif(N) * m), yb)))
  }
  d$Response <- as.integer(as.vector(t(Y)))
  attr(d, "truth") <- tr; attr(d, "re") <- UV
  attr(d, "shelter_cat") <- cs; attr(d, "panel") <- "B"
  d
}

## Panel A: conditionally UNDER-dispersed, no latent effect, no dynamics.
## The calibration VERIFIES the exact collapse rather than assuming it.
cub_sim_panelA <- function(N = 50L, Tn = 4L, m = 7L, cs = 7L,
                           seed_cov = 20240108L, nstart_pool = 3L,
                           target_R = 0.55) {
  d0 <- cub_sim_covariates(N, Tn, seed_cov)
  mu0 <- 4.30 + 0.90 * (as.integer(d0$Gender) - 1L) +
    0.35 * d0$Experience - 0.20 * d0$Inflation
  mkD <- function(y) {
    d <- d0; d$Response <- as.integer(y)
    cub_prepare(d, shelter_cat = cs)
  }
  cand <- list()
  for (se in c(0.80, 0.86, 0.92)) for (sd_ in 1:60) {
    set.seed(50000L + 100L * sd_ + round(se * 100))
    y <- pmin(pmax(round(mu0 + stats::rnorm(N * Tn, 0, se)), 1L), m)
    if (1L %in% y || !all(2:m %in% y)) next
    fd <- cub_floor_diag(mkD(y))
    if (!is.finite(fd[["R"]]) || fd[["R"]] > 0.85) next
    cand[[length(cand) + 1L]] <- list(
      y = y, se = se, seed = sd_, R = fd[["R"]],
      score = abs(fd[["R"]] - target_R))
  }
  if (!length(cand)) stop("cub_sim_panelA: no candidate met the criteria")
  ord <- order(vapply(cand, `[[`, 0, "score")); best <- NULL
  for (j in ord[seq_len(min(24L, length(ord)))]) {
    ca <- cand[[j]]; D <- mkD(ca$y)
    fsb <- try(cub_fit(D, "SB", nstart = 2L), silent = TRUE)
    f0  <- try(cub_fit(D, "0", nstart = nstart_pool), silent = TRUE)
    if (inherits(fsb, "try-error") || inherits(f0, "try-error")) next
    st <- f0$co$ev(f0$th, 1L, grad = FALSE, post = TRUE)$post
    ca$gap <- f0$ll - fsb$ll
    ca$max_ups <- max(st$ups); ca$max_del <- max(st$del)
    ca$pass <- (ca$gap < 1e-2 && ca$max_ups < 1e-2 && ca$max_del < 1e-2)
    if (is.null(best) || (ca$pass && !isTRUE(best$pass)) ||
        (isTRUE(ca$pass) == isTRUE(best$pass) && ca$gap < best$gap)) best <- ca
    if (isTRUE(best$pass)) break
  }
  if (is.null(best)) stop("cub_sim_panelA: no candidate could be fitted")
  d <- d0; d$Response <- as.integer(best$y)
  attr(d, "shelter_cat") <- cs; attr(d, "panel") <- "A"
  attr(d, "calibration") <- unlist(best[c("se", "seed", "R", "gap",
                                          "max_ups", "max_del", "pass")])
  d
}

## ---------------------------------------------------------------------
## 15. generating values, oracle theta, empirical Bayes calibration
## ---------------------------------------------------------------------

## The estimator targets the CENTRED intercepts.  Omitting this adjustment
## would report a spurious intercept bias.
cub_truth_vector <- function(D, tr) {
  c("xi:(Intercept)"       = tr$beta[[1]] + tr$rhof * D$ctr[["lag_f"]],
    "xi:GenderMale"        = tr$beta[[2]],
    "xi:Experience"        = tr$beta[[3]],
    "xi:Inflation"         = tr$beta[[4]],
    "xi:lag_f"             = tr$rhof,
    "ups:(Intercept)"      = tr$gam[[1]] + tr$rhou * D$ctr[["lag_u"]],
    "ups:EducationMedium"  = tr$gam[[2]],
    "ups:EducationHigh"    = tr$gam[[3]],
    "ups:Delivery_Speed"   = tr$gam[[4]],
    "ups:Customer_Service" = tr$gam[[5]],
    "ups:lag_u"            = tr$rhou,
    "del:(Intercept)"      = tr$alp[[1]],
    "del:Inflation"        = tr$alp[[2]],
    "sigma2_1" = tr$s1, "sigma2_2" = tr$s2,
    "sigma_12" = tr$rho12 * sqrt(tr$s1 * tr$s2), "rho_12" = tr$rho12)
}

cub_theta_oracle <- function(D, tr, spec = "RE-LAG") {
  co <- cub_core(D, spec, TRUE); tv <- cub_truth_vector(D, tr)
  th <- rep(NA_real_, co$p); names(th) <- co$nms
  k <- intersect(co$nms, names(tv)); th[k] <- tv[k]
  s1 <- sqrt(tr$s1); s2 <- sqrt(tr$s2); r <- tr$rho12
  th["log.L11"] <- log(s1); th["L21"] <- r * s2
  th["log.L22"] <- log(s2 * sqrt(1 - r^2))
  if (any(is.na(th)))
    stop("cub_theta_oracle: uncovered parameters: ",
         paste(co$nms[is.na(th)], collapse = ", "))
  list(th = as.numeric(th), co = co)
}

cub_eb_summary <- function(co, th, Qpred = 31L, re_true = NULL) {
  E <- co$ev(th, Qpred, grad = FALSE, pred = TRUE, quant = TRUE)
  P <- E$post; PR <- E$pred
  out <- c(sd_u = stats::sd(P$Up), sd_v = stats::sd(P$Vp),
           cor_uv = stats::cor(P$Up, P$Vp),
           cor_u_true = NA_real_, cor_v_true = NA_real_,
           cover_u = NA_real_, cover_v = NA_real_)
  if (!is.null(re_true)) {
    out[["cor_u_true"]] <- stats::cor(P$Up, re_true[, 1])
    out[["cor_v_true"]] <- stats::cor(P$Vp, re_true[, 2])
    out[["cover_u"]] <- mean(re_true[, 1] >= PR$Uq[, 1] &
                               re_true[, 1] <= PR$Uq[, 2])
    out[["cover_v"]] <- mean(re_true[, 2] >= PR$Vq[, 1] &
                               re_true[, 2] <= PR$Vq[, 2])
  }
  attr(out, "post") <- P; attr(out, "pred") <- PR
  out
}

## =====================================================================
## 16. the analysis pipeline
## =====================================================================

#' Control parameters for cubpanel_analyse
cubpanel_control <- function(
    Q = 15L, Qpred = 31L, Qcheck = c(9L, 15L, 21L, 27L), Qtol = 1e-2,
    gtol = 1e-3, nstart_pool = 3L, nstart_re = 1L,
    selfcheck = TRUE, refit_if_Q_inadequate = TRUE,
    rho_audit = TRUE, rho_grid = c(-0.3, 0, 0.3, 0.5, 0.7, 0.9),
    rho_budget = 160,
    profile = TRUE, profile_budget = 150,
    profile_grid = FALSE, profile_grid_budget = 180, profile_grid_n = 7L,
    chibar_B = 5000L, oos = TRUE, cold_start_audit = TRUE,
    nest_tol = NULL, progress = NULL) {
  as.list(environment())
}

#' Fit the specification ladder and every report block
#'
#' Returns a plain list of results.  All printing lives in
#' cubpanel-report.R, all plotting in cubpanel-figures.R, so that the
#' command-line and shiny paths produce identical numbers.
cubpanel_analyse <- function(D, control = cubpanel_control(),
                             truth = NULL, re_true = NULL) {
  ct <- control
  t0 <- Sys.time()
  say <- function(frac, msg) {
    if (is.function(ct$progress)) ct$progress(frac, msg)
    invisible(NULL)
  }
  an <- list(control = ct, D = D, started = t0)

  ## ---- block 1-2: data and dispersion -----------------------------
  say(0.02, "Data and dispersion diagnostic")
  an$decomp <- cub_decomp(D$y, D$gi)
  an$floor  <- cub_floor_diag(D)

  ## ---- block 3: self-checks ---------------------------------------
  if (isTRUE(ct$selfcheck)) {
    say(0.08, "Self-checks")
    an$checks <- cub_self_checks(D, Q = ct$Q, Qtol = ct$Qtol,
                                 Qcheck = ct$Qcheck)
  }
  nest_tol <- ct$nest_tol
  if (is.null(nest_tol))
    nest_tol <- max(1e-2,
                    if (!is.null(an$checks)) an$checks$implied_aggregate else 0)
  an$nest_tol <- nest_tol

  ## ---- blocks 4-5: ladder and chains ------------------------------
  say(0.18, "Specification ladder")
  lag_ok <- D$info$lag_ok
  spec_dyn <- if (lag_ok) "RE-LAG" else "RE"
  fits <- list()
  fits[["SB"]] <- cub_fit_rung(D, "SB", "Shifted binomial (boundary)", fits,
                               character(0), Q = ct$Q, nstart = 2L,
                               nest_tol = nest_tol, gtol = ct$gtol)
  fits[["0"]] <- cub_fit_rung(D, "0", "CUBpanel(0)", fits, "SB",
                              Q = ct$Q, nstart = ct$nstart_pool,
                              nest_tol = nest_tol, gtol = ct$gtol)
  if (lag_ok)
    fits[["LAG"]] <- cub_fit_rung(D, "LAG", "CUBpanel(LAG)", fits, "0",
                                  Q = ct$Q, nstart = ct$nstart_re,
                                  nest_tol = nest_tol, gtol = ct$gtol)
  fits[["RE"]] <- cub_fit_rung(D, "RE", "CUBpanel(RE)", fits, "0",
                               Q = ct$Q, nstart = ct$nstart_re,
                               nest_tol = nest_tol, gtol = ct$gtol)
  if (lag_ok) {
    fits[["RE-LAG"]] <- cub_fit_rung(
      D, "RE-LAG", "CUBpanel(RE-LAG)", fits, c("RE", "LAG"),
      Q = ct$Q, nstart = ct$nstart_re, nest_tol = nest_tol, gtol = ct$gtol)
    fits[["CRE-LAG"]] <- cub_fit_rung(
      D, "CRE-LAG", "CUBpanel(CRE-LAG)", fits, "RE-LAG",
      Q = ct$Q, nstart = ct$nstart_re, nest_tol = nest_tol, gtol = ct$gtol)
  }
  ## diagonal-Sigma variant of the richest dynamic rung
  fits[["DIAG"]] <- cub_fit(
    D, spec_dyn, rho_free = FALSE, Q = ct$Q, nstart = 1L,
    start = cub_pad_start(fits[[spec_dyn]], D, spec_dyn, FALSE),
    label = sprintf("CUBpanel(%s | Sigma diag)", spec_dyn), gtol = ct$gtol)
  fits[["DIAG"]]$nest_ok <- TRUE; fits[["DIAG"]]$nest_ok_lr <- TRUE

  an$fits <- fits
  an$spec_dyn <- spec_dyn
  say(0.42, "Fit assessment")
  an$stats <- lapply(fits, cub_fit_stats)
  an$diss  <- vapply(an$stats, `[[`, 0, "diss")

  ord <- c("SB", "0", if (lag_ok) "LAG", "RE", "DIAG",
           if (lag_ok) "RE-LAG", if (lag_ok) "CRE-LAG")
  ord <- ord[ord %in% names(fits)]
  an$ladder <- cub_compare(fits[ord], an$diss[ord])

  ch1 <- c("SB", "0", "RE", if (lag_ok) "RE-LAG", if (lag_ok) "CRE-LAG")
  ch1 <- ch1[ch1 %in% names(fits)]
  an$chain1 <- cub_chain(fits[ch1], an$diss[ch1], nest_tol)
  if (lag_ok) {
    ch2 <- c("SB", "0", "LAG", "RE-LAG", "CRE-LAG")
    ch2 <- ch2[ch2 %in% names(fits)]
    an$chain2 <- cub_chain(fits[ch2], an$diss[ch2], nest_tol)
  }
  an$selected <- list(bicN = an$ladder$best_bicN, aic = an$ladder$best_aic)

  ## ---- block 6: correlation audit ---------------------------------
  if (isTRUE(ct$rho_audit)) {
    say(0.48, "Correlation audit (chi^2_1 test and profile)")
    an$rho <- try(cub_rho_audit(fits[[spec_dyn]], fits[["DIAG"]],
                                grid = ct$rho_grid, budget = ct$rho_budget),
                  silent = TRUE)
    if (inherits(an$rho, "try-error")) an$rho <- NULL
  }

  ## ---- block 7: coefficients --------------------------------------
  say(0.58, "Coefficients, clustered standard errors, marginal effects")
  rep_rungs <- unique(c(spec_dyn, if (lag_ok) "CRE-LAG",
                        an$selected$bicN, an$selected$aic))
  rep_rungs <- rep_rungs[rep_rungs %in% names(fits)]
  an$coef <- list()
  for (nm in rep_rungs) {
    f <- fits[[nm]]
    V <- try(cub_vcov(f), silent = TRUE)
    if (inherits(V, "try-error")) next
    A <- try(cub_ame(f, V), silent = TRUE)
    if (inherits(A, "try-error")) A <- NULL
    tru <- if (!is.null(truth) && nm != "CRE-LAG") truth else NULL
    an$coef[[nm]] <- cub_coeftab(f, V, A, tru)
  }
  an$coef_rungs <- names(an$coef)
  an$primary_coef <- if (spec_dyn %in% names(an$coef)) spec_dyn
  else names(an$coef)[1]

  ## ---- block 8: profile interval ----------------------------------
  if (isTRUE(ct$profile)) {
    say(0.66, "Profile interval for sigma^2_1")
    an$profile_s1 <- try(cub_profile_s1(fits[[spec_dyn]],
                                        budget = ct$profile_budget),
                         silent = TRUE)
    if (inherits(an$profile_s1, "try-error")) an$profile_s1 <- NULL
    if (isTRUE(ct$profile_grid)) {
      an$profile_s1_grid <- cub_profile_s1_grid(
        fits[[spec_dyn]], ngrid = ct$profile_grid_n,
        budget = ct$profile_grid_budget)
    }
  }

  ## ---- block 9: test battery --------------------------------------
  say(0.72, "Test battery")
  an$tests <- cubpanel_tests(an)

  ## ---- blocks 10-11: fit of the selected rung, quadrature ---------
  say(0.80, "Quadrature certification on the selected rung")
  sel <- an$selected$bicN
  fsel <- fits[[sel]]
  Qseq <- sort(unique(as.integer(c(ct$Qcheck, ct$Q))))
  qq <- vapply(Qseq, function(q) fsel$co$ev(fsel$th, q, grad = FALSE)$ll, 0)
  quad <- data.frame(Q = Qseq, nodes = Qseq^2, loglik = qq,
                     deviation = qq - qq[length(qq)])
  dQ <- abs(qq[Qseq == ct$Q] - qq[length(qq)])
  an$quad <- list(table = quad, drift = dQ, tol = ct$Qtol,
                  adequate = dQ < ct$Qtol, refitted = FALSE,
                  ll_before = fsel$ll)
  if (dQ >= ct$Qtol && isTRUE(ct$refit_if_Q_inadequate)) {
    Qhi <- max(Qseq)
    fsel <- cub_fit(D, fsel$spec, fsel$rho_free, Q = Qhi, nstart = 1L,
                    start = fsel$th,
                    label = paste0(fits[[sel]]$label, sprintf(" (Q = %d)", Qhi)),
                    gtol = ct$gtol)
    an$quad$refitted <- TRUE; an$quad$Q_used <- Qhi
    an$quad$ll_after <- fsel$ll
    an$quad$ll_change <- fsel$ll - an$quad$ll_before
  }
  an$fit_selected <- fsel
  an$stats_selected <- cub_fit_stats(fsel, max(ct$Qpred, fsel$Q), quant = TRUE)

  ## ---- block 12: subject effects ----------------------------------
  say(0.86, "Subject effects")
  fdyn <- fits[[spec_dyn]]
  an$eb <- list(
    estimated = cub_eb_summary(fdyn$co, fdyn$th, ct$Qpred, re_true),
    sigma = { L <- fdyn$co$chol_of(fdyn$th)
      c(sigma_1 = L$L11, sigma_2 = sqrt(L$L21^2 + L$L22^2)) })
  if (!is.null(truth) && !is.null(re_true) && lag_ok) {
    tho <- try(cub_theta_oracle(D, attr(truth, "truth_list"), spec_dyn),
               silent = TRUE)
    if (!inherits(tho, "try-error")) {
      an$eb$oracle <- cub_eb_summary(tho$co, tho$th, ct$Qpred, re_true)
      an$eb$ll_oracle <- tho$co$ev(tho$th, ct$Q, grad = FALSE)$ll
      an$eb$ll_ml <- fdyn$ll
      an$eb$sigma_true <- c(
        sigma_1 = sqrt(attr(truth, "truth_list")$s1),
        sigma_2 = sqrt(attr(truth, "truth_list")$s2))
    }
  }

  ## ---- block 13: out-of-sample ------------------------------------
  if (isTRUE(ct$oos) && D$info$oos_ok && D$m >= 4L) {
    say(0.92, "Out-of-sample: held-out final wave")
    an$oos <- try(cubpanel_oos(D, fits[[spec_dyn]], ct), silent = TRUE)
    if (inherits(an$oos, "try-error")) {
      an$oos <- NULL
      an$oos_note <- "the out-of-sample block failed and is not reported"
    }
  } else {
    an$oos_note <- if (!D$info$oos_ok)
      "withheld: the held-out split requires min T_i >= 3" else
        "withheld"
  }

  ## ---- block 14: memory profile and collapse audit ----------------
  say(0.96, "Memory profile and collapse audit")
  an$memory <- cub_memory_profile(fdyn, max(ct$Qpred, fdyn$Q))
  an$verif  <- cubpanel_verify(an, ct)
  an$truth  <- truth
  an$secs   <- .el(t0)
  an$stamp  <- list(when = Sys.time(), R = as.character(getRversion()),
                    platform = R.version$platform)
  say(1, "Done")
  class(an) <- "cubpanel_analysis"
  an
}

## test battery ---------------------------------------------------------

cubpanel_tests <- function(an) {
  fits <- an$fits; ct <- an$control
  spec_dyn <- an$spec_dyn; D <- an$D
  rows <- list()
  add <- function(h, ref, stat, df, p, note = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      hypothesis = h, reference = ref, statistic = stat,
      df = df, p = p, note = note, stringsAsFactors = FALSE)
  }
  fd <- fits[[spec_dyn]]
  V <- if (spec_dyn %in% names(an$coef)) an$coef[[spec_dyn]]$V else NULL

  if (!is.null(V)) {
    drop_f <- c("(Intercept)", "lag_f")
    drop_u <- c("(Intercept)", "lag_u")
    wb <- cub_wald(fd, V, paste0("xi:", setdiff(colnames(fd$co$X), drop_f)))
    wg <- cub_wald(fd, V, paste0("ups:", setdiff(colnames(fd$co$W), drop_u)))
    wa <- cub_wald(fd, V, paste0("del:", setdiff(colnames(fd$co$Z),
                                                "(Intercept)")))
    if (!is.null(wb)) add("beta_{-0} = 0", "chi^2, Wald, clustered",
                          wb["W"], wb["df"], wb["p"])
    if (!is.null(wg)) add("gamma_{-0} = 0", "chi^2, Wald, clustered",
                          wg["W"], wg["df"], wg["p"])
    if (!is.null(wa)) add("omega_{-0} = 0", "chi^2, Wald, clustered",
                          wa["W"], wa["df"], wa["p"])
  }

  if (D$info$lag_ok && all(c("RE", "RE-LAG") %in% names(fits))) {
    lr <- max(2 * (fits[["RE-LAG"]]$ll - fits[["RE"]]$ll), 0)
    add("rho_f = rho_u = 0", "chi^2, LR", lr, 2,
        stats::pchisq(lr, 2, lower.tail = FALSE))
    jf <- which(fits[["RE-LAG"]]$co$nms == "xi:lag_f")
    ju <- which(fits[["RE-LAG"]]$co$nms == "ups:lag_u")
    if (length(jf) && length(ju)) {
      fe <- try(cub_fit(D, "RE-LAG", TRUE, Q = ct$Q, nstart = 1L,
                        tmap = cub_map_equal(fits[["RE-LAG"]], jf, ju),
                        label = "rho_f = rho_u"), silent = TRUE)
      if (!inherits(fe, "try-error")) {
        lre <- max(2 * (fits[["RE-LAG"]]$ll - fe$ll), 0)
        add("rho_f = rho_u", "chi^2, LR", lre, 1,
            stats::pchisq(lre, 1, lower.tail = FALSE),
            sprintf("rho_f = %+.3f, rho_u = %+.3f",
                    fits[["RE-LAG"]]$th["xi:lag_f"],
                    fits[["RE-LAG"]]$th["ups:lag_u"]))
      }
    }
  }

  lr_re <- max(2 * (fits[["RE"]]$ll - fits[["0"]]$ll), 0)
  if (!length(fits[["RE"]]$boundary_nat)) {
    VR <- try(cub_vcov(fits[["RE"]]), silent = TRUE)
    pcb <- if (inherits(VR, "try-error")) c(p = NA_real_, B = ct$chibar_B)
    else cub_chibar(cub_vc_natural(fits[["RE"]], VR)$V12, lr_re,
                    B = ct$chibar_B)
    add("Sigma = 0", sprintf("chi-bar^2, simulated, B = %d", ct$chibar_B),
        lr_re, NA, pcb[["p"]], "exact projection onto the psd cone")
  } else {
    add("Sigma = 0", "chi-bar^2 WITHHELD", lr_re, NA, NA_real_,
        sprintf("plug-in information at the frontier: %s",
                paste(fits[["RE"]]$boundary_nat, collapse = ", ")))
  }
  add("Sigma = 0", "chi^2_3, conservative companion", lr_re, 3,
      stats::pchisq(lr_re, 3, lower.tail = FALSE))

  if (!is.null(an$rho))
    add("rho_12 = 0", "chi^2_1, LR, interior point", an$rho$lr, 1, an$rho$p,
        if (an$rho$informative) "informative" else
          "uninformative corner: Sigma diagonal is implied")

  if ("CRE-LAG" %in% names(fits)) {
    fC <- fits[["CRE-LAG"]]
    lrm <- max(2 * (fC$ll - fits[["RE-LAG"]]$ll), 0)
    dfm <- fC$p - fits[["RE-LAG"]]$p
    add("CRE block = 0, jointly", "chi^2, LR", lrm, dfm,
        stats::pchisq(lrm, dfm, lower.tail = FALSE))
    VC <- if ("CRE-LAG" %in% names(an$coef)) an$coef[["CRE-LAG"]]$V else NULL
    if (!is.null(VC)) {
      wbar <- cub_wald(fC, VC, c(paste0("xi:", D$cre$xbar),
                                 paste0("ups:", D$cre$wbar)))
      wini <- cub_wald(fC, VC, c(paste0("xi:", D$cre$xinit),
                                 paste0("ups:", D$cre$winit)))
      if (!is.null(wbar))
        add("  (a) Mundlak means = 0", "chi^2, Wald, clustered",
            wbar["W"], wbar["df"], wbar["p"],
            "tests random-effects exogeneity")
      if (!is.null(wini))
        add("  (b) initial conditions = 0", "chi^2, Wald, clustered",
            wini["W"], wini["df"], wini["p"],
            "tests the conditioning on the first wave")
    }
  }

  lr_sb <- max(2 * (fits[["0"]]$ll - fits[["SB"]]$ll), 0)
  df_sb <- fits[["0"]]$p - fits[["SB"]]$p
  add("Shifted binomial vs CUBpanel(0)", "nominal only", lr_sb, df_sb,
      stats::pchisq(lr_sb, df_sb, lower.tail = FALSE),
      "both mixing intercepts diverge under H0; slopes unidentified")

  do.call(rbind, rows)
}

## out-of-sample block -------------------------------------------------

cubpanel_oos <- function(D, f_dyn, ct) {
  sp <- cub_oos_prep(D); Dtr <- sp$train; Dte <- sp$test
  stopifnot(identical(Dtr$subj, Dte$subj))
  ftr <- cub_fit(Dtr, f_dyn$spec, f_dyn$rho_free, Q = f_dyn$Q, nstart = 1L,
                 start = f_dyn$th, label = "training fit (CUBpanel)",
                 gtol = ct$gtol)
  Pc  <- cub_pred_cond(ftr, Dte, max(ct$Qpred, ftr$Q))
  fsb <- cub_fit(Dtr, "SB", nstart = 2L, Q = ct$Q,
                 label = "training fit (shifted binomial)", gtol = ct$gtol)
  Psb <- cub_core(Dte, "SB")$ev(fsb$th, 1L, grad = FALSE,
                                pred = TRUE)$pred$Pm
  Xtr <- cub_info_set(Dtr, f_dyn$spec); Xte <- cub_info_set(Dte, f_dyn$spec)
  stopifnot(identical(colnames(Xtr), colnames(Xte)))
  pof <- cub_po_fit(Dtr$y, Xtr, Dtr$gi, D$m)
  Pp  <- cub_po_pred(pof, Dtr$y, Xtr, Dtr$gi, Xte, Dte$gi)
  sc <- rbind(cub_score_pred(Pc, Dte$y),
              cub_score_pred(Pp, Dte$y),
              cub_score_pred(Psb, Dte$y))
  rownames(sc) <- c(sprintf("CUBpanel(%s)", f_dyn$spec),
                    "Random-effects proportional odds",
                    "Shifted binomial")
  direction <- c(logscore = -1, RPS = -1, hit = +1, diss = -1)
  winner <- vapply(names(direction), function(k)
    if (sign(sc[1, k] - sc[2, k]) == sign(direction[[k]]))
      "CUBpanel" else "prop-odds", "")
  nwin <- sum(winner == "CUBpanel")
  list(scores = sc, n_test = Dte$n, regressors = colnames(Xtr),
       ll = c(CUBpanel = ftr$ll, prop_odds = pof$ll, shifted_binomial = fsb$ll),
       p  = c(CUBpanel = ftr$p, prop_odds = length(pof$th),
              shifted_binomial = fsb$p),
       winner = winner, n_win = nwin,
       verdict = if (nwin >= 3L) "advantage CUBpanel"
       else if (nwin <= 1L) "advantage proportional odds"
       else "predictive parity",
       arbiter = if (sc[1, "RPS"] < sc[2, "RPS"]) "CUBpanel" else "prop-odds")
}

## collapse audit ------------------------------------------------------

cubpanel_verify <- function(an, ct) {
  fits <- an$fits; D <- an$D
  f0 <- fits[["0"]]
  p0 <- f0$co$ev(f0$th, 1L, grad = FALSE, post = TRUE)$post
  V1 <- f0$ll - an$floor[["ll_sb_glm"]]
  V2 <- max(p0$del); V3 <- mean((1 - p0$del) * p0$ups / D$m)
  gsel <- an$stats[[an$selected$bicN]]
  V4 <- (gsel$EV + gsel$VE) / gsel$var_obs
  V5 <- NA; cold <- NULL
  richest <- if ("CRE-LAG" %in% names(fits)) "CRE-LAG" else an$spec_dyn
  nested_by <- if (richest == "CRE-LAG") an$spec_dyn else "RE"
  if (isTRUE(ct$cold_start_audit) && richest %in% names(fits)) {
    cold <- try(cub_fit(D, fits[[richest]]$spec, fits[[richest]]$rho_free,
                        Q = ct$Q, nstart = 1L, start = NULL,
                        label = sprintf("%s (cold start)", richest),
                        gtol = ct$gtol), silent = TRUE)
    if (!inherits(cold, "try-error"))
      V5 <- cold$ll < fits[[nested_by]]$ll - 1e-6
  }
  tab <- data.frame(
    id = c("V1", "V2", "V3", "V4", "V5"),
    verification = c(
      "loglik CUBpanel(0) less shifted-binomial loglik from glm.fit",
      "largest fitted shelter weight over all rows",
      "mean uniform mass delivered per category",
      "model-implied variance over observed variance",
      sprintf("audit flagged a cold-start %s failure", richest)),
    value = c(sprintf("%.3e", V1), sprintf("%.3e", V2), sprintf("%.3e", V3),
              sprintf("%.4f", V4),
              if (isTRUE(V5)) "yes" else if (isFALSE(V5)) "no" else "n/a"),
    reference = c("0 if the collapse is exact", "0 if there is no shelter",
                  sprintf("1/m = %.4f if upsilon > 0", 1 / D$m),
                  "1 under correct dispersion (REPORTED, not tested)",
                  "detection required if present"),
    stringsAsFactors = FALSE)
  list(table = tab, V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5,
       sd_ratio = sqrt(V4),
       implied_total = gsel$EV + gsel$VE,
       mean_ups = mean(p0$ups), max_ups = max(p0$ups),
       mean_del = mean(p0$del), max_del = max(p0$del),
       post0 = p0, cold = if (inherits(cold, "try-error")) NULL else cold,
       collapsed = abs(V1) < 1e-2 && max(p0$ups) < 1e-2 && max(p0$del) < 1e-2)
}

## =====================================================================
## 17. replication driver (Monte Carlo)
## =====================================================================

CUB_MC_PARS <- c("xi:GenderMale", "xi:Experience", "xi:Inflation", "xi:lag_f",
                 "ups:EducationHigh", "ups:Delivery_Speed",
                 "ups:Customer_Service", "ups:lag_u", "del:Inflation",
                 "sigma2_1", "sigma2_2", "rho_12")

cub_mc_one <- function(r, N = 300L, Tn = 7L, specs = c("RE-LAG", "CRE-LAG"),
                       Q = 11L, truth = CUB_TRUTH_B, seed0 = 7e5,
                       seed_cov = 20240101L, fixed_design = TRUE,
                       tests = TRUE) {
  sc <- if (fixed_design) seed_cov else as.integer(seed0 + 5e5 + r)
  dat <- try(cub_sim_panelB(N = N, Tn = Tn, seed = as.integer(seed0 + r),
                            seed_cov = sc, truth = truth), silent = TRUE)
  if (inherits(dat, "try-error")) return(NULL)
  D <- try(cub_prepare(dat, shelter_cat = attr(dat, "shelter_cat")),
           silent = TRUE)
  if (inherits(D, "try-error")) return(NULL)
  tv <- cub_truth_vector(D, truth); fits <- list()
  for (sp in c("0", "RE", "RE-LAG", if ("CRE-LAG" %in% specs) "CRE-LAG")) {
    par <- switch(sp, "0" = character(0), "RE" = "0",
                  "RE-LAG" = "RE", "CRE-LAG" = "RE-LAG")
    f <- try(cub_fit_rung(D, sp, sp, fits, par, TRUE, Q = Q, nstart = 1L),
             silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    fits[[sp]] <- f
  }
  lr_sig <- max(2 * (fits[["RE"]]$ll - fits[["0"]]$ll), 0)
  lr_rho <- NA_real_
  if (tests) {
    fd <- try(cub_fit(D, "RE-LAG", FALSE, Q = Q, nstart = 1L,
                      start = cub_pad_start(fits[["RE-LAG"]], D,
                                            "RE-LAG", FALSE)), silent = TRUE)
    if (!inherits(fd, "try-error"))
      lr_rho <- max(2 * (fits[["RE-LAG"]]$ll - fd$ll), 0)
  }
  rows <- list()
  for (sp in specs) {
    f <- fits[[sp]]; if (is.null(f)) next
    V <- try(cub_vcov(f), silent = TRUE); if (inherits(V, "try-error")) next
    cc <- try(cub_coeftab(f, V), silent = TRUE)
    if (inherits(cc, "try-error")) next
    tb <- cc$tab; nb <- length(f$boundary_nat) > 0L
    for (pp in CUB_MC_PARS) {
      if (!(pp %in% rownames(tb))) next
      est <- tb[pp, "est"]; se <- tb[pp, "se"]
      ## in the CRE rung the variance components are CONDITIONAL on
      ## (y_i1, xbar) and are not comparable with the marginal truth
      cond <- sp == "CRE-LAG" &&
        pp %in% c("sigma2_1", "sigma2_2", "rho_12")
      true <- if (cond) NA_real_ else unname(tv[pp])
      rows[[length(rows) + 1L]] <- data.frame(
        rep = r, N = N, Tn = Tn, spec = sp, par = pp, est = est, se = se,
        true = true, comparable = !cond,
        cover = if (cond || !is.finite(se) || !is.finite(est)) NA_integer_
        else as.integer(abs(est - true) / se < 1.96),
        frontier = as.integer(nb), kkt = as.integer(f$ok),
        lr_sigma = lr_sig, lr_rho = lr_rho, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

cubpanel_mc <- function(R, N = 300L, Tn = 7L, specs = c("RE-LAG", "CRE-LAG"),
                        Q = 11L, ncores = 1L, budget = Inf, cache = NULL,
                        R_min = 100L, fixed_design = TRUE, tests = TRUE) {
  if (!is.null(cache) && file.exists(cache)) {
    old <- readRDS(cache); done <- max(old$rep, 0L)
    message(sprintf("resuming from cache: %d replications already complete",
                    length(unique(old$rep))))
  } else { old <- NULL; done <- 0L }
  par_ok <- ncores > 1L && .Platform$OS.type == "unix"
  t0 <- Sys.time()
  out <- if (is.null(old)) list() else list(old)
  nok <- if (is.null(old)) 0L else length(unique(old$rep))
  while (nok < R && .el(t0) < budget) {
    idx <- (done + 1L):min(done + ncores, R); done <- max(idx)
    res <- if (par_ok)
      parallel::mclapply(idx, cub_mc_one, N = N, Tn = Tn, specs = specs,
                         Q = Q, fixed_design = fixed_design, tests = tests,
                         mc.cores = ncores)
    else lapply(idx, cub_mc_one, N = N, Tn = Tn, specs = specs, Q = Q,
                fixed_design = fixed_design, tests = tests)
    res <- res[!vapply(res, function(z)
      is.null(z) || inherits(z, "try-error"), TRUE)]
    if (length(res)) {
      out[[length(out) + 1L]] <- do.call(rbind, res)
      nok <- nok + length(res)
      if (!is.null(cache)) saveRDS(do.call(rbind, out), cache)
    }
    if (done >= R) break
  }
  if (!nok) return(NULL)
  df <- do.call(rbind, out)
  attr(df, "secs") <- .el(t0); attr(df, "R_target") <- R
  attr(df, "R_min") <- R_min; attr(df, "ncores") <- ncores
  attr(df, "fixed_design") <- fixed_design
  df
}

cub_mc_summary <- function(df) {
  if (is.null(df)) return(NULL)
  key <- interaction(df$N, df$Tn, df$spec, df$par, drop = TRUE, sep = "|")
  agg <- do.call(rbind, lapply(split(df, key), function(z) {
    cmp <- z$comparable & is.finite(z$est)
    data.frame(N = z$N[1], Tn = z$Tn[1], spec = z$spec[1], par = z$par[1],
               R_ok = sum(is.finite(z$est)),
               true = if (any(cmp)) z$true[cmp][1] else NA_real_,
               mean_est = mean(z$est[is.finite(z$est)]),
               bias = if (any(cmp)) mean(z$est[cmp] - z$true[cmp]) else NA_real_,
               rmse = if (any(cmp))
                 sqrt(mean((z$est[cmp] - z$true[cmp])^2)) else NA_real_,
               sd_est = stats::sd(z$est[is.finite(z$est)]),
               mean_se = mean(z$se[is.finite(z$se)]),
               cover95 = if (any(!is.na(z$cover)))
                 mean(z$cover, na.rm = TRUE) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  agg <- agg[order(agg$N, agg$Tn, agg$spec, match(agg$par, CUB_MC_PARS)), ]
  rownames(agg) <- NULL
  u <- unique(df[, c("rep", "N", "Tn", "spec", "frontier", "kkt",
                     "lr_sigma", "lr_rho")])
  dg <- do.call(rbind, lapply(
    split(u, interaction(u$N, u$Tn, u$spec, drop = TRUE, sep = "|")),
    function(z) data.frame(
      N = z$N[1], Tn = z$Tn[1], spec = z$spec[1], R = nrow(z),
      KKT_rate = mean(z$kkt), frontier_rate = mean(z$frontier),
      reject_Sigma0_5pct =
        mean(stats::pchisq(z$lr_sigma, 3, lower.tail = FALSE) < 0.05),
      reject_rho12_5pct =
        mean(stats::pchisq(z$lr_rho, 1, lower.tail = FALSE) < 0.05,
             na.rm = TRUE),
      stringsAsFactors = FALSE)))
  rownames(dg) <- NULL
  list(param = agg, diag = dg)
}