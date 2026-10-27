# bandwidth_joint.R -- the bandwidth rules that target the RD-DID estimate itself:
# .bw_rot() is the rule-of-thumb bandwidth shared by every period;
# .bw_stage() and .bw_pilots() choose the pilot bandwidths d and b for the aggregate, as
# rdrobust's selector chooses them for one period;
# .bw_constants() estimates the variance constants and the bias design moments at h_rot and the
# curvature at the common b;
# .bw_joint() minimises the aggregate AMSE over one common h (bwselect = "joint", the default).
# Labels such as eq:common_h_opt, lem:agg-var or "B.4 P3" in the internal notes below point into the
# paper's Appendix B; the code-to-paper map is dev/appB_map.md (developer material, not shipped).

# Original file header:
# Bandwidth selection for the aggregate RD-DID estimator.
#
# Two modes, matching Section 5.4 (sec:est-bw) and Appendix B.4 (app:est-bw) of
# Leventer and Nevo; the code <-> equation map is dev/appB_map.md (section 2.5).
#   * "cct"   - per-period MSE-optimal bandwidth (Imbens-Kalyanaraman / CCT) applied to
#               each D-hat_t on its own (B.4 P1). Delegated to rdrobust::rdbwselect.
#               Leaves the cross-period bias cancellation unexploited.
#   * "joint" - one common h* minimising the asymptotic MSE of the AGGREGATE estimator
#               (B.4 P3, eq:common_h_opt):
#                 h* = ( ((p+1)!)^2/(2(p+1)) * V^S / B^2 )^{1/(2p+3)} n^{-1/(2p+3)},
#               with B = sum_tau w~_tau B_tau (signed: comparison biases can offset the
#               RD-period bias) and V^S the scheme-specific aggregate variance constant
#               (lem:agg-var at a common h).
# Both use the conventional variance (CCT's mserd convention); the regularization
# of B.4 ("Regularization") is applied in "joint". (The period-specific "iter" rule was
# removed after tag v0.4.0.9000-iter; it was dropped from the paper, T1633.)


#' Mass points of the running variable, as rdrobust's `masspoints = "adjust"`
#'
#' Per period and side: the number of rows `N` and of distinct values `M`. `mass` is TRUE when
#' some period has a side with `1 - M / N >= 0.2` (rdrobust then sets `bwcheck = 10`); `bw_min` is
#' the largest, over periods and sides, of the 10th smallest distinct distance to the cutoff
#' (fewer if the side has fewer distinct values), times 1 + sqrt(eps) (rdrobust 4.1); `m_bar` is the average number of
#' distinct values per period (rdrobust's `M`, which equals `N` without ties); `radius` the
#' largest distance to the cutoff (rdrobust's `bw_max`).
#' @keywords internal
#' @noRd
.bw_mass <- function(plist, c = 0) {
  cutoff <- c
  per <- lapply(plist, function(d) {
    x <- d$x[is.finite(d$x)]
    lapply(list(x[x < cutoff], x[x >= cutoff]), function(xs) {
      u <- sort(unique(abs(xs - cutoff)))
      list(N = length(xs), M = length(u), bw_min = u[min(10L, length(u))] * (1 + sqrt(.Machine$double.eps)))
    })
  })
  sides <- unlist(per, recursive = FALSE)
  list(mass   = any(vapply(sides, function(z) 1 - z$M / z$N >= 0.2, logical(1))),
       bw_min = max(vapply(sides, function(z) z$bw_min, numeric(1))),
       m_bar  = mean(vapply(per, function(z) z[[1]]$M + z[[2]]$M, numeric(1))),
       radius = max(vapply(plist, function(d) max(abs(d$x - cutoff), na.rm = TRUE), numeric(1))))
}

#' Rule-of-thumb preliminary bandwidth shared by every period (App. B.4, "Bandwidth estimation")
#'
#' rdrobust's rule of thumb for the preliminary bandwidth of `rdbwselect()`,
#' `C_K * min(sd, IQR / 1.349) * M^{-1/5}` (`C_K` = 2.576 triangular, 1.843 uniform, 2.34
#' Epanechnikov; the IQR from type-2 quantiles), computed once from the running variables of all
#' periods pooled, with `M` the average number of distinct running-variable values per period
#' (the number of rows when there are no ties; rdrobust's `masspoints = "adjust"`). As rdrobust,
#' it is capped at the largest distance to the cutoff and, when the running variable has mass
#' points ([.bw_mass()]), raised to at least the 10th distinct distance on each side. Every
#' period's variance constant is estimated at it, so the cross-period covariances are estimated
#' at a bandwidth ratio of one.
#' @keywords internal
#' @noRd
.bw_rot <- function(plist, kernel = "triangular", c = 0) {
  x <- unlist(lapply(plist, function(d) d$x), use.names = FALSE)
  x <- x[is.finite(x)]
  iqr <- unname(stats::quantile(x, 0.75, type = 2) - stats::quantile(x, 0.25, type = 2))
  s <- min(stats::sd(x), iqr / 1.349)
  if (!is.finite(s) || s <= 0) s <- stats::sd(x)
  ck <- switch(kernel, triangular = 2.576, uniform = 1.843, epanechnikov = 2.34)
  mp <- .bw_mass(plist, c = c)
  h  <- min(ck * s * mp$m_bar^(-1 / 5), mp$radius)
  if (mp$mass) h <- max(h, mp$bw_min)
  h
}

#' One side of one period for a stage of the selector
#'
#' The order-`o` local polynomial fit at bandwidth `hh` on one side: the coefficient on
#' `(x - c)^nu`, its per-unit influence times the HC1 residual of this fit (`g`), and rdrobust's
#' `BConst`, `hh^nu e_nu' (X'A X)^{-1} X'A ((x - c)/hh)^{o+1}` (`bconst`, the paper's design moment
#' `B_{t,(+-),nu,o}(hh)` divided by `nu!`).
#' @keywords internal
#' @noRd
.bw_side_poly <- function(y, x, id, cutoff, o, nu, hh, kernel) {
  w <- .rd_kweight((x - cutoff) / hh, kernel)
  k <- w > 0
  if (length(unique(x[k])) <= o + 1L)
    .abort("too few distinct running-variable values on one side of the cutoff for the ",
           "order-", o, " fit at bandwidth ", signif(hh, 4), " in the bandwidth selector; ",
           "use bwselect = \"cct\" or pass h.")
  y <- y[k]; x <- x[k]; id <- id[k]; w <- w[k]
  X  <- outer(x - cutoff, 0:o, `^`)
  Gi <- .qrXXinv(sqrt(w) * X)
  A  <- Gi %*% t(X * w)                                 # coefficients = A y
  beta <- drop(A %*% y)
  ns  <- length(y)
  res <- sqrt(ns / (ns - (o + 1L))) * (y - drop(X %*% beta))
  bconst <- hh^nu * drop(Gi %*% crossprod(X * w, ((x - cutoff) / hh)^(o + 1L)))[nu + 1L]
  list(coef = beta[nu + 1L], g = A[nu + 1L, ] * res, id = id, bconst = bconst)
}

#' Variance of sum_t coef_t (term_{t,+} - term_{t,-}) from per-side influence vectors
#'
#' `G[[t]][[side]]` holds `g` and `id`. The own-period blocks always enter; the cross-period
#' blocks by `cross`: `"none"`, `"same"` (same side, matched by id) or `"same_opp"` (same side
#' minus opposite side), the rules of [.bw_constants()].
#' @keywords internal
#' @noRd
.bw_agg_var <- function(G, coef, cross) {
  keys <- names(coef)
  ms <- function(a, sa, b, sb) .match_sum(G[[a]][[sa]]$id, G[[a]][[sa]]$g, G[[b]][[sb]]$id, G[[b]][[sb]]$g)
  v <- 0
  for (a in keys) for (b in keys) {
    cab <- if (a == b) {
      sum(G[[a]][["+"]]$g^2) + sum(G[[a]][["-"]]$g^2)
    } else if (cross == "none") {
      0
    } else {
      same <- ms(a, "+", b, "+") + ms(a, "-", b, "-")
      opp  <- if (cross == "same_opp") ms(a, "+", b, "-") + ms(a, "-", b, "+") else 0
      same - opp
    }
    v <- v + coef[[a]] * coef[[b]] * cab
  }
  v
}

#' One stage of the selector for the aggregate (App. B.4, "Bandwidth estimation")
#'
#' rdrobust's `rdrobust_bw()` stage, in its own units, for the aggregate
#' `sum_t coef_t (m^{(nu)}_{t,+} - m^{(nu)}_{t,-})` instead of one period's jump: every period at
#' the same `h_V` (variance, design moments) and `h_B` (the order-`oB` fit for the `(o+1)`-th
#' derivative),
#'   V  = (2 nu + 1) h_V^{2 nu + 1} Var(sum_t coef_t (beta_{nu,t,+} - beta_{nu,t,-})),
#'   B  = sqrt(2 (o + 1 - nu)) sum_t coef_t (BConst_{t,+} beta_{o+1,t,+} - BConst_{t,-} beta_{o+1,t,-}),
#'   R  = reg_const 2 (o + 1 - nu) Var(the sum in B, BConst held fixed),
#'   bw = (V / (B^2 + R))^{1 / (2 o + 3)},
#' the paper's `h^opt(nu, o, V-bar-hat, B-bar-hat)` with `B^2` replaced by `B^2 + R`. `h_B = "range"`
#' takes, in each period and on each side, the largest distance of the running variable to the
#' cutoff, times 1 + sqrt(eps) so that the farthest unit keeps a positive kernel weight (rdrobust
#' 4.1; rdrobust 2.2 used the distance itself, a difference of order 1e-8 to 1e-4 in d). The cross-period covariances follow the scheme as in [.bw_constants()]: in `V` the
#' same-side ones under `"pc"` only; in `R` none under `"cs"`, same side under `"pc"`, same minus
#' opposite side under `"pv"`.
#' @keywords internal
#' @noRd
.bw_stage <- function(plist, coef, scheme, o, nu, oB, h_V, h_B, reg_const = 0, c = 0,
                      kernel = "triangular") {
  cutoff <- c
  keys <- names(coef)
  GV <- list(); GB <- list(); Bsum <- 0
  for (k in keys) {
    d <- plist[[k]]; GV[[k]] <- list(); GB[[k]] <- list()
    for (sd in c("+", "-")) {
      s  <- if (sd == "+") d$x >= cutoff else d$x < cutoff
      hb <- if (identical(h_B, "range")) max(abs(d$x[s] - cutoff)) * (1 + sqrt(.Machine$double.eps)) else h_B
      fv <- .bw_side_poly(d$y[s], d$x[s], d$id[s], cutoff, o, nu, h_V, kernel)
      fb <- .bw_side_poly(d$y[s], d$x[s], d$id[s], cutoff, oB, o + 1L, hb, kernel)
      sg <- if (sd == "+") 1 else -1
      Bsum <- Bsum + coef[[k]] * sg * fv$bconst * fb$coef
      GV[[k]][[sd]] <- list(g = fv$g, id = fv$id)
      GB[[k]][[sd]] <- list(g = fv$bconst * fb$g, id = fb$id)
    }
  }
  V <- (2 * nu + 1) * h_V^(2 * nu + 1) *
    .bw_agg_var(GV, coef, if (scheme == "pc") "same" else "none")
  B <- sqrt(2 * (o + 1 - nu)) * Bsum
  R <- if (reg_const > 0)
    reg_const * 2 * (o + 1 - nu) * .bw_agg_var(GB, coef, c(cs = "none", pc = "same", pv = "same_opp")[[scheme]])
  else 0
  if (!is.finite(V) || V <= 0)
    .abort("the bandwidth selector's variance scale is not positive; use bwselect = \"cct\" or pass h.")
  list(V = V, B = B, R = R, bw = (V / (B^2 + R))^(1 / (2 * o + 3)))
}

#' Pilot bandwidths d and b for the aggregate (App. B.4, "Bandwidth estimation")
#'
#' rdrobust's first two stages (`bwselect = "mserd"`) for the aggregate, every period at the
#' shared `h_rot`:
#'   d = h^opt(q+1, q+1, V-bar-hat(q+1, q+1, h_rot), B-bar-hat(q+1, q+1, q+2, h_rot, range)),
#'       not regularized;
#'   b = h^opt(p+1, q,   V-bar-hat(p+1, q, h_rot),   B-bar-hat(p+1, q, q+1, h_rot, d)),
#'       regularized when `regularize`.
#' Each is capped at the largest distance of the running variable to the cutoff over all
#' periods, as rdrobust's `bwrestrict = TRUE`, and, when the running variable has mass points,
#' `d` is raised to at least the 10th distinct distance on each side ([.bw_mass()]), as
#' rdrobust's `bwcheck`. With one period and `coef = 1` these are
#' rdbwselect's `d` (internal) and `b`; under a balanced PC panel, its `d` and `b` on the composite
#' outcome `sum_t coef_t Y_t` (at the same `h_rot`).
#' @return list with `d`, `b` and the two stages (`stage_d`, `stage_b`: `V`, `B`, `R`, `bw`).
#' @keywords internal
#' @noRd
.bw_pilots <- function(plist, coef, scheme = "cs", h_rot = NULL, c = 0, p = 1L, q = 2L,
                       kernel = "triangular", regularize = TRUE, reg_const = 3) {
  cutoff <- c
  if (is.null(h_rot)) h_rot <- .bw_rot(plist[names(coef)], kernel = kernel, c = cutoff)
  mp <- .bw_mass(plist[names(coef)], c = cutoff)
  radius <- mp$radius
  st_d <- .bw_stage(plist, coef, scheme, o = q + 1L, nu = q + 1L, oB = q + 2L,
                    h_V = h_rot, h_B = "range", reg_const = 0, c = cutoff, kernel = kernel)
  d <- min(st_d$bw, radius)
  if (mp$mass) d <- max(d, mp$bw_min)
  st_b <- .bw_stage(plist, coef, scheme, o = q, nu = p + 1L, oB = q + 1L,
                    h_V = h_rot, h_B = d, reg_const = if (regularize) reg_const else 0,
                    c = cutoff, kernel = kernel)
  b <- min(st_b$bw, radius)
  list(d = d, b = b, stage_d = st_d, stage_b = st_b)
}

#' Per-period asymptotic constants for the aggregate AMSE (App. B.4)
#'
#' Used by the joint selector; the paper's "Bandwidth estimation" paragraph of App. B.4.
#' Two sets of fits, as in CCT's own selector:
#'   * curvature: every period at the common pilot bandwidth `b`; gives the order-q
#'     coefficient on `(x - c)^{p+1}` at `b` on each side (`curv`);
#'   * design side: every period at ONE rule-of-thumb bandwidth `h_rot` ([.bw_rot()]); gives
#'     the bias design moments `bmom` (rdrobust's `BConst`), so that
#'     `B-hat_t = (p+1)! [bmom_+(h_rot) curv_+(b) - bmom_-(h_rot) curv_-(b)]`, the
#'     covariance matrix of the `B-hat_t` for the regularization, with the cross-period
#'     covariances the scheme allows (none under `"cs"`, same side under `"pc"`, same minus
#'     opposite side under `"pv"`, as for the jumps), and
#'     the variance constant `V-hat_t = n_t h_rot V-hat(D_t(h_rot))` (`v_const`) and, under
#'     `"pc"`, the same-side cross-period term `h_rot Cov-hat(D_t(h_rot), D_s(h_rot))`, estimated
#'     directly at a bandwidth ratio of one (the ratio the common bandwidth needs). Under
#'     `"cs"` the covariance is zero and under `"pv"` it is O(1/n) = o(1/(nh)) (`lem:cov-pv`),
#'     so the limit has no cross term and none is used.
#' Nothing here depends on which period is the RD period, so the selector built on it is
#' invariant to that label; this matters when the estimator aggregates several RD periods.
#'
#' @param plist named per-period data list; `coef` named period coefficients.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`.
#' @param b the common pilot bandwidth of the curvature fits (from [.bw_pilots()]).
#' @param h_rot optional bandwidth for the variance side; defaults to [.bw_rot()].
#' @return list with `keys`, `coef` (in `keys` order), `bias_const`, `var_const`, `n_obs`,
#'   `var_bias` (K x K covariance of the `B-hat_t`), `var_cross` (K x K same-side cross term,
#'   zero off PC), `h_rot`, `b`, `pilot_fits`, `var_fits`.
#' @keywords internal
#' @noRd
.bw_constants <- function(plist, coef, scheme = "cs", b, h_rot = NULL,
                          c = 0, p = 1L, q = 2L, kernel = "triangular") {
  cutoff <- c
  keys <- names(coef)
  K    <- length(keys)
  fp1  <- factorial(p + 1)
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= 0)
    .abort("`b` must be one positive pilot bandwidth, common to all periods.")
  if (is.null(h_rot)) h_rot <- .bw_rot(plist[keys], kernel = kernel, c = cutoff)

  fit_at <- function(k, h, b)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = h, b = b, id = plist[[k]]$id,
              c = cutoff, p = p, q = q, kernel = kernel)

  # The estimated bias of D_t on a side is h^{p+1} * bmom * curv: a design moment at h and the
  # curvature (the order-q coefficient on (x - c)^{p+1}) from the fit at b. As in rdrobust's
  # selector, the curvature has its own bandwidth b, but the design moment, like the variance,
  # uses only the running variable, so both are taken at the shared h_rot.
  # curvature fits: every period at the common b (main bandwidth b too, so that the active set,
  # and with it the HC1 factor of the curvature's influence, is the b window alone)
  pilot_fits <- stats::setNames(lapply(keys, function(k) fit_at(k, b, b)), keys)
  # variance-side fits: every period at the common rule-of-thumb bandwidth (the variance
  # constants and the design moments)
  var_fits <- stats::setNames(lapply(keys, function(k) fit_at(k, h_rot, h_rot)), keys)
  period_coef <- unname(vapply(keys, function(k) coef[[k]], numeric(1)))

  # B-hat_t = (p+1)! [bmom_+(h_rot) curv_+(b) - bmom_-(h_rot) curv_-(b)]. psi is the
  # per-unit influence of one side's term, (p+1)! bmom(h_rot) g_curv(b), unsigned: the
  # minus sign of the side below is applied through "same - opp" in the covariances.
  psi <- function(k, side)
    fp1 * var_fits[[k]]$sides[[side]]$bmom * pilot_fits[[k]]$sides[[side]]$g_curv
  bias_const <- vapply(keys, function(k)
    fp1 * (var_fits[[k]]$sides[["+"]]$bmom * pilot_fits[[k]]$sides[["+"]]$curv -
           var_fits[[k]]$sides[["-"]]$bmom * pilot_fits[[k]]$sides[["-"]]$curv), numeric(1))

  # Cov(B-hat_t, B-hat_s) from the psi vectors, matched by id across periods like the jumps'
  # g: this is what the regularization term adds to the squared bias of the aggregate
  ids <- function(k, side) pilot_fits[[k]]$sides[[side]]$id
  ms  <- function(a, sa, b, sb) .match_sum(ids(a, sa), psi(a, sa), ids(b, sb), psi(b, sb))
  var_bias <- matrix(0, K, K, dimnames = list(keys, keys))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    a <- keys[i]; bb <- keys[j]
    var_bias[i, j] <- if (i == j) {
      sum(psi(a, "+")^2) + sum(psi(a, "-")^2)
    } else if (scheme == "cs") {
      0
    } else {
      same <- ms(a, "+", bb, "+") + ms(a, "-", bb, "-")
      opp  <- if (scheme == "pv") ms(a, "+", bb, "-") + ms(a, "-", bb, "+") else 0
      same - opp
    }
  }

  var_const <- vapply(keys, function(k) var_fits[[k]]$v_const, numeric(1))
  n_obs     <- vapply(keys, function(k) var_fits[[k]]$n,       numeric(1))
  var_cross <- matrix(0, K, K, dimnames = list(keys, keys))
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K)) for (j in seq_len(K)) if (i != j)
      var_cross[i, j] <- h_rot * .cross_cov(var_fits[[keys[i]]], var_fits[[keys[j]]], bc = FALSE)$pc
    if (all(var_cross == 0))
      .warn("scheme = \"pc\" but no unit is active in two periods' windows: ",
            "the cross-period term is identically zero (same as scheme = \"cs\").")
  }
  list(keys = keys, coef = period_coef, bias_const = bias_const, var_const = var_const,
       n_obs = n_obs, var_bias = var_bias, var_cross = var_cross, h_rot = h_rot, b = b,
       pilot_fits = pilot_fits, var_fits = var_fits)
}

#' Joint AMSE-optimal common bandwidth for the aggregate estimator
#'
#' Feasible plug-in for eq:common_h_opt (App. B.4, "Bandwidth estimation"): the last line of
#' the paper's display, after `h_rot` ([.bw_rot()]) and the pilots `d` and `b` for the aggregate
#' ([.bw_pilots()]). With the constants of [.bw_constants()] (curvature at the common `b`;
#' variance and bias design moments at the shared rule-of-thumb `h_rot`), the aggregate objective
#' at a common `h` is
#'   AMSE^S(h) = (h^{p+1}/(p+1)!)^2 (B^2 + reg) + Veff / h,
#'   B    = sum_tau w~_tau B-hat_tau                       (signed; can cancel),
#'   Veff = sum_tau w~_tau^2 V-hat_tau / n_tau
#'          + 1{S = PC} sum_{t != s} w~_t w~_s h_rot Cov-hat(D_t(h_rot), D_s(h_rot)),
#'   reg  = lambda Var-hat(sum_tau w~_tau B-hat_tau)   (with the cross-period covariances),
#' (`Veff` is `V^S(t_RD) / n` of `lem:agg-var` at a common h, so the `n^{-1/(2p+3)}`
#' factor is absorbed), whose minimizer is
#'   h* = ( ((p+1)!)^2/(2(p+1)) * Veff / (B^2 + reg) )^{1/(2p+3)}.
#' Every period's bias correction uses the same pilot bandwidth `b` (B.4 after
#' eq:common_h_opt; consistent with Assumption R(g), as h*/b -> 0).
#'
#' @param plist named per-period data list; `coef` named period coefficients
#'   (RD period = +1, comparisons = -w_t); `t_rd` the RD-period label (kept for
#'   call compatibility; the selector does not use it).
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`; sets which variance scales the
#'   bandwidth. Under `"pc"` the same-side cross-period covariances enter (they are
#'   O(1/(nh)), `lem:cov-pc`); under `"pv"` they are O(1/n) = o(1/(nh)) (`lem:cov-pv`)
#'   and drop, so the covariance-free CS form is used (`lem:agg-var`).
#' @param b optional common pilot bandwidth; defaults to the joint `b` of [.bw_pilots()].
#' @param regularize add the App. B.4 regularization term to the squared bias so a
#'   small, noisy estimated bias constant cannot inflate `h*` (mirrors the
#'   `regularize = TRUE` default of [rdrobust::rdbwselect]); also regularizes `b`. Default
#'   `TRUE`.
#' @param reg_const multiple of the bias-estimate variance used as the
#'   regularization term (the paper's `lambda`; default 3, the CCT convention).
#' @param constants optional precomputed [.bw_constants()] output (reuses the
#'   constants without refitting; its `b` is then the pilot).
#' @return list with the common `h`, the common pilot `b`, the pilot `d` of `b` (`NA` when
#'   `b` or `constants` is given), the bias constant `B`, the variance scale `Veff`, the
#'   regularization term `reg`, the rule-of-thumb `h_rot`, and the `constants`.
#' @keywords internal
#' @noRd
.bw_joint <- function(plist, coef, t_rd = NULL, scheme = "cs", b = NULL,
                      c = 0, p = 1L, q = 2L, kernel = "triangular",
                      regularize = TRUE, reg_const = 3, constants = NULL) {
  cutoff <- c
  d <- NA_real_
  if (is.null(constants)) {
    h_rot <- .bw_rot(plist[names(coef)], kernel = kernel, c = cutoff)
    if (is.null(b)) {
      pil <- .bw_pilots(plist, coef, scheme = scheme, h_rot = h_rot, c = cutoff, p = p, q = q,
                        kernel = kernel, regularize = regularize, reg_const = reg_const)
      d <- pil$d
      b <- pil$b
    }
    constants <- .bw_constants(plist, coef, scheme = scheme, b = b, h_rot = h_rot,
                               c = cutoff, p = p, q = q, kernel = kernel)
  }
  const <- constants
  fp1   <- factorial(p + 1)
  K     <- length(const$keys)

  # aggregate bias constant: sum_t coef_t B-hat_t
  B_agg <- sum(const$coef * const$bias_const)
  # variance scale: own-period terms sum_t coef_t^2 V-hat_t / n_t, plus, under "pc", the
  # same-side cross-period terms (estimated at h_rot in every period, so at rho = 1)
  cf   <- const$coef
  Veff <- sum(cf^2 * const$var_const / const$n_obs) + drop(crossprod(cf, const$var_cross %*% cf))
  # regularization: reg_const x Var(sum_t coef_t B-hat_t), with the cross-period covariances of
  # the B-hat_t. At a common h it carries the same h^{2(p+1)} factor as B_agg^2, so it simply
  # joins B_agg^2 in the denominator of h*.
  reg   <- if (regularize) reg_const * drop(crossprod(cf, const$var_bias %*% cf)) else 0
  denom <- B_agg^2 + reg

  if (!is.finite(denom) || denom <= 0)
    .abort("aggregate bias constant B is ~0 and regularization is off: the ",
           "AMSE-optimal bandwidth diverges. Use bwselect = \"cct\" or pass h.")
  if (!is.finite(Veff) || Veff <= 0)
    .abort("aggregate variance scale is not positive; check the sampling scheme and ",
           "the per-period fits.")
  # the closed-form minimiser of the AMSE in h: for p = 1 the leading constant is 1 and the
  # exponent 1/5
  h_star <- (fp1^2 / (2 * (p + 1)) * Veff / denom)^(1 / (2 * p + 3))
  # capped at the largest distance to the cutoff, as rdrobust's bwrestrict
  radius <- max(vapply(plist, function(d) max(abs(d$x - cutoff), na.rm = TRUE), numeric(1)))
  if (is.finite(radius) && h_star > radius) {
    .warn("joint AMSE-optimal bandwidth h* = ", signif(h_star, 4),
          " exceeds the running-variable radius ", signif(radius, 4),
          " (the period biases nearly cancel) and is set to the radius; consider ",
          "bwselect = \"cct\" or a fixed h.")
    h_star <- radius
  }
  list(h = h_star, b = const$b, d = d, B = B_agg, Veff = Veff, reg = reg, h_rot = const$h_rot,
       constants = const)
}
