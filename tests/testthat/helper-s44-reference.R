# ---------------------------------------------------------------------------
# helper-s44-reference.R
#
# Independent reference implementation of the four validation tests of
# Section 4.4 of Leventer and Nevo:
#
#   ass:type-cont   (A8)  pi_{t,(+)}(v) = pi_{t,(-)}(v)      -> ref_typecont()
#   ass:comp-stable (A9)  pi_{tRD,(+)}(v) = pi_{t0,(+)}(v)   -> ref_compstable()
#   ass:homog       (A10) alpha_{t0,0}(0) = alpha_{t0,0}(1)  -> ref_homog()
#   ass:trend-cell  (A7)  alpha_{tRD,0}(v) = alpha_{t0,0}(v) -> ref_trendcell()
#
# Written *only* from the Section 4.4 prose and from the Appendix B machinery
# already in helper-appB-reference.R (ref_period_fit / ref_sigma /
# ref_cov_beta0 / ref_hc1 / ref_V_D).  R/test_typecont.R, R/test_compstable.R,
# R/test_homog.R, R/test_trendcell.R and R/test_helpers.R -- the code under
# test -- were NOT read.  Every local-linear jump below is a ref_period_fit(),
# and every variance is an Appendix B.2 sandwich.
#
# ---------------------------------------------------------------------------
# Assumptions made where the Section 4.4 prose leaves a choice (all reported).
# ---------------------------------------------------------------------------
#
# S1. Type.  "A unit's type in period t is its side of the cutoff in the other
#     period s": V_{i,s} = 1{R_{i,s} >= c}.  Only two-period panels are used
#     for ass:type-cont, so "the other period" is unambiguous.  For ass:homog
#     and ass:trend-cell the type/cell is the tRD side V_{i,tRD}, which the
#     package calls `type_by = "rd_side"`.  Types are labelled "+" (V = 1) and
#     "-" (V = 0).
#
# S2. The tested statistic.  The paper's test regresses the *indicator*
#     1{V_{i,s} = 1} on R_{i,t} by local-linear RD.  With a binary type the
#     "-" indicator is 1 minus the "+" indicator, so its jump is the exact
#     negative; a per-period Wald over both type indicators therefore has rank
#     1, and the reference uses the single "+" contrast per period (df = 1 per
#     period, df = #periods jointly).
#
# S3. Sides.  Inherited from Appendix B (helper A4): "+" is R >= c, so a unit
#     exactly at the cutoff counts as above.  In the reflected dataset of
#     ass:comp-stable this would put a period-t0 unit with R_{i,t0} == c on the
#     "+" side rather than "just below" as the prose says; the DGP builders
#     below draw R from continuous distributions, and `s44_no_ties()` asserts
#     that no unit sits exactly at the cutoff, so the case never arises.
#
# S4. Balanced panels.  A unit's type is defined only if it is observed in the
#     period that defines the type.  The builders return balanced panels
#     (every unit in every period) except `s44_dgp_split()`, which splits the
#     *comparison* periods across disjoint unit sets while keeping every unit
#     in tRD, so every type is still defined.
#
# S5. Variance of a jump.  V(D-hat) = V(beta0_+) + V(beta0_-) - 2 Cov(+, -).
#     For an ordinary regression the two sides use disjoint rows and the cross
#     term is exactly zero, recovering the Appendix B ref_V_D().  It is written
#     in full because the reflected regression of ass:comp-stable is the one
#     case where the two sides share *units* (a unit above the cutoff in both
#     periods contributes a row to each side), and the prose explicitly
#     requires that dependence to enter the variance.
#
# S6. Covariance of two jumps.  Cov(D_t, D_s) = Cov(+,+) + Cov(-,-)
#     - Cov(+,-) - Cov(-,+), each term an Appendix B.2 sandwich with
#     Sigma_{t,s} supported on shared unit ids.  The sampling scheme selects
#     which terms are kept: "cs" none (cross-period covariance zero), "pc"
#     the same-side pair only, "pv" all four.  This mirrors eq:var-cs /
#     eq:var-pc / eq:var-pv, whose reference implementation ref_aggregate()
#     makes exactly the same three choices.  The *diagonal* V(D_t) never
#     depends on the scheme.
#
# S7. Sandwiches as influence sums.  Each sandwich is evaluated as
#     sum_u g_t(u) g_s(u) with g(u) the unit's total intercept influence times
#     its residual, rather than by forming Psi.  This is algebraically the
#     same object as ref_cov_beta0() -- `s44_selfcheck_cov()` asserts the two
#     agree to 1e-12 on a live design -- but it also expresses the case
#     ref_cov_beta0() cannot, namely a Sigma with off-diagonal entries linking
#     two *rows of one regression* that belong to the same unit (S5).
#
# S8. Finite-sample factor.  The package's HC1 convention (helper-appB
#     assumption A8 and ref_hc1()) is applied inside the reference here, since
#     Section 4.4 quotes variances from Appendix B and the package's Appendix B
#     implementation carries the factor.  `hc1 = FALSE` switches it off.
#     Residual variant follows the package: p-fit at h for conventional, q-fit
#     at b for bias-corrected.  All tests below set b = h.
#
# S9. Bandwidths.  Every regression a test runs -- per period, per type, per
#     cell, and the reflected one -- uses the same fixed h, with b = h and
#     p = 1, q = 2.  Fixed h is passed explicitly so nothing is data-driven.
#
# S10. Wald.  theta' Sigma^{-1} theta ~ chi^2(K) with K = length(theta); no
#     small-sample correction.  Sigma is solved with solve().
#
# S11. Trend contrasts.  For `trend = "constant"` the null is that a cell's
#     jump is the same in every comparison period, so the contrasts are
#     D_t(v) - D_{t1}(v) for each cell v and each comparison period t after the
#     first.  Contrasts for different cells use disjoint unit sets, so their
#     covariance is exactly zero.
#
# S12. Rounding of the package's period labels.  Period values are compared as
#     characters (the package keys its per-period lists by the printed value),
#     and the builders use character period labels so the two orders agree.
#
# S13. Contrast orientation and names.  Section 4.4 states null hypotheses
#     ("test whether the within-type jumps are equal"), not an orientation, so
#     both the sign of a contrast and its name are free; the Wald statistic is
#     invariant to either.  The reference adopts the package's *reporting*
#     convention so that contrast vectors can be compared elementwise:
#     ass:homog uses type "+" as the baseline, contrast (-) - (+), named
#     "<period>::-"; ass:trend-cell uses the first comparison period as the
#     baseline, contrast D_t(v) - D_{t1}(v), named "<cell>::<t>-<t1>".  Both
#     were read off the `reference` flag of the package's returned
#     `$period_type_jumps` / `$cell_period_jumps` tables -- an interface fact,
#     not a formula -- and the conformance tests assert that flag explicitly so
#     the convention is visible where it is used.
# ---------------------------------------------------------------------------


## --- small utilities --------------------------------------------------------

# Assert that no observation sits exactly at the cutoff (S3).
s44_no_ties <- function(x, c = 0) {
  stopifnot(!any(x == c))
  invisible(TRUE)
}

# Pad one (unit, x, y) triple onto a common unit universe, in the universe's
# order (helper-appB A2/A5: units outside N_t get y = 0, R = c, inN = FALSE).
s44_pad <- function(unit, x, y, ids, c = 0) {
  m   <- match(ids, unit)
  inN <- !is.na(m)
  yy <- rep(0, length(ids)); rr <- rep(c, length(ids))
  yy[inN] <- y[m[inN]]
  rr[inN] <- x[m[inN]]
  list(y = yy, r = rr, inN = inN)
}


## --- one local-linear regression, carrying its rows' unit ids ---------------

# A single local-linear RD fit.  `unit` gives the unit id of each row; rows are
# kept in the given order and no padding is done, so the same unit may appear
# on both sides (the reflected design of ass:comp-stable).  `ids` optionally
# pads the fit onto a wider unit universe first, which is what the panel
# regressions do so that two periods' residual vectors line up.
s44_fit <- function(y, x, unit, h, c = 0, kernel = "triangular",
                    ids = NULL, p = 1L, q = 2L) {
  if (is.null(ids)) {
    f <- ref_period_fit(y, x, h = h, b = h, c = c, p = p, q = q,
                        kernel = kernel, inN = rep(TRUE, length(y)),
                        n = length(y))
    f$unit <- unit
  } else {
    pd <- s44_pad(unit, x, y, ids, c)
    f <- ref_period_fit(pd$y, pd$r, h = h, b = h, c = c, p = p, q = q,
                        kernel = kernel, inN = pd$inN, n = length(ids))
    f$unit <- ids
  }
  f
}

# Row-wise intercept influence a[i] = (1/n) e_0' Gamma^{-1} X' A  (conventional)
# or (1/n) e_0' Gamma^{-1} Q  (bias-corrected), i.e. the linear map from Y to
# beta0-hat on that side.  Rows off the side (or outside N_t) get exactly 0.
s44_influence <- function(f, side, bc = FALSE) {
  s <- f$sides[[side]]
  if (bc) as.vector(s$Gpi[1L, ] %*% s$Q) / f$n
  else    as.vector(s$Gpi[1L, ] %*% t(s$a_h * f$X_p)) / f$n
}

# g[i] = a[i] * eps_i, the Appendix B.2 summand.  Residual variant follows the
# package (S8): p-fit at h conventionally, q-fit at b for the BC sandwich.
s44_g <- function(f, side, bc = FALSE, hc1 = TRUE) {
  e <- if (bc) f$eps_q else f$eps_p
  if (hc1) e <- e * sqrt(ref_hc1(f, side, bc))
  s44_influence(f, side, bc) * e
}

# g aggregated to the unit level (a unit contributing two rows contributes the
# sum of their influences).  Returns a vector named by unit id.
s44_gu <- function(f, side, bc = FALSE, hc1 = TRUE) {
  g <- s44_g(f, side, bc = bc, hc1 = hc1)
  u <- as.character(f$unit)
  if (!anyDuplicated(u)) return(stats::setNames(g, u))
  vapply(split(g, u), sum, numeric(1))
}

# Cov(beta0_{f1,side1}, beta0_{f2,side2}) = sum over shared units of g1 * g2
# (S7).  Units in only one of the two fits contribute nothing.
s44_cov <- function(f1, side1, f2, side2, bc = FALSE, hc1 = TRUE) {
  g1 <- s44_gu(f1, side1, bc, hc1)
  g2 <- s44_gu(f2, side2, bc, hc1)
  sh <- intersect(names(g1), names(g2))
  if (!length(sh)) return(0)
  sum(g1[sh] * g2[sh])
}

# V(D-hat) for one regression, including the shared-unit cross term (S5).
s44_V_jump <- function(f, bc = FALSE, hc1 = TRUE) {
  s44_cov(f, "+", f, "+", bc, hc1) + s44_cov(f, "-", f, "-", bc, hc1) -
    2 * s44_cov(f, "+", f, "-", bc, hc1)
}

# Cov(D-hat_1, D-hat_2) under a sampling scheme (S6).
s44_cov_jump <- function(f1, f2, bc = FALSE, hc1 = TRUE,
                         scheme = c("pv", "pc", "cs")) {
  scheme <- match.arg(scheme)
  if (scheme == "cs") return(0)
  out <- s44_cov(f1, "+", f2, "+", bc, hc1) + s44_cov(f1, "-", f2, "-", bc, hc1)
  if (scheme == "pc") return(out)
  out - s44_cov(f1, "+", f2, "-", bc, hc1) - s44_cov(f1, "-", f2, "+", bc, hc1)
}

# Wald statistic theta' Sigma^{-1} theta (S10).
s44_wald <- function(theta, Sigma) {
  theta <- as.vector(theta)
  Sigma <- as.matrix(Sigma)
  stat <- as.vector(crossprod(theta, solve(Sigma, theta)))
  list(stat = stat, df = length(theta),
       p = stats::pchisq(stat, df = length(theta), lower.tail = FALSE))
}

# S7 self-check: on two ordinary padded fits the influence-sum form and the
# Appendix B.2 Psi sandwich of ref_cov_beta0() are the same number.
s44_selfcheck_cov <- function(f1, f2, bc = FALSE, tol = 1e-12) {
  wh <- if (bc) "q" else "p"
  sg <- ref_sigma(f1, f2, which = wh)
  for (s1 in c("+", "-")) for (s2 in c("+", "-")) {
    a <- s44_cov(f1, s1, f2, s2, bc = bc, hc1 = FALSE)
    b <- ref_cov_beta0(f1, s1, f2, s2, sg, bc = bc)
    stopifnot(isTRUE(all.equal(a, b, tolerance = tol)))
  }
  invisible(TRUE)
}


## --- ass:type-cont (A8) -----------------------------------------------------

# For each period t: local-linear RD of the type indicator 1{V_{i,s} = 1} on
# R_{i,t}, where s is the other period (S1/S2).  Returns the per-period jumps
# and Wald statistics, and the joint Wald over periods with the cross-period
# covariance of S6.
ref_typecont <- function(data, x, time, id, c = 0, h, kernel = "triangular",
                         bc = FALSE, scheme = c("pv", "pc", "cs"),
                         hc1 = TRUE) {
  scheme <- match.arg(scheme)
  tv <- as.character(data[[time]])
  periods <- sort(unique(tv))
  stopifnot(length(periods) == 2L)          # S1: "the other period"
  s44_no_ties(data[[x]], c)

  ids <- sort(unique(data[[id]]))
  per <- lapply(periods, function(p) data[tv == p, , drop = FALSE])
  names(per) <- periods

  # V_{i,s} = 1{R_{i,s} >= c}, looked up on the unit universe
  side_of <- function(d) {
    v <- rep(NA_real_, length(ids))
    v[match(d[[id]], ids)] <- as.numeric(d[[x]] >= c)
    v
  }
  V <- lapply(per, side_of)

  fits <- list()
  for (k in seq_along(periods)) {
    p <- periods[k]; s <- periods[3L - k]
    d <- per[[p]]
    w <- V[[s]][match(d[[id]], ids)]        # type indicator 1{V_{i,s} = 1}
    stopifnot(!any(is.na(w)))               # S4: balanced panel
    fits[[p]] <- s44_fit(w, d[[x]], d[[id]], h = h, c = c, kernel = kernel,
                         ids = ids)
  }

  Dh <- vapply(periods, function(p) if (bc) fits[[p]]$D_bc else fits[[p]]$D,
               numeric(1))
  names(Dh) <- periods

  K <- length(periods)
  S <- matrix(0, K, K, dimnames = list(periods, periods))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    S[i, j] <- if (i == j) s44_V_jump(fits[[periods[i]]], bc, hc1)
               else s44_cov_jump(fits[[periods[i]]], fits[[periods[j]]],
                                 bc, hc1, scheme)
  }

  pp <- lapply(periods, function(p)
    s44_wald(Dh[[p]], S[p, p, drop = FALSE]))
  names(pp) <- periods

  list(jumps = Dh, cov = S, per_period = pp, joint = s44_wald(Dh, S),
       fits = fits)
}


## --- ass:comp-stable (A9) ---------------------------------------------------

# The reflected dataset.  Units above the cutoff in t0 keep -(R_{i,t0} - c) and
# carry the outcome 1{V_{i,tRD} = 1}; units above the cutoff in tRD keep
# R_{i,tRD} - c and carry 1{V_{i,t0} = 1}.  The local-linear jump at the
# artificial cutoff 0 is pi-hat_{tRD,(+)}(1) - pi-hat_{t0,(+)}(1).
ref_compstable_data <- function(data, x, time, id, t_rd, t0, c = 0) {
  tv <- as.character(data[[time]])
  d0 <- data[tv == as.character(t0),   , drop = FALSE]
  d1 <- data[tv == as.character(t_rd), , drop = FALSE]
  s44_no_ties(c(d0[[x]], d1[[x]]), c)

  V0 <- stats::setNames(as.numeric(d0[[x]] >= c), as.character(d0[[id]]))
  V1 <- stats::setNames(as.numeric(d1[[x]] >= c), as.character(d1[[id]]))

  a0 <- d0[d0[[x]] >= c, , drop = FALSE]         # t0 units above the cutoff
  a1 <- d1[d1[[x]] >= c, , drop = FALSE]         # tRD units above the cutoff
  stopifnot(all(as.character(a0[[id]]) %in% names(V1)),
            all(as.character(a1[[id]]) %in% names(V0)))   # S4

  rbind(
    data.frame(unit = as.character(a0[[id]]),
               xr   = -(a0[[x]] - c),
               w    = as.numeric(V1[as.character(a0[[id]])]),
               grp  = "t0", stringsAsFactors = FALSE),
    data.frame(unit = as.character(a1[[id]]),
               xr   = a1[[x]] - c,
               w    = as.numeric(V0[as.character(a1[[id]])]),
               grp  = "trd", stringsAsFactors = FALSE))
}

# `shared = FALSE` drops the cross term, i.e. pretends the two groups are
# independent; the tests use it to show the dependence correction bites.
ref_compstable <- function(data, x, time, id, t_rd, t0, c = 0, h,
                           kernel = "triangular", bc = FALSE, hc1 = TRUE,
                           shared = TRUE) {
  cd <- ref_compstable_data(data, x, time, id, t_rd, t0, c)
  f  <- s44_fit(cd$w, cd$xr, cd$unit, h = h, c = 0, kernel = kernel)

  jump <- if (bc) f$D_bc else f$D
  V <- if (shared) s44_V_jump(f, bc, hc1)
       else s44_cov(f, "+", f, "+", bc, hc1) + s44_cov(f, "-", f, "-", bc, hc1)
  list(jump = jump, V = V, cross = s44_cov(f, "+", f, "-", bc, hc1),
       wald = s44_wald(jump, matrix(V, 1, 1)), fit = f, cdata = cd)
}


## --- within-type / within-cell jumps (shared by A10 and A7) -----------------

# One local-linear RD of `y` on `x` per (comparison period, type), on the
# subsample of units of that type.  Type is the tRD side (S1).
s44_type_fits <- function(data, y, x, time, id, comparisons, t_rd, c = 0, h,
                          kernel = "triangular") {
  tv <- as.character(data[[time]])
  comparisons <- as.character(comparisons)
  s44_no_ties(data[[x]], c)

  drd <- data[tv == as.character(t_rd), , drop = FALSE]
  typ <- stats::setNames(ifelse(drd[[x]] >= c, "+", "-"), as.character(drd[[id]]))

  out <- list()
  for (p in comparisons) {
    d <- data[tv == p, , drop = FALSE]
    tt <- typ[as.character(d[[id]])]
    stopifnot(!any(is.na(tt)))              # S4
    for (v in c("+", "-")) {
      dv  <- d[tt == v, , drop = FALSE]
      ids <- sort(unique(as.character(dv[[id]])))
      out[[paste(p, v, sep = "::")]] <-
        s44_fit(dv[[y]], dv[[x]], as.character(dv[[id]]), h = h, c = c,
                kernel = kernel, ids = ids)
    }
  }
  out
}

# Covariance of two within-type jumps.  Different types use disjoint unit sets,
# so the covariance is exactly zero there; s44_cov() returns 0 on its own
# (empty intersection), and it is short-circuited here for clarity.
s44_cov_typejump <- function(fits, k1, k2, bc, hc1, scheme) {
  if (identical(k1, k2)) return(s44_V_jump(fits[[k1]], bc, hc1))
  v1 <- sub("^.*::", "", k1); v2 <- sub("^.*::", "", k2)
  if (v1 != v2) return(0)
  s44_cov_jump(fits[[k1]], fits[[k2]], bc, hc1, scheme)
}


## --- ass:homog (A10) --------------------------------------------------------

# Within each comparison period, the difference of the two within-type jumps,
# jointly across comparison periods.  `ref_type` is the type held fixed as the
# baseline; the contrast is (other type) - (reference type), labelled
# "<period>::<other type>".  The default "+" and that labelling are the
# package's *reporting* convention, read off the `reference` flag of its
# `$period_type_jumps` table -- an interface fact, not a formula.  The Wald
# statistic is invariant to the choice.
ref_homog <- function(data, y, x, time, id, comparisons, t_rd, c = 0, h,
                      kernel = "triangular", bc = FALSE,
                      scheme = c("pv", "pc", "cs"), hc1 = TRUE,
                      ref_type = "-") {
  # Baseline = the all-below type ("-", V_{i,tRD} = 0), so the contrast is
  # D-hat(1) - D-hat(0) as in the paper's statement; the package uses the same
  # baseline (locale-independent radix order, all-below first).
  scheme <- match.arg(scheme)
  comparisons <- as.character(comparisons)
  oth <- setdiff(c("+", "-"), ref_type)
  fits <- s44_type_fits(data, y, x, time, id, comparisons, t_rd, c, h, kernel)

  Djt <- vapply(names(fits), function(k) if (bc) fits[[k]]$D_bc else fits[[k]]$D,
                numeric(1))

  kp <- paste(comparisons, oth,      sep = "::")   # non-reference type
  km <- paste(comparisons, ref_type, sep = "::")   # reference type
  theta <- stats::setNames(Djt[kp] - Djt[km], paste(comparisons, oth, sep = "::"))

  K <- length(comparisons)
  S <- matrix(0, K, K, dimnames = list(names(theta), names(theta)))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    # Cov(D_i(v') - D_i(v), D_j(v') - D_j(v)); cross-type terms are zero
    S[i, j] <- s44_cov_typejump(fits, kp[i], kp[j], bc, hc1, scheme) +
               s44_cov_typejump(fits, km[i], km[j], bc, hc1, scheme) -
               s44_cov_typejump(fits, kp[i], km[j], bc, hc1, scheme) -
               s44_cov_typejump(fits, km[i], kp[j], bc, hc1, scheme)
  }

  w <- s44_wald(theta, S)
  list(contrasts = theta, cov = S, statistic = w$stat, df = w$df,
       p_value = w$p, jumps = Djt, fits = fits)
}


## --- ass:trend-cell (A7) ----------------------------------------------------

# Within each cell, D-hat_t(v) - D-hat_{t1}(v) for every comparison period t
# after the first (S11), jointly over cells and periods.
ref_trendcell <- function(data, y, x, time, id, comparisons, t_rd, c = 0, h,
                          kernel = "triangular", bc = FALSE,
                          scheme = c("pv", "pc", "cs"), hc1 = TRUE) {
  scheme <- match.arg(scheme)
  comparisons <- as.character(comparisons)
  stopifnot(length(comparisons) >= 2L)
  fits <- s44_type_fits(data, y, x, time, id, comparisons, t_rd, c, h, kernel)

  Djt <- vapply(names(fits), function(k) if (bc) fits[[k]]$D_bc else fits[[k]]$D,
                numeric(1))

  base <- comparisons[1L]
  rest <- comparisons[-1L]
  cells <- c("+", "-")
  # one contrast per (cell, later period), labelled "<cell>::<period>-<base>"
  # (the package's reporting convention, read off its `$cell_period_jumps`)
  lab <- as.vector(outer(cells, rest,
                         function(v, p) paste0(v, "::", p, "-", base)))
  key <- function(p, v) paste(p, v, sep = "::")

  theta <- stats::setNames(numeric(length(lab)), lab)
  cur <- stats::setNames(as.list(lab), lab)
  bas <- stats::setNames(as.list(lab), lab)
  for (p in rest) for (v in cells) {
    l <- paste0(v, "::", p, "-", base)
    theta[[l]] <- Djt[[key(p, v)]] - Djt[[key(base, v)]]
    cur[[l]] <- key(p, v); bas[[l]] <- key(base, v)
  }
  theta <- theta[lab]

  K <- length(lab)
  S <- matrix(0, K, K, dimnames = list(lab, lab))
  for (a in lab) for (bq in lab) {
    S[a, bq] <-
      s44_cov_typejump(fits, cur[[a]], cur[[bq]], bc, hc1, scheme) -
      s44_cov_typejump(fits, cur[[a]], bas[[bq]], bc, hc1, scheme) -
      s44_cov_typejump(fits, bas[[a]], cur[[bq]], bc, hc1, scheme) +
      s44_cov_typejump(fits, bas[[a]], bas[[bq]], bc, hc1, scheme)
  }

  w <- s44_wald(theta, S)
  list(contrasts = theta, cov = S, statistic = w$stat, df = w$df,
       p_value = w$p, jumps = Djt, fits = fits)
}


## --- DGP builders -----------------------------------------------------------

# (i) Two-period panel, time-varying running variable.
#
#   eta_i ~ U(-1, 1);  R_{i,t} = eta_i + e_{i,t},  e ~ N(0, sd_e)
#
# so units switch sides between periods.  `share_jump != 0` shifts the tRD
# running variable by that amount for units above the cutoff in t0, which makes
# P(V_{i,tRD} = 1 | R_{i,t0} = r) jump at r = c: it turns the ass:type-cont and
# ass:comp-stable violations on.  `homog_break != 0` makes the period-t0
# confounding jump differ by tRD type, turning the ass:homog violation on.
# A unit-level outcome component u_i is shared across periods, so cross-period
# covariances are non-zero.
s44_dgp_pv2 <- function(seed, n = 2000, c = 0, share_jump = 0, homog_break = 0,
                        sd_e = 0.35, sd_u = 0.4, sd_y = 0.3,
                        jump_t0 = 0.4, jump_rd = 0.9) {
  set.seed(seed)
  eta <- stats::runif(n, -1, 1)
  u   <- stats::rnorm(n, 0, sd_u)
  R0  <- eta + stats::rnorm(n, 0, sd_e)
  V0  <- as.numeric(R0 >= c)
  R1  <- eta + stats::rnorm(n, 0, sd_e) + share_jump * V0
  V1  <- as.numeric(R1 >= c)
  Y0 <- ref_mfun(R0, c, 0) + (jump_t0 + homog_break * V1) * (R0 >= c) +
    u + stats::rnorm(n, 0, sd_y)
  Y1 <- ref_mfun(R1, c, 0) + jump_rd * (R1 >= c) + u + stats::rnorm(n, 0, sd_y)
  d <- rbind(
    data.frame(id = seq_len(n), time = "t0",  x = R0, y = Y0,
               stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = "trd", x = R1, y = Y1,
               stringsAsFactors = FALSE))
  s44_no_ties(d$x, c)
  d
}

# (ii) Three-period panel: comparison periods t1, t2 with a census-like
# (time-constant) running variable, and an RD period whose running variable
# moves, so cells are not collinear with the comparison-period side.
#
#   R_{i,t1} = R_{i,t2} = eta_i,   R_{i,tRD} = eta_i + e_i
#   cell v_i = 1{R_{i,tRD} >= c}
#   Y_{i,t} = m(R_{i,t}) + alpha_t(v_i) 1{R_{i,t} >= c} + u_i + eps_{i,t}
#
# with alpha_{t1}(v) = a0 + a1 v and alpha_{t2}(v) = alpha_{t1}(v)
# + break_t2 (0.5 + v).  `break_t2 = 0` satisfies ass:trend-cell's testable
# implication in both cells; `a1 != 0` breaks ass:homog without breaking it.
s44_dgp_cell3 <- function(seed, n = 2000, c = 0, break_t2 = 0,
                          a0 = 0.4, a1 = 0, sd_e = 0.35, sd_u = 0.4,
                          sd_y = 0.3, jump_rd = 0.9) {
  set.seed(seed)
  eta <- stats::runif(n, -1, 1)
  u   <- stats::rnorm(n, 0, sd_u)
  Rrd <- eta + stats::rnorm(n, 0, sd_e)
  v   <- as.numeric(Rrd >= c)
  a_t1 <- a0 + a1 * v
  a_t2 <- a_t1 + break_t2 * (0.5 + v)
  Y1 <- ref_mfun(eta, c, 0) + a_t1 * (eta >= c) + u + stats::rnorm(n, 0, sd_y)
  Y2 <- ref_mfun(eta, c, 0) + a_t2 * (eta >= c) + u + stats::rnorm(n, 0, sd_y)
  Yr <- ref_mfun(Rrd, c, 0) + jump_rd * (Rrd >= c) + u + stats::rnorm(n, 0, sd_y)
  d <- rbind(
    data.frame(id = seq_len(n), time = "t1",  x = eta, y = Y1,
               stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = "t2",  x = eta, y = Y2,
               stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = "trd", x = Rrd, y = Yr,
               stringsAsFactors = FALSE))
  s44_no_ties(d$x, c)
  d
}

# (iii) "Disjoint" designs: balanced panels in which no unit is inside the
# bandwidth in more than one period, so every cross-period sandwich is exactly
# zero and every joint statistic must collapse to a sum of per-period ones.
#
# Simply deleting rows will not do it: the package needs a unit observed in
# every period to give it a type at all (S4), so half the units are pushed a
# distance `far` (>> h) away from the cutoff in the period where they are meant
# not to count, on a randomly signed side so that their *type* still varies.
# Those rows carry exactly zero kernel weight and hence zero influence, but
# they keep the panel balanced.
#
# Two periods, for ass:type-cont.  Group A is near the cutoff in t0 and far in
# tRD; group B the other way round.
s44_dgp_disjoint2 <- function(seed, n = 3000, c = 0, far = 2,
                              sd_u = 0.4, sd_y = 0.3,
                              jump_t0 = 0.4, jump_rd = 0.9) {
  set.seed(seed)
  eta <- stats::runif(n, -1, 1)
  u   <- stats::rnorm(n, 0, sd_u)
  A   <- rep(c(TRUE, FALSE), length.out = n)[sample.int(n)]
  s0  <- sample(c(-1, 1), n, replace = TRUE)
  s1  <- sample(c(-1, 1), n, replace = TRUE)
  R0  <- ifelse(A, eta, eta + far * s0) + c
  R1  <- ifelse(A, eta + far * s1, eta) + c
  Y0 <- ref_mfun(R0, c, 0) + jump_t0 * (R0 >= c) + u + stats::rnorm(n, 0, sd_y)
  Y1 <- ref_mfun(R1, c, 0) + jump_rd * (R1 >= c) + u + stats::rnorm(n, 0, sd_y)
  d <- rbind(
    data.frame(id = seq_len(n), time = "t0",  x = R0, y = Y0,
               stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = "trd", x = R1, y = Y1,
               stringsAsFactors = FALSE))
  s44_no_ties(d$x, c)
  d
}

# Three periods, for ass:homog and ass:trend-cell.  Every unit is near the
# cutoff in tRD (so every cell is defined and non-degenerate); group A is near
# the cutoff in t1 only, group B in t2 only.
s44_dgp_disjoint3 <- function(seed, n = 3000, c = 0, far = 2, break_t2 = 0,
                              a0 = 0.4, a1 = 0, sd_e = 0.35, sd_u = 0.4,
                              sd_y = 0.3, jump_rd = 0.9) {
  set.seed(seed)
  eta <- stats::runif(n, -1, 1)
  u   <- stats::rnorm(n, 0, sd_u)
  A   <- rep(c(TRUE, FALSE), length.out = n)[sample.int(n)]
  s1  <- sample(c(-1, 1), n, replace = TRUE)
  s2  <- sample(c(-1, 1), n, replace = TRUE)
  R1  <- ifelse(A, eta, eta + far * s1) + c
  R2  <- ifelse(A, eta + far * s2, eta) + c
  Rrd <- eta + stats::rnorm(n, 0, sd_e) + c
  v   <- as.numeric(Rrd >= c)
  a_t1 <- a0 + a1 * v
  a_t2 <- a_t1 + break_t2 * (0.5 + v)
  Y1 <- ref_mfun(R1, c, 0) + a_t1 * (R1 >= c) + u + stats::rnorm(n, 0, sd_y)
  Y2 <- ref_mfun(R2, c, 0) + a_t2 * (R2 >= c) + u + stats::rnorm(n, 0, sd_y)
  Yr <- ref_mfun(Rrd, c, 0) + jump_rd * (Rrd >= c) + u + stats::rnorm(n, 0, sd_y)
  d <- rbind(
    data.frame(id = seq_len(n), time = "t1",  x = R1,  y = Y1,
               stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = "t2",  x = R2,  y = Y2,
               stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = "trd", x = Rrd, y = Yr,
               stringsAsFactors = FALSE))
  s44_no_ties(d$x, c)
  d
}
