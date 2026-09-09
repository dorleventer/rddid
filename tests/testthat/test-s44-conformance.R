# ---------------------------------------------------------------------------
# Conformance of rddid's four validation tests against Section 4.4 of
# Leventer and Nevo.
#
# Every expectation compares a package quantity to the independent reference in
# helper-s44-reference.R, which was written from the Section 4.4 prose plus the
# Appendix B machinery of helper-appB-reference.R, without reading
# R/test_typecont.R, R/test_compstable.R, R/test_homog.R, R/test_trendcell.R or
# R/test_helpers.R.  Test titles start with the paper label they check.
#
# Conventions that Section 4.4 does not pin down are listed as S1-S13 in the
# reference's header.  The two that are visible in the numbers are re-asserted
# here rather than left implicit:
#
#  (i) Orientation and naming of contrasts (S13).  The paper states equality
#      nulls, not orientations.  The tests read the package's own `reference`
#      flag out of `$period_type_jumps` / `$cell_period_jumps` and check that
#      it is the type/period the reference uses as baseline, so a silent flip
#      of the package's convention shows up as a failure here and not as an
#      unexplained sign difference.
#
# (ii) Sampling scheme (S6).  "cs" drops every cross-period covariance, "pc"
#      keeps the same-side pair, "pv" keeps all four -- the same three choices
#      eq:var-cs / eq:var-pc / eq:var-pv make for the aggregate variance.  The
#      scheme is passed explicitly wherever it can change a number.
#
# Bandwidths are fixed (`h = s44_h`, hence b = h) so nothing is data-driven.
# ---------------------------------------------------------------------------

s44_tol <- 1e-10   # relative, expect_equal semantics
s44_h   <- 0.30

# Labelled scalar comparison, as in test-appB-conformance.R: `info=` is
# silently ignored by expect_equal() under testthat's 3rd edition, so the design
# label goes into the names instead.
s44_eq <- function(actual, expected, lab, tolerance = s44_tol) {
  a <- as.vector(actual); e <- as.vector(expected)
  nm <- if (length(a) == 1L) lab else paste0(lab, "[", seq_along(a), "]")
  expect_equal(stats::setNames(a, nm), stats::setNames(e, nm),
               tolerance = tolerance)
}

# Appendix B variance of one local-linear jump, with the package's HC1 factor
# (helper-appB A8 / ref_hc1) applied here in the test source so the departure
# from a literal reading of B.2 stays visible.
s44_V_appB <- function(f, bc = FALSE) {
  sg <- ref_sigma(f, f, which = if (bc) "q" else "p")
  ref_cov_beta0(f, "+", f, "+", sg, bc) * ref_hc1(f, "+", bc) +
    ref_cov_beta0(f, "-", f, "-", sg, bc) * ref_hc1(f, "-", bc)
}

# One-sided local-linear intercept at the cutoff by plain weighted least
# squares -- a route through lm() rather than through any matrix code in either
# the package or the Appendix B reference.
s44_one_sided <- function(y, x, h, c = 0, above = TRUE) {
  keep <- if (above) x >= c else x < c
  yy <- y[keep]; xx <- x[keep] - c
  w  <- ref_kernel(xx / h, "triangular")
  ok <- w > 0
  unname(stats::coef(stats::lm(yy[ok] ~ xx[ok], weights = w[ok]))[1L])
}

s44_pkg_typecont <- function(d, seed = 1L, ...) {
  set.seed(seed)
  rd_typecont(d, x = "x", time = "time", id = "id", c = 0, h = s44_h,
              kernel = "triangular", ...)
}

s44_pkg_compstable <- function(d, seed = 1L, ...) {
  set.seed(seed)
  rd_compstable(d, x = "x", time = "time", id = "id", t_rd = "trd",
                comparisons = "t0", c = 0, h = s44_h,
                kernel = "triangular", ...)
}

s44_pkg_homog <- function(d, comparisons, ...) {
  rd_homog(d, y = "y", x = "x", time = "time", id = "id",
           comparisons = comparisons, t_rd = "trd", c = 0, h = s44_h,
           kernel = "triangular", min_n = 10L, type_by = "rd_side", ...)
}

s44_pkg_trendcell <- function(d, comparisons, ...) {
  rd_trendcell(d, y = "y", x = "x", time = "time", id = "id",
               comparisons = comparisons, t_rd = "trd", c = 0, h = s44_h,
               kernel = "triangular", min_n = 10L, type_by = "rd_side",
               trend = "constant", ...)
}


# ---------------------------------------------------------------------------
test_that("ass:type-cont -- per-period side-share jumps and Wald match", {
  # A8 tested "in each period t, [by] a local-linear RD of the side indicator
  # 1{V_is = 1} on R_it", period by period and jointly.
  designs <- list(hold = s44_dgp_pv2(101, share_jump = 0),
                  fail = s44_dgp_pv2(404, share_jump = 0.30))

  for (nm in names(designs)) {
    d <- designs[[nm]]
    for (bc in c(FALSE, TRUE)) {
      lab <- paste0(nm, "/bc=", bc)
      r <- ref_typecont(d, "x", "time", "id", c = 0, h = s44_h, bc = bc,
                        scheme = "pv")
      p <- s44_pkg_typecont(d, scheme = "auto", bc = bc)

      # the running variable moves across periods, so "auto" must land on pv
      expect_identical(p$meta$scheme, "pv")

      s44_eq(p$ll_wald$stat, r$joint$stat, paste0("joint stat ", lab))
      s44_eq(p$ll_wald$df,   r$joint$df,   paste0("joint df ", lab))
      s44_eq(p$ll_wald$p,    r$joint$p,    paste0("joint p ", lab))
      expect_identical(p$ll_wald$df, 2L)          # two periods, binary type

      for (per in names(r$per_period)) {
        s44_eq(p$per_period[[per]]$ll_wald$stat, r$per_period[[per]]$stat,
               paste0("stat ", per, " ", lab))
        s44_eq(p$per_period[[per]]$ll_wald$df, r$per_period[[per]]$df,
               paste0("df ", per, " ", lab))
        expect_identical(p$per_period[[per]]$ll_wald$df, 1L)
      }
    }

    # The jump the reference feeds the Wald is the local-linear RD of the type
    # indicator: rebuild the indicator here and run ref_period_fit on it.
    r <- ref_typecont(d, "x", "time", "id", c = 0, h = s44_h, scheme = "pv")
    d0 <- d[d$time == "t0", ];  d1 <- d[d$time == "trd", ]
    V0 <- stats::setNames(as.numeric(d0$x >= 0), as.character(d0$id))
    V1 <- stats::setNames(as.numeric(d1$x >= 0), as.character(d1$id))
    w_t0  <- as.numeric(V1[as.character(d0$id)])   # 1{V_{i,tRD} = 1}
    w_trd <- as.numeric(V0[as.character(d1$id)])   # 1{V_{i,t0}  = 1}

    f0 <- ref_period_fit(w_t0,  d0$x, h = s44_h, b = s44_h, c = 0)
    f1 <- ref_period_fit(w_trd, d1$x, h = s44_h, b = s44_h, c = 0)
    s44_eq(f0$D, r$jumps[["t0"]],  paste0("appB jump t0 ", nm))
    s44_eq(f1$D, r$jumps[["trd"]], paste0("appB jump trd ", nm))

    # and the same two jumps by weighted lm(), i.e. pi-hat_(+) - pi-hat_(-)
    s44_eq(s44_one_sided(w_t0, d0$x, s44_h, 0, TRUE) -
             s44_one_sided(w_t0, d0$x, s44_h, 0, FALSE),
           r$jumps[["t0"]], paste0("lm jump t0 ", nm))
    s44_eq(s44_one_sided(w_trd, d1$x, s44_h, 0, TRUE) -
             s44_one_sided(w_trd, d1$x, s44_h, 0, FALSE),
           r$jumps[["trd"]], paste0("lm jump trd ", nm))

    # S7: the influence-sum form of every sandwich is ref_cov_beta0()'s Psi
    expect_silent(s44_selfcheck_cov(r$fits[["t0"]], r$fits[["trd"]], bc = FALSE))
    expect_silent(s44_selfcheck_cov(r$fits[["t0"]], r$fits[["trd"]], bc = TRUE))
  }

  # The DGP switch has to move the statistic, or the comparisons above are
  # vacuous.  Sanity bounds only (not part of the conformance check): the
  # null design must not reject at 5%, the violated one must reject at 1%.
  quiet <- ref_typecont(designs$hold, "x", "time", "id", c = 0, h = s44_h,
                        scheme = "pv")$joint$stat
  loud  <- ref_typecont(designs$fail, "x", "time", "id", c = 0, h = s44_h,
                        scheme = "pv")$joint$stat
  expect_lt(quiet, stats::qchisq(0.95, 2))
  expect_gt(loud,  stats::qchisq(0.99, 2))
})


# ---------------------------------------------------------------------------
test_that("ass:type-cont -- joint Wald uses the cross-period covariance", {
  # "the off-diagonal element ... is non-zero in panel data, as the same units
  # enter both periods' regressions."
  d <- s44_dgp_pv2(101)

  r_cs <- ref_typecont(d, "x", "time", "id", c = 0, h = s44_h, scheme = "cs")
  r_pv <- ref_typecont(d, "x", "time", "id", c = 0, h = s44_h, scheme = "pv")

  # under "cs" the off-diagonal is dropped, so the joint is the sum of the
  # per-period chi-squares; under "pv" it is not
  expect_identical(r_cs$cov[1, 2], 0)
  expect_gt(abs(r_pv$cov[1, 2]), 0)

  sum_per <- sum(vapply(r_cs$per_period, function(z) z$stat, numeric(1)))
  s44_eq(r_cs$joint$stat, sum_per, "ref cs joint == sum of per-period")

  p_cs <- s44_pkg_typecont(d, seed = 2L, scheme = "cs", bc = FALSE)
  p_pv <- s44_pkg_typecont(d, seed = 3L, scheme = "pv", bc = FALSE)
  p_au <- s44_pkg_typecont(d, seed = 4L, scheme = "auto", bc = FALSE)

  s44_eq(p_cs$ll_wald$stat, r_cs$joint$stat, "pkg cs joint")
  s44_eq(p_pv$ll_wald$stat, r_pv$joint$stat, "pkg pv joint")
  s44_eq(p_au$ll_wald$stat, r_pv$joint$stat, "pkg auto joint == pv")

  # the package's own per-period statistics do not depend on the scheme, and
  # under "cs" the joint is exactly their sum
  pkg_per <- vapply(p_cs$per_period, function(z) z$ll_wald$stat, numeric(1))
  s44_eq(p_cs$ll_wald$stat, sum(pkg_per), "pkg cs joint == sum of per-period")
  s44_eq(vapply(p_pv$per_period, function(z) z$ll_wald$stat, numeric(1)),
         pkg_per, "per-period stats are scheme-free")

  # and with the cross term the joint genuinely differs from that sum
  expect_gt(abs(p_pv$ll_wald$stat / sum(pkg_per) - 1), 1e-3)
})


# ---------------------------------------------------------------------------
test_that("ass:comp-stable -- reflected-cutoff jump and its dependence-corrected variance match", {
  # A9: reflect the t0-above units to -(R - c), append the tRD-above units at
  # R - c, cross-assign the other period's side indicator, and test the jump.
  designs <- list(hold = s44_dgp_pv2(101, share_jump = 0),
                  fail = s44_dgp_pv2(404, share_jump = 0.30))

  for (nm in names(designs)) {
    d <- designs[[nm]]
    for (bc in c(FALSE, TRUE)) {
      lab <- paste0(nm, "/bc=", bc)
      r <- ref_compstable(d, "x", "time", "id", t_rd = "trd", t0 = "t0",
                          c = 0, h = s44_h, bc = bc)
      p <- s44_pkg_compstable(d, scheme = "auto", bc = bc)
      pr <- p$pairs[["trd::t0"]]

      s44_eq(pr$ll_wald$stat, r$wald$stat, paste0("pair stat ", lab))
      s44_eq(pr$ll_wald$df,   r$wald$df,   paste0("pair df ", lab))
      s44_eq(pr$ll_wald$p,    r$wald$p,    paste0("pair p ", lab))
      expect_identical(pr$ll_wald$df, 1L)

      # a single (t_rd, t0) pair: the joint equals the pair
      s44_eq(p$joint$ll_wald$stat, r$wald$stat, paste0("joint == pair ", lab))
    }

    # the jump is pi-hat_{tRD,(+)}(1) - pi-hat_{t0,(+)}(1): two *one-sided*
    # local-linear intercepts, both taken from above the cutoff in the original
    # coordinates, computed here by weighted lm()
    r  <- ref_compstable(d, "x", "time", "id", "trd", "t0", c = 0, h = s44_h)
    d0 <- d[d$time == "t0", ];  d1 <- d[d$time == "trd", ]
    V0 <- stats::setNames(as.numeric(d0$x >= 0), as.character(d0$id))
    V1 <- stats::setNames(as.numeric(d1$x >= 0), as.character(d1$id))
    pi_rd <- s44_one_sided(as.numeric(V0[as.character(d1$id)]), d1$x, s44_h, 0, TRUE)
    pi_t0 <- s44_one_sided(as.numeric(V1[as.character(d0$id)]), d0$x, s44_h, 0, TRUE)
    s44_eq(pi_rd - pi_t0, r$jump, paste0("reflected jump ", nm))

    # the two sides of the reflected regression share units, so the cross term
    # is non-zero and dropping it changes the variance by a visible amount
    rn <- ref_compstable(d, "x", "time", "id", "trd", "t0", c = 0, h = s44_h,
                         shared = FALSE)
    expect_gt(abs(r$cross), 0)
    s44_eq(rn$V - r$V, 2 * r$cross, paste0("cross term accounting ", nm))
    expect_gt(abs(rn$V / r$V - 1), 1e-3)

    # the package matches the dependence-corrected variance, not the naive one
    p <- s44_pkg_compstable(d, seed = 5L, scheme = "auto", bc = FALSE)
    s44_eq(p$pairs[["trd::t0"]]$ll_wald$stat, r$wald$stat,
           paste0("pkg == shared-unit variance ", nm))
    expect_gt(abs(p$pairs[["trd::t0"]]$ll_wald$stat / rn$wald$stat - 1), 1e-3)
    expect_gt(p$pairs[["trd::t0"]]$n_both, 0)   # units above the cutoff twice
  }

  # sanity bounds only: the null design must not reject at 5%, the violated
  # one must reject at 1%
  quiet <- ref_compstable(designs$hold, "x", "time", "id", "trd", "t0",
                          c = 0, h = s44_h)$wald$stat
  loud  <- ref_compstable(designs$fail, "x", "time", "id", "trd", "t0",
                          c = 0, h = s44_h)$wald$stat
  expect_lt(quiet, stats::qchisq(0.95, 1))
  expect_gt(loud,  stats::qchisq(0.99, 1))
})


# ---------------------------------------------------------------------------
test_that("ass:homog -- within-type jumps on disjoint subsamples, equality Wald matches", {
  # A10: D-hat_{t0}(0) on {V_{i,tRD} = 0} and D-hat_{t0}(1) on {V_{i,tRD} = 1},
  # then test that the two within-type jumps are equal.

  ## --- one comparison period: df = 1, disjoint subsamples ------------------
  for (nm in c("hold", "fail")) {
    d <- s44_dgp_pv2(66, homog_break = if (nm == "fail") 0.6 else 0)
    for (bc in c(FALSE, TRUE)) {
      lab <- paste0(nm, "/bc=", bc)
      r <- ref_homog(d, "y", "x", "time", "id", comparisons = "t0",
                     t_rd = "trd", c = 0, h = s44_h, bc = bc, scheme = "pv")
      p <- s44_pkg_homog(d, "t0", scheme = "auto", bc = bc)

      s44_eq(p$statistic, r$statistic, paste0("stat ", lab))
      s44_eq(p$df,        r$df,        paste0("df ", lab))
      s44_eq(p$p_value,   r$p_value,   paste0("p ", lab))
      expect_identical(p$df, 1L)
      s44_eq(p$contrasts, r$contrasts, paste0("contrast ", lab))
      expect_identical(names(p$contrasts), names(r$contrasts))
      s44_eq(as.vector(p$cov_matrix), as.vector(r$cov), paste0("cov ", lab))

      # S13: the package's baseline is the all-below type, the reference's baseline too
      ptj <- p$period_type_jumps
      expect_identical(ptj$type[ptj$reference], "-")

      # stat = (D-hat(1) - D-hat(0))^2 / (V(1) + V(0)), each piece a separate
      # local-linear RD on its own subsample and its own Appendix B variance
      d0  <- d[d$time == "t0", ];  drd <- d[d$time == "trd", ]
      typ <- stats::setNames(drd$x >= 0, as.character(drd$id))
      v1  <- typ[as.character(d0$id)]
      f1  <- ref_period_fit(d0$y[v1],  d0$x[v1],  h = s44_h, b = s44_h, c = 0)
      f0  <- ref_period_fit(d0$y[!v1], d0$x[!v1], h = s44_h, b = s44_h, c = 0)
      D1  <- if (bc) f1$D_bc else f1$D
      D0  <- if (bc) f0$D_bc else f0$D
      stat <- (D1 - D0)^2 / (s44_V_appB(f1, bc) + s44_V_appB(f0, bc))
      s44_eq(stat, p$statistic, paste0("closed-form stat ", lab))
      s44_eq(abs(D1 - D0), abs(p$contrasts[[1L]]),
             paste0("closed-form contrast ", lab))

      # with one comparison period there is no cross-period term to model, so
      # the answer cannot depend on the scheme
      for (sch in c("cs", "pc", "pv"))
        s44_eq(s44_pkg_homog(d, "t0", scheme = sch, bc = bc)$statistic,
               p$statistic, paste0("scheme-free ", sch, " ", lab))
    }
  }

  # sanity bounds only: the DGP switch has to move the statistic across the
  # 5% / 1% critical values
  expect_lt(ref_homog(s44_dgp_pv2(66, homog_break = 0), "y", "x", "time", "id",
                      "t0", "trd", c = 0, h = s44_h, scheme = "pv")$statistic,
            stats::qchisq(0.95, 1))
  expect_gt(ref_homog(s44_dgp_pv2(66, homog_break = 0.6), "y", "x", "time", "id",
                      "t0", "trd", c = 0, h = s44_h, scheme = "pv")$statistic,
            stats::qchisq(0.99, 1))

  ## --- two comparison periods, shared units, "pc" data ---------------------
  # R is time-constant across t1 and t2 (census-like) and moves only in tRD, so
  # the cells are not collinear with the comparison-period side.
  d3 <- s44_dgp_cell3(202, a1 = 0.25, break_t2 = 0)
  for (bc in c(FALSE, TRUE)) {
    lab <- paste0("pc/bc=", bc)
    r <- ref_homog(d3, "y", "x", "time", "id", comparisons = c("t1", "t2"),
                   t_rd = "trd", c = 0, h = s44_h, bc = bc, scheme = "pc")
    p <- s44_pkg_homog(d3, c("t1", "t2"), scheme = "pc", bc = bc)

    expect_identical(p$df, 2L)
    s44_eq(p$statistic, r$statistic, paste0("joint stat ", lab))
    s44_eq(p$df,        r$df,        paste0("joint df ", lab))
    s44_eq(p$contrasts, r$contrasts, paste0("contrasts ", lab))
    expect_identical(names(p$contrasts), names(r$contrasts))
    s44_eq(as.vector(p$cov_matrix), as.vector(r$cov), paste0("cov ", lab))

    # the units are shared across the two comparison periods, so the joint is
    # not the sum of the two per-period statistics
    expect_gt(abs(r$cov[1, 2]), 0)
    per <- vapply(c("t1", "t2"), function(pp)
      s44_pkg_homog(d3, pp, scheme = "pc", bc = bc)$statistic, numeric(1))
    expect_gt(abs(p$statistic / sum(per) - 1), 1e-3)

    # x is time-constant across t1 and t2, so side membership is too: the
    # opposite-side cross terms vanish and "pv" collapses onto "pc"
    s44_eq(ref_homog(d3, "y", "x", "time", "id", c("t1", "t2"), "trd",
                     c = 0, h = s44_h, bc = bc, scheme = "pv")$statistic,
           r$statistic, paste0("pv == pc on time-constant x ", lab))
  }
})


# ---------------------------------------------------------------------------
test_that("ass:trend-cell -- within-cell cross-period contrasts and Wald match", {
  # A7's suggestive test: with T0 = {t1, t2}, D_{t1}(v) = D_{t2}(v) for each
  # cell v, the cell being the period-tRD side.
  designs <- list(hold = s44_dgp_cell3(606, a1 = 0.2, break_t2 = 0),
                  fail = s44_dgp_cell3(606, a1 = 0.2, break_t2 = 0.4))

  for (nm in names(designs)) {
    d3 <- designs[[nm]]
    for (bc in c(FALSE, TRUE)) {
      lab <- paste0(nm, "/bc=", bc)
      r <- ref_trendcell(d3, "y", "x", "time", "id", comparisons = c("t1", "t2"),
                         t_rd = "trd", c = 0, h = s44_h, bc = bc, scheme = "pc")
      p <- s44_pkg_trendcell(d3, c("t1", "t2"), scheme = "pc", bc = bc)

      # two cells x one cross-period contrast each
      expect_identical(p$df, 2L)
      s44_eq(p$statistic, r$statistic, paste0("stat ", lab))
      s44_eq(p$df,        r$df,        paste0("df ", lab))
      s44_eq(p$p_value,   r$p_value,   paste0("p ", lab))
      s44_eq(p$contrasts, r$contrasts, paste0("contrasts ", lab))
      expect_identical(names(p$contrasts), names(r$contrasts))
      s44_eq(as.vector(p$cov_matrix), as.vector(r$cov), paste0("cov ", lab))

      # S13: the package's baseline period is the reference's baseline period
      cpj <- p$cell_period_jumps
      expect_true(all(cpj$period[cpj$reference] == "t1"))

      # per-cell jumps: D-hat_t(v) from a local-linear RD on the cell subsample
      drd <- d3[d3$time == "trd", ]
      typ <- stats::setNames(ifelse(drd$x >= 0, "+", "-"), as.character(drd$id))
      for (per in c("t1", "t2")) for (v in c("+", "-")) {
        dv <- d3[d3$time == per, ]
        dv <- dv[typ[as.character(dv$id)] == v, ]
        f  <- ref_period_fit(dv$y, dv$x, h = s44_h, b = s44_h, c = 0)
        s44_eq(if (bc) f$D_bc else f$D, r$jumps[[paste(per, v, sep = "::")]],
               paste0("cell jump ", per, "/", v, " ", lab))
      }

      # each contrast is (later period) - (first comparison period) within cell
      for (v in c("+", "-"))
        s44_eq(r$jumps[[paste("t2", v, sep = "::")]] -
                 r$jumps[[paste("t1", v, sep = "::")]],
               r$contrasts[[paste0(v, "::t2-t1")]],
               paste0("contrast build ", v, " ", lab))

      # the two cells are disjoint unit sets, so their contrasts are
      # uncorrelated; the within-cell cross-period term is not zero
      expect_identical(r$cov[1, 2], 0)
      expect_identical(p$cov_matrix[1, 2], 0)
      s44_eq(p$statistic, sum(r$contrasts^2 / diag(r$cov)),
             paste0("diagonal Wald ", lab))
    }
  }

  # sanity bounds only, as above
  expect_lt(ref_trendcell(designs$hold, "y", "x", "time", "id", c("t1", "t2"),
                          "trd", c = 0, h = s44_h, scheme = "pc")$statistic,
            stats::qchisq(0.95, 2))
  expect_gt(ref_trendcell(designs$fail, "y", "x", "time", "id", c("t1", "t2"),
                          "trd", c = 0, h = s44_h, scheme = "pc")$statistic,
            stats::qchisq(0.99, 2))
})


# ---------------------------------------------------------------------------
test_that("Section 4.4 -- with no shared units across periods the joint statistics are sums of per-period ones", {
  # Structural check.  The package needs a unit in every period to give it a
  # type at all, so "no shared units" is realised by keeping the panel balanced
  # and pushing half the units far outside the bandwidth in each period: their
  # rows carry exactly zero kernel weight, hence zero influence, so no unit
  # contributes to two periods' regressions and every Sigma_{t,s} (t != s)
  # vanishes.  Each joint statistic must then be the sum of its parts.

  ## --- ass:type-cont -------------------------------------------------------
  d2 <- s44_dgp_disjoint2(707)
  for (bc in c(FALSE, TRUE)) {
    lab <- paste0("typecont/bc=", bc)
    r <- ref_typecont(d2, "x", "time", "id", c = 0, h = s44_h, bc = bc,
                      scheme = "pv")
    expect_identical(r$cov[1, 2], 0)      # exactly, not to tolerance

    p <- s44_pkg_typecont(d2, seed = 6L, scheme = "pv", bc = bc)
    per <- vapply(p$per_period, function(z) z$ll_wald$stat, numeric(1))
    s44_eq(p$ll_wald$stat, sum(per),      paste0("joint == sum ", lab))
    s44_eq(p$ll_wald$stat, r$joint$stat,  paste0("joint == ref ", lab))

    # and the scheme becomes irrelevant once the off-diagonal is zero
    s44_eq(s44_pkg_typecont(d2, seed = 7L, scheme = "cs", bc = bc)$ll_wald$stat,
           p$ll_wald$stat, paste0("cs == pv ", lab))
  }

  ## --- ass:homog and ass:trend-cell ---------------------------------------
  d3 <- s44_dgp_disjoint3(808, a1 = 0.2, break_t2 = 0.15)
  for (bc in c(FALSE, TRUE)) {
    lab <- paste0("homog/bc=", bc)
    rh <- ref_homog(d3, "y", "x", "time", "id", c("t1", "t2"), "trd",
                    c = 0, h = s44_h, bc = bc, scheme = "pv")
    ph <- s44_pkg_homog(d3, c("t1", "t2"), scheme = "pv", bc = bc)
    expect_identical(rh$cov[1, 2], 0)
    expect_identical(ph$cov_matrix[1, 2], 0)
    s44_eq(ph$statistic, rh$statistic, paste0("joint == ref ", lab))
    per <- vapply(c("t1", "t2"), function(pp)
      s44_pkg_homog(d3, pp, scheme = "pv", bc = bc)$statistic, numeric(1))
    s44_eq(ph$statistic, sum(per), paste0("joint == sum of per-period ", lab))

    lab <- paste0("trendcell/bc=", bc)
    rt <- ref_trendcell(d3, "y", "x", "time", "id", c("t1", "t2"), "trd",
                        c = 0, h = s44_h, bc = bc, scheme = "pv")
    pt <- s44_pkg_trendcell(d3, c("t1", "t2"), scheme = "pv", bc = bc)
    s44_eq(pt$statistic, rt$statistic, paste0("joint == ref ", lab))
    s44_eq(as.vector(pt$cov_matrix), as.vector(rt$cov), paste0("cov ", lab))
    # cells are disjoint and periods no longer share units, so the covariance
    # of the two contrasts is diagonal and each variance is a plain sum
    expect_identical(rt$cov[1, 2], 0)
    for (v in c("+", "-")) {
      k <- paste0(v, "::t2-t1")
      vv <- s44_V_appB(rt$fits[[paste("t1", v, sep = "::")]], bc) +
            s44_V_appB(rt$fits[[paste("t2", v, sep = "::")]], bc)
      s44_eq(rt$cov[k, k], vv, paste0("variance is a plain sum ", v, " ", lab))
    }
  }
})
