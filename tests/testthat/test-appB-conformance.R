# ---------------------------------------------------------------------------
# Conformance of rddid against the Appendix B / Section 5.2 equations.
#
# Every expectation compares a package quantity to the independent reference in
# helper-appB-reference.R, which was written from the paper's equations alone
# (literal matrix form, base R) without reading R/rd_period.R or R/aggregate.R.
# Test titles start with the paper label they check.
#
# Two package conventions are NOT in the paper.  They are applied here, in the
# test source, so that a reader sees exactly where the package departs from a
# literal reading of Appendix B:
#
#  (i) BC residuals.  The paper's B.2 uses eps-hat (own-side p-order fit at
#      h_t) in Sigma-hat for every sandwich.  The package uses the own-side
#      q-order fit at b_t for the bias-corrected sandwich (rdrobust
#      convention).  Both are computed below; the package is asserted equal to
#      the q-at-b variant and the gap to the p-at-h variant is reported.
#
# (ii) HC1 factor.  The package multiplies each residual by
#      sqrt(n_side / (n_side - (p+1))) for the conventional fit and
#      sqrt(n_side / (n_side - (q+1))) for the bias-corrected fit, where
#      n_side counts the units on that side with positive *pilot*-window kernel
#      weight.  The reference has no such factor, so every variance comparison
#      below multiplies the reference by f (or by sqrt(f_t f_s) for a
#      cross-period, cross-side term).  Each such test also checks that the
#      unfactored reference is within 2% of the package, i.e. that the
#      correction is a finite-sample detail and not a structural difference.
# ---------------------------------------------------------------------------

appB_tol  <- 1e-10   # relative, expect_equal semantics
appB_kern <- c("triangular", "epanechnikov", "uniform")
appB_pq   <- list(c(1L, 2L), c(2L, 3L))
appB_cut  <- c(0, 0.3)                  # second cutoff catches a hard-coded c
appB_h    <- c(0.30, 0.22, 0.26)        # period-specific bandwidths
appB_b    <- 1.4 * appB_h               # pilot band, b_t >= h_t throughout

# Labelled scalar comparison.  `info=` is deliberately NOT used: it is silently
# ignored by expect_equal() under testthat's 3rd edition and, worse, suppresses
# the supplied `tolerance` if these tests are ever run under the 2nd edition.
# Naming both scalars puts the design label in the failure output instead, and
# comparing one scalar at a time keeps a single bad design from being averaged
# away against the others.
appB_eq <- function(actual, expected, lab, tolerance = appB_tol) {
  expect_equal(stats::setNames(actual, lab), stats::setNames(expected, lab),
               tolerance = tolerance)
}


test_that("App. B preamble -- kernel shapes agree and the 1/h, 1/n scalings cancel", {
  u <- seq(-1.6, 1.6, length.out = 201)
  for (k in appB_kern) {
    # the reference's K() is the same shape the package weights with
    expect_equal(ref_kernel(u, k), rddid:::.rd_kweight(u, k), tolerance = appB_tol)
  }

  d <- ref_dgp_single(11, n = 1500, c = 0.3)
  for (pq in appB_pq) {
    # A3: K_h(u) = K(u/h)/h vs the un-normalised K(u/h)
    r1 <- ref_period_fit(d$y, d$x, h = 0.3, b = 0.42, c = 0.3,
                         p = pq[1], q = pq[2], divide_by_h = TRUE)
    r0 <- ref_period_fit(d$y, d$x, h = 0.3, b = 0.42, c = 0.3,
                         p = pq[1], q = pq[2], divide_by_h = FALSE)
    expect_equal(r0$D, r1$D, tolerance = appB_tol)
    expect_equal(r0$D_bc, r1$D_bc, tolerance = appB_tol)
    expect_equal(ref_V_D(r0, bc = FALSE), ref_V_D(r1, bc = FALSE), tolerance = appB_tol)
    expect_equal(ref_V_D(r0, bc = TRUE),  ref_V_D(r1, bc = TRUE),  tolerance = appB_tol)
  }

  # A2: the literal 1/n factors cancel, so padding a period into a larger unit
  # universe (as the panel formulas do) leaves every quantity unchanged
  dl   <- ref_dgp_panel(12, n = 1200, c = 0, design = "pv", drop_frac = 0.1)
  rfs  <- ref_panel_fits(dl, appB_h, appB_b, c = 0)
  solo <- ref_period_fit(dl[["2"]]$y, dl[["2"]]$x, h = appB_h[2], b = appB_b[2], c = 0)
  expect_equal(solo$D,    rfs[["2"]]$D,    tolerance = appB_tol)
  expect_equal(solo$D_bc, rfs[["2"]]$D_bc, tolerance = appB_tol)
  expect_equal(ref_V_D(solo), ref_V_D(rfs[["2"]]), tolerance = appB_tol)
  expect_equal(ref_V_D(solo, bc = TRUE), ref_V_D(rfs[["2"]], bc = TRUE),
               tolerance = appB_tol)
})


test_that("eq:wls_rd -- conventional intercepts and D-hat match the matrix form", {
  for (cut in appB_cut) for (pq in appB_pq) for (k in appB_kern) {
    d  <- ref_dgp_single(101, n = 1500, c = cut)
    h  <- 0.30; b <- 1.4 * h
    pk <- rd_period(d$y, d$x, h = h, b = b, id = d$id, c = cut,
                    p = pq[1], q = pq[2], kernel = k)
    rf <- ref_period_fit(d$y, d$x, h = h, b = b, c = cut,
                         p = pq[1], q = pq[2], kernel = k)
    lab <- sprintf("c=%g p=%d q=%d kernel=%s", cut, pq[1], pq[2], k)
    appB_eq(pk$sides[["+"]]$beta0, rf$sides[["+"]]$beta0, lab)
    appB_eq(pk$sides[["-"]]$beta0, rf$sides[["-"]]$beta0, lab)
    appB_eq(pk$D, rf$D, lab)
    # beta-hat^{(1)} = 1! e_1' beta-hat is the reported local slope
    appB_eq(pk$sides[["+"]]$slope, rf$sides[["+"]]$beta_p[2], lab)
    appB_eq(pk$sides[["-"]]$slope, rf$sides[["-"]]$beta_p[2], lab)
  }
})


test_that("eq:bc_intercept == eq:bc_Q -- two forms agree; package BC intercepts match", {
  for (cut in appB_cut) for (pq in appB_pq) for (k in appB_kern) {
    d  <- ref_dgp_single(102, n = 1600, c = cut)
    h  <- 0.28; b <- 1.5 * h
    pk <- rd_period(d$y, d$x, h = h, b = b, id = d$id, c = cut,
                    p = pq[1], q = pq[2], kernel = k)
    rf <- ref_period_fit(d$y, d$x, h = h, b = b, c = cut,
                         p = pq[1], q = pq[2], kernel = k)
    lab <- sprintf("c=%g p=%d q=%d kernel=%s", cut, pq[1], pq[2], k)
    for (sd in c("+", "-")) {
      # eq:bc_intercept and eq:bc_Q are the same estimator (a reference-internal
      # identity: two floating-point paths through an ill-conditioned inverse,
      # so 1e-8 rather than the 1e-10 used for package-vs-reference checks)
      appB_eq(rf$sides[[sd]]$beta0_bc, rf$sides[[sd]]$beta0_bc_Q, lab,
              tolerance = 1e-8)
      appB_eq(pk$sides[[sd]]$beta0_bc, rf$sides[[sd]]$beta0_bc, lab)
    }
    appB_eq(pk$D_bc, rf$D_bc, lab)
  }
})


test_that("B.1 bias constant B_{t,(s),p}(h) is finite and D - D^BC equals the estimated bias", {
  for (cut in appB_cut) for (pq in appB_pq) {
    d  <- ref_dgp_single(103, n = 1700, c = cut)
    h  <- 0.26; b <- 1.6 * h
    p  <- pq[1]
    pk <- rd_period(d$y, d$x, h = h, b = b, id = d$id, c = cut,
                    p = pq[1], q = pq[2])
    rf <- ref_period_fit(d$y, d$x, h = h, b = b, c = cut, p = pq[1], q = pq[2])
    lab <- sprintf("c=%g p=%d q=%d", cut, pq[1], pq[2])

    Bp <- rf$sides[["+"]]$B; Bm <- rf$sides[["-"]]$B
    expect_true(is.finite(Bp) && is.finite(Bm))
    # estimated bias of D-hat: (h^{p+1}/(p+1)!) [beta^{(p+1)}_+ B_+ - beta^{(p+1)}_- B_-]
    bias <- (h^(p + 1) / factorial(p + 1)) *
      (rf$sides[["+"]]$beta_deriv_p1 * Bp - rf$sides[["-"]]$beta_deriv_p1 * Bm)
    appB_eq(rf$D - rf$D_bc, bias, lab)
    appB_eq(pk$D - pk$D_bc, bias, lab)
    # the correction is not numerically trivial in this design
    expect_true(abs(bias) > 1e-4 * abs(rf$D))
  }
})


test_that("B.2 single-period sandwich -- V(D-hat) matches (with HC1 factors made explicit)", {
  bc_resid_gap <- c()
  for (cut in appB_cut) for (pq in appB_pq) for (k in appB_kern) {
    d  <- ref_dgp_single(104, n = 1600, c = cut)
    h  <- 0.30; b <- 1.4 * h
    pk <- rd_period(d$y, d$x, h = h, b = b, id = d$id, c = cut,
                    p = pq[1], q = pq[2], kernel = k)
    rf <- ref_period_fit(d$y, d$x, h = h, b = b, c = cut,
                         p = pq[1], q = pq[2], kernel = k)
    lab <- sprintf("c=%g p=%d q=%d kernel=%s", cut, pq[1], pq[2], k)

    ## ---- conventional: Sigma-hat from the p-fit at h (paper) --------------
    sg   <- ref_sigma(rf, rf, which = "p")
    Vp   <- ref_cov_beta0(rf, "+", rf, "+", sg, bc = FALSE)
    Vm   <- ref_cov_beta0(rf, "-", rf, "-", sg, bc = FALSE)
    f_p  <- ref_hc1(rf, "+", bc = FALSE)   # n_side / (n_side - (p+1))
    f_m  <- ref_hc1(rf, "-", bc = FALSE)
    appB_eq(pk$V_D, f_p * Vp + f_m * Vm, lab)
    # V(beta0) side by side, and V(D) = V(+) + V(-) since the sides are disjoint
    appB_eq(sum(pk$sides[["+"]]$g^2), f_p * Vp, lab)
    appB_eq(sum(pk$sides[["-"]]$g^2), f_m * Vm, lab)
    appB_eq(Vp + Vm, ref_V_D(rf, bc = FALSE), lab)
    # the HC1 factor is a small finite-sample detail
    expect_true(abs(pk$V_D - (Vp + Vm)) / pk$V_D < 0.02)

    ## ---- bias corrected: Psi^BC = (1/n) Q Sigma Q' ------------------------
    f_pb <- ref_hc1(rf, "+", bc = TRUE)    # n_side / (n_side - (q+1))
    f_mb <- ref_hc1(rf, "-", bc = TRUE)
    sg_q <- ref_sigma(rf, rf, which = "q") # package convention: q-fit at b
    Vbp  <- ref_cov_beta0(rf, "+", rf, "+", sg_q, bc = TRUE)
    Vbm  <- ref_cov_beta0(rf, "-", rf, "-", sg_q, bc = TRUE)
    appB_eq(pk$V_D_bc, f_pb * Vbp + f_mb * Vbm, lab)
    appB_eq(sum(pk$sides[["+"]]$g_bc^2), f_pb * Vbp, lab)
    appB_eq(sum(pk$sides[["-"]]$g_bc^2), f_mb * Vbm, lab)
    expect_true(abs(pk$V_D_bc - (Vbp + Vbm)) / pk$V_D_bc < 0.02)

    # informational: the paper's text uses eps-hat (p-fit at h) everywhere.
    # Record, do not assert, how far that variant is from the package.
    sg_p <- ref_sigma(rf, rf, which = "p")
    Vbp2 <- ref_cov_beta0(rf, "+", rf, "+", sg_p, bc = TRUE)
    Vbm2 <- ref_cov_beta0(rf, "-", rf, "-", sg_p, bc = TRUE)
    gap  <- abs(pk$V_D_bc - (f_pb * Vbp2 + f_mb * Vbm2)) / pk$V_D_bc
    expect_true(is.finite(gap))
    bc_resid_gap <- c(bc_resid_gap, gap)
  }
  message(sprintf(
    "[informational] V_D_bc under the paper's eps-hat (p-fit at h) vs the package's q-fit at b: relative gap %.3g to %.3g over %d designs",
    min(bc_resid_gap), max(bc_resid_gap), length(bc_resid_gap)))
})


test_that("eq:cross-decomp / C^same, C^opp -- cross-period intercept covariances match g-vector products", {
  for (cut in appB_cut) for (des in c("pc", "pv")) {
    dl  <- ref_dgp_panel(205, n = 1500, c = cut, design = des,
                         drop_frac = if (des == "pv") 0.1 else 0)
    rfs <- ref_panel_fits(dl, appB_h, appB_b, c = cut)
    pfs <- ref_pkg_fits(dl,  appB_h, appB_b, c = cut)
    nm  <- names(dl)

    for (bc in c(FALSE, TRUE)) {
      wh <- if (bc) "q" else "p"
      for (i in seq_along(nm)) for (j in seq_along(nm))
        for (st in c("+", "-")) for (ss in c("+", "-")) {
          lab <- sprintf("c=%g %s bc=%s (%s,%s)-(%s,%s)", cut, des, bc,
                         nm[i], st, nm[j], ss)
          sg  <- ref_sigma(rfs[[i]], rfs[[j]], which = wh)
          ref <- ref_cov_beta0(rfs[[i]], st, rfs[[j]], ss, sg, bc = bc)
          fac <- sqrt(ref_hc1(rfs[[i]], st, bc) * ref_hc1(rfs[[j]], ss, bc))
          pkg <- ref_pkg_gcov(pfs[[i]], st, pfs[[j]], ss, bc = bc)
          appB_eq(pkg, fac * ref, lab)
          if (abs(pkg) > 0) expect_true(abs(pkg - ref) / abs(pkg) < 0.02)
        }

      # Cov(b-_t, b+_s) = Cov(b+_s, b-_t): the reference builds both, so the
      # 2P x 2P covariance is symmetric
      Cpm <- outer(seq_along(nm), seq_along(nm), Vectorize(function(i, j)
        ref_cov_beta0(rfs[[i]], "+", rfs[[j]], "-",
                      ref_sigma(rfs[[i]], rfs[[j]], which = wh), bc = bc)))
      Cmp <- outer(seq_along(nm), seq_along(nm), Vectorize(function(i, j)
        ref_cov_beta0(rfs[[i]], "-", rfs[[j]], "+",
                      ref_sigma(rfs[[i]], rfs[[j]], which = wh), bc = bc)))
      expect_true(max(abs(Cmp - t(Cpm))) <= appB_tol * max(abs(Cpm)))
      # within a period the two sides are disjoint: C^opp_{t,t} = 0 exactly
      expect_true(all(diag(Cpm) == 0), info = sprintf("c=%g %s bc=%s", cut, des, bc))
    }
  }
})


test_that("eq:var-cs, eq:var-pc, eq:var-pv -- .aggregate_fits matches under CS / PC / PV data", {
  weights <- list(c("3" = 1, "1" = -0.5,  "2" = -0.5),
                  c("3" = 1, "1" = -0.75, "2" = -0.25))
  grid <- list(list(des = "cs", cut = 0,   k = "triangular",   pq = c(1L, 2L)),
               list(des = "pc", cut = 0,   k = "triangular",   pq = c(1L, 2L)),
               list(des = "pv", cut = 0,   k = "triangular",   pq = c(1L, 2L)),
               list(des = "cs", cut = 0.3, k = "triangular",   pq = c(2L, 3L)),
               list(des = "pc", cut = 0.3, k = "epanechnikov", pq = c(1L, 2L)),
               list(des = "pv", cut = 0.3, k = "uniform",      pq = c(2L, 3L)))

  for (g in grid) {
    dl  <- ref_dgp_panel(206, n = 1500, c = g$cut, design = g$des,
                         drop_frac = if (g$des == "pv") 0.1 else 0)
    rfs <- ref_panel_fits(dl, appB_h, appB_b, c = g$cut,
                          p = g$pq[1], q = g$pq[2], kernel = g$k)
    pfs <- ref_pkg_fits(dl,  appB_h, appB_b, c = g$cut,
                        p = g$pq[1], q = g$pq[2], kernel = g$k)
    for (cf in weights) for (bc in c(FALSE, TRUE)) {
      lab <- sprintf("%s c=%g kernel=%s p=%d bc=%s w=(%g,%g)", g$des, g$cut, g$k,
                     g$pq[1], bc, -cf[["1"]], -cf[["2"]])
      a <- rddid:::.aggregate_fits(pfs, cf, bc = bc)
      r <- ref_agg_from_fits(rfs, cf, rd = "3", bc = bc, hc1 = TRUE)
      appB_eq(unname(a[["est"]]),  unname(r[["est"]]), lab)
      appB_eq(unname(a[["V_cs"]]), unname(r[["V_cs"]]), lab)
      appB_eq(unname(a[["V_pc"]]), unname(r[["V_pc"]]), lab)
      appB_eq(unname(a[["V_pv"]]), unname(r[["V_pv"]]), lab)
      # HC1-free reference: same numbers to within the finite-sample factor
      r0 <- ref_agg_from_fits(rfs, cf, rd = "3", bc = bc, hc1 = FALSE)
      for (v in c("V_cs", "V_pc", "V_pv"))
        expect_true(abs(a[[v]] - r0[[v]]) / abs(a[[v]]) < 0.02, info = paste(lab, v))
    }
  }
})


test_that("eq:var-cs/pc/pv -- PC data: C^opp = 0 and V^PV == V^PC; CS data: V^CS == V^PC == V^PV", {
  cf <- c("3" = 1, "1" = -0.6, "2" = -0.4)

  ## panel, time-constant running variable: no unit ever switches side, so
  ## every C^opp entry vanishes and eq:var-pv collapses onto eq:var-pc
  for (cut in appB_cut) {
    dl  <- ref_dgp_panel(207, n = 1500, c = cut, design = "pc")
    rfs <- ref_panel_fits(dl, appB_h, appB_b, c = cut)
    pfs <- ref_pkg_fits(dl,  appB_h, appB_b, c = cut)
    for (bc in c(FALSE, TRUE)) {
      r <- ref_agg_from_fits(rfs, cf, rd = "3", bc = bc)
      Copp <- attr(r, "Cpm") + attr(r, "Cmp")
      expect_true(all(Copp == 0), info = sprintf("PC c=%g bc=%s", cut, bc))
      expect_equal(unname(r[["V_pv"]]), unname(r[["V_pc"]]), tolerance = appB_tol)
      a <- rddid:::.aggregate_fits(pfs, cf, bc = bc)
      expect_equal(unname(a[["V_pv"]]), unname(a[["V_pc"]]), tolerance = appB_tol)
      # and PC genuinely differs from CS here (the cross terms are not empty)
      expect_true(abs(a[["V_pc"]] - a[["V_cs"]]) / a[["V_cs"]] > 1e-3)
    }
  }

  ## repeated cross-section: no shared ids, so every cross-period covariance
  ## vanishes and the three variances coincide
  for (cut in appB_cut) {
    dl  <- ref_dgp_panel(208, n = 1500, c = cut, design = "cs")
    rfs <- ref_panel_fits(dl, appB_h, appB_b, c = cut)
    pfs <- ref_pkg_fits(dl,  appB_h, appB_b, c = cut)
    for (bc in c(FALSE, TRUE)) {
      r <- ref_agg_from_fits(rfs, cf, rd = "3", bc = bc)
      Csame <- attr(r, "Cpp") + attr(r, "Cmm")
      expect_true(all(Csame[upper.tri(Csame)] == 0), info = sprintf("CS c=%g bc=%s", cut, bc))
      expect_equal(unname(r[["V_pc"]]), unname(r[["V_cs"]]), tolerance = appB_tol)
      expect_equal(unname(r[["V_pv"]]), unname(r[["V_cs"]]), tolerance = appB_tol)
      a <- rddid:::.aggregate_fits(pfs, cf, bc = bc)
      expect_equal(unname(a[["V_pc"]]), unname(a[["V_cs"]]), tolerance = appB_tol)
      expect_equal(unname(a[["V_pv"]]), unname(a[["V_cs"]]), tolerance = appB_tol)
    }
  }

  ## panel, time-varying running variable: units do switch side, so the PV
  ## correction is live and the three variances separate
  dl  <- ref_dgp_panel(209, n = 1500, c = 0, design = "pv", drop_frac = 0.1)
  rfs <- ref_panel_fits(dl, appB_h, appB_b, c = 0)
  r   <- ref_agg_from_fits(rfs, cf, rd = "3", bc = FALSE)
  Copp <- attr(r, "Cpm") + attr(r, "Cmp")
  expect_true(any(abs(Copp) > 0))
  expect_true(abs(r[["V_pv"]] - r[["V_pc"]]) / r[["V_pc"]] > 1e-6)
})
