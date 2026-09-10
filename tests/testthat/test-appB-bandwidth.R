# Appendix B.4 (app:est-bw) conformance: bandwidth objectives and selectors.
# Each test title starts with the paper object it checks; see dev/appB_map.md §2.5.
# The hand-coded formulas below are written from the paper, not from R/bandwidth.R.

.b4_panel <- function(scheme = c("pc", "cs", "pv"), n = 1500, seed = 11) {
  scheme <- match.arg(scheme)
  set.seed(seed)
  Tn <- 3; t_rd <- 3
  R0 <- runif(n, -1, 1)
  ue <- rnorm(n, 0, 0.4)
  rows <- lapply(seq_len(Tn), function(t) {
    R <- switch(scheme,
      cs = runif(n, -1, 1),
      pc = R0,
      pv = pmin(1, pmax(-1, R0 + rnorm(n, 0, 0.15))))
    id <- if (scheme == "cs") (t - 1) * n + seq_len(n) else seq_len(n)
    m <- 0.3 * t + 0.5 * R + (0.4 + 0.3 * t) * R^2 + 0.4 * (R >= 0) + (t == t_rd) * (R >= 0)
    u <- if (scheme == "cs") rnorm(n, 0, 0.5) else 0.35 * rnorm(n) + ue[id]
    data.frame(id = id, year = t, R = R, Y = m + u)
  })
  do.call(rbind, rows)
}

.b4_plist <- function(d, c = 0) {
  periods <- c(3, 1, 2)
  stats::setNames(lapply(periods, function(tv) {
    rows <- d$year == tv
    data.frame(y = d$Y[rows], x = d$R[rows], id = d$id[rows])
  }), as.character(periods))
}
.b4_coef <- c(`3` = 1, `1` = -0.5, `2` = -0.5)

test_that("B.4 constants general p — rd_period b_const = (p+1)! (D - D_bc) / h^{p+1}", {
  d <- .b4_panel("pc"); d3 <- d[d$year == 3, ]
  for (p in 1:3) {
    f <- rd_period(d3$Y, d3$R, h = 0.35, b = 0.5, id = d3$id, p = p, q = p + 1L)
    expect_equal(f$b_const, factorial(p + 1) * (f$D - f$D_bc) / 0.35^(p + 1), tolerance = 1e-12)
    expect_equal(f$v_const, f$n * 0.35 * f$V_D, tolerance = 1e-12)   # eq:per-period-orders
  }
})

test_that("B.4 common h (eq:common_h_opt) — .bw_joint returns the closed form in its own constants, per-period pilots, p = 1, 2", {
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  # period-specific pilots (the selector fits each period at its own pilot)
  pb <- list(`3` = c(h = 0.40, b = 0.55), `1` = c(h = 0.30, b = 0.45), `2` = c(h = 0.50, b = 0.62))
  for (p in 1:2) {
    q <- p + 1L
    jb <- .bw_joint(pl, .b4_coef, "3", scheme = "cs", pilot_bws = pb, p = p, q = q)
    fp1 <- factorial(p + 1)
    # closed form from B.4 P3 with the regularization in the denominator
    h_closed <- (fp1^2 / (2 * (p + 1)) * jb$Veff / (jb$B^2 + jb$reg))^(1 / (2 * p + 3))
    expect_equal(jb$h, h_closed, tolerance = 1e-12)
    # b_t = h* b_t^pilot / h_t^pilot, per period (B.4 after eq:common_h_opt)
    for (k in names(.b4_coef))
      expect_equal(unname(jb$b[[k]]), jb$h * pb[[k]][["b"]] / pb[[k]][["h"]], tolerance = 1e-12)
    # hand-built constants, each period at ITS OWN pilot
    fits <- lapply(names(.b4_coef), function(k)
      rd_period(pl[[k]]$y, pl[[k]]$x, h = pb[[k]][["h"]], b = pb[[k]][["b"]], id = pl[[k]]$id, p = p, q = q))
    names(fits) <- names(.b4_coef)
    # B = sum_tau coef_tau B-hat_tau (B.4 P3), from rd_period's b_const
    B_hand <- sum(.b4_coef * vapply(fits, function(f) f$b_const, numeric(1)))
    expect_equal(jb$B, B_hand, tolerance = 1e-12)
    # Veff = sum_tau coef_tau^2 V-hat_tau / n_tau = sum coef^2 h_tau V_D(h_tau)
    # (eq:per-period-orders; lem:agg-var at a common h)
    V_hand <- sum(.b4_coef^2 * vapply(names(.b4_coef), function(k) pb[[k]][["h"]] * fits[[k]]$V_D, numeric(1)))
    expect_equal(jb$Veff, V_hand, tolerance = 1e-12)
    # regularization: lambda * sum coef^2 ((p+1)!/h_tau^{p+1})^2 Var(D - D_bc), each at its own pilot
    vd <- sum(.b4_coef^2 * vapply(names(.b4_coef), function(k) {
      f <- fits[[k]]
      (fp1 / pb[[k]][["h"]]^(p + 1))^2 * (sum(f$sides[["+"]]$g_diff^2) + sum(f$sides[["-"]]$g_diff^2)) }, numeric(1)))
    expect_equal(jb$reg, 3 * vd, tolerance = 1e-12)
    if (p == 1L) expect_equal(jb$h, (jb$Veff / (jb$B^2 + jb$reg))^(1 / 5), tolerance = 1e-12)
  }
})

test_that("B.4 common h — h* is the scalar minimizer of the iterative rule's own objective (cs and pc)", {
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  pb <- list(`3` = c(h = 0.40, b = 0.55), `1` = c(h = 0.30, b = 0.45), `2` = c(h = 0.50, b = 0.62))
  for (sc in c("cs", "pc")) {
    jb <- .bw_joint(pl, .b4_coef, "3", scheme = sc, pilot_bws = pb)
    ib <- .bw_joint_iter(pl, .b4_coef, "3", scheme = sc, pilot_bws = pb, start = "cct")
    f1 <- function(h) ib$amse_fun(rep(h, 3))
    h_opt <- stats::optimize(f1, interval = c(0.02, 1), tol = 1e-10)$minimum
    expect_equal(jb$h, h_opt, tolerance = 1e-6, info = sc)
    # the objective at h* equals the closed-form objective of .bw_joint
    expect_equal(f1(jb$h), (jb$h^2 / 2)^2 * (jb$B^2 + jb$reg) + jb$Veff / jb$h, tolerance = 1e-10, info = sc)
  }
})

test_that("B.4 joint selectors do not depend on which period carries the t_rd label", {
  # an aggregate over two RD periods (3 and 2) with comparison period 1, coef (1, 1, -2):
  # the same linear combination whichever RD period is labelled t_rd
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  coef_a <- c(`3` = 1, `2` = 1, `1` = -2)
  coef_b <- c(`2` = 1, `3` = 1, `1` = -2)
  for (sc in c("cs", "pc")) {
    ja <- .bw_joint(pl, coef_a, "3", scheme = sc)
    jb <- .bw_joint(pl, coef_b, "2", scheme = sc)
    expect_equal(ja$h, jb$h, tolerance = 1e-12, info = sc)
    expect_equal(ja$b[names(coef_a)], jb$b[names(coef_a)], tolerance = 1e-12, info = sc)
    expect_equal(c(ja$B, ja$Veff, ja$reg), c(jb$B, jb$Veff, jb$reg), tolerance = 1e-12, info = sc)
    # the descent cycles periods in names(coef) order, so tighten its stopping rule
    # to compare the two fixed points rather than two stopping points
    ia <- .bw_joint_iter(pl, coef_a, "3", scheme = sc, tol = 1e-10, maxit = 500L)
    ib <- .bw_joint_iter(pl, coef_b, "2", scheme = sc, tol = 1e-10, maxit = 500L)
    expect_equal(unlist(ia$bws[names(coef_a)]), unlist(ib$bws[names(coef_a)]), tolerance = 1e-7, info = sc)
    expect_equal(ia$objective, ib$objective, tolerance = 1e-10, info = sc)
  }
  # through rddid(): relabelling t_rd with the matching comparison weights
  ra <- rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, comparisons = c(2, 1),
              weights = c(-1, 2), bwselect = "joint", scheme = "pc")
  rb <- rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 2, comparisons = c(3, 1),
              weights = c(-1, 2), bwselect = "joint", scheme = "pc")
  cols <- c("est", "se_cs", "se_pc", "se_pv")
  expect_equal(ra$estimates[, cols], rb$estimates[, cols], tolerance = 1e-10)
  expect_equal(ra$bandwidth$h, rb$bandwidth$h, tolerance = 1e-12)
  for (k in c("1", "2", "3")) expect_equal(ra$fits[[k]]$b, rb$fits[[k]]$b, tolerance = 1e-12)
})

test_that("B.4 scheme switch — PV drops the cross-period term (lem:cov-pv): joint and iter agree with CS", {
  d <- .b4_panel("pv"); pl <- .b4_plist(d)
  pilot <- c(h = 0.4, b = 0.55)
  pb <- stats::setNames(rep(list(pilot), 3), names(.b4_coef))
  j_cs <- .bw_joint(pl, .b4_coef, "3", scheme = "cs", pilot_bws = pb)
  j_pv <- .bw_joint(pl, .b4_coef, "3", scheme = "pv", pilot_bws = pb)
  expect_equal(j_cs$h, j_pv$h, tolerance = 1e-12)
  i_cs <- .bw_joint_iter(pl, .b4_coef, "3", scheme = "cs", pilot_bws = pb, start = "cct")
  i_pv <- .bw_joint_iter(pl, .b4_coef, "3", scheme = "pv", pilot_bws = pb, start = "cct")
  expect_equal(unlist(i_cs$bws), unlist(i_pv$bws), tolerance = 1e-10)
  # and PC differs from CS on the same data (the cross term is live)
  i_pc <- .bw_joint_iter(pl, .b4_coef, "3", scheme = "pc", pilot_bws = pb, start = "cct")
  expect_false(isTRUE(all.equal(unlist(i_cs$bws), unlist(i_pc$bws), tolerance = 1e-6)))
})

test_that("B.4 objective (eq:amse-att, lem:agg-var) — amse_fun equals a hand-coded objective, PC, p = 1 and 2", {
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  for (p in 1:2) {
    q <- p + 1L
    pb <- list(`3` = c(h = 0.40, b = 0.55), `1` = c(h = 0.30, b = 0.45), `2` = c(h = 0.36, b = 0.50))
    ib <- .bw_joint_iter(pl, .b4_coef, "3", scheme = "pc", pilot_bws = pb, start = "cct",
                         p = p, q = q, maxit = 3L)
    keys <- names(.b4_coef); cf <- unname(.b4_coef)
    fitp <- lapply(keys, function(k) rd_period(pl[[k]]$y, pl[[k]]$x, h = pb[[k]][["h"]], b = pb[[k]][["b"]],
                                              id = pl[[k]]$id, p = p, q = q))
    names(fitp) <- keys
    fp1 <- factorial(p + 1)
    bt <- vapply(fitp, function(f) f$b_const, numeric(1))
    vt <- vapply(fitp, function(f) f$v_const, numeric(1))
    nt <- vapply(fitp, function(f) f$n, numeric(1))
    h0 <- vapply(keys, function(k) pb[[k]][["h"]], numeric(1))
    varb <- vapply(fitp, function(f) sum(f$sides[["+"]]$g_diff^2) + sum(f$sides[["-"]]$g_diff^2), numeric(1)) *
      (fp1 / h0^(p + 1))^2
    # per-side pilot same-side covariances and their h-free scales (lem:cov-pc)
    kap <- function(i, j, side) {
      gi <- fitp[[i]]$sides[[side]]; gj <- fitp[[j]]$sides[[side]]
      m <- match(gi$id, gj$id); ok <- !is.na(m)
      P <- sum(gi$g[ok] * gj$g[m[ok]])
      P * h0[j] / .kc_c(p, side, h0[i] / h0[j])
    }
    hand <- function(hv) {
      Bbar <- sum(cf * hv^(p + 1) * bt) / fp1
      pen  <- 3 * sum(cf^2 * (hv^(p + 1) / fp1)^2 * varb)
      vv   <- sum(cf^2 * vt / (nt * hv))
      cov  <- 0
      for (i in 1:2) for (j in (i + 1):3) {
        rho <- hv[i] / hv[j]
        cov <- cov + 2 * cf[i] * cf[j] *
          (kap(i, j, "+") * .kc_c(p, "+", rho) + kap(i, j, "-") * .kc_c(p, "-", rho)) / hv[j]
      }
      unname(Bbar^2 + pen + vv + cov)
    }
    for (hv in list(c(0.3, 0.3, 0.3), c(0.45, 0.25, 0.33), c(0.2, 0.5, 0.4)))
      expect_equal(ib$amse_fun(hv), hand(hv), tolerance = 1e-10)
    # returned bandwidths are a coordinate-wise minimizer of the objective (up to optimize's tolerance)
    h <- vapply(keys, function(k) ib$bws[[k]][["h"]], numeric(1))
    expect_equal(ib$objective, hand(h), tolerance = 1e-10)
    if (ib$niter < 3L) {   # converged: each coordinate is a local minimizer
      for (j in 1:3) for (d in c(-0.02, 0.02)) {
        hv <- h; hv[j] <- hv[j] * (1 + d)
        expect_gte(hand(hv), hand(h) - 1e-9)
      }
    }
  }
})

test_that("B.4 PC cross term (lem:cov-pc) — symmetric in (t, s) and equals the pilot covariance at rho = 1", {
  for (p in 1:2) for (kern in c("triangular", "uniform", "epanechnikov")) {
    for (rho in c(0.5, 0.8, 1, 1.6, 2.5)) {
      # c(1/rho) = rho c(rho): the term kappa c(h_t/h_s)/h_s is invariant to the ordering
      expect_equal(.kc_c(p, "+", 1 / rho, kern), rho * .kc_c(p, "+", rho, kern), tolerance = 1e-8)
    }
    expect_equal(.kc_c(p, "+", 1, kern), .kc_v(p, "+", kern), tolerance = 1e-10)
  }
  # at a common bandwidth the new term reduces to the plug-in covariance itself:
  # kappa * c(1) / h0 = P (the old max(h) form and the new form coincide at rho = 1)
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  f1 <- rd_period(pl[["1"]]$y, pl[["1"]]$x, h = 0.4, b = 0.55, id = pl[["1"]]$id)
  f2 <- rd_period(pl[["2"]]$y, pl[["2"]]$x, h = 0.4, b = 0.55, id = pl[["2"]]$id)
  cc <- .cross_cov(f1, f2)
  expect_equal(cc$pc, cc$pc_p + cc$pc_m, tolerance = 1e-12)
  kap_p <- cc$pc_p * 0.4 / .kc_c(1, "+", 1)
  expect_equal(kap_p * .kc_c(1, "+", 1) / 0.4, cc$pc_p, tolerance = 1e-12)
})

test_that("B.4 uniform p = 0 sanity — c(rho) = min(1, 1/rho) so the PC term is kappa / max(h_t, h_s)", {
  for (rho in c(0.25, 0.5, 1, 2, 4))
    expect_equal(.kc_c(0, "+", rho, "uniform"), min(1, 1 / rho), tolerance = 1e-8)
})
