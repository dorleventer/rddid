# Appendix B.4 (app:est-bw) conformance: bandwidth objectives and selectors.
# Each test title starts with the paper object it checks; see dev/appB_map.md §2.5.
# The hand-coded formulas below are written from the paper, not from R/bandwidth_joint.R.

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

# B-hat_t from scratch (App. B.4, "Bandwidth estimation"): on each side, the design moment
# e_0'(X'A(h)X)^{-1} X'A(h) ((x - c)/h)^{p+1} at h (an lm intercept) times the coefficient on
# (x - c)^{p+1} of the order-q fit at b (an lm coefficient), scaled by (p+1)!; triangular kernel
.b4_Bhat <- function(y, x, h, b, p, q) {
  tri <- function(u) pmax(0, 1 - abs(u))
  side <- function(s) {
    ys <- y[s]; xs <- x[s]; wh <- tri(xs / h); wb <- tri(xs / b)
    Xh <- outer(xs, 0:p, `^`); Xb <- outer(xs, 0:q, `^`)
    bm <- coef(lm((xs / h)^(p + 1) ~ Xh - 1, weights = wh, subset = wh > 0))[[1]]
    cu <- coef(lm(ys ~ Xb - 1, weights = wb, subset = wb > 0))[[p + 2]]
    bm * cu
  }
  factorial(p + 1) * (side(x >= 0) - side(x < 0))
}

test_that("B.4 constants general p — rd_period b_const = (p+1)! (D - D_bc) / h^{p+1}", {
  d <- .b4_panel("pc"); d3 <- d[d$year == 3, ]
  for (p in 1:3) {
    f <- rd_period(d3$Y, d3$R, h = 0.35, b = 0.5, id = d3$id, p = p, q = p + 1L)
    expect_equal(f$b_const, factorial(p + 1) * (f$D - f$D_bc) / 0.35^(p + 1), tolerance = 1e-12)
    expect_equal(f$v_const, f$n * 0.35 * f$V_D, tolerance = 1e-12)   # eq:per-period-orders
    # the two factors of the estimated bias on each side: D - D_bc = h^{p+1} bmom curv, and
    # the influence of D - D_bc is h^{p+1} bmom times that of curv
    for (sd in c("+", "-")) {
      S <- f$sides[[sd]]
      # (a difference of two intercepts: at p = 3 the cancellation leaves ~1e-10 relative error)
      expect_equal(S$beta0 - S$beta0_bc, 0.35^(p + 1) * S$bmom * S$curv, tolerance = 1e-8)
      expect_equal(S$g_diff, 0.35^(p + 1) * S$bmom * S$g_curv, tolerance = 1e-10)
    }
  }
})

test_that("B.4 bias design moment — rd_period's bmom and curv, from scratch", {
  # bmom is the intercept of the kernel-weighted order-p regression of ((x - c)/h)^{p+1} on the
  # design (rdrobust's BConst in rdbwselect; checked against rdrobust 2.2's internal
  # rdrobust_bw on 2026-10-27, not tested here because that internal is not a stable API);
  # curv is the coefficient on (x - c)^{p+1} of the order-q fit at b
  d <- .b4_panel("pc"); d3 <- d[d$year == 3, ]
  tri <- function(u) pmax(0, 1 - abs(u))
  for (p in 1:2) {
    q <- p + 1L
    f <- rd_period(d3$Y, d3$R, h = 0.35, b = 0.5, id = d3$id, p = p, q = q)
    for (sd in c("+", "-")) {
      s <- if (sd == "+") d3$R >= 0 else d3$R < 0
      xs <- d3$R[s]; ys <- d3$Y[s]; wh <- tri(xs / 0.35); wb <- tri(xs / 0.5)
      Xh <- outer(xs, 0:p, `^`); Xb <- outer(xs, 0:q, `^`)
      bm <- coef(lm((xs / 0.35)^(p + 1) ~ Xh - 1, weights = wh, subset = wh > 0))[[1]]
      cu <- coef(lm(ys ~ Xb - 1, weights = wb, subset = wb > 0))[[p + 2]]
      expect_equal(f$sides[[sd]]$bmom, bm, tolerance = 1e-10)
      expect_equal(f$sides[[sd]]$curv, cu, tolerance = 1e-8)
    }
  }
})

test_that("B.4 common h (eq:common_h_opt) — .bw_joint returns the closed form in its own constants: curvature at per-period pilots, variance and design moments at h_rot, p = 1, 2", {
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
    # h_rot: rdrobust's rule of thumb on the pooled running variables, N = mean rows per period
    xa <- unlist(lapply(pl, `[[`, "x"))
    h_rot <- 2.576 * min(sd(xa), IQR(xa) / 1.349) * mean(vapply(pl, nrow, 0))^(-1 / 5)
    expect_equal(jb$h_rot, h_rot, tolerance = 1e-12)
    # B = sum_tau coef_tau B-hat_tau (B.4 P3): design moments at h_rot, curvature at each
    # period's own pilot b_tau
    B_hand <- sum(.b4_coef * vapply(names(.b4_coef), function(k)
      .b4_Bhat(pl[[k]]$y, pl[[k]]$x, h = h_rot, b = pb[[k]][["b"]], p = p, q = q), numeric(1)))
    expect_equal(jb$B, B_hand, tolerance = 1e-10)
    # Veff = sum_tau coef_tau^2 V-hat_tau / n_tau with V-hat_tau = n_tau h_rot V_D(h_rot): every
    # period at the SAME h_rot (App. B.4 "Bandwidth estimation"; eq:per-period-orders)
    V_hand <- sum(.b4_coef^2 * vapply(names(.b4_coef), function(k)
      h_rot * rd_period(pl[[k]]$y, pl[[k]]$x, h = h_rot, b = h_rot, id = pl[[k]]$id, p = p, q = q)$V_D,
      numeric(1)))
    expect_equal(jb$Veff, V_hand, tolerance = 1e-12)
    # regularization under CS (no cross-period covariances):
    # lambda * sum coef^2 Var(B-hat_tau), B-hat_tau's side terms (p+1)! bmom(h_rot) curv(b_tau)
    vd <- sum(.b4_coef^2 * vapply(names(.b4_coef), function(k) {
      f <- fits[[k]]
      fr <- rd_period(pl[[k]]$y, pl[[k]]$x, h = h_rot, b = h_rot, id = pl[[k]]$id, p = p, q = q)
      fp1^2 * sum(vapply(c("+", "-"), function(sd)
        fr$sides[[sd]]$bmom^2 * sum(f$sides[[sd]]$g_curv^2), numeric(1))) }, numeric(1)))
    expect_equal(jb$reg, 3 * vd, tolerance = 1e-12)
    if (p == 1L) expect_equal(jb$h, (jb$Veff / (jb$B^2 + jb$reg))^(1 / 5), tolerance = 1e-12)
  }
})

test_that("B.4 the joint selector does not depend on which period carries the t_rd label", {
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
  }
  # through didisc(): relabelling t_rd with the matching comparison weights
  ra <- didisc(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, comparisons = c(2, 1),
              trend = c(-1, 2), bwselect = "joint", scheme = "pc")
  rb <- didisc(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 2, comparisons = c(3, 1),
              trend = c(-1, 2), bwselect = "joint", scheme = "pc")
  cols <- c("est", "se_cs", "se_pc", "se_pv")
  expect_equal(ra$estimates[, cols], rb$estimates[, cols], tolerance = 1e-10)
  expect_equal(ra$bandwidth$h, rb$bandwidth$h, tolerance = 1e-12)
  for (k in c("1", "2", "3")) expect_equal(ra$fits[[k]]$b, rb$fits[[k]]$b, tolerance = 1e-12)
})

test_that("B.4 scheme switch — PV drops the cross-period variance term (lem:cov-pv) but keeps the bias estimates' covariances in the regularization", {
  d <- .b4_panel("pv"); pl <- .b4_plist(d)
  pilot <- c(h = 0.4, b = 0.55)
  pb <- stats::setNames(rep(list(pilot), 3), names(.b4_coef))
  j_cs <- .bw_joint(pl, .b4_coef, "3", scheme = "cs", pilot_bws = pb)
  j_pv <- .bw_joint(pl, .b4_coef, "3", scheme = "pv", pilot_bws = pb)
  expect_equal(j_cs$Veff, j_pv$Veff, tolerance = 1e-12)       # the limit has no cross term
  expect_equal(j_cs$B, j_pv$B, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(j_cs$reg, j_pv$reg, tolerance = 1e-6)))   # covariances of B-hat
  # without regularization the two coincide
  expect_equal(.bw_joint(pl, .b4_coef, "3", scheme = "cs", pilot_bws = pb, regularize = FALSE)$h,
               .bw_joint(pl, .b4_coef, "3", scheme = "pv", pilot_bws = pb, regularize = FALSE)$h,
               tolerance = 1e-12)
  # and PC differs from CS on the same data (the cross term is live)
  j_pc <- .bw_joint(pl, .b4_coef, "3", scheme = "pc", pilot_bws = pb)
  expect_false(isTRUE(all.equal(j_cs$Veff, j_pc$Veff, tolerance = 1e-6)))
})

test_that("B.4 objective (eq:amse-att, lem:agg-var) — the joint rule's objective equals a hand-coded one at a common h, PC, p = 1 and 2", {
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  for (p in 1:2) {
    q <- p + 1L
    pb <- list(`3` = c(h = 0.40, b = 0.55), `1` = c(h = 0.30, b = 0.45), `2` = c(h = 0.36, b = 0.50))
    jb <- .bw_joint(pl, .b4_coef, "3", scheme = "pc", pilot_bws = pb, p = p, q = q)
    keys <- names(.b4_coef); cf <- unname(.b4_coef); hr <- jb$h_rot
    fit <- function(k, h, b) rd_period(pl[[k]]$y, pl[[k]]$x, h = h, b = b, id = pl[[k]]$id, p = p, q = q)
    fitp <- stats::setNames(lapply(keys, function(k) fit(k, pb[[k]][["h"]], pb[[k]][["b"]])), keys)
    fitv <- stats::setNames(lapply(keys, function(k) fit(k, hr, hr)), keys)
    fp1 <- factorial(p + 1)
    h0 <- vapply(keys, function(k) pb[[k]][["h"]], numeric(1))
    bt <- vapply(keys, function(k)
      .b4_Bhat(pl[[k]]$y, pl[[k]]$x, h = hr, b = pb[[k]][["b"]], p = p, q = q), numeric(1))
    # same-side products of two periods' influence vectors, matched by id
    same <- function(fa, fb, el) sum(vapply(c("+", "-"), function(sd) {
      a <- fa$sides[[sd]]; b <- fb$sides[[sd]]; m <- match(a$id, b$id); ok <- !is.na(m)
      sum(a[[el]][ok] * b[[el]][m[ok]]) }, numeric(1)))
    # Cov(B-hat_t, B-hat_s): per side, (p+1)! bmom(h_rot) times the curvature's influence
    # g_curv at each period's pilot
    psi <- lapply(keys, function(k) { f <- fitp[[k]]
      for (sd in c("+", "-")) f$sides[[sd]]$psi <- fp1 * fitv[[k]]$sides[[sd]]$bmom * f$sides[[sd]]$g_curv
      f })
    Vb <- outer(seq_along(keys), seq_along(keys), Vectorize(function(i, j)
      same(psi[[i]], psi[[j]], "psi")))
    # V(D_t(h_rot)) on the diagonal, the same-side covariances at (h_rot, h_rot) off it
    S <- outer(seq_along(keys), seq_along(keys), Vectorize(function(i, j)
      if (i == j) fitv[[i]]$V_D else same(fitv[[i]], fitv[[j]], "g")))
    hand <- function(h) {
      Bbar <- sum(cf * bt) * h^(p + 1) / fp1
      pen  <- 3 * (h^(p + 1) / fp1)^2 * drop(t(cf) %*% Vb %*% cf)
      vv   <- hr * drop(t(cf) %*% S %*% cf) / h        # variance scales as 1/h at a common h
      Bbar^2 + pen + vv
    }
    obj <- function(h) (h^(p + 1) / fp1)^2 * (jb$B^2 + jb$reg) + jb$Veff / h
    for (h in c(0.2, 0.3, 0.45)) expect_equal(obj(h), hand(h), tolerance = 1e-10)
    h_opt <- stats::optimize(hand, interval = c(0.02, 1), tol = 1e-10)$minimum
    expect_equal(jb$h, h_opt, tolerance = 1e-6)
    expect_true(abs(Vb[1, 2]) > 0)                     # the bias estimates' covariance is live
  }
})

test_that("B.4 PC kernel constant (lem:cov-pc) — c(1/rho) = rho c(rho) and c(1) = v; same-side covariance = sum of the two sides", {
  for (p in 1:2) for (kern in c("triangular", "uniform", "epanechnikov")) {
    for (rho in c(0.5, 0.8, 1, 1.6, 2.5)) {
      # c(1/rho) = rho c(rho): the term kappa c(h_t/h_s)/h_s is invariant to the ordering
      expect_equal(.kc_c(p, "+", 1 / rho, kern), rho * .kc_c(p, "+", rho, kern), tolerance = 1e-8)
    }
    expect_equal(.kc_c(p, "+", 1, kern), .kc_v(p, "+", kern), tolerance = 1e-10)
  }
  d <- .b4_panel("pc"); pl <- .b4_plist(d)
  f1 <- rd_period(pl[["1"]]$y, pl[["1"]]$x, h = 0.4, b = 0.55, id = pl[["1"]]$id)
  f2 <- rd_period(pl[["2"]]$y, pl[["2"]]$x, h = 0.4, b = 0.55, id = pl[["2"]]$id)
  cc <- .cross_cov(f1, f2)
  expect_equal(cc$pc, cc$pc_p + cc$pc_m, tolerance = 1e-12)
})

test_that("B.4 uniform p = 0 sanity — c(rho) = min(1, 1/rho) so the PC term is kappa / max(h_t, h_s)", {
  for (rho in c(0.25, 0.5, 1, 2, 4))
    expect_equal(.kc_c(0, "+", rho, "uniform"), min(1, 1 / rho), tolerance = 1e-8)
})
