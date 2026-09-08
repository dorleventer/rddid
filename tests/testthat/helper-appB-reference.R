# ---------------------------------------------------------------------------
# helper-appB-reference.R
#
# Independent reference implementation of the RD-DID per-period estimator,
# bias correction, and variance formulas of Leventer & Nevo:
#
#   Appendix B preamble : eq:wls_rd  (weighted LS in matrix form)
#   Appendix B.1        : eq:bc_intercept, eq:bc_Q (bias correction)
#   Appendix B.2        : the Gamma/Psi sandwiches, eq:cross-decomp
#   Section 5.2         : eq:var-cs, eq:var-pc, eq:var-pv (aggregate variance)
#
# Written *only* from the equations, in literal matrix form, base R.  It is
# deliberately not a translation of R/rd_period.R or R/aggregate.R (neither
# file was read); it exists so that tests/testthat/test-appB-conformance.R can
# check the package against the paper rather than against itself.
#
# ---------------------------------------------------------------------------
# Assumptions made where the equations leave a choice (all reported):
#
# A1. Diagonal matrices (A_{t,(+/-)}, Sigma_t, Sigma_{t,s}) are stored as their
#     diagonal vectors and applied by elementwise scaling: `t(X) %*% (a * X)`
#     *is* `t(X) %*% A %*% X`.  The omitted products are exact zeros, so the
#     two agree bit-for-bit, not merely to tolerance.
#
# A2. `n` is the size of the unit universe.  Every 1/n in the paper is written
#     out literally.  They cancel in every estimator and sandwich, so a
#     stand-alone period fit (n = n_t) and the same period padded into a panel
#     (n = |union of ids|) return identical numbers; this is asserted in the
#     conformance tests ("padding invariance").
#
# A3. Kernel: K_h(u) = K(u/h)/h, with the standard shapes on [-1, 1]
#     (triangular 1-|u|, Epanechnikov 0.75(1-u^2), uniform 0.5).  The 1/h is a
#     free multiplicative constant of A: `divide_by_h = FALSE` drops it and
#     every returned quantity is unchanged (also asserted in the tests).
#
# A4. Sides: "+" is R >= c and "-" is R < c, per the indicators in A_{t,(+/-)}.
#     A unit exactly at the cutoff is on the "+" side.
#
# A5. For i not in N_t the paper sets Y_{i,t} = 0 with zero weight; we also set
#     R_{i,t} = c so that the design row is well defined.  Both choices are
#     irrelevant: the unit carries zero weight in A_{t,(+/-)}(h_t) and in
#     A_{t,(+/-)}(b_t), its residual is set to 0, and Sigma_{t,s} is supported
#     on N_{t,s} only.
#
# A6. Residuals: eps-hat_{i,t} is the own-side p-order fit at h_t (paper).  For
#     the bias-corrected sandwich the reference also exposes the own-side
#     q-order fit at b_t ("q_b"), the rdrobust convention; which one the paper
#     intends for Psi^BC is a known open item, so both are computed.
#
# A7. e_{v,p} is the (v+1)-th standard basis vector, so the v-th derivative
#     estimate reads coefficient [v + 1] of the fitted vector.
#
# A8. No finite-sample (HC1-type) correction is applied anywhere in this file:
#     it is the literal paper.  The conformance tests apply the package's
#     documented factors explicitly, so the correction is visible in the test
#     source rather than buried in the reference.
# ---------------------------------------------------------------------------


## --- linear algebra --------------------------------------------------------

# Inverse of a symmetric positive-definite Gram matrix after symmetric diagonal
# equilibration: G = D G* D with D = diag(sqrt(diag G)), so G^{-1} = D^{-1}
# G*^{-1} D^{-1}. Exact algebra; it brings the condition number of the
# unscaled polynomial Gram matrix (~1e6 at order 3) down to a few hundred, so
# the reference's numbers do not depend on the BLAS at the 1e-10 level.
ref_solve_equil <- function(G) {
  s <- sqrt(diag(G))
  Gi <- solve(G / outer(s, s))
  Gi / outer(s, s)
}

## --- kernels ---------------------------------------------------------------

# K(u): shape only, zero outside [-1, 1].
ref_kernel <- function(u, kernel = c("triangular", "epanechnikov", "uniform")) {
  kernel <- match.arg(kernel)
  inwin <- abs(u) <= 1
  switch(kernel,
    triangular   = (1 - abs(u)) * inwin,
    epanechnikov = 0.75 * (1 - u^2) * inwin,
    uniform      = 0.5 * inwin)
}

# K_h(u) = K(u / h) / h   (A3: the 1/h is optional, results are invariant).
ref_Kh <- function(u, h, kernel = "triangular", divide_by_h = TRUE) {
  k <- ref_kernel(u / h, kernel)
  if (divide_by_h) k / h else k
}


## --- design matrix ---------------------------------------------------------

# X_{t,p}: row i is X_p(R_{i,t})' = [1, (R-c), ..., (R-c)^p].
ref_Xp <- function(r, c, p) {
  X <- outer(r - c, seq.int(0L, as.integer(p)), "^")
  dimnames(X) <- NULL
  X
}


## --- one period ------------------------------------------------------------

# Fits both sides of one period: conventional (order p at h) and pilot
# (order q at b), the bias-correction pieces of B.1, and the residual vectors
# of B.2.  `inN` marks the units observed in this period (A5).
ref_period_fit <- function(y, r, h, b = h, c = 0, p = 1L, q = 2L,
                           kernel = "triangular",
                           inN = rep(TRUE, length(y)),
                           n = length(y),
                           divide_by_h = TRUE) {
  p <- as.integer(p); q <- as.integer(q)
  stopifnot(q > p, h > 0, b > 0,
            length(y) == length(r), length(inN) == length(y))

  y <- ifelse(inN, y, 0)   # A5
  r <- ifelse(inN, r, c)   # A5

  X_p <- ref_Xp(r, c, p)
  X_q <- ref_Xp(r, c, q)

  Kh <- ref_Kh(r - c, h, kernel, divide_by_h)
  Kb <- ref_Kh(r - c, b, kernel, divide_by_h)

  # 1{i in N_t} 1{R_{i,t} >= c} and 1{i in N_t} 1{R_{i,t} < c}
  ind <- list("+" = as.numeric(inN & r >= c),
              "-" = as.numeric(inN & r <  c))

  z <- (r - c) / h              # z_{i,t}(h_t)
  Z <- z^(p + 1L)               # Z_t(h_t)

  fit_side <- function(sd) {
    a_h <- ind[[sd]] * Kh       # diag of A_{t,(sd)}(h_t)
    a_b <- ind[[sd]] * Kb       # diag of A_{t,(sd)}(b_t)

    # Gamma_{t,(sd),p}(h_t) = (1/n) X_p' A(h) X_p ; same at order q, band b
    G_p <- (1 / n) * (t(X_p) %*% (a_h * X_p))
    G_q <- (1 / n) * (t(X_q) %*% (a_b * X_q))
    Gpi <- ref_solve_equil(G_p)
    Gqi <- ref_solve_equil(G_q)

    # eq:wls_rd : beta-hat = (1/n) Gamma^{-1} X' A Y
    beta_p <- (1 / n) * as.vector(Gpi %*% (t(X_p) %*% (a_h * y)))
    beta_q <- (1 / n) * as.vector(Gqi %*% (t(X_q) %*% (a_b * y)))

    # B.1: vartheta = (1/n) X_p' A(h) Z ; B = e_0' Gamma_p(h)^{-1} vartheta
    theta <- (1 / n) * as.vector(t(X_p) %*% (a_h * Z))
    Bc    <- as.vector(Gpi[1L, ] %*% theta)

    # beta-hat^{(v)} = v! e_{v}' beta-hat   (A7)
    beta0    <- beta_p[1L]
    beta_p1  <- factorial(p + 1L) * beta_q[p + 2L]   # (p+1)-th deriv, order q at b

    # eq:bc_intercept
    beta0_bc <- beta0 - (h^(p + 1L) / factorial(p + 1L)) * beta_p1 * Bc

    # eq:bc_Q :  Q = X_p' A(h) - h^{p+1} vartheta e_{p+1,q}' Gamma_q(b)^{-1} X_q' A(b)
    XtA_p <- t(a_h * X_p)                                   # (p+1) x n
    XtA_q <- t(a_b * X_q)                                   # (q+1) x n
    Q <- XtA_p - h^(p + 1L) *
      outer(theta, as.vector(Gqi[p + 2L, ] %*% XtA_q))      # (p+1) x n
    beta0_bcQ <- as.vector(Gpi[1L, ] %*% ((1 / n) * (Q %*% y)))

    # units on this side inside the *pilot* window (used only by the tests,
    # for the package's finite-sample factor; the reference itself never
    # uses it -- see A8)
    n_pilot <- sum(ind[[sd]] > 0 & Kb > 0)

    list(a_h = a_h, a_b = a_b,
         Gamma_p = G_p, Gamma_q = G_q, Gpi = Gpi, Gqi = Gqi,
         beta_p = beta_p, beta_q = beta_q,
         theta = theta, B = Bc,
         beta0 = beta0, beta_deriv_p1 = beta_p1,
         beta0_bc = beta0_bc, beta0_bc_Q = beta0_bcQ,
         Q = Q, n_pilot = n_pilot)
  }

  sides <- list("+" = fit_side("+"), "-" = fit_side("-"))

  # eq:bc_intercept and eq:bc_Q are algebraically the same estimator. The two
  # forms take different floating-point paths through Gamma_q(b)^{-1}, whose
  # unscaled condition number is ~1e6 at q = 3 (columns (r-c)^k with |r-c| <=
  # b), so the agreement to demand is that of an identity computed on an
  # ill-conditioned matrix (~1e-10 on some BLAS, run-to-run on threaded
  # OpenBLAS), not machine precision. A formula error would be O(1).
  stopifnot(isTRUE(all.equal(sides[["+"]]$beta0_bc, sides[["+"]]$beta0_bc_Q,
                             tolerance = 1e-8)),
            isTRUE(all.equal(sides[["-"]]$beta0_bc, sides[["-"]]$beta0_bc_Q,
                             tolerance = 1e-8)))

  # residuals (A6): own-side fit, evaluated for every i in N_t
  own_fit <- function(X, bp, bm) {
    as.vector(X %*% bp) * (r >= c) + as.vector(X %*% bm) * (r < c)
  }
  eps_p <- (y - own_fit(X_p, sides[["+"]]$beta_p, sides[["-"]]$beta_p)) * inN
  eps_q <- (y - own_fit(X_q, sides[["+"]]$beta_q, sides[["-"]]$beta_q)) * inN

  list(n = n, c = c, p = p, q = q, h = h, b = b, kernel = kernel,
       inN = inN, y = y, r = r, X_p = X_p, X_q = X_q,
       sides = sides,
       eps_p = eps_p, eps_q = eps_q,
       D    = sides[["+"]]$beta0    - sides[["-"]]$beta0,
       D_bc = sides[["+"]]$beta0_bc - sides[["-"]]$beta0_bc)
}


## --- B.2 sandwiches --------------------------------------------------------

# Diagonal of Sigma-hat_{t,s}: eps_{i,t} eps_{i,s} on N_{t,s}, 0 elsewhere.
# `which` picks the residual variant (A6): "p" = p-fit at h, "q" = q-fit at b.
# The N_{t,s} mask is belt-and-braces: ref_period_fit already zeroes eps off
# N_t, and A_t / Q_t are zero there too.  It is kept because the paper states
# the restriction, not because any number depends on it.
ref_sigma <- function(ft, fs = ft, which = c("p", "q")) {
  which <- match.arg(which)
  et <- if (which == "p") ft$eps_p else ft$eps_q
  es <- if (which == "p") fs$eps_p else fs$eps_q
  et * es * (ft$inN & fs$inN)
}

# (0,0) entry of
#   Cov(beta_{t,(st),p}(h_t), beta_{s,(ss),p}(h_s))
#     = (1/n) Gamma_t^{-1} Psi_{t,s} Gamma_s^{-1},
# with Psi_{t,s} = (1/n) X_t' A_t Sigma_{t,s} A_s X_s   (conventional), or
#      Psi^BC    = (1/n) Q_t Sigma_{t,s} Q_s'           (bias-corrected;
#                                                        Gamma unchanged).
ref_cov_beta0 <- function(ft, side_t, fs, side_s, sigma, bc = FALSE) {
  stopifnot(ft$n == fs$n, length(sigma) == ft$n)
  n  <- ft$n
  st <- ft$sides[[side_t]]
  ss <- fs$sides[[side_s]]
  Psi <- if (bc) {
    (1 / n) * (st$Q %*% (sigma * t(ss$Q)))
  } else {
    (1 / n) * (t(ft$X_p) %*% ((st$a_h * sigma * ss$a_h) * fs$X_p))
  }
  M <- (1 / n) * (st$Gpi %*% Psi %*% ss$Gpi)
  M[1L, 1L]
}

# V(D-hat_t) = V(beta0_+) + V(beta0_-)  (disjoint sides).
ref_V_D <- function(ft, bc = FALSE, which = if (bc) "q" else "p") {
  sg <- ref_sigma(ft, ft, which = which)
  ref_cov_beta0(ft, "+", ft, "+", sg, bc) + ref_cov_beta0(ft, "-", ft, "-", sg, bc)
}


## --- the package's finite-sample factor (NOT in the paper) ------------------

# Variance-scale HC1 factor n_side / (n_side - k), k = p+1 (conventional) or
# q+1 (bias-corrected); n_side counts units on that side with positive *pilot*
# kernel weight.  Residuals are scaled by its square root, so a variance term
# built from residuals of (t, side_t) and (s, side_s) is scaled by
# sqrt(f_t * f_s).
ref_hc1 <- function(f, side, bc = FALSE) {
  ns <- f$sides[[side]]$n_pilot
  k  <- if (bc) f$q + 1L else f$p + 1L
  ns / (ns - k)
}


## --- Section 5.2 aggregation ------------------------------------------------

# Cpp[t,s] = Cov(b+_t, b+_s), Cmm = Cov(b-_t, b-_s),
# Cpm[t,s] = Cov(b+_t, b-_s), Cmp[t,s] = Cov(b-_t, b+_s).
# `coefs` are the signed weights (RD period +1, comparison periods -w);
# `rd` names the RD period.  Returns c(est, V_cs, V_pc, V_pv).
ref_aggregate <- function(D, coefs, Cpp, Cmm, Cpm, Cmp, rd) {
  nm <- names(coefs)
  stopifnot(!is.null(nm), all(nm %in% names(D)), rd %in% nm,
            identical(dimnames(Cpp)[[1L]], nm))
  D  <- D[nm]
  Cpp <- Cpp[nm, nm, drop = FALSE]; Cmm <- Cmm[nm, nm, drop = FALSE]
  Cpm <- Cpm[nm, nm, drop = FALSE]; Cmp <- Cmp[nm, nm, drop = FALSE]

  # Cov(b-_t, b+_s) is the transpose of Cov(b+_s, b-_t): self-check.
  stopifnot(isTRUE(all.equal(Cmp, t(Cpm), tolerance = 1e-12,
                             check.attributes = FALSE)))

  Csame <- Cpp + Cmm                  # C^same
  Copp  <- Cpm + Cmp                  # C^opp
  # sides within a period are disjoint => C^opp_{t,t} = 0 exactly
  stopifnot(all(diag(Copp) == 0))

  t0 <- setdiff(nm, rd)
  w  <- -coefs[t0]                    # coef_{t0} = -w_{t0}
  Vd <- diag(Csame)                   # V(D-hat_t)

  # eq:var-cs / eq:var-pc / eq:var-pv, written out as in the paper
  V_cs <- Vd[[rd]] + sum(w^2 * Vd[t0])

  cross_same_00 <- 0
  cross_opp_00  <- 0
  for (a in t0) for (bq in t0) if (a != bq) {
    cross_same_00 <- cross_same_00 + w[[a]] * w[[bq]] * Csame[a, bq]
    cross_opp_00  <- cross_opp_00  + w[[a]] * w[[bq]] * Copp[a, bq]
  }
  V_pc <- V_cs - 2 * sum(w * Csame[rd, t0]) + cross_same_00
  V_pv <- V_pc + 2 * sum(w * Copp[rd, t0]) - cross_opp_00

  # self-check: the same three objects as quadratic forms in the signed
  # weights (V^PV full, V^PC same-side only, V^CS diagonal only)
  cf   <- as.vector(coefs)
  CovD <- Csame - Copp
  q_pv <- as.vector(cf %*% CovD  %*% cf)
  q_pc <- as.vector(cf %*% Csame %*% cf)
  q_cs <- sum(cf^2 * Vd)
  stopifnot(isTRUE(all.equal(V_pv, q_pv, tolerance = 1e-10)),
            isTRUE(all.equal(V_pc, q_pc, tolerance = 1e-10)),
            isTRUE(all.equal(V_cs, q_cs, tolerance = 1e-10)))

  c(est = sum(coefs * D), V_cs = V_cs, V_pc = V_pc, V_pv = V_pv)
}

# Convenience: build the four covariance matrices from a named list of
# ref_period_fit objects (all padded to a common unit universe) and aggregate.
# `hc1 = TRUE` applies the package's finite-sample factor sqrt(f_t f_s) to each
# covariance entry, so the result is comparable to .aggregate_fits().
ref_agg_from_fits <- function(rfits, coefs, rd, bc = FALSE, hc1 = TRUE,
                              bc_resid = c("q_b", "p_h")) {
  bc_resid <- match.arg(bc_resid)
  nm <- names(coefs)
  wh <- if (bc && bc_resid == "q_b") "q" else "p"
  P  <- length(nm)
  mk <- function() matrix(0, P, P, dimnames = list(nm, nm))
  Cpp <- mk(); Cmm <- mk(); Cpm <- mk(); Cmp <- mk()
  for (i in seq_len(P)) for (j in seq_len(P)) {
    ft <- rfits[[nm[i]]]; fs <- rfits[[nm[j]]]
    sg <- ref_sigma(ft, fs, which = wh)
    one <- function(sd_t, sd_s) {
      v <- ref_cov_beta0(ft, sd_t, fs, sd_s, sg, bc = bc)
      if (hc1) v <- v * sqrt(ref_hc1(ft, sd_t, bc) * ref_hc1(fs, sd_s, bc))
      v
    }
    Cpp[i, j] <- one("+", "+"); Cmm[i, j] <- one("-", "-")
    Cpm[i, j] <- one("+", "-"); Cmp[i, j] <- one("-", "+")
  }
  D <- vapply(rfits[nm], function(f) if (bc) f$D_bc else f$D, numeric(1))
  names(D) <- nm
  out <- ref_aggregate(D, coefs, Cpp, Cmm, Cpm, Cmp, rd = rd)
  attr(out, "Cpp") <- Cpp; attr(out, "Cmm") <- Cmm
  attr(out, "Cpm") <- Cpm; attr(out, "Cmp") <- Cmp
  out
}


## --- test fixtures ----------------------------------------------------------

# outcome mean: curvature through order 3 (so the (p+1)-th derivative that the
# bias correction estimates is non-zero for p = 1 and p = 2) plus a jump.
ref_mfun <- function(x, c = 0, jump = 0) {
  u <- x - c
  0.6 * u + 0.9 * u^2 - 0.4 * u^3 + jump * (x >= c)
}

ref_dgp_single <- function(seed, n = 1500, c = 0, jump = 0.5, sd = 0.3) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  data.frame(id = seq_len(n), x = x,
             y = ref_mfun(x, c, jump) + rnorm(n, 0, sd))
}

# Three-period designs.  "cs": fresh units each period.  "pc": same units, same
# running variable every period.  "pv": same units, running variable moves, so
# some units switch side (and period 2 drops `drop_frac` of them, exercising
# N_{t,s}).  Panel designs carry a unit-level error component, so Sigma_{t,s}
# is non-zero.
ref_dgp_panel <- function(seed, n = 1500, c = 0,
                          design = c("cs", "pc", "pv"),
                          jumps = c(0.5, 0.2, 0.35),
                          sd = 0.3, sd_u = 0.4, drop_frac = 0) {
  design <- match.arg(design)
  set.seed(seed)
  P <- length(jumps)
  u <- rnorm(n, 0, sd_u)
  x0 <- runif(n, -1, 1)
  out <- vector("list", P)
  for (t in seq_len(P)) {
    if (design == "cs") {
      xt  <- runif(n, -1, 1)
      idt <- (t - 1L) * 1000000L + seq_len(n)
      yt  <- ref_mfun(xt, c, jumps[t]) + rnorm(n, 0, sd)
    } else {
      idt <- seq_len(n)
      xt  <- if (design == "pc") x0 else x0 + rnorm(n, 0, 0.35)
      yt  <- ref_mfun(xt, c, jumps[t]) + u + rnorm(n, 0, sd)
    }
    d <- data.frame(id = idt, x = xt, y = yt)
    if (drop_frac > 0 && t == 2L) {
      keep <- sample.int(nrow(d), size = round((1 - drop_frac) * nrow(d)))
      d <- d[sort(keep), , drop = FALSE]
    }
    out[[t]] <- d
  }
  names(out) <- as.character(seq_len(P))
  out
}

# Pad one period's data onto a common unit universe (A2/A5).
ref_pad <- function(dat, ids, c = 0) {
  m   <- match(ids, dat$id)
  inN <- !is.na(m)
  y <- rep(0, length(ids)); r <- rep(c, length(ids))
  y[inN] <- dat$y[m[inN]]
  r[inN] <- dat$x[m[inN]]
  list(y = y, r = r, inN = inN)
}

# Reference fits for a whole panel, all on the union-of-ids universe.
ref_panel_fits <- function(dlist, h, b, c = 0, p = 1L, q = 2L,
                           kernel = "triangular") {
  ids <- sort(unique(unlist(lapply(dlist, function(d) d$id))))
  nn  <- length(ids)
  nm  <- names(dlist)
  out <- vector("list", length(nm)); names(out) <- nm
  for (i in seq_along(nm)) {
    pd <- ref_pad(dlist[[i]], ids, c)
    out[[i]] <- ref_period_fit(pd$y, pd$r, h = h[i], b = b[i], c = c,
                               p = p, q = q, kernel = kernel,
                               inN = pd$inN, n = nn)
  }
  out
}

# Package fits for the same panel.
ref_pkg_fits <- function(dlist, h, b, c = 0, p = 1L, q = 2L,
                         kernel = "triangular") {
  nm  <- names(dlist)
  out <- vector("list", length(nm)); names(out) <- nm
  for (i in seq_along(nm)) {
    d <- dlist[[i]]
    out[[i]] <- rd_period(d$y, d$x, h = h[i], b = b[i], id = d$id, c = c,
                          p = as.integer(p), q = as.integer(q), kernel = kernel)
  }
  out
}

# Cov(beta0_{t,(side_t)}, beta0_{s,(side_s)}) implied by the package's returned
# influence vectors: the sum of g_t * g_s over shared unit ids.
ref_pkg_gcov <- function(ft, side_t, fs, side_s, bc = FALSE) {
  st <- ft$sides[[side_t]]; ss <- fs$sides[[side_s]]
  gt <- if (bc) st$g_bc else st$g
  gs <- if (bc) ss$g_bc else ss$g
  m  <- match(st$id, ss$id)
  ok <- !is.na(m)
  if (!any(ok)) return(0)
  sum(gt[ok] * gs[m[ok]])
}
