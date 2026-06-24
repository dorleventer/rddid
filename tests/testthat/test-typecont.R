# Tests for rd_typecont() / test_helpers.R
# Three objectives:
#  (a) package functions reproduce standalone ll_cov_p / ck_perm_p / mccrary_p
#      p-values on the same simulated cross section (within tolerance)
#  (b) under the null scenario p-values are not systematically tiny
#  (c) under a sorting scenario the LL-Wald rejects

# ---- shared DGP (mirrors dgp_s3.R, self-contained) -------------------------
dgp_s3_local <- function(n, scenario = "null", p = 0.5, w = 0.5, sigma_R = 1,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  probs <- switch(scenario,
    null           = c(pu1 = 0,   pu0 = 0, pd1 = 0, pd0 = 0),
    a7_pooled      = c(pu1 = p,   pu0 = 0, pd1 = 0, pd0 = 0),
    within_notnec  = c(pu1 = p,   pu0 = p, pd1 = 0, pd0 = 0),
    pooled_notsuff = c(pu1 = p,   pu0 = 0, pd1 = 0, pd0 = p),
    stop("unknown scenario: ", scenario))

  eta <- rnorm(n)
  R1  <- eta + rnorm(n, 0, sigma_R)
  R2  <- eta + rnorm(n, 0, sigma_R)
  V1  <- as.integer(R1 >= 0)

  below <- R2 > -w & R2 < 0
  above <- R2 < w & R2 > 0
  pmove <- numeric(n)
  pmove[below] <- ifelse(V1[below] == 1, probs["pu1"], probs["pu0"])
  pmove[above] <- ifelse(V1[above] == 1, probs["pd1"], probs["pd0"])
  mover <- runif(n) < pmove
  R2[mover] <- -R2[mover]

  data.frame(id = seq_len(n), R = R2, type = V1)
}

# Convert single-period cross section (with type) to 2-period long panel
# consistent with rd_typecont() expectations.
xsec_to_panel <- function(df) {
  # Period 2 is the RD period: R = df$R, side = 1{R >= 0}
  # Period 1 is the comparison: R1 such that type (= V1) is the side in period 1.
  # Here we set R1 = type (0 or 1) shifted so it's above/below 0.
  # Actually, the dgp already produced R1 separately; for the cross-section test
  # we only have R2 and type. We RECONSTRUCT a minimal consistent panel:
  #   period 1 running variable: if type==1, set R1 = +0.1; if type==0, R1 = -0.1
  #   so 1{R1 >= 0} == type exactly, and the period-2 running variable is R2.
  n <- nrow(df)
  R1_fake <- ifelse(df$type == 1L, 0.1, -0.1)
  rbind(
    data.frame(id = df$id, time = 1L, R = R1_fake, type_check = df$type),
    data.frame(id = df$id, time = 2L, R = df$R,    type_check = df$type)
  )
}

# ---- standalone reference implementations (from s3_tests_labtest.R) ---------
ll_cov_p_ref <- function(x, g, h) {
  s  <- as.numeric(x >= 0)
  w  <- pmax(0, 1 - abs(x) / h)
  inw <- w > 0
  Z  <- cbind(1, s, x, s * x)
  Zw <- Z * w
  XtWXi <- tryCatch(solve(crossprod(Z, Zw)), error = function(e) NULL)
  if (is.null(XtWXi)) return(NA_real_)
  beta <- XtWXi %*% crossprod(Zw, g)
  tau  <- beta[2L]
  u    <- as.vector(g - Z %*% beta)
  a    <- as.vector(Zw %*% XtWXi[2L, ])
  infl <- a * u
  infl[!inw] <- 0
  se <- sqrt(sum(infl^2))
  2 * pnorm(-abs(tau / se))
}

ck_perm_p_ref <- function(x, g, q, S = 499L, seed = 42L) {
  set.seed(seed)
  rt   <- order(x[x >= 0])[seq_len(min(q, sum(x >= 0)))]
  lt   <- order(-x[x < 0])[seq_len(min(q, sum(x < 0)))]
  gr   <- g[x >= 0][rt]
  gl   <- g[x < 0][lt]
  pool <- c(gr, gl)
  nr   <- length(gr)
  obs  <- abs(mean(gr) - mean(gl))
  perm <- replicate(S, {
    idx <- sample.int(length(pool), nr)
    abs(mean(pool[idx]) - mean(pool[-idx]))
  })
  (1 + sum(perm >= obs)) / (S + 1)
}

mccrary_p_ref <- function(x, h) {
  n <- length(x)
  if (n < 30L) return(NA_real_)
  binw <- 2 * sd(x) * n^(-1 / 2)
  lo   <- floor(min(x) / binw) * binw
  hi   <- ceiling(max(x) / binw) * binw
  brks <- seq(lo, hi, by = binw)
  mids <- brks[-length(brks)] + binw / 2
  cnt  <- tabulate(findInterval(x, brks, rightmost.closed = TRUE),
                   nbins = length(mids))
  Y    <- cnt / (n * binw)
  fit_side <- function(keep) {
    g <- mids[keep]; y <- Y[keep]; w <- pmax(0, 1 - abs(g) / h); use <- w > 0
    g <- g[use]; y <- y[use]; w <- w[use]
    if (length(g) < 3L) return(NULL)
    Z <- cbind(1, g); Zw <- Z * w
    XtWXi <- tryCatch(solve(crossprod(Z, Zw)), error = function(e) NULL)
    if (is.null(XtWXi)) return(NULL)
    beta <- XtWXi %*% crossprod(Zw, y)
    vj   <- pmax(y, 1e-8) / (n * binw)
    Vb   <- XtWXi %*% crossprod(Zw, vj * Zw) %*% XtWXi
    c(f0 = beta[1L], var = Vb[1L, 1L])
  }
  Rr <- fit_side(mids > 0); Ll <- fit_side(mids < 0)
  if (is.null(Rr) || is.null(Ll) || Rr["f0"] <= 0 || Ll["f0"] <= 0)
    return(NA_real_)
  theta <- log(Rr["f0"]) - log(Ll["f0"])
  se    <- sqrt(Rr["var"] / Rr["f0"]^2 + Ll["var"] / Ll["f0"]^2)
  unname(2 * pnorm(-abs(theta / se)))
}

# ============================================================================
# (a) Reproduction: package matches standalone on the same cross section
# ============================================================================

test_that(".mccrary matches mccrary_p_ref on same inputs", {
  set.seed(7)
  n <- 3000
  x <- rnorm(n, 0, 0.8)   # centred at 0
  h <- 0.5

  expect_equal(.mccrary(x, h), mccrary_p_ref(x, h), tolerance = 1e-12)
})

test_that(".ck_perm matches ck_perm_p_ref on same inputs", {
  set.seed(99)
  n <- 2000
  x <- rnorm(n, 0, 0.8)
  g <- as.integer(rnorm(n) > 0)   # random binary covariate
  q <- 75L; S <- 199L
  # Both use the same seed
  ref_p <- ck_perm_p_ref(x, g, q, S, seed = 123L)

  set.seed(123L)
  pkg_p <- .ck_perm(x, g, q, S)

  expect_equal(pkg_p, ref_p, tolerance = 1e-12)
})

test_that("rd_typecont LL jump matches ll_cov_p_ref on the same data (HC0 vs HC1 noted)", {
  # ll_cov_p_ref uses raw residuals (HC0); rd_period uses HC1 df correction.
  # The ESTIMATE (theta) must match; the SE and hence p-value will differ
  # slightly. We test the jump itself.
  set.seed(2025)
  n <- 2000; h <- 0.5
  d  <- dgp_s3_local(n, "null", seed = 2025)
  panel <- xsec_to_panel(d)

  # Package estimate of the jump for type==1 in period 2
  fit_pkg <- rd_period(
    y  = as.numeric(panel$type_check[panel$time == 2L] == 1L),
    x  = panel$R[panel$time == 2L],
    h  = h, b = h,
    id = panel$id[panel$time == 2L],
    c  = 0, p = 1L, q = 2L, kernel = "triangular"
  )
  # Standalone: uses the same running variable and type indicator
  ref_p <- ll_cov_p_ref(d$R, as.numeric(d$type == 1L), h)

  # Jump estimates must agree to machine precision
  # (both run the same LL regression; ref uses column 2 of the same design matrix)
  # The ref beta[2] corresponds to the side indicator = rd_period $D
  # Verify they match:
  expect_equal(fit_pkg$D, local({
    x_cs <- d$R; g_cs <- as.numeric(d$type == 1L)
    s  <- as.numeric(x_cs >= 0); w  <- pmax(0, 1 - abs(x_cs) / h)
    Z  <- cbind(1, s, x_cs, s * x_cs); Zw <- Z * w
    XtWXi <- solve(crossprod(Z, Zw))
    beta <- XtWXi %*% crossprod(Zw, g_cs)
    as.numeric(beta[2L])
  }), tolerance = 1e-10)
})

test_that("rd_typecont McCrary pooled matches mccrary_p_ref on same cross section", {
  set.seed(11)
  n <- 3000; h <- 0.5
  d     <- dgp_s3_local(n, "null", seed = 11)
  panel <- xsec_to_panel(d)

  out <- rd_typecont(panel, x = "R", time = "time", id = "id",
                     c = 0, h = h, q = 75L, S = 99L, kernel = "triangular")

  # Pooled McCrary for period 2 should match the standalone on d$R
  pkg_p <- out$mccrary_pooled$p[out$mccrary_pooled$period == "2"]
  ref_p <- mccrary_p_ref(d$R, h)

  # They use the same centred x (d$R is already at cutoff 0)
  expect_equal(pkg_p, ref_p, tolerance = 1e-10)
})

# ============================================================================
# (b) Under the null, p-values are not systematically tiny
# ============================================================================

test_that("rd_typecont null scenario p-values are not all < 0.05", {
  set.seed(42)
  n <- 3000
  d     <- dgp_s3_local(n, "null", seed = 42)
  panel <- xsec_to_panel(d)

  out <- rd_typecont(panel, x = "R", time = "time", id = "id",
                     c = 0, h = 0.5, q = 75L, S = 499L, kernel = "triangular")

  # Under the null, we do NOT expect rejection.  Use 0.001 as a very
  # conservative threshold so this test does not generate false failures.
  expect_gt(out$ll_wald$p,  0.001)
  expect_gt(out$ck_perm$p,  0.001)
  # McCrary pooled p-values per period
  expect_true(all(out$mccrary_pooled$p > 0.001, na.rm = TRUE))
})

# ============================================================================
# (c) Under a sorting scenario the LL-Wald rejects
# ============================================================================

test_that("rd_typecont rejects LL-Wald under a7_pooled (type-1 mass moves up)", {
  set.seed(123)
  n <- 5000
  d     <- dgp_s3_local(n, "a7_pooled", p = 0.7, seed = 123)
  panel <- xsec_to_panel(d)

  out <- rd_typecont(panel, x = "R", time = "time", id = "id",
                     c = 0, h = 0.5, q = 75L, S = 99L, kernel = "triangular")

  # Strong violation: expect rejection at 1% level
  expect_lt(out$ll_wald$p, 0.01)
  expect_lt(out$ck_perm$p, 0.05)
})

# ============================================================================
# Structural / output checks
# ============================================================================

test_that("rd_typecont returns an rd_typecont object with expected fields", {
  set.seed(5)
  n <- 1000
  d     <- dgp_s3_local(n, "null", seed = 5)
  panel <- xsec_to_panel(d)

  out <- rd_typecont(panel, x = "R", time = "time", id = "id",
                     c = 0, h = 0.5, S = 49L)

  expect_s3_class(out, "rd_typecont")
  expect_named(out, c("ll_wald", "ck_perm", "mccrary_within", "mccrary_pooled", "meta"))
  expect_named(out$ll_wald, c("stat", "df", "p"))
  expect_named(out$ck_perm, c("stat", "p"))
  expect_true(is.data.frame(out$mccrary_within))
  expect_true(is.data.frame(out$mccrary_pooled))
  # p-values in [0,1]
  expect_true(out$ll_wald$p >= 0 && out$ll_wald$p <= 1)
  expect_true(out$ck_perm$p >= 0 && out$ck_perm$p <= 1)
  expect_true(all(out$mccrary_pooled$p >= 0 & out$mccrary_pooled$p <= 1,
                  na.rm = TRUE))
})

test_that("print.rd_typecont runs without error", {
  set.seed(6)
  n <- 500
  d     <- dgp_s3_local(n, "null", seed = 6)
  panel <- xsec_to_panel(d)
  out   <- rd_typecont(panel, x = "R", time = "time", id = "id",
                       c = 0, h = 0.5, S = 49L)
  expect_output(print(out), "LL-Wald")
})

test_that("rd_typecont works with 3 periods", {
  # 3-period panel: type in each period = 2-bit integer (2^2 = 4 possible types)
  set.seed(77)
  n <- 800
  eta <- rnorm(n)
  panel3 <- do.call(rbind, lapply(1:3, function(t) {
    data.frame(id = seq_len(n), time = t,
               R  = eta + rnorm(n, 0, 1))
  }))

  out <- rd_typecont(panel3, x = "R", time = "time", id = "id",
                     c = 0, h = 0.5, S = 49L)
  expect_s3_class(out, "rd_typecont")
  expect_equal(length(out$meta$periods), 3L)
  # With P=3 periods, n_types = 2^(3-1) = 4 possible type values
  expect_lte(length(out$meta$type_values), 4L)  # at most 4
  expect_true(out$ll_wald$df >= 0)
})

# ============================================================================
# P>=3 CK permutation: correctness under no-sorting and sorting DGPs
# ============================================================================

# Helper: 3-period no-sorting panel (running variables independent of each
# other; no unit is nudged across the cutoff based on type).
make_panel3_null <- function(n, seed) {
  set.seed(seed)
  eta <- rnorm(n)
  do.call(rbind, lapply(1:3, function(t) {
    data.frame(id = seq_len(n), time = t,
               R  = eta + rnorm(n, 0, 1))
  }))
}

# Helper: 3-period sorting panel. In each period, units with all-other-periods
# above the cutoff (type == 3 in binary encoding) are nudged from just below to
# just above the cutoff with probability p_sort.
make_panel3_sort <- function(n, seed, p_sort = 0.6) {
  set.seed(seed)
  eta <- rnorm(n)
  R   <- matrix(eta + matrix(rnorm(n * 3, 0, 1), n, 3), n, 3)
  # In each period, check if all OTHER periods are above cutoff (type = 3)
  for (t in 1:3) {
    other <- setdiff(1:3, t)
    all_above <- R[, other[1]] >= 0 & R[, other[2]] >= 0
    # Nudge: move from (-0.3, 0) to (0, 0.3) for units just below cutoff
    near_below <- R[, t] >= -0.3 & R[, t] < 0
    movers <- all_above & near_below & (runif(n) < p_sort)
    R[movers, t] <- -R[movers, t]  # reflect to just above 0
  }
  do.call(rbind, lapply(1:3, function(t) {
    data.frame(id = seq_len(n), time = t, R = R[, t])
  }))
}

test_that("CK permutation (P=3) does not over-reject under no-sorting DGP", {
  # Under the null, the Canay-Kamat p-value should not be systematically tiny.
  # We run 5 independent draws and require at least 3 of them to exceed 0.01,
  # which gives <1% probability of failure if the true size is <=5%.
  seeds  <- c(201L, 202L, 203L, 204L, 205L)
  p_vals <- vapply(seeds, function(s) {
    panel3 <- make_panel3_null(1500L, s)
    out    <- rd_typecont(panel3, x = "R", time = "time", id = "id",
                          c = 0, h = 0.5, q = 75L, S = 399L)
    out$ck_perm$p
  }, numeric(1L))
  # At least 3 out of 5 must not reject at 0.01 level
  expect_gte(sum(p_vals > 0.01), 3L)
})

test_that("CK permutation (P=3) rejects under sorting DGP", {
  # Under strong type-sorting, p-value should be small.
  panel3 <- make_panel3_sort(3000L, seed = 301L, p_sort = 0.7)
  out    <- rd_typecont(panel3, x = "R", time = "time", id = "id",
                        c = 0, h = 0.5, q = 75L, S = 499L)
  expect_lt(out$ck_perm$p, 0.05)
})

test_that("CK permutation (P=3): observed stat is deterministic and p in [0,1]", {
  # The observed statistic is computed from the data (no randomness); it must
  # be identical across two calls.  p-value must be a valid probability.
  # This also indirectly verifies the shared-draw structure: the observed stat
  # is the sum of |mean_right - mean_left| across (period, type) cells, which
  # is the same formula used for permuted stats — so the test is internally
  # consistent only if both sides use the same grouping.
  panel3 <- make_panel3_null(800L, seed = 401L)

  out1 <- rd_typecont(panel3, x = "R", time = "time", id = "id",
                      c = 0, h = 0.5, q = 50L, S = 99L)
  out2 <- rd_typecont(panel3, x = "R", time = "time", id = "id",
                      c = 0, h = 0.5, q = 50L, S = 99L)

  # Observed stat is purely deterministic — must match exactly
  expect_equal(out1$ck_perm$stat, out2$ck_perm$stat, tolerance = 1e-12)
  # p-value is a valid probability
  expect_gte(out1$ck_perm$p, 0)
  expect_lte(out1$ck_perm$p, 1)
})
