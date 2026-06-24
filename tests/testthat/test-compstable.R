# Tests for rd_compstable() — Assumption A7 composition-stability reflection test.
#
# Three objectives:
#  (a) Under the null (composition stable across periods), p-values are not
#      systematically tiny.
#  (b) Under the alternative (composition differs across periods), the test rejects.
#  (c) Structural / output checks + print method.
#
# DGPs are 2-period panels built inline.  Some units are above the cutoff in
# BOTH periods, exercising the cross-side / unit-level handling.

# ============================================================================
# Shared DGP helpers
# ============================================================================

# Null DGP: above-cutoff type mix is the SAME in both periods.
# eta drives both running variables → units cluster together, so many are
# above (or below) in both periods (exercising the n_both path).
#
# Period 1 (comparison t0=1), Period 2 (RD t_rd=2):
#   R_{i,1} = eta_i + eps1_i
#   R_{i,2} = eta_i + eps2_i
# No sorting / manipulation, so type shares (= side of other period) are
# continuous at 0 in both periods.
make_panel_null <- function(n, sigma_eta = 0.6, sigma_eps = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  eta  <- rnorm(n, 0, sigma_eta)
  R1   <- eta + rnorm(n, 0, sigma_eps)
  R2   <- eta + rnorm(n, 0, sigma_eps)
  rbind(
    data.frame(id = seq_len(n), time = 1L, R = R1, stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = 2L, R = R2, stringsAsFactors = FALSE)
  )
}

# Alternative DGP: units that are above the cutoff in period 2 (t_rd=2) tend
# to have a DIFFERENT mix of "period-1 sides" (b) than the above-cutoff units
# in period 1 (t0=1).
# We achieve this by sorting a fraction `p_sort` of period-2 above-cutoff
# units toward the positive side: their period-1 running variable is shifted
# upward, changing their "type" (= 1{R1 >= 0}) systematically.
make_panel_alt <- function(n, p_sort = 0.6, sigma_eta = 0.6, sigma_eps = 0.8,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  eta  <- rnorm(n, 0, sigma_eta)
  R1   <- eta + rnorm(n, 0, sigma_eps)
  R2   <- eta + rnorm(n, 0, sigma_eps)

  # For units above the cutoff in period 2, shift their period-1 running
  # variable upward with probability p_sort (so their period-1 side becomes 1).
  above2 <- R2 >= 0
  movers  <- above2 & (runif(n) < p_sort)
  # Shift: if currently below 0 in period 1, flip to above; if already above,
  # shift further right (doesn't change type but is fine)
  R1[movers & R1 < 0] <- abs(R1[movers & R1 < 0]) + 0.05

  rbind(
    data.frame(id = seq_len(n), time = 1L, R = R1, stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = 2L, R = R2, stringsAsFactors = FALSE)
  )
}

# ============================================================================
# (a) Null scenario: p-values not systematically tiny
# ============================================================================

test_that("rd_compstable null scenario: p-values not tiny", {
  dat <- make_panel_null(n = 3000, seed = 42L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, q = 75L, S = 499L,
                       kernel = "triangular")

  expect_s3_class(out, "rd_compstable")
  pr <- out$pairs[["2::1"]]
  expect_false(is.null(pr))

  # Under the null, do NOT expect systematic rejection.
  # Use a very conservative 0.001 threshold to avoid false CI failures.
  expect_gt(pr$ll_wald$p, 0.001)
  expect_gt(pr$ck_perm$p, 0.001)
})

# ============================================================================
# (b) Alternative scenario: test rejects
# ============================================================================

test_that("rd_compstable alternative: LL-Wald rejects when composition shifts", {
  dat <- make_panel_alt(n = 5000, p_sort = 0.8, seed = 7L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, q = 75L, S = 99L,
                       kernel = "triangular")

  expect_s3_class(out, "rd_compstable")
  pr <- out$pairs[["2::1"]]
  expect_false(is.null(pr))

  # Strong violation: expect LL-Wald rejection at 5% (should be much smaller).
  expect_lt(pr$ll_wald$p, 0.05)
})

# ============================================================================
# (c) Structural / output checks
# ============================================================================

test_that("rd_compstable returns expected structure", {
  dat <- make_panel_null(n = 800, seed = 5L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, S = 49L)

  # Class
  expect_s3_class(out, "rd_compstable")

  # Top-level names
  expect_named(out, c("pairs", "joint", "meta"))

  # Pair result
  pr <- out$pairs[["2::1"]]
  expect_false(is.null(pr))
  expect_named(pr, c("ll_wald", "ck_perm", "type_values",
                     "scheme", "q", "n_trd", "n_t0", "n_both"))

  # ll_wald
  expect_named(pr$ll_wald, c("stat", "df", "p"))
  expect_true(is.numeric(pr$ll_wald$stat))
  expect_true(pr$ll_wald$df >= 0)
  expect_true(pr$ll_wald$p  >= 0 && pr$ll_wald$p <= 1)

  # ck_perm
  expect_named(pr$ck_perm, c("stat", "p"))
  expect_true(is.numeric(pr$ck_perm$stat))
  expect_true(pr$ck_perm$p >= 0 && pr$ck_perm$p <= 1)

  # meta
  expect_equal(out$meta$t_rd, 2L)
  expect_equal(out$meta$comparisons, 1L)
})

test_that("rd_compstable n_both > 0 when eta drives both periods", {
  # With a strong eta (low eps noise), many units are above in both periods.
  dat <- make_panel_null(n = 1000, sigma_eta = 0.9, sigma_eps = 0.3, seed = 9L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, S = 49L)

  pr <- out$pairs[["2::1"]]
  # With strong eta correlation there should be units above in both periods.
  expect_gt(pr$n_both, 0L)
  # And scheme should be "pv" (auto-detected).
  expect_equal(pr$scheme, "pv")
})

test_that("rd_compstable default q follows the Canay-Kamat rule of thumb", {
  dat <- make_panel_null(n = 3000, seed = 42L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, S = 99L)

  # Default q = NULL -> per-pair rule of thumb, flagged in meta.
  expect_identical(out$meta$q, "rot")
  ub <- ceiling(3000^0.9 / log(3000))
  q_used <- unlist(out$meta$q_used)
  expect_true(all(q_used >= 10L & q_used <= ub))
  # Per-pair q is echoed into the pair output and matches meta$q_used.
  expect_equal(out$pairs[["2::1"]]$q, out$meta$q_used[["2::1"]])

  # A supplied q overrides the rule of thumb on every pair.
  out_fix <- rd_compstable(dat, x = "R", time = "time", id = "id",
                           t_rd = 2L, comparisons = 1L,
                           c = 0, h = 0.5, q = 30L, S = 99L)
  expect_true(all(unlist(out_fix$meta$q_used) == 30L))
})

test_that("rd_compstable multiple comparison periods works", {
  # 3-period panel: t_rd=3, comparisons = c(1, 2)
  set.seed(55L)
  n   <- 600
  eta <- rnorm(n)
  dat <- do.call(rbind, lapply(1:3, function(t) {
    data.frame(id = seq_len(n), time = t,
               R  = eta + rnorm(n, 0, 0.8),
               stringsAsFactors = FALSE)
  }))

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 3L, comparisons = c(1L, 2L),
                       c = 0, h = 0.5, S = 49L)

  expect_s3_class(out, "rd_compstable")
  expect_length(out$pairs, 2L)
  expect_true("3::1" %in% names(out$pairs))
  expect_true("3::2" %in% names(out$pairs))

  # Joint result aggregates both pairs
  expect_gte(out$joint$ll_wald$df, out$pairs[["3::1"]]$ll_wald$df)
})

test_that("print.rd_compstable runs without error", {
  dat <- make_panel_null(n = 500, seed = 6L)
  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, S = 49L)
  expect_output(print(out), "Composition-stability")
  expect_output(print(out), "nec & suff")
})

test_that("rd_compstable scheme = 'cs' overrides auto", {
  dat <- make_panel_null(n = 600, seed = 11L)
  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, S = 49L, scheme = "cs")
  pr <- out$pairs[["2::1"]]
  expect_equal(pr$scheme, "cs")
  expect_true(pr$ll_wald$p >= 0 && pr$ll_wald$p <= 1)
})

test_that("rd_compstable: NULL comparisons defaults to all non-t_rd periods", {
  set.seed(22L)
  n   <- 400
  dat <- do.call(rbind, lapply(1:3, function(t) {
    data.frame(id = seq_len(n), time = t, R = rnorm(n),
               stringsAsFactors = FALSE)
  }))
  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 3L, comparisons = NULL,
                       c = 0, h = 0.5, S = 29L)
  # Should have pairs for t0 = 1 and t0 = 2
  expect_length(out$pairs, 2L)
})

# ============================================================================
# (d) High-dual-sided-unit regime: Sigma diagonal cross-side fix + permutation
#     consistency under high cross-period correlation.
#
# DGP: sigma_eta >> sigma_eps so the majority of above-cutoff units are above
# in BOTH periods (high n_both).  This specifically exercises the pv diagonal
# fix (Defect 1) and the side-specific permutation aggregation (Defect 2).
# ============================================================================

# Helper: high-correlation null (many dual-sided units, no composition shift)
make_panel_highcorr_null <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  eta <- rnorm(n, 0, 1.2)         # dominant component → high cross-period correlation
  R1  <- eta + rnorm(n, 0, 0.2)
  R2  <- eta + rnorm(n, 0, 0.2)
  rbind(
    data.frame(id = seq_len(n), time = 1L, R = R1, stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = 2L, R = R2, stringsAsFactors = FALSE)
  )
}

# Helper: high-correlation alternative (composition shift between the two periods
# among above-cutoff units despite high cross-period correlation).
make_panel_highcorr_alt <- function(n, p_sort = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  eta <- rnorm(n, 0, 1.2)
  R1  <- eta + rnorm(n, 0, 0.2)
  R2  <- eta + rnorm(n, 0, 0.2)
  above2  <- R2 >= 0
  movers  <- above2 & (runif(n) < p_sort)
  R1[movers & R1 < 0] <- abs(R1[movers & R1 < 0]) + 0.05
  rbind(
    data.frame(id = seq_len(n), time = 1L, R = R1, stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = 2L, R = R2, stringsAsFactors = FALSE)
  )
}

test_that("rd_compstable high-dual-sided: size controlled under the null", {
  # With sigma_eta = 1.2, sigma_eps = 0.2 the vast majority of above-cutoff
  # units are above in BOTH periods (n_both is large).  Under the null the
  # pv-diagonal fix should keep Sigma well-conditioned and p-values reasonable.
  dat <- make_panel_highcorr_null(n = 2000, seed = 101L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, q = 75L, S = 499L,
                       kernel = "triangular")

  pr <- out$pairs[["2::1"]]
  expect_false(is.null(pr))
  # Many dual-sided units → pv scheme auto-detected.
  expect_equal(pr$scheme, "pv")
  expect_gt(pr$n_both, 50L)      # sanity check that we're in the right regime

  # Size control: neither test should systematically reject under the null.
  expect_gt(pr$ll_wald$p, 0.001)
  expect_gt(pr$ck_perm$p, 0.001)
})

test_that("rd_compstable high-dual-sided: rejects under composition shift", {
  # Strong composition shift + high cross-period correlation.
  dat <- make_panel_highcorr_alt(n = 4000, p_sort = 0.85, seed = 202L)

  out <- rd_compstable(dat, x = "R", time = "time", id = "id",
                       t_rd = 2L, comparisons = 1L,
                       c = 0, h = 0.5, q = 75L, S = 199L,
                       kernel = "triangular")

  pr <- out$pairs[["2::1"]]
  expect_false(is.null(pr))
  # Scheme should be pv (dual-sided units present).
  expect_equal(pr$scheme, "pv")

  # Strong violation: LL-Wald should reject at 5%.
  expect_lt(pr$ll_wald$p, 0.05)
})
