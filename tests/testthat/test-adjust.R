# Tests for rd_adjust(): composition-adjusted RD-DID estimator (Theorem 3).
# DGP mirrors code/simulations/s35_comp_adjust_sim.R: P=2, period 1 comparison,
# period 2 RD; period-1 jump depends on partner side v2, period-2 jump on v1.

.adj_draw <- function(n, delta, a1_lo, a1_hi, a2_lo, a2_hi,
                      rho = 0.6, s = 0.5, sigma = 1, tau = 5) {
  eta <- rnorm(n); e1 <- rnorm(n, 0, sqrt(1 - rho^2)); e2 <- rnorm(n, 0, sqrt(1 - rho^2))
  R1 <- rho * eta + e1 + delta; R2 <- rho * eta + e2
  v1 <- as.numeric(R1 >= 0); v2 <- as.numeric(R2 >= 0)
  a1 <- ifelse(v2 == 1, a1_hi, a1_lo); a2 <- ifelse(v1 == 1, a2_hi, a2_lo)
  Y1 <- s * R1 + a1 * v1 + rnorm(n, 0, sigma)
  Y2 <- s * R2 + a2 * v2 + tau * v2 + rnorm(n, 0, sigma)
  id <- seq_len(n)
  rbind(data.frame(id, time = 1L, x = R1, y = Y1),
        data.frame(id, time = 2L, x = R2, y = Y2))
}

test_that("adjusted estimator recovers the ATT when A8 fails but A10 holds", {
  set.seed(101)
  # composition shift (delta=2 -> A8 violated) + stationary confounding (A10 holds)
  d <- .adj_draw(8000, delta = 2, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_adjust(d, y = "y", x = "x", time = "time", id = "id",
                 t_rd = 2, comparisons = 1, h = 0.5, se = "none")
  expect_s3_class(r, "rd_adjust")
  expect_equal(unname(r$att_adj["conv"]), 5, tolerance = 0.4)      # recovers tau
  expect_gt(abs(r$att_unadj["conv"] - 5), 0.8)                     # unadjusted biased
  # within-type and reweighting forms agree (A7 adding-up gap ~ 0)
  expect_lt(abs(r$addingup_gap["conv"]), 0.3)
})

test_that("adjusted nests the unadjusted when composition is stable", {
  set.seed(102)
  d <- .adj_draw(8000, delta = 0, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_adjust(d, y = "y", x = "x", time = "time", id = "id",
                 t_rd = 2, comparisons = 1, h = 0.5, se = "none")
  expect_equal(unname(r$att_adj["conv"]), unname(r$att_unadj["conv"]), tolerance = 0.3)
  expect_equal(unname(r$att_adj["conv"]), 5, tolerance = 0.4)
})

test_that("adjusted estimator fails when the within-cell trend (A10) is violated", {
  set.seed(103)
  # composition shift + NON-stationary confounding (a2 = a1 + 2) -> A10 violated
  d <- .adj_draw(8000, delta = 2, a1_lo = 1, a1_hi = 3, a2_lo = 3, a2_hi = 5)
  r <- rd_adjust(d, y = "y", x = "x", time = "time", id = "id",
                 t_rd = 2, comparisons = 1, h = 0.5, se = "none")
  expect_gt(abs(r$att_adj["conv"] - 5), 1)   # does not recover tau
})

test_that("cluster bootstrap returns positive finite SEs and the documented shape", {
  set.seed(104)
  d <- .adj_draw(3000, delta = 2, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_adjust(d, y = "y", x = "x", time = "time", id = "id",
                 t_rd = 2, comparisons = 1, h = 0.6, se = "bootstrap", B = 100)
  expect_true(is.finite(r$att_adj_se["conv"]) && r$att_adj_se["conv"] > 0)
  expect_true(all(c("label", "pi_plus", "att_by_type", "att_by_type_se") %in%
                    names(r$within_type)))
  expect_equal(nrow(r$within_type), 2L)
  expect_true(all(is.finite(r$within_type$att_by_type_se)))
})
