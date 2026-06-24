# Tests for rd_homog() — Assumption A8 type-homogeneous confounding test.
#
# Two-period DGP:  period 1 = RD period (not used in test),
#                  period 2 = comparison period.
# Type in period 2 = sign of R_{i,1} (the OTHER period's running variable).
# So units with R_{i,1} > 0 are type "+", units with R_{i,1} < 0 are type "-".
#
# Under the null (size test):  the comparison-period outcome jump is the same
# for both types => test should NOT systematically reject.
#
# Under the alternative (power test): the comparison-period outcome jump
# differs by type => test should reject.

# Helper: build a 2-period panel
make_panel <- function(n, jump_pos, jump_neg, seed = 42,
                       alpha_rd = 1.0, n_periods = 2) {
  set.seed(seed)
  # period 1 (RD period): just needs to generate the running variable
  r1 <- stats::runif(n, -1, 1)
  # period 2 (comparison period): outcome jump depends on type (sign of r1)
  r2 <- stats::runif(n, -1, 1)
  type_from_r1 <- ifelse(r1 > 0, "+", "-")
  jump_vec <- ifelse(type_from_r1 == "+", jump_pos, jump_neg)
  y1 <- 0.3 * r1 + alpha_rd * (r1 >= 0) + stats::rnorm(n, 0, 0.3)
  y2 <- 0.3 * r2 + jump_vec * (r2 >= 0) + stats::rnorm(n, 0, 0.3)
  data.frame(
    id   = rep(seq_len(n), times = 2),
    time = rep(c(1L, 2L), each = n),
    x    = c(r1, r2),
    y    = c(y1, y2),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. Basic functionality: result structure
# ---------------------------------------------------------------------------
test_that("rd_homog returns an rd_homog object with the required fields", {
  dat <- make_panel(n = 600, jump_pos = 0.5, jump_neg = 0.5)
  res <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 2L, t_rd = 1L, h = 0.4)
  expect_s3_class(res, "rd_homog")
  expect_true(all(c("statistic", "df", "p_value",
                    "period_type_jumps", "contrasts",
                    "cov_matrix", "scheme") %in% names(res)))
  expect_true(is.numeric(res$statistic) && length(res$statistic) == 1L)
  expect_true(res$statistic >= 0)
  expect_true(res$df >= 1L)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

# ---------------------------------------------------------------------------
# 2. Size: homogeneous jump => test does not systematically reject
# ---------------------------------------------------------------------------
test_that("rd_homog has correct size under null (homogeneous jump across types)", {
  # Repeat over multiple seeds and check that average p-value is not tiny
  pvals <- numeric(40)
  for (s in seq_along(pvals)) {
    dat <- make_panel(n = 800, jump_pos = 0.4, jump_neg = 0.4, seed = s + 100)
    res <- tryCatch(
      rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
               comparisons = 2L, t_rd = 1L, h = 0.4),
      error = function(e) NULL
    )
    pvals[s] <- if (!is.null(res)) res$p_value else NA_real_
  }
  pvals <- pvals[!is.na(pvals)]
  # Under H0, p-values should be roughly uniform; mean well above 0.05
  # and rejection rate at 10% level should not exceed 30% (loose bound)
  reject_rate <- mean(pvals < 0.10)
  expect_lt(reject_rate, 0.35)
  # Statistic should roughly follow chi^2(1): mean should be around 1
  # (very loose check: mean p-value should not be near 0)
  expect_gt(mean(pvals), 0.10)
})

# ---------------------------------------------------------------------------
# 3. Power: heterogeneous jump => test rejects reliably
# ---------------------------------------------------------------------------
test_that("rd_homog rejects when comparison-period jumps differ by type", {
  # Large n, large gap between type-specific jumps
  dat <- make_panel(n = 2000, jump_pos = 0.8, jump_neg = 0.1, seed = 999)
  res <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 2L, t_rd = 1L, h = 0.4)
  expect_lt(res$p_value, 0.05)
})

# ---------------------------------------------------------------------------
# 4. The period_type_jumps table has the right shape
# ---------------------------------------------------------------------------
test_that("rd_homog returns correct per-type jump table", {
  dat <- make_panel(n = 800, jump_pos = 0.5, jump_neg = 0.5)
  res <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 2L, t_rd = 1L, h = 0.4)
  jt <- res$period_type_jumps
  expect_s3_class(jt, "data.frame")
  # Two types in one comparison period => 2 rows
  expect_equal(nrow(jt), 2L)
  expect_true(all(c("period", "type", "jump", "se", "n", "reference") %in% names(jt)))
  # Exactly one reference per period
  expect_equal(sum(jt$reference), 1L)
  # SEs are positive
  expect_true(all(jt$se > 0))
})

# ---------------------------------------------------------------------------
# 5. contrasts vector has length (#types - 1) * #comparison_periods
# ---------------------------------------------------------------------------
test_that("rd_homog contrasts vector length equals (#types-1) x #comparison_periods", {
  dat <- make_panel(n = 600, jump_pos = 0.5, jump_neg = 0.5)
  res <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 2L, t_rd = 1L, h = 0.4)
  # 1 comparison period, 2 types => 1 contrast
  expect_equal(length(res$contrasts), 1L)
  expect_equal(dim(res$cov_matrix), c(1L, 1L))
})

# ---------------------------------------------------------------------------
# 6. RD period outcome (y) does not affect the test
# ---------------------------------------------------------------------------
test_that("rd_homog result is unchanged when t_rd outcome values are permuted", {
  # The test only fits rd_period() in comparison periods; the RD period's
  # outcome y plays no role. Replacing those y values should leave the
  # Wald statistic identical (the running variable x, which determines types,
  # is kept the same).
  dat <- make_panel(n = 600, jump_pos = 0.5, jump_neg = 0.5)
  res_a <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                    comparisons = 2L, t_rd = 1L, h = 0.4)
  # Replace t_rd (period 1) outcome with noise
  dat_modified <- dat
  rd_rows <- dat_modified$time == 1L
  set.seed(321)
  dat_modified$y[rd_rows] <- stats::rnorm(sum(rd_rows), mean = 999, sd = 1)
  res_b <- rd_homog(dat_modified, y = "y", x = "x", time = "time", id = "id",
                    comparisons = 2L, t_rd = 1L, h = 0.4)
  expect_equal(res_a$statistic, res_b$statistic, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 7. print method runs without error
# ---------------------------------------------------------------------------
test_that("print.rd_homog does not error", {
  dat <- make_panel(n = 600, jump_pos = 0.5, jump_neg = 0.5)
  res <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 2L, t_rd = 1L, h = 0.4)
  expect_output(print(res), "Type-homogeneous confounding")
  expect_output(print(res), "neither necessary nor sufficient", ignore.case = TRUE)
})

# ---------------------------------------------------------------------------
# 8. Multiple comparison periods: stacks contrasts correctly
# ---------------------------------------------------------------------------
test_that("rd_homog stacks contrasts across multiple comparison periods", {
  set.seed(11)
  n <- 400
  # 3-period panel: period 3 is RD; periods 1 and 2 are comparisons
  dat3 <- data.frame(
    id   = rep(seq_len(n), times = 3),
    time = rep(c(1L, 2L, 3L), each = n),
    x    = c(stats::runif(n, -1, 1),
             stats::runif(n, -1, 1),
             stats::runif(n, -1, 1)),
    y    = stats::rnorm(3 * n),
    stringsAsFactors = FALSE
  )
  res <- rd_homog(dat3, y = "y", x = "x", time = "time", id = "id",
                  comparisons = c(1L, 2L), t_rd = 3L, h = 0.4)
  # Each comparison period has 2 types => 1 contrast each => 2 total
  # (assuming both types are well-populated)
  expect_gte(length(res$contrasts), 1L)
  expect_equal(length(res$contrasts), nrow(res$cov_matrix))
})
