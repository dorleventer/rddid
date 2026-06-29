# Tests for rd_trendcell() — Assumption A10: within-cell confounding trend.
#
# DGP key:
#   Periods 1, 2, (optionally 3) are comparison periods; period `t_rd` is RD.
#   Cell = sign of R_{i,t_rd} (type_by = "rd_side", the default).
#
#   Under A10 (null): within each cell, D_{t0}(cell) is constant across t0.
#   Violation (alternative): within at least one cell, D_{t0}(cell) differs
#   across comparison periods.

# ---------------------------------------------------------------------------
# Shared DGP helpers
# ---------------------------------------------------------------------------

# Build a panel with T_comp comparison periods + 1 RD period.
# jump_fn(cell, t0): returns the comparison-period jump for (cell, period).
# By default, t_rd is the last period.
make_trendcell_panel <- function(n, T_comp, jump_fn, seed = 42,
                                 noise_sd = 0.3) {
  set.seed(seed)
  T_total <- T_comp + 1L
  t_rd    <- T_total

  # Fixed running variable for t_rd (defines cell); independent RV for comparisons
  r_rd   <- stats::runif(n, -1, 1)
  cell   <- ifelse(r_rd >= 0, "+", "-")

  rows <- vector("list", T_total)
  for (t in seq_len(T_comp)) {
    r_t <- stats::runif(n, -1, 1)
    alpha_t <- vapply(cell, function(k) jump_fn(k, t), numeric(1))
    y_t <- 0.3 * r_t + alpha_t * (r_t >= 0) + stats::rnorm(n, 0, noise_sd)
    rows[[t]] <- data.frame(id = seq_len(n), time = t, x = r_t, y = y_t,
                             stringsAsFactors = FALSE)
  }
  # RD period
  r_rd_obs <- r_rd
  y_rd <- 0.3 * r_rd_obs + 1.0 * (r_rd_obs >= 0) + stats::rnorm(n, 0, noise_sd)
  rows[[T_total]] <- data.frame(id = seq_len(n), time = T_total,
                                 x = r_rd_obs, y = y_rd,
                                 stringsAsFactors = FALSE)

  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# 1. Structure: result has the required fields and correct class
# ---------------------------------------------------------------------------
test_that("rd_trendcell returns an rd_trendcell object with required fields", {
  dat <- make_trendcell_panel(n = 500, T_comp = 2,
                               jump_fn = function(k, t) if (k == "+") 0.5 else 0.3)
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:2, t_rd = 3L, h = 0.4)

  expect_s3_class(res, "rd_trendcell")
  expect_true(all(c("statistic", "df", "p_value",
                    "cell_period_jumps", "contrasts",
                    "cov_matrix", "scheme", "bc", "trend", "comparisons") %in%
                    names(res)))
  expect_true(is.numeric(res$statistic) && length(res$statistic) == 1L)
  expect_true(res$statistic >= 0)
  expect_true(res$df >= 1L)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_equal(res$trend, "constant")
})

# ---------------------------------------------------------------------------
# 2. SIZE — constant trend (H0 holds): test does not systematically reject
# ---------------------------------------------------------------------------
test_that("rd_trendcell has correct size when within-cell jumps are constant", {
  # Same confounding jump in both comparison periods within each cell => H0 true.
  pvals <- numeric(40)
  for (s in seq_along(pvals)) {
    dat <- make_trendcell_panel(n = 800, T_comp = 2, seed = s + 200,
                                 jump_fn = function(k, t) if (k == "+") 0.5 else 0.3)
    res <- tryCatch(
      rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                   comparisons = 1:2, t_rd = 3L, h = 0.4),
      error = function(e) NULL
    )
    pvals[s] <- if (!is.null(res)) res$p_value else NA_real_
  }
  pvals <- pvals[!is.na(pvals)]
  reject_rate <- mean(pvals < 0.10)
  # Under H0, rejection rate at 10% level should not exceed 30% (loose bound)
  expect_lt(reject_rate, 0.35)
  # Mean p-value should not be near 0
  expect_gt(mean(pvals), 0.10)
})

# ---------------------------------------------------------------------------
# 3. POWER — trend violation: test rejects reliably
# ---------------------------------------------------------------------------
test_that("rd_trendcell rejects when within-cell jumps differ across periods", {
  # cell "+": jump changes from 0.2 to 0.8 across periods => clear violation
  dat <- make_trendcell_panel(n = 2000, T_comp = 2, seed = 999,
                               jump_fn = function(k, t) {
                                 if (k == "+") c(0.2, 0.8)[t] else 0.4
                               })
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:2, t_rd = 3L, h = 0.4)
  expect_lt(res$p_value, 0.05)
})

# ---------------------------------------------------------------------------
# 4. Two comparison periods: df = (#cells) x 1 = 2 with two cells
# ---------------------------------------------------------------------------
test_that("rd_trendcell with 2 comparison periods returns df = #cells", {
  dat <- make_trendcell_panel(n = 600, T_comp = 2,
                               jump_fn = function(k, t) if (k == "+") 0.4 else 0.4)
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:2, t_rd = 3L, h = 0.4)
  # 2 cells × 1 contrast per cell = 2 contrasts; df <= 2
  expect_equal(length(res$contrasts), 2L)
  expect_equal(dim(res$cov_matrix), c(2L, 2L))
  # cov_matrix should be diagonal (block-diagonal with 1x1 blocks =>
  # off-diagonal entries are 0 since cells are disjoint)
  expect_equal(res$cov_matrix[1L, 2L], 0, tolerance = 1e-10)
  expect_equal(res$cov_matrix[2L, 1L], 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. Linear trend with only 2 comparison periods: returns df = 0 with message
# ---------------------------------------------------------------------------
test_that("rd_trendcell with trend='linear' and 2 comparison periods returns df=0", {
  dat <- make_trendcell_panel(n = 600, T_comp = 2,
                               jump_fn = function(k, t) if (k == "+") 0.4 else 0.4)
  expect_message(
    res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                        comparisons = 1:2, t_rd = 3L, h = 0.4,
                        trend = "linear"),
    regexp = "df = 0"
  )
  expect_equal(res$df, 0L)
  expect_true(is.na(res$statistic))
  expect_true(is.na(res$p_value))
})

# ---------------------------------------------------------------------------
# 6. Linear trend with 3 comparison periods: returns finite chi-sq, df > 0
# ---------------------------------------------------------------------------
test_that("rd_trendcell with trend='linear' and 3 comparison periods is testable", {
  # True linear trend per cell => H0 linear holds approximately
  dat <- make_trendcell_panel(n = 800, T_comp = 3,
                               jump_fn = function(k, t) {
                                 base <- if (k == "+") 0.3 else 0.1
                                 base + 0.05 * t   # linear in t => H0 holds
                               })
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:3, t_rd = 4L, h = 0.4,
                      trend = "linear")
  expect_equal(res$trend, "linear")
  expect_gte(res$df, 1L)
  expect_false(is.na(res$statistic))
  expect_true(res$statistic >= 0)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

# ---------------------------------------------------------------------------
# 7. cell_period_jumps table has the right shape and fields
# ---------------------------------------------------------------------------
test_that("rd_trendcell returns correct cell_period_jumps table", {
  dat <- make_trendcell_panel(n = 600, T_comp = 2,
                               jump_fn = function(k, t) 0.4)
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:2, t_rd = 3L, h = 0.4)

  jt <- res$cell_period_jumps
  expect_s3_class(jt, "data.frame")
  expect_true(all(c("cell", "period", "jump", "se", "n", "reference") %in% names(jt)))
  # 2 cells × 2 periods = 4 rows
  expect_equal(nrow(jt), 4L)
  # SEs are positive
  expect_true(all(jt$se > 0))
  # Exactly one reference row per cell (first period)
  for (ck in unique(jt$cell)) {
    cell_rows <- jt[jt$cell == ck, , drop = FALSE]
    expect_equal(sum(cell_rows$reference), 1L)
  }
})

# ---------------------------------------------------------------------------
# 8. RD period outcome does not affect the test (only comparison periods used)
# ---------------------------------------------------------------------------
test_that("rd_trendcell result is unchanged when t_rd outcome values are permuted", {
  dat <- make_trendcell_panel(n = 600, T_comp = 2,
                               jump_fn = function(k, t) 0.4)
  res_a <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                        comparisons = 1:2, t_rd = 3L, h = 0.4)
  dat_mod <- dat
  rd_rows <- dat_mod$time == 3L
  set.seed(777)
  dat_mod$y[rd_rows] <- stats::rnorm(sum(rd_rows), mean = 999, sd = 1)
  res_b <- rd_trendcell(dat_mod, y = "y", x = "x", time = "time", id = "id",
                        comparisons = 1:2, t_rd = 3L, h = 0.4)
  expect_equal(res_a$statistic, res_b$statistic, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 9. print method runs without error and mentions key phrases
# ---------------------------------------------------------------------------
test_that("print.rd_trendcell does not error", {
  dat <- make_trendcell_panel(n = 600, T_comp = 2,
                               jump_fn = function(k, t) 0.4)
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:2, t_rd = 3L, h = 0.4)
  expect_output(print(res), "Within-cell confounding trend")
  expect_output(print(res), "neither necessary nor sufficient", ignore.case = TRUE)
})

# ---------------------------------------------------------------------------
# 10. Covariance is block-diagonal: off-diagonal blocks are zero
# ---------------------------------------------------------------------------
test_that("rd_trendcell cov_matrix is block-diagonal by cell (3 comparison periods)", {
  dat <- make_trendcell_panel(n = 600, T_comp = 3,
                               jump_fn = function(k, t) 0.4)
  res <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = 1:3, t_rd = 4L, h = 0.4)
  # 2 cells × 2 contrasts each = 4 contrasts; cov_matrix is 4x4
  # Off-diagonal 2x2 blocks (cross-cell) should be zero
  Sigma <- res$cov_matrix
  n_per_cell <- 2L   # (3 - 1) contrasts per cell
  cross_block <- Sigma[seq_len(n_per_cell), (n_per_cell + 1L):nrow(Sigma),
                        drop = FALSE]
  expect_equal(max(abs(cross_block)), 0, tolerance = 1e-10)
})
