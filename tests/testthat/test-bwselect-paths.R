# Tests strengthening CCT and rot bandwidth path coverage.
#
# (1) bwselect = "cct" vs "rot" produce DIFFERENT statistics (confirms the
#     bandwidth selection path actually engages, not just a code path exercise).
# (2) Per-cell rot branch (bwselect = "rot", h = NULL) for rd_homog and
#     rd_trendcell — previously untested; every other test uses explicit h
#     or bwselect = "cct".

# ---------------------------------------------------------------------------
# Shared minimal DGP: 3-period panel, period 3 = RD, periods 1 and 2 = comparisons.
# Running variable is i.i.d. Uniform(-1, 1) so CCT and rot bandwidths differ.
# ---------------------------------------------------------------------------
make_bwtest_panel <- function(n, seed = 1L) {
  set.seed(seed)
  T_total <- 3L
  t_rd    <- 3L
  r_rd    <- stats::runif(n, -1, 1)
  cell    <- ifelse(r_rd >= 0, "+", "-")
  alpha   <- ifelse(cell == "+", 0.4, 0.3)

  rows <- vector("list", T_total)
  for (t in seq_len(2L)) {
    r_t <- stats::runif(n, -1, 1)
    y_t <- 0.3 * r_t + alpha * (r_t >= 0) + stats::rnorm(n, 0, 0.3)
    rows[[t]] <- data.frame(id = seq_len(n), time = t, x = r_t, y = y_t,
                             stringsAsFactors = FALSE)
  }
  y_rd <- 0.3 * r_rd + 1.0 * (r_rd >= 0) + stats::rnorm(n, 0, 0.3)
  rows[[3L]] <- data.frame(id = seq_len(n), time = 3L, x = r_rd, y = y_rd,
                             stringsAsFactors = FALSE)
  do.call(rbind, rows)
}

# 2-period panel for rd_typecont
make_typecont_panel <- function(n, seed = 1L) {
  set.seed(seed)
  eta <- stats::rnorm(n)
  R1  <- eta + stats::rnorm(n)
  R2  <- eta + stats::rnorm(n)
  rbind(
    data.frame(id = seq_len(n), time = 1L, R = R1, stringsAsFactors = FALSE),
    data.frame(id = seq_len(n), time = 2L, R = R2, stringsAsFactors = FALSE)
  )
}

# ---------------------------------------------------------------------------
# (1) bwselect = "cct" and "rot" produce different statistics
# ---------------------------------------------------------------------------
test_that("bwselect='cct' and 'rot' produce different statistics (rd_typecont)", {
  skip_if_not_installed("rdrobust")
  dat <- make_typecont_panel(n = 800L, seed = 7L)

  out_cct <- rd_typecont(dat, x = "R", time = "time", id = "id",
                         c = 0, h = NULL, bwselect = "cct",
                         q = 30L, S = 49L, kernel = "triangular")
  out_rot <- rd_typecont(dat, x = "R", time = "time", id = "id",
                         c = 0, h = NULL, bwselect = "rot",
                         q = 30L, S = 49L, kernel = "triangular")

  # Different bandwidth paths should yield different LL-Wald statistics.
  expect_false(isTRUE(all.equal(out_cct$ll_wald$stat, out_rot$ll_wald$stat,
                                tolerance = 1e-6)))
  # Both must be finite and produce valid p-values.
  expect_true(is.finite(out_cct$ll_wald$stat))
  expect_true(is.finite(out_rot$ll_wald$stat))
  expect_true(out_cct$ll_wald$p >= 0 && out_cct$ll_wald$p <= 1)
  expect_true(out_rot$ll_wald$p  >= 0 && out_rot$ll_wald$p  <= 1)
  # bwselect is stored correctly in meta.
  expect_equal(out_cct$meta$bwselect, "cct")
  expect_equal(out_rot$meta$bwselect, "rot")
  # h is NA for cct (per-cell), numeric for rot (full-sample IQR).
  expect_true(is.na(out_cct$meta$h))
  expect_true(is.finite(out_rot$meta$h))
})

test_that("bwselect='cct' and 'rot' produce different statistics (rd_homog)", {
  skip_if_not_installed("rdrobust")
  dat <- make_bwtest_panel(n = 800L, seed = 11L)

  out_cct <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = c(1L, 2L), t_rd = 3L,
                      c = 0, h = NULL, bwselect = "cct",
                      type_by = "rd_side", bc = FALSE)
  out_rot <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = c(1L, 2L), t_rd = 3L,
                      c = 0, h = NULL, bwselect = "rot",
                      type_by = "rd_side", bc = FALSE)

  # Different bandwidth paths should yield different Wald statistics.
  expect_false(isTRUE(all.equal(out_cct$statistic, out_rot$statistic,
                                tolerance = 1e-6)))
  # Both must return valid results.
  expect_true(is.finite(out_cct$statistic))
  expect_true(is.finite(out_rot$statistic))
  expect_true(out_cct$p_value >= 0 && out_cct$p_value <= 1)
  expect_true(out_rot$p_value >= 0 && out_rot$p_value <= 1)
})

# ---------------------------------------------------------------------------
# (2) Per-cell rot branch for rd_homog and rd_trendcell
# ---------------------------------------------------------------------------
test_that("rd_homog bwselect='rot' h=NULL (per-cell 0.2*range) runs end-to-end", {
  dat <- make_bwtest_panel(n = 600L, seed = 3L)

  out <- rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
                  comparisons = c(1L, 2L), t_rd = 3L,
                  c = 0, h = NULL, bwselect = "rot",
                  type_by = "rd_side", bc = FALSE)

  expect_s3_class(out, "rd_homog")
  expect_true(is.finite(out$statistic))
  expect_gte(out$df, 1L)
  expect_true(out$p_value >= 0 && out$p_value <= 1)
  # Jump table must be populated (rot still fits rd_period per cell)
  expect_gt(nrow(out$period_type_jumps), 0L)
})

test_that("rd_trendcell bwselect='rot' h=NULL (per-cell 0.2*range) runs end-to-end", {
  dat <- make_bwtest_panel(n = 600L, seed = 5L)

  out <- rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
                      comparisons = c(1L, 2L), t_rd = 3L,
                      c = 0, h = NULL, bwselect = "rot",
                      type_by = "rd_side", trend = "constant", bc = FALSE)

  expect_s3_class(out, "rd_trendcell")
  expect_true(is.finite(out$statistic))
  expect_gte(out$df, 1L)
  expect_true(out$p_value >= 0 && out$p_value <= 1)
  expect_gt(nrow(out$cell_period_jumps), 0L)
})
