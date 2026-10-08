# Regression tests for the five behaviour bugs found during the 2026-10-08 UX sweep (T1648).

make_panel <- function(n, n_periods, step_sd, seed, R_lo = -1) {
  set.seed(seed)
  R0 <- runif(n, R_lo, 1)
  do.call(rbind, lapply(seq_len(n_periods), function(t) {
    R <- if (step_sd > 0) R0 + (t - 1) * rnorm(n, 0, step_sd) else R0
    data.frame(id = seq_len(n), year = t, R = R,
               Y = R + 0.5 * (R >= 0) + rnorm(n, 0, 0.3))
  }))
}

test_that("bug 1: the per-cell CCT bandwidth is chosen for the fit's own p", {
  d3 <- rddid_sim_pv[rddid_sim_pv$year == 1, ]
  cb1 <- rddid:::.cell_bandwidth(d3$Y, d3$R, 0, "triangular", NULL, "cct", p = 1L)
  cb2 <- rddid:::.cell_bandwidth(d3$Y, d3$R, 0, "triangular", NULL, "cct", p = 2L)
  expect_equal(cb1, rd_bw_cct(d3$Y, d3$R, c = 0, p = 1L))
  expect_equal(cb2, rd_bw_cct(d3$Y, d3$R, c = 0, p = 2L))
  expect_false(isTRUE(all.equal(cb1, cb2)))
  # end to end: rd_homog(p = 2) now asks rd_bw_cct for p = 2 (and p = 1 by default)
  seen <- new.env(); seen$p <- integer(0)
  tracer <- bquote(assign("p", c(get("p", envir = .(seen)), p), envir = .(seen)))
  suppressMessages(trace("rd_bw_cct", tracer, where = asNamespace("rddid"), print = FALSE))
  on.exit(suppressMessages(untrace("rd_bw_cct", where = asNamespace("rddid"))), add = TRUE)
  invisible(rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                     p = 2L, q = 3L))
  expect_true(length(seen$p) > 0 && all(seen$p == 2L))
  seen$p <- integer(0)
  invisible(rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3))
  expect_true(length(seen$p) > 0 && all(seen$p == 1L))
})

test_that("bug 2: linear within-type contrasts annihilate a linear trend with unequal spacing", {
  # the contrast rows for periods (1, 2, 4): weights (t2 - t1, -(t2 - t0), t1 - t0) x 2/(t2 - t0)
  all_meta <- list(a = list(period = "1"), b = list(period = "2"), c = list(period = "4"))
  cm <- rddid:::.trendcell_contrast_matrix("linear", "+", c("a", "b", "c"), all_meta, 1L, 3L)
  tv <- c(1, 2, 4)
  expect_equal(as.numeric(cm$C %*% (3 + 0.7 * tv)), 0)          # linear jump -> zero contrast
  expect_equal(as.numeric(cm$C), 2 * c(2, -3, 1) / 3)
  # equally spaced periods give exactly (1, -2, 1), so earlier results are unchanged
  all_meta <- list(a = list(period = "1"), b = list(period = "2"), c = list(period = "3"))
  cm <- rddid:::.trendcell_contrast_matrix("linear", "+", c("a", "b", "c"), all_meta, 1L, 3L)
  expect_identical(as.numeric(cm$C), c(1, -2, 1))
  # end to end: with comparison periods 1, 2, 4 the test does not reject a linear trend
  # systematically (one draw: the statistic is finite and the contrast is not the slope)
  d <- make_panel(1500, 5, 0.15, seed = 11)
  d$Y <- d$Y + 0.4 * d$year * (d$R >= 0)     # confounding jump linear in the year
  r <- suppressMessages(rd_trendcell(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 5,
                                     comparisons = c(1, 2, 4), trend = "linear"))
  expect_true(is.finite(r$statistic))
  expect_gt(r$p_value, 0.001)
})

test_that("bug 3: a failed type fit is skipped, not a 'subscript out of bounds' error", {
  d <- make_panel(30, 2, 0.3, seed = 3)
  expect_error(r <- rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 2, h = 0.01), NA)
  expect_s3_class(r, "rd_compstable")
  expect_length(r$pairs, 1L)
})

test_that("bug 4: a constant outcome gives a clear error, not a chi-squared from rounding noise", {
  dc <- transform(rddid_sim_pv, Y = 1)
  expect_error(rd_homog(dc, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, h = 0.4),
               "zero to working precision")
  expect_error(rd_trendcell(dc, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, h = 0.4),
               "zero to working precision")
  # a normal outcome is far from the threshold
  r <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_true(is.finite(r$statistic) && r$df > 0)
})

test_that("bug 5: the reference flag marks the type actually used as reference", {
  # type_by = "pattern": in period 1 a unit's type is (side in 2, side in 3). Period 2 is
  # shifted far above the cutoff, so the types "-+" and "--" are rare there and dropped by
  # min_n; "+-" must then be the reference of period 1 (radix-decreasing order), while period 2
  # (types = sides in 1 and 3, all common) keeps the all-below "--" as its reference.
  set.seed(7); n <- 800
  R0 <- runif(n, -1, 1)
  d <- rbind(data.frame(id = 1:n, year = 1, R = R0),
             data.frame(id = 1:n, year = 2, R = R0 + 0.9),
             data.frame(id = 1:n, year = 3, R = R0 + rnorm(n, 0, 0.3)))
  d$Y <- d$R + 0.5 * (d$R >= 0) + rnorm(nrow(d), 0, 0.3)
  # rare cells trigger rd_bw_cct's fallback warnings before min_n drops them
  r <- suppressWarnings(suppressMessages(
    rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
             type_by = "pattern", min_n = 25)))
  tab <- r$period_type_jumps
  per <- split(tab, tab$period)
  for (tp in names(per)) {
    if (tp %in% r$comparisons) {                 # the period contributed a contrast
      expect_equal(sum(per[[tp]]$reference), 1L)
      first_fitted <- sort(per[[tp]]$type, method = "radix", decreasing = TRUE)[1L]
      expect_equal(per[[tp]]$type[per[[tp]]$reference], first_fitted)
    } else {                                     # a single fitted type: no contrast, no reference
      expect_equal(sum(per[[tp]]$reference), 0L)
    }
  }
  expect_false("--" %in% per[["1"]]$type)                   # the all-below cell was dropped
  expect_equal(per[["1"]]$type[per[["1"]]$reference], "+-") # so the next type stands in
})

# ---- second bug hunt (2026-10-08, estimator surface) ----------------------------------------

test_that("an NA running variable does not read as a side switch", {
  d <- rddid_sim; d$R[1] <- NA
  expect_message(f <- rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3), "dropped")
  expect_equal(f$scheme, "pc")
  f0 <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_equal(f$bandwidth$h, f0$bandwidth$h, tolerance = 1e-3)   # one row cannot move h much
  expect_equal(nobs(f), 2999L)                                     # rows actually used
})

test_that("duplicated unit-period rows and duplicated comparisons are errors", {
  d <- rbind(rddid_sim, rddid_sim[1, ])
  expect_error(rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3), "once per period")
  expect_error(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                     comparisons = c(1, 1)), "more than once")
  dp <- rbind(rddid_sim_pv, rddid_sim_pv[1, ])
  expect_error(rd_typecont(dp, x = "R", time = "year", id = "id"), "once per period")
  expect_error(rd_compstable(dp, x = "R", time = "year", id = "id", t_rd = 3), "once per period")
})

test_that("named numeric trend weights are matched by period name", {
  a <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct",
             comparisons = c(1, 2), trend = c(`2` = 0.3, `1` = 0.7))
  b <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct",
             comparisons = c(1, 2), trend = c(0.7, 0.3))
  expect_equal(a$estimates$est, b$estimates$est)
  expect_error(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                     trend = c(`5` = 0.5, `1` = 0.5)), "names of the numeric")
})

test_that("q reaches the CCT bandwidth selector", {
  d3 <- rddid_sim[rddid_sim$year == 3, ]
  bw <- rd_bw_cct(d3$Y, d3$R, c = 0, p = 1L, q = 3L)
  ref <- rdrobust::rdbwselect(y = d3$Y, x = d3$R, c = 0, p = 1, q = 3, kernel = "triangular",
                              bwselect = "mserd")$bws
  expect_equal(unname(bw[["h"]]), unname(ref[1, 1]))
  expect_equal(unname(bw[["b"]]), unname(ref[1, 3]))
  # the default (q = p + 1) is unchanged
  expect_equal(rd_bw_cct(d3$Y, d3$R), rd_bw_cct(d3$Y, d3$R, q = 2L))
})

test_that("input checks: b without h, kernel case, factor t_rd", {
  expect_error(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, b = 0.4),
               "supply `h` too")
  a <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct",
             kernel = "Uniform")
  expect_equal(a$kernel, "uniform")
  df <- rddid_sim; df$year <- factor(df$year)
  f <- rddid(df, y = "Y", x = "R", time = "year", id = "id", t_rd = factor(3, levels = 1:3),
             bwselect = "cct")
  expect_equal(f$t_rd, "3")
})

# ---- second bug hunt (2026-10-08, tests and plots) --------------------------------------------

test_that("plot() of a test works when t_rd was passed as a variable", {
  skip_if_not_installed("ggplot2")
  rd <- 3
  hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = rd)
  tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = rd)
  expect_s3_class(plot(hg), "ggplot")
  expect_s3_class(plot(tr), "ggplot")
  expect_equal(hg$t_rd, 3)
})

test_that("untestable results are NA, not chi-squared(0) = 0 with p = 1", {
  expect_true(is.na(rddid:::.joint_wald(numeric(0), matrix(numeric(0), 0, 0))$stat))
  expect_true(is.na(rddid:::.typecont_wald(1, matrix(1), integer(0))$p))
  # a share whose jump has a zero standard error is undefined, not a rejection
  r <- rddid:::.joint_wald(c(1, 0.5), diag(c(1e-30, 0.01)))
  expect_true(is.na(r$stat) && r$df == 0L)
})

test_that("comparisons may not include t_rd, and must exist, in every test", {
  expect_error(rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3,
                             comparisons = c(1, 3)), "must not be one of")
  expect_error(rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                        comparisons = c(1, 3)), "must not be one of")
  expect_error(rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                            comparisons = c(1, 9)), "not in `data`")
  expect_error(rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 7),
               "not a period")
})

test_that("rd_trendcell refuses type_by = 'pattern' with a clear message", {
  expect_error(rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                            type_by = "pattern"), "not defined here")
})

test_that("the reference flag is set only in periods that contributed a contrast", {
  hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  tab <- hg$period_type_jumps
  per <- split(tab, tab$period)
  for (tp in names(per))
    expect_equal(sum(per[[tp]]$reference), as.integer(tp %in% hg$comparisons))
})

test_that("plot guards: bad pair index, cross-section switchers", {
  skip_if_not_installed("ggplot2")
  cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
  expect_error(plot(cs, pair = 5), "index 1..2")
  xs <- rddid_sim_pv; xs$id <- seq_len(nrow(xs))
  expect_error(plot_switchers(xs, x = "R", time = "year", id = "id", periods = c(1, 3)),
               "never repeat")
})
