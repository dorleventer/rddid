# Tests for `estimand = c("att", "atu")` (dev/atu_estimand_plan.md S2, S4
# items 4-6). Data defined inline per the task spec; no dependence on
# helper-mirror.R (owned by another agent/task).
#
# On any FAILURE at the stated tolerance, do not loosen the tolerance or skip
# -- report the configuration and the max abs discrepancy.

set.seed(11); n <- 600
R1 <- rnorm(n); R2 <- R1 + 0.3 + rnorm(n, 0, 0.5)     # drift: composition stability fails
mk <- function(t, R) data.frame(id = seq_len(n), t = t, R = R,
                                 Y = R + 1 * (R >= 0) + 0.5 * (t == 2) * (R >= 0) + rnorm(n, 0, 0.5))
d  <- rbind(mk(1, R1), mk(2, R2))
dm <- d; dm$R <- -dm$R                                  # hand-mirrored (cutoff 0)

d3 <- rbind(mk(1, R1), mk(2, R1 + rnorm(n, 0, 0.2)), mk(3, R2))

# ============================================================================
# 1. Default and validation
# ============================================================================

test_that("estimand accepts att/atu, errors on other values, and stores the field", {
  # rddid(): $estimand
  r_att <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv")
  r_atu <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "atu")
  r_def <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv")
  expect_identical(r_att$estimand, "att")
  expect_identical(r_atu$estimand, "atu")
  expect_identical(r_def$estimand, "att")
  expect_error(rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                      comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "x"))

  # rd_typecont(): $meta$estimand
  tc_att <- rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct")
  tc_atu <- rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct", estimand = "atu")
  expect_identical(tc_att$meta$estimand, "att")
  expect_identical(tc_atu$meta$estimand, "atu")
  expect_error(rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct", estimand = "x"))

  # rd_compstable(): $meta$estimand
  cs_att <- rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                           comparisons = 1, bwselect = "cct")
  cs_atu <- rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                           comparisons = 1, bwselect = "cct", estimand = "atu")
  expect_identical(cs_att$meta$estimand, "att")
  expect_identical(cs_atu$meta$estimand, "atu")
  expect_error(rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                              comparisons = 1, bwselect = "cct", estimand = "x"))

  # rd_homog(): $estimand
  h_att <- rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                     comparisons = 1, bwselect = "cct", scheme = "pv")
  h_atu <- rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                     comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "atu")
  expect_identical(h_att$estimand, "att")
  expect_identical(h_atu$estimand, "atu")
  expect_error(rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                         comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "x"))

  # rd_trendcell(): $estimand
  tr_att <- rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                          comparisons = c(1, 2), bwselect = "cct", scheme = "pv")
  tr_atu <- rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                          comparisons = c(1, 2), bwselect = "cct", scheme = "pv", estimand = "atu")
  expect_identical(tr_att$estimand, "att")
  expect_identical(tr_atu$estimand, "atu")
  expect_error(rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                             comparisons = c(1, 2), bwselect = "cct", scheme = "pv", estimand = "x"))
})

# ============================================================================
# 2. Label only: rddid, rd_typecont, rd_homog, rd_trendcell
# ============================================================================

test_that("rddid(): atu is numerically identical to att (label only)", {
  r_att <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv")
  r_atu <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "atu")
  expect_equal(r_att$estimates, r_atu$estimates, tolerance = 1e-12)
})

test_that("rd_typecont(): atu is numerically identical to att (label only)", {
  tc_att <- rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct")
  tc_atu <- rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct", estimand = "atu")
  expect_equal(tc_att$ll_wald$stat, tc_atu$ll_wald$stat, tolerance = 1e-12)
  expect_equal(tc_att$ll_wald$df,   tc_atu$ll_wald$df,   tolerance = 1e-12)
  expect_equal(tc_att$ll_wald$p,    tc_atu$ll_wald$p,    tolerance = 1e-12)
  # per-period components (the jump tables typecont carries)
  expect_identical(names(tc_att$per_period), names(tc_atu$per_period))
  for (pk in names(tc_att$per_period)) {
    expect_equal(tc_att$per_period[[pk]]$ll_wald$stat,
                 tc_atu$per_period[[pk]]$ll_wald$stat, tolerance = 1e-12, info = pk)
    expect_equal(tc_att$per_period[[pk]]$ll_wald$df,
                 tc_atu$per_period[[pk]]$ll_wald$df, tolerance = 1e-12, info = pk)
    expect_equal(tc_att$per_period[[pk]]$ll_wald$p,
                 tc_atu$per_period[[pk]]$ll_wald$p, tolerance = 1e-12, info = pk)
  }
})

test_that("rd_homog(): atu is numerically identical to att (label only)", {
  h_att <- rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                     comparisons = 1, bwselect = "cct", scheme = "pv")
  h_atu <- rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                     comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "atu")
  expect_equal(h_att$statistic, h_atu$statistic, tolerance = 1e-12)
  expect_equal(h_att$df,        h_atu$df,        tolerance = 1e-12)
  expect_equal(h_att$p_value,   h_atu$p_value,   tolerance = 1e-12)
  expect_equal(h_att$period_type_jumps, h_atu$period_type_jumps, tolerance = 1e-12)
})

test_that("rd_trendcell(): atu is numerically identical to att (label only)", {
  tr_att <- rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                          comparisons = c(1, 2), bwselect = "cct", scheme = "pv")
  tr_atu <- rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                          comparisons = c(1, 2), bwselect = "cct", scheme = "pv", estimand = "atu")
  expect_equal(tr_att$statistic, tr_atu$statistic, tolerance = 1e-12)
  expect_equal(tr_att$df,        tr_atu$df,        tolerance = 1e-12)
  expect_equal(tr_att$p_value,   tr_atu$p_value,   tolerance = 1e-12)
  expect_equal(tr_att$cell_period_jumps, tr_atu$cell_period_jumps, tolerance = 1e-12)
})

# ============================================================================
# 3. compstable mirrors
# ============================================================================

test_that("rd_compstable(atu) on d equals rd_compstable(att) on hand-mirrored dm, and differs from att on d", {
  cs_atu_d  <- rd_compstable(d,  x = "R", time = "t", id = "id", t_rd = 2,
                              comparisons = 1, bwselect = "cct", estimand = "atu")
  cs_att_dm <- rd_compstable(dm, x = "R", time = "t", id = "id", t_rd = 2,
                              comparisons = 1, bwselect = "cct", estimand = "att")
  cs_att_d  <- rd_compstable(d,  x = "R", time = "t", id = "id", t_rd = 2,
                              comparisons = 1, bwselect = "cct", estimand = "att")

  pk <- names(cs_atu_d$pairs)[1]
  expect_identical(pk, names(cs_att_dm$pairs)[1])
  pa <- cs_atu_d$pairs[[pk]]
  pb <- cs_att_dm$pairs[[pk]]

  expect_equal(pa$ll_wald$stat, pb$ll_wald$stat, tolerance = 1e-12)
  expect_equal(pa$jumps,        pb$jumps,        tolerance = 1e-12)
  expect_equal(pa$jump_se,      pb$jump_se,      tolerance = 1e-12)
  expect_identical(pa$n_trd, pb$n_trd)
  expect_identical(pa$n_t0,  pb$n_t0)
  expect_identical(pa$n_both, pb$n_both)

  pc <- cs_att_d$pairs[[names(cs_att_d$pairs)[1]]]
  expect_false(isTRUE(all.equal(pa$ll_wald$stat, pc$ll_wald$stat)))
  expect_false(isTRUE(all.equal(pa$n_trd, pc$n_trd)))

  # meta$c: the user's cutoff (default 0)
  expect_identical(cs_atu_d$meta$c, 0)

  # meta$c with a non-zero cutoff, and statistic invariance to the shift
  d_shift <- d
  d_shift$R <- d_shift$R + 0.1
  cs_atu_c0   <- rd_compstable(d,       x = "R", time = "t", id = "id", t_rd = 2,
                                comparisons = 1, bwselect = "cct", estimand = "atu", c = 0)
  cs_atu_c0.1 <- rd_compstable(d_shift, x = "R", time = "t", id = "id", t_rd = 2,
                                comparisons = 1, bwselect = "cct", estimand = "atu", c = 0.1)
  expect_identical(cs_atu_c0.1$meta$c, 0.1)
  pk2 <- names(cs_atu_c0.1$pairs)[1]
  expect_equal(cs_atu_c0.1$pairs[[pk2]]$ll_wald$stat,
               cs_atu_c0$pairs[[names(cs_atu_c0$pairs)[1]]]$ll_wald$stat,
               tolerance = 1e-10)
})

# ============================================================================
# 4. From the equations: below-side reflected fit built by hand
# ============================================================================

test_that("rd_compstable(atu) jump matches a from-the-equations below-side reflected fit", {
  # period-1 units with R1 < 0: x = R1 (negative side), ind = 1{R2 < 0}
  keep1 <- R1 < 0
  x1    <- R1[keep1]
  ind1  <- as.integer(R2[keep1] < 0)
  # period-2 units with R2 < 0: x = -R2 (positive side), ind = 1{R1 < 0}
  keep2 <- R2 < 0
  x2    <- -R2[keep2]
  ind2  <- as.integer(R1[keep2] < 0)

  x   <- c(x1, x2)
  ind <- c(ind1, ind2)

  bw <- rd_bw_cct(ind, x, c = 0)
  f  <- rd_period(ind, x, h = bw[["h"]], b = bw[["b"]], c = 0)

  cs_bc <- rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                          comparisons = 1, bwselect = "cct", estimand = "atu", bc = TRUE)
  cs_conv <- rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                            comparisons = 1, bwselect = "cct", estimand = "atu", bc = FALSE)

  pk_bc   <- names(cs_bc$pairs)[1]
  pk_conv <- names(cs_conv$pairs)[1]

  expect_equal(unname(f$D_bc), unname(cs_bc$pairs[[pk_bc]]$jumps[1]), tolerance = 1e-8)
  expect_equal(unname(f$D),    unname(cs_conv$pairs[[pk_conv]]$jumps[1]), tolerance = 1e-8)
})

# ============================================================================
# 5. Ties at the cutoff
# ============================================================================

test_that("estimand = 'atu' errors on ties at the cutoff; 'att' does not", {
  dt <- d
  dt$R[1] <- 0
  expect_error(
    rd_compstable(dt, x = "R", time = "t", id = "id", t_rd = 2,
                  comparisons = 1, bwselect = "cct", estimand = "atu"),
    regexp = "x == c"
  )
  expect_error(
    rd_compstable(dt, x = "R", time = "t", id = "id", t_rd = 2,
                  comparisons = 1, bwselect = "cct", estimand = "atu"),
    regexp = "between support points"
  )
  expect_error(
    rd_compstable(dt, x = "R", time = "t", id = "id", t_rd = 2,
                  comparisons = 1, bwselect = "cct", estimand = "att"),
    NA
  )
})

# ============================================================================
# 6. Print methods
# ============================================================================

test_that("print methods label ATU designs and stay silent on 'estimand' for att", {
  r_att <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv")
  r_atu <- rddid(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                 comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "atu")
  expect_output(print(r_atu), "ATU\\(t_RD\\)")
  expect_output(print(r_atu), "estimand: ATU")
  expect_false(grepl("estimand", paste(capture.output(print(r_att)), collapse = "\n")))

  cs_att <- rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                           comparisons = 1, bwselect = "cct")
  cs_atu <- rd_compstable(d, x = "R", time = "t", id = "id", t_rd = 2,
                           comparisons = 1, bwselect = "cct", estimand = "atu")
  expect_output(print(cs_atu), "below-cutoff shares")
  expect_false(grepl("estimand", paste(capture.output(print(cs_att)), collapse = "\n")))

  tc_att <- rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct")
  tc_atu <- rd_typecont(d, x = "R", time = "t", id = "id", bwselect = "cct", estimand = "atu")
  print(tc_atu)  # exercise the atu print path
  expect_false(grepl("estimand", paste(capture.output(print(tc_att)), collapse = "\n")))

  h_att <- rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                     comparisons = 1, bwselect = "cct", scheme = "pv")
  h_atu <- rd_homog(d, y = "Y", x = "R", time = "t", id = "id", t_rd = 2,
                     comparisons = 1, bwselect = "cct", scheme = "pv", estimand = "atu")
  print(h_atu)
  expect_false(grepl("estimand", paste(capture.output(print(h_att)), collapse = "\n")))

  tr_att <- rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                          comparisons = c(1, 2), bwselect = "cct", scheme = "pv")
  tr_atu <- rd_trendcell(d3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                          comparisons = c(1, 2), bwselect = "cct", scheme = "pv", estimand = "atu")
  print(tr_atu)
  expect_false(grepl("estimand", paste(capture.output(print(tr_att)), collapse = "\n")))
})
