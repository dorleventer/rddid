# Tests for rd_c(): cross-period composition term C (Proposition cdisc).
# DGP mirrors test-adjust.R's .adj_draw (P=2, binary partner-side types): period
# 1 is the comparison role, period 2 the RD role; the period-1 jump depends on
# the partner side v2, the period-2 jump on v1.

.c_draw <- function(n, delta, a1_lo, a1_hi, a2_lo, a2_hi,
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

test_that("rd_c returns the documented object structure", {
  set.seed(301)
  d <- .c_draw(4000, delta = 1, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_c(d, y = "y", x = "x", time = "time", id = "id", t = 1, s = 2,
            c = 0, h = 0.5, se = "none")
  expect_s3_class(r, "rd_c")
  expect_true(all(c("C", "C_bc", "C_se", "C_bc_se", "t", "s", "cells",
                    "bwselect", "kernel", "se", "B") %in% names(r)))
  expect_true(all(c("label", "dpi", "dpi_bc", "D", "D_bc", "pi_s_plus",
                    "pi_t_plus", "contrib", "contrib_bc") %in% names(r$cells)))
  expect_equal(nrow(r$cells), 2L)        # binary sides
  expect_identical(c(r$t, r$s), c(1, 2))
  expect_true(is.finite(r$C) && is.finite(r$C_bc))
  # C is the sum of the per-cell contributions
  expect_equal(r$C, sum(r$cells$contrib), tolerance = 1e-10)
})

test_that("C reconciles rd_adjust: D_s - D_t - C == att_adj when A7 adding-up ~0", {
  set.seed(101)
  d <- .c_draw(20000, delta = 2, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  H <- 0.5
  adj <- rd_adjust(d, y = "y", x = "x", time = "time", id = "id",
                   t_rd = 2, comparisons = 1, c = 0, h = H, se = "none")
  Dper <- function(t) { dd <- d[d$time == t, ]
    rd_period(dd$y, dd$x, h = H, b = H, id = dd$id, c = 0)$D }
  cc <- rd_c(d, y = "y", x = "x", time = "time", id = "id", t = 1, s = 2,
             c = 0, h = H, se = "none")
  # the correction form (D_s - D_t - C) reproduces rd_adjust's composition-
  # adjusted ATT; on this clean DGP the adding-up gap is negligible
  expect_lt(abs(adj$addingup_gap["conv"]), 0.05)
  expect_equal(Dper(2) - Dper(1) - cc$C, unname(adj$att_adj["conv"]),
               tolerance = 0.02)
  # exact algebraic identity: att_unadj - C == att_reweight - gap_t0
  cl <- cc$cells
  gap_t0 <- Dper(1) - sum(cl$pi_t_plus * cl$D)
  att_rw <- Dper(2) - sum(cl$pi_s_plus * cl$D)
  expect_equal(Dper(2) - Dper(1) - cc$C, att_rw - gap_t0, tolerance = 1e-8)
})

test_that("C is near zero when composition is stable (delta = 0)", {
  set.seed(102)
  d <- .c_draw(20000, delta = 0, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  cc <- rd_c(d, y = "y", x = "x", time = "time", id = "id", t = 1, s = 2,
             c = 0, h = 0.5, se = "none")
  expect_lt(abs(cc$C), 0.6)             # no composition shift -> C ~ 0
})

test_that("C is not symmetric in (t, s)", {
  set.seed(303)
  d <- .c_draw(8000, delta = 2, a1_lo = 1, a1_hi = 3, a2_lo = 2, a2_hi = 6)
  c_ts <- rd_c(d, y = "y", x = "x", time = "time", id = "id", t = 1, s = 2,
               c = 0, h = 0.5, se = "none")$C
  c_st <- rd_c(d, y = "y", x = "x", time = "time", id = "id", t = 2, s = 1,
               c = 0, h = 0.5, se = "none")$C
  expect_false(isTRUE(all.equal(c_ts, c_st, tolerance = 1e-3)))
})

test_that("cluster bootstrap returns positive finite SEs", {
  set.seed(304)
  d <- .c_draw(3000, delta = 2, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_c(d, y = "y", x = "x", time = "time", id = "id", t = 1, s = 2,
            c = 0, h = 0.6, se = "bootstrap", B = 80)
  expect_true(is.finite(r$C_se) && r$C_se > 0)
  expect_true(is.finite(r$C_bc_se) && r$C_bc_se > 0)
  expect_equal(r$B, 80L)
  expect_true(r$n_boot_ok > 0L)
})

test_that("errors on missing columns and identical periods", {
  set.seed(305)
  d <- .c_draw(1000, delta = 1, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  expect_error(rd_c(d, y = "nope", x = "x", time = "time", id = "id",
                    t = 1, s = 2, h = 0.5), "not found")
  expect_error(rd_c(d, y = "y", x = "x", time = "time", id = "id",
                    t = 1, s = 1, h = 0.5), "different periods")
})
