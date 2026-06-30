# Tests for rd_att(): the correction-form composition-adjusted ATT family.

# small synthetic two-period panel with cross-cutoff switching (so S and C are
# non-degenerate). Running variable varies across periods => units switch sides.
make_panel <- function(n = 400, seed = 1) {
  set.seed(seed)
  x0 <- stats::runif(n, -1, 1)
  x1 <- x0 + stats::rnorm(n, 0, 0.35)          # some switch sides
  mk <- function(xx, t, eff)
    data.frame(id = seq_len(n), time = t, x = xx,
               y = 1 + 0.5 * xx + eff * (xx >= 0) + stats::rnorm(n, 0, 0.5))
  rbind(mk(x0, 0L, 1.0), mk(x1, 1L, 1.8))
}

test_that("rd_att returns the documented structure", {
  d <- make_panel()
  fit <- rd_att(d, "y", "x", "time", "id", comparisons = 0L, t_rd = 1L,
                h = 0.6, se = "none")
  expect_s3_class(fit, "rd_att")
  expect_named(fit, c("att", "per_period", "cterm", "cagg", "t_rd",
                      "comparisons", "weights", "trend", "c", "h", "bwselect",
                      "kernel", "se", "B", "n_boot_ok", "call"))
  expect_setequal(fit$att$estimator, c("unadj", "s", "c", "sc"))
  expect_named(fit$att, c("estimator", "conv", "bc", "se", "bc_se"))
  expect_named(fit$per_period, c("period", "role", "D", "D_bc", "D_se",
                                 "D_bc_se", "S", "S_bc", "S_se", "S_bc_se",
                                 "adj_s_D", "adj_s_D_bc", "adj_s_D_se",
                                 "adj_s_D_bc_se"))
  expect_named(fit$cterm, c("comparison", "C", "C_bc", "C_se", "C_bc_se"))
  expect_equal(nrow(fit$per_period), 2L)
  expect_equal(fit$per_period$adj_s_D, fit$per_period$D - fit$per_period$S)
})

test_that("two-period point ATTs reduce to the hand-assembled primitives", {
  d <- make_panel()
  hh <- 0.6
  fit <- rd_att(d, "y", "x", "time", "id", comparisons = 0L, t_rd = 1L,
                h = hh, trend = "constant", se = "none")

  # hand assembly from the public primitives at the same fixed bandwidth
  d0 <- d[d$time == 0L, ]; d1 <- d[d$time == 1L, ]
  D_t0  <- rd_period(d0$y, d0$x, h = hh, c = 0)$D
  D_trd <- rd_period(d1$y, d1$x, h = hh, c = 0)$D
  sa <- rd_sadjust(d, "y", "x", "time", "id", comparisons = 0L, t_rd = 1L,
                   h = hh, se = "none")
  S_t0  <- sa$table$S[sa$table$role == "comparison"]
  S_trd <- sa$table$S[sa$table$role == "rd"]
  cc <- rd_c(d, "y", "x", "time", "id", t = 0L, s = 1L, h = hh, se = "none")$C

  unadj <- D_trd - D_t0
  s_adj <- (D_trd - S_trd) - (D_t0 - S_t0)
  c_adj <- D_trd - D_t0 - cc
  sc    <- (D_trd - S_trd) - (D_t0 - S_t0) - cc

  pick <- function(e) fit$att$conv[fit$att$estimator == e]
  expect_equal(pick("unadj"), unadj, tolerance = 1e-6)
  expect_equal(pick("s"),     s_adj, tolerance = 1e-6)
  expect_equal(pick("c"),     c_adj, tolerance = 1e-6)
  expect_equal(pick("sc"),    sc,    tolerance = 1e-6)

  # the reduction identities the docstring promises
  expect_equal(pick("c"),  D_trd - D_t0 - cc, tolerance = 1e-6)
  expect_equal(pick("sc"), (D_trd - S_trd) - (D_t0 - S_t0) - cc, tolerance = 1e-6)
})

test_that("ATT_sc point equals a manual assembly of the same primitives", {
  d <- make_panel(seed = 7)
  hh <- 0.55
  fit <- rd_att(d, "y", "x", "time", "id", comparisons = 0L, t_rd = 1L,
                h = hh, se = "none")
  # rebuild from per_period + cterm tables (internal consistency of the object)
  pp <- fit$per_period
  rd <- pp[pp$role == "rd", ]; co <- pp[pp$role == "comparison", ]
  manual_sc <- (rd$D - rd$S) - (co$D - co$S) - fit$cterm$C
  expect_equal(fit$att$conv[fit$att$estimator == "sc"], manual_sc,
               tolerance = 1e-6)
})

test_that("shared-bootstrap SEs are finite and positive", {
  d <- make_panel()
  fit <- rd_att(d, "y", "x", "time", "id", comparisons = 0L, t_rd = 1L,
                h = 0.6, se = "bootstrap", B = 60L)
  expect_true(all(is.finite(fit$att$se)))
  expect_true(all(is.finite(fit$att$bc_se)))
  expect_true(all(fit$att$se > 0))
  expect_true(all(is.finite(fit$per_period$D_se)))
  expect_true(all(is.finite(fit$cterm$C_se)))
  expect_true(is.finite(fit$cagg[["se"]]))
  expect_gt(fit$n_boot_ok, 0L)
})
