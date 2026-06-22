make_panel <- function(n = 4000, alpha = 0.3, att = 0.5, curv = 0,
                        varyR = FALSE, seed = 11) {
  set.seed(seed)
  u <- rnorm(n)
  x0 <- runif(n, -1, 1)
  mk <- function(jump, t) {
    x <- if (varyR) x0 + rnorm(n, 0, 0.5) else x0
    data.frame(time = t, id = seq_len(n), x = x,
               y = 0.4 * x + curv * x^2 + ifelse(x >= 0, jump, 0) +
                 0.5 * u + rnorm(n, 0, 0.3))
  }
  rbind(mk(alpha, 1), mk(alpha, 2), mk(att + alpha, 3))
}

test_that("rddid recovers the ATT within sampling error", {
  dat <- make_panel(curv = 0)                 # linear mean -> unbiased
  fit <- rddid(dat, y = "y", x = "x", time = "time", id = "id", t_rd = 3,
               weights = "constant", bwselect = "joint")
  est <- fit$estimates["Robust", "est"]
  se  <- fit$estimates["Robust", "se"]
  expect_lt(abs(est - 0.5), 3 * se)
})

test_that("sampling scheme is auto-detected", {
  # panel, time-constant R -> PC
  d_pc <- make_panel(varyR = FALSE)
  expect_equal(rddid(d_pc, "y", "x", "time", "id", t_rd = 3)$scheme, "pc")
  # panel, time-varying R -> PV
  d_pv <- make_panel(varyR = TRUE)
  expect_equal(rddid(d_pv, "y", "x", "time", "id", t_rd = 3)$scheme, "pv")
  # repeated cross-section (no id) -> CS
  expect_equal(rddid(d_pc, "y", "x", "time", t_rd = 3)$scheme, "cs")
})

test_that("constant and linear weights are admissible and correct", {
  dat <- make_panel()
  fc <- rddid(dat, "y", "x", "time", "id", t_rd = 3, comparisons = c(1, 2),
              weights = "constant")
  expect_equal(unname(fc$weights), c(0.5, 0.5))
  expect_equal(sum(fc$weights), 1)
  fl <- rddid(dat, "y", "x", "time", "id", t_rd = 3, comparisons = c(1, 2),
              weights = "linear")
  expect_equal(unname(fl$weights), c(-1, 2))      # line through t=1,2 to t=3
  expect_equal(sum(fl$weights), 1)
})

test_that("fixed bandwidth is honoured", {
  dat <- make_panel()
  fit <- rddid(dat, "y", "x", "time", "id", t_rd = 3, h = 0.25)
  expect_equal(fit$bandwidth$method, "fixed")
  expect_equal(fit$bandwidth$h, 0.25)
})

test_that("repeated cross-section: PC and PV SEs equal CS", {
  dat <- make_panel()
  fit <- rddid(dat, "y", "x", "time", t_rd = 3)   # no id -> CS
  e <- fit$estimates["Robust", ]
  expect_equal(e$se_pc, e$se_cs, tolerance = 1e-10)
  expect_equal(e$se_pv, e$se_cs, tolerance = 1e-10)
})
