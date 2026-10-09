# Fixes from the 2026-10-09 verification of the UX sweep (T1653): one test per fix.

test_that("rd_period takes one main and one pilot bandwidth, not a vector", {
  set.seed(1)
  x <- stats::runif(400, -1, 1)
  y <- 0.5 * x + (x >= 0) + stats::rnorm(400, sd = 0.2)
  expect_error(rd_period(y, x, h = c(0.4, 0.5)), "`h` must be a single positive number")
  expect_error(rd_period(y, x, h = 0.4, b = c(0.5, 0.6, 0.7)),
               "`b` must be a single positive number")
  expect_error(rd_period(y, x, h = 0.4, b = NA_real_), "`b` must be a single positive number")
  expect_error(rd_period(y, x, h = -1), "`h` must be a single positive number")
  # the internal callers pass named length-one vectors (bw["h"]): still accepted, same fit
  fit_named <- rd_period(y, x, h = c(h = 0.4), b = c(b = 0.6))
  fit_plain <- rd_period(y, x, h = 0.4, b = 0.6)
  expect_equal(fit_named$D, fit_plain$D)
  expect_equal(fit_named$D_bc, fit_plain$D_bc)
  expect_equal(fit_named$V_D_bc, fit_plain$V_D_bc)
})

test_that("a share with no variation and no jump near the cutoff is left out of the Wald test", {
  # the no-information share (jump 0, se 0) is dropped exactly: same test as without it
  full <- rddid:::.joint_wald(c(0, 0.5), diag(c(0, 0.01)))
  rest <- rddid:::.joint_wald(0.5, matrix(0.01))
  expect_equal(full$stat, rest$stat)
  expect_equal(full$p, rest$p)
  expect_identical(full$df, 1L)
  expect_identical(full$dropped, 1L)
  expect_null(rest$dropped)
  # its zero covariances leave with it
  S  <- matrix(c(0, 0,     0,
                 0, 0.02,  0.005,
                 0, 0.005, 0.03), 3, 3)
  r3 <- rddid:::.joint_wald(c(0, 0.3, -0.2), S)
  r2 <- rddid:::.joint_wald(c(0.3, -0.2), S[2:3, 2:3])
  expect_equal(r3$stat, r2$stat)
  expect_identical(r3$df, 2L)
  # every share without information: nothing testable
  r0 <- rddid:::.joint_wald(c(0, 0), diag(c(0, 0)))
  expect_true(is.na(r0$stat))
  expect_identical(r0$df, 0L)
  expect_identical(r0$dropped, 2L)
  # a deterministic jump (se 0, jump not 0) still makes the whole test undefined
  rdet <- rddid:::.joint_wald(c(1, 0, 0.5), diag(c(0, 0, 0.01)))
  expect_true(is.na(rdet$stat))
  expect_true(isTRUE(rdet$degenerate))
})
