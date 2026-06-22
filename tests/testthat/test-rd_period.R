test_that("rd_period matches rdrobust to machine precision", {
  skip_if_not_installed("rdrobust")
  set.seed(1)
  n <- 4000
  x <- runif(n, -1, 1)
  m <- 0.4 * x + 0.8 * x^2 - 0.3 * x^3 + ifelse(x >= 0, 0.5, 0)
  y <- m + rnorm(n, 0, 0.3)
  h <- 0.25; b <- 0.40

  fit <- rd_period(y, x, h = h, b = b, kernel = "triangular")
  rr <- rdrobust::rdrobust(y, x, c = 0, h = h, b = b, p = 1, q = 2,
                           kernel = "triangular", vce = "hc1")

  expect_equal(fit$D,           unname(rr$coef["Conventional", 1]),   tolerance = 1e-10)
  expect_equal(fit$D_bc,        unname(rr$coef["Bias-Corrected", 1]), tolerance = 1e-10)
  expect_equal(sqrt(fit$V_D),   unname(rr$se["Conventional", 1]),     tolerance = 1e-10)
  expect_equal(sqrt(fit$V_D_bc),unname(rr$se["Robust", 1]),           tolerance = 1e-10)
})

test_that("rd_period recovers a known jump and is reproducible", {
  set.seed(42)
  n <- 8000
  x <- runif(n, -1, 1)
  # linear mean: local-linear RD is unbiased, so D concentrates on the true jump
  y <- 0.3 * x + ifelse(x >= 0, 0.7, 0) + rnorm(n, 0, 0.2)
  fit <- rd_period(y, x, h = 0.25, b = 0.4)
  expect_lt(abs(fit$D - 0.7), 3 * sqrt(fit$V_D))   # within 3 SE of the truth
  expect_true(fit$V_D > 0 && fit$V_D_bc > 0)
})
