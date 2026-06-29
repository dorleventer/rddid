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

test_that("rd_period sides include a finite slope on both sides", {
  set.seed(99)
  n <- 3000
  x <- runif(n, -1, 1)
  y <- 0.4 * x + ifelse(x >= 0, 0.5, 0) + rnorm(n, 0, 0.25)
  fit <- rd_period(y, x, h = 0.3, b = 0.45)
  expect_true(is.finite(fit$sides[["+"]]$slope))
  expect_true(is.finite(fit$sides[["-"]]$slope))
})

test_that("rd_period slope matches independent lm on active window", {
  set.seed(11)
  n <- 4000
  x <- runif(n, -1, 1)
  c0 <- 0; h <- 0.30; b <- 0.30    # b = h so pilot window = main window

  # true DGP: different slope on each side
  true_slope_R <- 0.6
  true_slope_L <- -0.3
  y <- ifelse(x >= c0, true_slope_R * (x - c0), true_slope_L * (x - c0)) +
    ifelse(x >= c0, 0.5, 0) + rnorm(n, 0, 0.2)

  fit <- rd_period(y, x, h = h, b = b, c = c0, kernel = "triangular")

  # --- right side ("+") ---
  keep_R   <- x >= c0
  u_R      <- (x[keep_R] - c0) / h
  w_R      <- pmax(0, 1 - abs(u_R))          # triangular kernel
  active_R <- w_R > 0
  xs_R     <- x[keep_R][active_R]
  ys_R     <- y[keep_R][active_R]
  ws_R     <- w_R[active_R]
  lm_R     <- lm(ys_R ~ I(xs_R - c0), weights = ws_R)
  expect_equal(fit$sides[["+"]]$beta0, unname(coef(lm_R)[1L]), tolerance = 1e-6)
  expect_equal(fit$sides[["+"]]$slope, unname(coef(lm_R)[2L]), tolerance = 1e-6)

  # --- left side ("-") ---
  keep_L   <- x < c0
  u_L      <- (x[keep_L] - c0) / h
  w_L      <- pmax(0, 1 - abs(u_L))
  active_L <- w_L > 0
  xs_L     <- x[keep_L][active_L]
  ys_L     <- y[keep_L][active_L]
  ws_L     <- w_L[active_L]
  lm_L     <- lm(ys_L ~ I(xs_L - c0), weights = ws_L)
  expect_equal(fit$sides[["-"]]$beta0, unname(coef(lm_L)[1L]), tolerance = 1e-6)
  expect_equal(fit$sides[["-"]]$slope, unname(coef(lm_L)[2L]), tolerance = 1e-6)
})
