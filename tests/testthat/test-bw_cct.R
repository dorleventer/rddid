test_that("rd_bw_cct returns finite positive h on standard continuous DGP", {
  skip_if_not_installed("rdrobust")
  set.seed(7)
  x <- runif(2000, -1, 1)
  y <- 0.3 * x + 0.5 * (x >= 0) + rnorm(2000, 0, 0.3)
  bw <- rd_bw_cct(y, x)
  expect_true(is.finite(bw["h"]) && bw["h"] > 0)
  expect_true(is.finite(bw["b"]) && bw["b"] > 0)
})

test_that("rd_bw_cct fallback fires and returns 0.5*IQR when rdbwselect fails", {
  # Force fallback by supplying data entirely on one side of the cutoff.
  # rdbwselect errors internally ("missing value where TRUE/FALSE needed")
  # because there are no right-side observations, so tryCatch catches the error.
  set.seed(13)
  x_onesided <- runif(50, -1, -0.01)   # all strictly below c = 0
  y_onesided  <- runif(50)
  h0_expected <- 0.5 * stats::IQR(x_onesided)
  if (!is.finite(h0_expected) || h0_expected <= 0)
    h0_expected <- stats::sd(x_onesided)
  expect_message(
    bw <- rd_bw_cct(y_onesided, x_onesided),
    regexp = "rd_bw_cct: CCT unavailable"
  )
  expect_equal(unname(bw["h"]), h0_expected, tolerance = 1e-10)
  expect_equal(unname(bw["b"]), h0_expected, tolerance = 1e-10)
})
