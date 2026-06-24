# Regression tests for the internal helpers shared across the package:
# the cutoff-as-treated convention, sampling-scheme detection, and the
# scheme-combined cross-period covariance.

test_that(".build_types treats the cutoff as treated and drops unobserved units", {
  d <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3, 4),   # unit 4 is observed only in period 1
    time = c(1, 2, 1, 2, 1, 2, 1),
    x    = c(0.0, -0.3, 0.5, 0.4, -0.2, 0.1, 0.2)
  )
  bt <- .build_types(d, "x", "time", "id", c = 0)

  # Unit 1 sits exactly at the cutoff in period 1 -> treated ("+").
  expect_equal(bt$wide$side_1[bt$wide$id == 1], "+")

  # Period-2 type = sign pattern of the period-1 side; unit 1 -> "+".
  pt2 <- bt$period_types[["2"]]
  expect_equal(pt2$type[pt2$id == 1], "+")

  # Unit 4 has no period-2 observation, so its type is undefined in BOTH
  # periods' frames and it is dropped.
  expect_false(4 %in% bt$period_types[["1"]]$id)
  expect_false(4 %in% bt$period_types[["2"]]$id)
})

test_that(".scheme_from_long classifies cs / pc / pv", {
  cs <- data.frame(period = c(1, 1, 2, 2), id = c(1, 2, 3, 4), side = c(1, 0, 1, 0))
  expect_equal(.scheme_from_long(cs), "cs")          # no repeated unit

  pc <- data.frame(period = c(1, 2, 1, 2), id = c(1, 1, 2, 2), side = c(1, 1, 0, 0))
  expect_equal(.scheme_from_long(pc), "pc")          # repeated, never switches

  pv <- data.frame(period = c(1, 2, 1, 2), id = c(1, 1, 2, 2), side = c(0, 1, 0, 0))
  expect_equal(.scheme_from_long(pv), "pv")          # unit 1 switches side
})

test_that(".detect_scheme treats x == c as treated (a cutoff unit is not a switch)", {
  # Unit 1: x = 0 in period 1, x = 0.3 in period 2 -> both treated -> no switch.
  plist <- list("1" = list(id = c(1, 2), x = c(0.0, 0.5)),
                "2" = list(id = c(1, 2), x = c(0.3, 0.6)))
  expect_equal(.detect_scheme(plist, c = 0), "pc")   # would be "pv" under sign(x - c)
})

test_that(".cov_scheme reproduces the .cross_cov scheme combinations", {
  set.seed(1)
  n <- 300
  x <- runif(n, -1, 1)
  y <- 0.3 * x + 0.5 * (x >= 0) + rnorm(n, 0, 0.3)
  f <- rd_period(y = y, x = x, h = 0.5, b = 0.5, id = seq_len(n),
                 c = 0, p = 1L, q = 2L)
  cc <- .cross_cov(f, f)
  expect_equal(.cov_scheme(f, f, "cs"), 0)
  expect_equal(.cov_scheme(f, f, "pc"), cc$pc)
  expect_equal(.cov_scheme(f, f, "pv"), cc$pc - cc$pv)
})
