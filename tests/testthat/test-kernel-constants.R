# Tests for the internal kernel-constants module (R/kernel_constants.R).
# Ports the self-test of the validated source module
# (rd-did/code/simulations/appb_verify/R/kernel_constants.R) as testthat
# expectations. See R/kernel_constants.R for the objects under test and their
# citations (Appendix B.3 of the paper).

# Compare by MAX ABSOLUTE difference, exactly mirroring the source module's
# own .kc_check(). testthat::expect_equal()'s `tolerance` is magnitude-scaled
# (effectively relative once |value| > 1, e.g. for v_+ ~ 10), which would
# silently loosen the tolerances ported below, so they are checked directly.
expect_kc_close <- function(got, want, tol, label = NULL) {
  got  <- as.numeric(got)
  want <- as.numeric(want)
  if (length(want) == 1L && length(got) > 1L) want <- rep(want, length(got))
  expect_equal(length(got), length(want), info = label)
  expect_true(max(abs(got - want)) <= tol, info = label)
}

.kc_kernels <- c("triangular", "epanechnikov", "uniform")
.kc_ps      <- 1:2
.kc_sides   <- c("+", "-")
.kc_rhos    <- c(0.5, 0.6, 1, 1.4, 2, 7 / 3)

test_that("Simpson agrees with Gauss-Legendre: Gamma~, Psi~, vartheta~, b, v (all kernels)", {
  for (kernel in .kc_kernels) {
    for (p in .kc_ps) {
      for (side in .kc_sides) {
        lbl <- sprintf("kernel=%s p=%d side=%s", kernel, p, side)
        expect_kc_close(.kc_gamma(p, side, kernel, rule = "simpson"),
                        .kc_gamma(p, side, kernel, rule = "gl"),
                        1e-8, paste("Gamma~", lbl))
        expect_kc_close(.kc_psi(p, side, kernel, rule = "simpson"),
                        .kc_psi(p, side, kernel, rule = "gl"),
                        1e-8, paste("Psi~", lbl))
        expect_kc_close(.kc_theta(p, side, kernel = kernel, rule = "simpson"),
                        .kc_theta(p, side, kernel = kernel, rule = "gl"),
                        1e-8, paste("vartheta~", lbl))
        expect_kc_close(.kc_b(p, side, kernel, rule = "simpson"),
                        .kc_b(p, side, kernel, rule = "gl"),
                        1e-8, paste("b", lbl))
        expect_kc_close(.kc_v(p, side, kernel, rule = "simpson"),
                        .kc_v(p, side, kernel, rule = "gl"),
                        1e-8, paste("v", lbl))
      }
    }
  }
})

test_that("Simpson agrees with Gauss-Legendre: Omega~(rho), c(rho) (all kernels)", {
  for (kernel in .kc_kernels) {
    for (p in .kc_ps) {
      for (side in .kc_sides) {
        for (rho in .kc_rhos) {
          lbl <- sprintf("kernel=%s p=%d side=%s rho=%.4f", kernel, p, side, rho)
          expect_kc_close(.kc_omega(p, side, rho, kernel, rule = "simpson"),
                          .kc_omega(p, side, rho, kernel, rule = "gl"),
                          1e-8, paste("Omega~", lbl))
          expect_kc_close(.kc_c(p, side, rho, kernel, rule = "simpson"),
                          .kc_c(p, side, rho, kernel, rule = "gl"),
                          1e-8, paste("c", lbl))
        }
      }
    }
  }
})

test_that("Omega~(1) == Psi~; c(1) == v; c(1/rho) == rho * c(rho) (all kernels)", {
  for (kernel in .kc_kernels) {
    for (p in .kc_ps) {
      for (side in .kc_sides) {
        lbl <- sprintf("kernel=%s p=%d side=%s", kernel, p, side)
        expect_kc_close(.kc_omega(p, side, 1, kernel), .kc_psi(p, side, kernel),
                        1e-10, paste("Omega~(1)==Psi~", lbl))
        expect_kc_close(.kc_c(p, side, 1, kernel), .kc_v(p, side, kernel),
                        1e-10, paste("c(1)==v", lbl))
        for (rho in .kc_rhos) {
          expect_kc_close(.kc_c(p, side, 1 / rho, kernel),
                          rho * .kc_c(p, side, rho, kernel),
                          1e-8, sprintf("c(1/rho)==rho*c(rho) %s rho=%.4f", lbl, rho))
        }
      }
    }
  }
})

test_that("kernel constants are symmetric across sides (all kernels)", {
  for (kernel in .kc_kernels) {
    for (p in .kc_ps) {
      lbl <- sprintf("kernel=%s p=%d", kernel, p)
      expect_kc_close(.kc_b(p, "-", kernel), (-1)^(p + 1) * .kc_b(p, "+", kernel),
                      1e-10, paste("b_- == (-1)^(p+1) b_+", lbl))
      expect_kc_close(.kc_v(p, "-", kernel), .kc_v(p, "+", kernel),
                      1e-10, paste("v_- == v_+", lbl))
    }
  }
})

test_that("closed-form kernel constants, triangular kernel", {
  expect_kc_close(.kc_b(1, "+", "triangular"), -0.1, 1e-6)
  expect_kc_close(.kc_v(1, "+", "triangular"), 4.8, 1e-6)
  expect_kc_close(.kc_c(1, "+", 0.5, "triangular"), 5.7, 1e-6)
  expect_kc_close(.kc_c(1, "+", 2, "triangular"), 2.85, 1e-6)
  expect_kc_close(.kc_b(2, "+", "triangular"), 0.0285714285714286, 1e-6)
  expect_kc_close(.kc_v(2, "+", "triangular"), 10.2857142857143, 1e-4)
})

test_that("closed-form kernel constants, uniform kernel, p = 0", {
  # K = 0.5 on [-1, 1]; Omega~_+(rho) = int_0^min(1,1/rho) 0.25 dv, so
  # Gamma~_+ = 0.5, Psi~_+ = 0.25, v_+ = 1, c_+(rho) = min(1, 1/rho).
  expect_kc_close(.kc_gamma(0, "+", "uniform"), 0.5, 1e-8)
  expect_kc_close(.kc_psi(0, "+", "uniform"), 0.25, 1e-8)
  expect_kc_close(.kc_v(0, "+", "uniform"), 1, 1e-8)
  for (rho in c(0.5, 1, 2, 3)) {
    expect_kc_close(.kc_c(0, "+", rho, "uniform"), min(1, 1 / rho), 1e-8,
                    sprintf("c_+(rho=%.4f)", rho))
  }
})

test_that("memoisation caches results and repeated calls are identical", {
  rm(list = ls(.kc_cache), envir = .kc_cache)
  expect_equal(length(ls(.kc_cache)), 0L)

  g1 <- .kc_gamma(2, "+", "epanechnikov")
  expect_gt(length(ls(.kc_cache)), 0L)            # first call populated the cache
  n_after_gamma <- length(ls(.kc_cache))
  g2 <- .kc_gamma(2, "+", "epanechnikov")          # same key -> cache hit
  expect_identical(g1, g2)
  expect_equal(length(ls(.kc_cache)), n_after_gamma)

  o1 <- .kc_omega(2, "-", 1.4, "uniform")
  o2 <- .kc_omega(2, "-", 1.4, "uniform")
  expect_identical(o1, o2)

  # a different rho is a genuinely different key -> a new cache entry
  o3 <- .kc_omega(2, "-", 2.5, "uniform")
  expect_gt(length(ls(.kc_cache)), n_after_gamma + 1L)
})
