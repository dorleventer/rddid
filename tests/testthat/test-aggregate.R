mfun <- function(x, jump) 0.5 * x + 0.6 * x^2 + ifelse(x >= 0, jump, 0)

test_that("repeated cross-section: disjoint ids collapse PC and PV to CS", {
  set.seed(7)
  N <- 3000
  mk <- function(jump, idoff) {
    x <- runif(N, -1, 1)
    data.frame(id = idoff + seq_len(N), x = x, y = mfun(x, jump) + rnorm(N, 0, 0.3))
  }
  d0 <- mk(0.5, 0); d1 <- mk(0.2, 1e6); d2 <- mk(0.1, 2e6)
  fits <- list("0" = rd_period(d0$y, d0$x, h = .25, b = .4, id = d0$id),
               "1" = rd_period(d1$y, d1$x, h = .25, b = .4, id = d1$id),
               "2" = rd_period(d2$y, d2$x, h = .25, b = .4, id = d2$id))
  coef <- c("0" = 1, "1" = -0.5, "2" = -0.5)
  a <- rddid:::.aggregate_fits(fits, coef, bc = FALSE)
  expect_equal(unname(a["V_pc"]), unname(a["V_cs"]), tolerance = 1e-12)
  expect_equal(unname(a["V_pv"]), unname(a["V_cs"]), tolerance = 1e-12)
})

test_that("panel time-varying R: PV variance equals per-unit clustered score SS", {
  set.seed(7)
  n <- 2500
  u <- rnorm(n)
  xa <- runif(n, -1, 1); xb <- xa + rnorm(n, 0, .5); xc <- runif(n, -1, 1)
  pdat <- function(x, jump) data.frame(id = seq_len(n), x = x,
                                       y = mfun(x, jump) + 0.5 * u + rnorm(n, 0, .3))
  pf <- list("0" = pdat(xa, 0.5), "1" = pdat(xb, 0.2), "2" = pdat(xc, 0.1))
  fits <- lapply(pf, function(d) rd_period(d$y, d$x, h = .3, b = .45, id = d$id))
  coef <- c("0" = 1, "1" = -0.5, "2" = -0.5)
  a <- rddid:::.aggregate_fits(fits, coef, bc = FALSE)

  # independent: cluster the influence x residual scores by unit across periods
  score <- setNames(numeric(n), as.character(seq_len(n)))
  for (key in names(coef)) {
    f <- fits[[key]]; cf <- coef[[key]]
    sp <- f$sides[["+"]]; sm <- f$sides[["-"]]
    score[as.character(sp$id)] <- score[as.character(sp$id)] + cf * (+1) * sp$g
    score[as.character(sm$id)] <- score[as.character(sm$id)] + cf * (-1) * sm$g
  }
  expect_equal(unname(a["V_pv"]), sum(score^2), tolerance = 1e-12)
})
