# weighting = "min_variance": the minimum-variance comparison weights (paper: weighting choice and its
# appendix). One test per item of the implementation plan (rd-did docs/rddid_mv_option_plan.md).

# independent solver for checks: the KKT system [2 Xi, L'; L, 0] (w, mu) = (2 xi, l)
kkt_weights <- function(Xi, xi, L, l) {
  m <- nrow(Xi); k <- nrow(L)
  A <- rbind(cbind(2 * Xi, t(L)), cbind(L, matrix(0, k, k)))
  solve(A, c(2 * xi, l))[seq_len(m)]
}
# the plug-in blocks, built from a fit's per-period fits as the standard errors use them
blocks_of <- function(fit) {
  comps <- as.character(fit$comparisons); rd <- as.character(fit$t_rd); m <- length(comps)
  Xi <- matrix(0, m, m); xi <- numeric(m)
  for (i in seq_len(m)) {
    Xi[i, i] <- fit$fits[[comps[i]]]$V_D
    xi[i] <- rddid:::.cov_scheme(fit$fits[[comps[i]]], fit$fits[[rd]], fit$scheme)
    for (j in seq_len(m)[-i]) Xi[i, j] <- rddid:::.cov_scheme(fit$fits[[comps[i]]], fit$fits[[comps[j]]], fit$scheme)
  }
  list(Xi = Xi, xi = xi)
}
args_pv <- list(data = rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
run <- function(...) suppressMessages(do.call(rddid, utils::modifyList(args_pv, list(...))))

test_that("the closed form equals the KKT solution, constant and linear constraints", {
  set.seed(3)
  for (k in 1:2) {
    G <- matrix(stats::rnorm(25), 5); S <- G %*% t(G) + diag(5)
    Xi <- S[2:5, 2:5]; xi <- S[2:5, 1]
    L <- if (k == 1) matrix(1, 1, 4) else rbind(1, c(-4, -3, -2, -1)); l <- if (k == 1) 1 else c(1, 0)
    w <- rddid:::.mv_solve(Xi, xi, L, l)
    expect_equal(w, kkt_weights(Xi, xi, L, l), tolerance = 1e-10)
    expect_lt(max(abs(L %*% w - l)), 1e-10)
  }
})

test_that("repeated cross-section, constant trend: inverse-variance weights", {
  fit <- suppressMessages(rddid(rddid_sim, y = "Y", x = "R", time = "year", t_rd = 3, h = 0.3,
                                weighting = "min_variance"))           # no id: repeated cross-section
  expect_identical(fit$scheme, "cs")
  v <- vapply(c("1", "2"), function(k) fit$fits[[k]]$V_D, numeric(1))
  expect_equal(unname(fit$weights), unname((1 / v) / sum(1 / v)), tolerance = 1e-12)
})

test_that("fixed h: weights, estimate and standard error match a computation by hand", {
  ols <- run(h = 0.3); mv <- run(h = 0.3, weighting = "min_variance")
  b <- blocks_of(ols); w <- kkt_weights(b$Xi, b$xi, matrix(1, 1, 2), 1)
  expect_equal(unname(mv$weights), w, tolerance = 1e-10)
  D <- vapply(c("3", "1", "2"), function(k) ols$fits[[k]]$D, numeric(1))
  expect_equal(mv$estimates["Conventional", "est"], D[["3"]] - sum(w * D[c("1", "2")]), tolerance = 1e-12)
  v <- ols$fits[["3"]]$V_D - 2 * sum(w * b$xi) + c(t(w) %*% b$Xi %*% w)
  expect_equal(mv$estimates["Conventional", "se"]^2, v, tolerance = 1e-10)
  expect_false(mv$weights_detail$pinned)
  expect_equal(unname(mv$weights_detail$pilot), c(0.5, 0.5))
})

test_that("joint: the three steps reproduce by hand", {
  mv <- run(weighting = "min_variance")
  step1 <- run()                                                     # common h with the "ols" weights
  b <- blocks_of(step1); w <- kkt_weights(b$Xi, b$xi, matrix(1, 1, 2), 1)
  step3 <- run(trend = stats::setNames(w, c("1", "2")))              # common h with those weights
  expect_equal(unname(mv$weights), w, tolerance = 1e-10)
  expect_equal(mv$estimates, step3$estimates, tolerance = 1e-10)
  expect_equal(mv$bandwidth$h, step3$bandwidth$h, tolerance = 1e-10)
  expect_equal(mv$weights_detail$h, step1$bandwidth$h_by_period)
})

test_that("cct: weights computed once at the per-period CCT bandwidths", {
  mv <- run(bwselect = "cct", weighting = "min_variance"); ols <- run(bwselect = "cct")
  b <- blocks_of(ols); w <- kkt_weights(b$Xi, b$xi, matrix(1, 1, 2), 1)
  manual <- run(bwselect = "cct", trend = stats::setNames(w, c("1", "2")))
  expect_equal(mv$estimates, manual$estimates, tolerance = 1e-10)
  expect_equal(mv$bandwidth$h_by_period, ols$bandwidth$h_by_period)
})

test_that("pinned weights: as many comparison periods as restrictions gives the ols result", {
  lin_ols <- run(trend = "linear"); lin_mv <- run(trend = "linear", weighting = "min_variance")
  expect_true(lin_mv$weights_detail$pinned)
  expect_identical(lin_mv$estimates, lin_ols$estimates)
  expect_identical(lin_mv$weights, lin_ols$weights)
  one_ols <- run(comparisons = 2); one_mv <- run(comparisons = 2, weighting = "min_variance")
  expect_true(one_mv$weights_detail$pinned)
  expect_identical(one_mv$estimates, one_ols$estimates)
})

test_that("clear errors: numeric trend, iter, singular covariance", {
  expect_error(run(trend = c(0.3, 0.7), weighting = "min_variance"), "fixes them")
  expect_error(run(bwselect = "iter", weighting = "min_variance"), "not with `bwselect = \"iter\"`")
  expect_error(run(trend = "min_variance"), "is a value of `weighting`")
  dup <- rbind(rddid_sim_pv, transform(rddid_sim_pv[rddid_sim_pv$year == 1, ], year = 0))
  expect_error(suppressMessages(rddid(dup, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                                      h = 0.3, weighting = "min_variance")), "not positive")
})

test_that("the default is unchanged", {
  a <- run(); b <- run(weighting = "ols")
  expect_identical(a$estimates, b$estimates)
  expect_identical(a$coef, b$coef)
  expect_identical(a$bandwidth, b$bandwidth)
  expect_identical(a$weighting, "ols")
  expect_null(a$weights_detail)
})

test_that("print and glance report the weighting", {
  mv <- run(weighting = "min_variance")
  expect_output(print(mv), "minimum-variance weights")
  expect_identical(glance(mv)$weighting, "min_variance")
  expect_identical(glance(run())$weighting, "ols")
})

test_that("the ATU label gives the same weights and estimates", {
  att <- run(weighting = "min_variance"); atu <- run(weighting = "min_variance", estimand = "atu")
  expect_identical(att$weights, atu$weights)
  expect_identical(att$estimates, atu$estimates)
})
