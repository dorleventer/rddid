# Size corrections of the validation tests (T1678, from the T1652 size study): sup-t / sWald joint
# statistics, the full-covariance joint test of rd_compstable(), and the common bandwidth of
# rd_homog() / rd_trendcell().

test_that(".gauss_joint: one jump is the exact two-sided normal / t test", {
  g <- .gauss_joint(1.7, matrix(0.25, 1, 1), nu = 40)
  z <- 1.7 / 0.5
  expect_equal(g$supt_normal$stat, z)
  expect_equal(g$swald_normal$stat, z^2)
  expect_equal(g$supt_normal$p, 2 * pnorm(-z))
  expect_equal(g$swald_normal$p, pchisq(z^2, 1, lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(g$supt_t$p, 2 * pt(-z, df = 40))
  expect_equal(g$swald_t$p, 2 * pt(-z, df = 40))
})

test_that(".gauss_joint with independent jumps matches the closed forms", {
  K <- 6
  theta <- c(1.2, -0.4, 2.1, 0.3, -1.6, 0.9)
  nu <- c(30, 50, 80, 40, 60, 100)
  g <- .gauss_joint(theta, diag(K), nu = nu)
  # normal: P(max |Z_k| >= m) = 1 - (1 - 2 pnorm(-m))^K and sum z^2 ~ chi-squared(K);
  # studentized: P(max |T_k| >= m) = 1 - prod_k (1 - 2 pt(-m, nu_k)). The Monte Carlo standard
  # error with 1e5 draws is about 0.0013 at these p-values; allow four of them.
  m <- max(abs(theta))
  expect_lt(abs(g$supt_normal$p - (1 - (1 - 2 * pnorm(-m))^K)), 0.005)
  expect_lt(abs(g$swald_normal$p - pchisq(sum(theta^2), K, lower.tail = FALSE)), 0.005)
  expect_lt(abs(g$supt_t$p - (1 - prod(1 - 2 * pt(-m, df = nu)))), 0.005)
  expect_gt(g$supt_t$p, g$supt_normal$p)   # t references have heavier tails
  expect_identical(g$supt_t$cell, 3L)
  expect_identical(g$supt_t$df, 6L)
})

test_that(".gauss_joint is reproducible and leaves the caller's RNG state alone", {
  Sigma <- matrix(c(1, 0.6, 0.6, 1), 2)
  set.seed(99); a <- runif(1)
  g1 <- .gauss_joint(c(1.5, 2), Sigma, nu = c(50, 60))
  b <- runif(1)
  set.seed(99); a2 <- runif(1); b2 <- runif(1)
  expect_identical(b, b2)
  g2 <- .gauss_joint(c(1.5, 2), Sigma, nu = c(50, 60))
  expect_identical(g1, g2)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  .gauss_joint(c(1.5, 2), Sigma, nu = c(50, 60))
  expect_false(exists(".Random.seed", envir = globalenv()))
})

test_that(".gauss_joint drops zero-SE jumps without information and flags deterministic ones", {
  Sigma <- diag(c(1, 0, 1))
  g <- .gauss_joint(c(1, 0, 2), Sigma, nu = c(50, 50, 50))
  expect_identical(g$supt_t$df, 2L)
  expect_identical(g$supt_t$dropped, 1L)
  expect_identical(g$supt_t$cell, 3L)   # the index refers to the input vector
  g2 <- .gauss_joint(c(1, 0.5, 2), Sigma, nu = c(50, 50, 50))
  expect_true(isTRUE(g2$supt_t$degenerate))
  expect_true(is.na(g2$swald_normal$p))
})

test_that(".jump_df matches the Bell-McCaffrey degrees of freedom of the robust variance", {
  set.seed(11)
  x <- runif(1500, -1, 1); y <- as.numeric(runif(1500) < 0.4)
  f <- rd_period(y = y, x = x, h = 0.5, b = 0.8, id = seq_along(x))
  # full Bell-McCaffrey: V-hat = sum_s e_s' M_s' A_s M_s e_s, nu = (sum tr B)^2 / sum tr(B^2),
  # with M_s the residual maker of the order-q pilot fit (the bias-corrected residuals)
  bm_side <- function(xs) {
    w <- pmax(0, 1 - abs(xs) / 0.8); act <- w > 0; xs <- xs[act]; w <- w[act]
    side <- rddid:::.rd_side_fit(rep(0, length(xs)), xs, seq_along(xs), 0, 0.5, 0.8, 1L, 2L,
                                 "triangular")
    Rq <- outer(xs, 0:2, `^`); P <- Rq %*% solve(crossprod(Rq * w, Rq)) %*% t(Rq * w)
    M <- diag(length(xs)) - P; B <- crossprod(M, side$a_bc^2 * M)
    c(sum(diag(B)), sum(B * B))
  }
  tr <- bm_side(x[x >= 0]) + bm_side(x[x < 0])
  expect_lt(abs(.jump_df(f, bc = TRUE) / (tr[1]^2 / tr[2]) - 1), 0.03)
})

test_that("rd_typecont: joint = 'wald' reproduces the Wald, and all three are reported", {
  tc_s <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
  tc_w <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, joint = "wald")
  expect_identical(tc_w$statistic, tc_w$ll_wald$stat)
  expect_identical(tc_w$p_value, tc_w$ll_wald$p)
  expect_identical(tc_s$ll_wald, tc_w$ll_wald)
  expect_equal(tc_s$joint_tests, tc_w$joint_tests)
  expect_identical(tc_s$joint_method, "supt")
  expect_true(tc_s$studentize)
  expect_equal(tc_s$statistic, max(abs(tc_s$z)))
  expect_equal(tc_s$p_value, tc_s$joint_tests$p_value[1])
  tc_n <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, studentize = FALSE)
  expect_equal(tc_n$p_value, tc_s$joint_tests$p_value[2])
  sw <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, joint = "swald")
  expect_equal(sw$statistic, sum(sw$z^2))
})

test_that("rd_compstable: the joint test stacks the pairs with their cross-pair covariance", {
  cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
  th <- cs$joint$theta
  S  <- cs$joint$Sigma
  expect_length(th, sum(vapply(cs$pairs, function(p) length(p$jumps), 1L)))
  own <- rep(seq_along(cs$pairs), vapply(cs$pairs, function(p) length(p$jumps), 1L))
  # each pair's own block is the pair's covariance
  for (a in seq_along(cs$pairs))
    expect_equal(S[own == a, own == a], unname(cs$pairs[[a]]$cov_matrix))
  # the blocks across pairs are not zero: both pairs use the RD-period units
  expect_gt(max(abs(S[own == 1, own == 2])), 0)
  expect_equal(S, t(S))
  # the joint Wald uses the full covariance (no longer the sum of the pair statistics)
  summed <- sum(vapply(cs$pairs, function(p) p$ll_wald$stat, 1))
  expect_false(isTRUE(all.equal(cs$joint$ll_wald$stat, summed)))
  expect_equal(cs$joint$ll_wald$stat, .joint_wald(th, S)$stat)
  # under "cs" only the RD-period term enters across pairs
  fa <- cs$pairs[[1]]$fits[[cs$pairs[[1]]$kept[1]]]
  fb <- cs$pairs[[2]]$fits[[cs$pairs[[2]]$kept[1]]]
  cc <- .cross_cov(fa, fb, bc = TRUE)
  expect_identical(.compstable_cross_cov(fa, fb, "cs", TRUE), cc$pc_p)
  expect_identical(.compstable_cross_cov(fa, fb, "pc", TRUE), cc$pc)
  expect_identical(.compstable_cross_cov(fa, fb, "auto", TRUE), cc$pc - cc$pv)
})

test_that("rd_compstable with one comparison period: the joint test is the pair's test", {
  cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3,
                      comparisons = 1)
  pr <- cs$pairs[[1]]
  expect_identical(cs$joint$ll_wald$stat, pr$ll_wald$stat)
  expect_identical(cs$p_value, pr$test$p)
})

test_that("rd_homog / rd_trendcell: bwselect = 'common' is the pooled comparison-period CCT h", {
  d <- rddid_sim_pv
  pooled <- d[d$year %in% c(1, 2), ]
  h_ref <- unname(rd_bw_cct(pooled$Y, pooled$R, c = 0, p = 1L, kernel = "triangular")[["h"]])
  hg <- rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_identical(hg$bandwidth, list(method = "common", h = h_ref))
  # the same as passing that h
  hf <- rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, h = h_ref)
  expect_identical(hf$statistic, hg$statistic)
  expect_identical(hf$bandwidth$method, "fixed")
  tr <- rd_trendcell(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_identical(tr$bandwidth$h, h_ref)
  # per-cell CCT is still available
  hc <- rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct")
  expect_identical(hc$bandwidth, list(method = "cct", h = NA_real_))
  expect_false(identical(hc$statistic, hg$statistic))
})

test_that("rd_checks() equals the four separate calls and passes options through", {
  ch <- rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
  hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_identical(ch$table$p_value[1], tc$p_value)
  expect_identical(ch$table$p_value[3], hg$p_value)
  expect_identical(ch$typecont$p_value, tc$p_value)
  ch2 <- rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                   joint = "wald", outcome_bwselect = "cct", homog = list(bwselect = "common"))
  expect_identical(ch2$typecont$joint_method, "wald")
  expect_identical(ch2$homog$bandwidth$method, "common")      # the per-test list wins
  expect_identical(ch2$trendcell$bandwidth$method, "cct")
  expect_named(tidy(ch), c("test", "assumption", "statistic", "p.value", "method", "bandwidth",
                           "error"))
  ch3 <- rd_checks(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_true(all(!is.na(ch3$table$error)))                    # no switchers: every test says why
  expect_output(print(ch), "Reading at the 5% level")
})
