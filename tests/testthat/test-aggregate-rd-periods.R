# Aggregate over several RD periods (paper: Corollary cor:agg, Section 5.6 sec:est-agg).
# ATT^phi = sum_t phi_t ATT(t) is estimated by one linear combination of the per-period jumps:
# coefficient phi_t on each RD period and -w^phi_t0 = -sum_t phi_t w_t0(t) on each comparison
# period. Every expected value below is computed outside rddid()'s aggregation code: from the
# single-period estimates, from rd_period() fits, or from the per-unit influence vectors.

# five periods: comparisons 1, 2; RD periods 3, 4, 5. Confounding jump linear in t, effect
# 1, 1.2, 1.4 in the RD periods. "pc": running variable fixed; "pv": it drifts (units switch
# sides); "cs": new units every period.
.agg_data <- function(scheme = c("pc", "pv", "cs"), n = 1200, seed = 11) {
  scheme <- match.arg(scheme)
  set.seed(seed)
  base <- stats::runif(n, -1, 1)
  u    <- stats::rnorm(n, 0, 0.5)
  rows <- lapply(1:5, function(t) {
    R <- switch(scheme,
                pc = base,
                pv = pmin(pmax(base + stats::rnorm(n, 0, 0.08 * (t - 1)), -1), 1),
                cs = stats::runif(n, -1, 1))
    id <- if (scheme == "cs") paste0(t, "_", seq_len(n)) else seq_len(n)
    above <- as.numeric(R >= 0)
    tau <- if (t >= 3) 1 + 0.2 * (t - 3) else 0
    Y <- 0.8 * R + 0.6 * R^2 + (0.5 + 0.1 * t) * above + tau * above +
      (if (scheme == "cs") 0 else u) + stats::rnorm(n, 0, 0.3)
    data.frame(id = id, year = t, R = R, Y = Y)
  })
  do.call(rbind, rows)
}
.fit <- function(d, ...) suppressMessages(rddid(d, y = "Y", x = "R", time = "year", id = "id", ...))
.rd <- 3:5
.comps <- c(1, 2)

test_that("aggregate = phi-weighted sum of the single-period estimates (fixed h and cct; constant and linear; equal, unequal and non-summing phi)", {
  d <- .agg_data("pc")
  for (tr in c("constant", "linear")) for (phi in list(NULL, c(0.2, 0.3, 0.5), c(1, -0.5, 2))) {
    ph <- if (is.null(phi)) rep(1 / 3, 3) else phi
    for (bw in list(list(h = 0.4), list(bwselect = "cct"))) {
      agg <- do.call(.fit, c(list(d, t_rd = .rd, comparisons = .comps, trend = tr,
                                  rd_weights = phi), bw))
      one <- lapply(.rd, function(t)
        do.call(.fit, c(list(d, t_rd = t, comparisons = .comps, trend = tr), bw)))
      for (row in c("Conventional", "Robust"))
        expect_equal(agg$estimates[row, "est"],
                     sum(ph * vapply(one, function(f) f$estimates[row, "est"], numeric(1))),
                     tolerance = 1e-10, info = paste(tr, row, names(bw)))
    }
  }
})

test_that("CS standard error = sqrt(sum coef^2 V(D-hat_t)) from rd_period() fits", {
  d <- .agg_data("cs")
  phi <- c(0.2, 0.3, 0.5)
  agg <- .fit(d, t_rd = .rd, comparisons = .comps, rd_weights = phi, h = 0.4)
  expect_identical(agg$scheme, "cs")
  coef <- c(phi, -rep(sum(phi) / 2, 2))                      # constant trend: w_t0(t) = 1/2
  per  <- c(.rd, .comps)
  fits <- lapply(per, function(t) { s <- d[d$year == t, ]
    rd_period(s$Y, s$R, h = 0.4, b = 0.4, id = s$id) })
  expect_equal(agg$estimates["Conventional", "se_cs"],
               sqrt(sum(coef^2 * vapply(fits, function(f) f$V_D, numeric(1)))), tolerance = 1e-10)
  expect_equal(agg$estimates["Robust", "se_cs"],
               sqrt(sum(coef^2 * vapply(fits, function(f) f$V_D_bc, numeric(1)))), tolerance = 1e-10)
})

test_that("PC and PV standard errors = sqrt(c' Sigma c), Sigma built pair by pair from the influence vectors (RD-RD pairs included)", {
  d <- .agg_data("pv")
  phi <- c(0.5, 0.25, 0.25)
  agg <- .fit(d, t_rd = .rd, comparisons = .comps, rd_weights = phi, trend = "linear", h = 0.45)
  expect_identical(agg$scheme, "pv")
  per  <- c(.rd, .comps)
  fits <- lapply(per, function(t) { s <- d[d$year == t, ]
    rd_period(s$Y, s$R, h = 0.45, b = 0.45, id = s$id) })
  # comparison weights for target t (linear, comparisons 1 and 2): (2 - t, t - 1)
  wphi <- c(sum(phi * (2 - .rd)), sum(phi * (.rd - 1)))
  coef <- c(phi, -wphi)
  # Cov(b_side1_t, b_side2_s) = sum over units in both windows of g_t g_s
  pair <- function(ft, fs, s1, s2, g) {
    a <- ft$sides[[s1]]; b <- fs$sides[[s2]]
    m <- match(a$id, b$id); ok <- !is.na(m)
    sum(a[[g]][ok] * b[[g]][m[ok]])
  }
  for (g in c("g", "g_bc")) {
    K <- length(per); S_pc <- S_pv <- matrix(0, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      if (i == j) {
        v <- sum(fits[[i]]$sides[["+"]][[g]]^2) + sum(fits[[i]]$sides[["-"]][[g]]^2)
        S_pc[i, i] <- S_pv[i, i] <- v
      } else {
        same <- pair(fits[[i]], fits[[j]], "+", "+", g) + pair(fits[[i]], fits[[j]], "-", "-", g)
        opp  <- pair(fits[[i]], fits[[j]], "+", "-", g) + pair(fits[[i]], fits[[j]], "-", "+", g)
        S_pc[i, j] <- same
        S_pv[i, j] <- same - opp
      }
    }
    row <- if (g == "g") "Conventional" else "Robust"
    expect_equal(agg$estimates[row, "se_pc"], sqrt(drop(t(coef) %*% S_pc %*% coef)), tolerance = 1e-10)
    expect_equal(agg$estimates[row, "se_pv"], sqrt(drop(t(coef) %*% S_pv %*% coef)), tolerance = 1e-10)
    expect_true(any(abs(S_pv[1:3, 1:3][upper.tri(diag(3))]) > 0))   # RD-RD covariances are live
  }
})

test_that("linear trend: aggregated comparison weights are those of Corollary 1 at the mean RD date", {
  d <- .agg_data("pc")
  f <- .fit(d, t_rd = .rd, comparisons = .comps, trend = "linear", h = 0.4)
  expect_equal(unname(f$weights_by_target), rbind(c(-1, 2), c(-2, 3), c(-3, 4)), tolerance = 1e-12)
  expect_equal(f$weights, c(`1` = -2, `2` = 3), tolerance = 1e-12)      # t-bar = 4
  expect_equal(unname(f$coef), c(1, 1, 1, 6, -9) / 3, tolerance = 1e-12)
  # shifting the period labels changes nothing (the weights are computed with time centred at
  # the RD period, so calendar years lose no precision)
  d2 <- transform(d, year = year + 1998)
  f2 <- .fit(d2, t_rd = 2001:2003, comparisons = c(1999, 2000), trend = "linear", h = 0.4)
  expect_equal(unname(f2$weights), c(-2, 3), tolerance = 1e-12)
  expect_equal(f2$estimates, f$estimates, tolerance = 1e-10)
})

test_that("a single RD period is unchanged: rd_weights = 1 gives identical results", {
  d <- .agg_data("pc")
  a <- .fit(d, t_rd = 3, comparisons = .comps)
  b <- .fit(d, t_rd = 3, comparisons = .comps, rd_weights = 1)
  expect_identical(a$estimates, b$estimates)
  expect_identical(a$coef, b$coef)
  expect_identical(a$bandwidth, b$bandwidth)
  expect_equal(a$rd_weights, c(`3` = 1))
})

test_that("estimand = \"atu\" returns the same numbers for the aggregate", {
  d <- .agg_data("pv")
  a <- .fit(d, t_rd = .rd, comparisons = .comps, estimand = "att")
  u <- .fit(d, t_rd = .rd, comparisons = .comps, estimand = "atu")
  expect_identical(u$estimates, a$estimates)
  expect_identical(u$bandwidth, a$bandwidth)
})

test_that("joint rule: one h for the aggregate, the closed form in B, Veff and reg, with B from pilot fits", {
  d <- .agg_data("pc")
  phi <- c(0.2, 0.3, 0.5)
  f <- .fit(d, t_rd = .rd, comparisons = .comps, rd_weights = phi, bwselect = "joint")
  bw <- f$bandwidth
  expect_equal(bw$h, (bw$Veff / (bw$B^2 + bw$reg))^(1 / 5), tolerance = 1e-12)
  # B-bar^phi = sum_t coef_t B-hat_t, each B-hat_t from a pilot fit at that period's pilot pair
  per <- names(f$coef)
  bconst <- vapply(per, function(k) { s <- d[d$year == as.numeric(k), ]
    pb <- bw$pilot_bws[[k]]
    rd_period(s$Y, s$R, h = pb[["h"]], b = pb[["b"]], id = s$id)$b_const }, numeric(1))
  expect_equal(bw$B, sum(f$coef[per] * bconst), tolerance = 1e-10)
  # every period is fitted at the common h
  expect_true(all(abs(bw$h_by_period - bw$h) < 1e-12))
  # it is the aggregate's bandwidth, not any single period's
  for (t in .rd) {
    one <- .fit(d, t_rd = t, comparisons = .comps, bwselect = "joint")
    expect_false(isTRUE(all.equal(one$bandwidth$h, bw$h, tolerance = 1e-6)))
  }
})

test_that("same result as the one-RD-period device (other RD periods as comparisons with negative weights)", {
  d <- .agg_data("pv")
  phi <- c(0.2, 0.3, 0.5)
  wphi <- c(sum(phi * (2 - .rd)), sum(phi * (.rd - 1)))
  for (bw in list(list(h = 0.4), list(bwselect = "cct"), list(bwselect = "joint"))) {
    agg <- do.call(.fit, c(list(d, t_rd = .rd, comparisons = .comps, trend = "linear",
                                rd_weights = phi), bw))
    dev <- do.call(.fit, c(list(d, t_rd = 3, comparisons = c(4, 5, 1, 2),
                                trend = c(-phi[2:3], wphi) / phi[1]), bw))
    cols <- c("est", "se", "se_cs", "se_pc", "se_pv")
    expect_equal(as.matrix(agg$estimates[, cols]), phi[1] * as.matrix(dev$estimates[, cols]),
                 tolerance = 1e-10, info = names(bw))
    expect_equal(agg$bandwidth$h_by_period[names(dev$bandwidth$h_by_period)],
                 dev$bandwidth$h_by_period, tolerance = 1e-10, info = names(bw))
  }
})

test_that("a numeric matrix trend gives per-RD-period comparison weights", {
  d <- .agg_data("pc")
  M <- rbind(c(0.5, 0.5), c(-1, 2), c(0.3, 0.7))
  f <- .fit(d, t_rd = .rd, comparisons = .comps, trend = M, h = 0.4)
  expect_equal(unname(f$weights_by_target), M)
  expect_equal(unname(f$weights), colMeans(M), tolerance = 1e-12)
  expect_identical(f$weights_type, "custom")
  # named rows are matched to t_rd
  Mn <- M[3:1, ]; rownames(Mn) <- c("5", "4", "3"); colnames(Mn) <- c("1", "2")
  fn <- .fit(d, t_rd = .rd, comparisons = .comps, trend = Mn, h = 0.4)
  expect_equal(fn$estimates, f$estimates, tolerance = 1e-12)
})

test_that("any finite rd_weights are allowed, without a warning; named weights are matched to t_rd", {
  d <- .agg_data("pc")
  expect_no_warning(f <- .fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(1, -0.5, 2), h = 0.4))
  expect_equal(unname(f$weights), rep(2.5 / 2, 2), tolerance = 1e-12)
  fn <- .fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(`5` = 2, `3` = 1, `4` = -0.5), h = 0.4)
  expect_equal(fn$estimates, f$estimates, tolerance = 1e-12)
})

test_that("input errors", {
  d <- .agg_data("pc")
  expect_error(.fit(d, t_rd = c(3, 3), comparisons = .comps), "more than once")
  expect_error(.fit(d, t_rd = c(2, 3), comparisons = .comps), "must not be one of the `comparisons`")
  expect_error(.fit(d, t_rd = c(3, 9), comparisons = .comps), "`t_rd` = 9 is not a period in `year`", fixed = TRUE)
  expect_error(.fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(1, 2)), "one entry per RD period")
  expect_error(.fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(1, NA, 2)), "finite")
  expect_error(.fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(`3` = 1, `4` = 1, `9` = 1)),
               "names of `rd_weights`")
  expect_error(.fit(d, t_rd = .rd, comparisons = .comps, trend = matrix(0.5, 2, 2)),
               "one row per RD period")
  expect_error(.fit(d, t_rd = .rd, comparisons = .comps, weighting = "min_variance"),
               "single RD period")
})

test_that("methods: summary roles, glance, print and plot with several RD periods", {
  d <- .agg_data("pc")
  f <- .fit(d, t_rd = .rd, comparisons = .comps, h = 0.4)
  tab <- summary(f)$per_period
  expect_identical(tab$role, c("comparison", "comparison", "RD", "RD", "RD"))
  expect_equal(tab$coef, c(-0.5, -0.5, 1 / 3, 1 / 3, 1 / 3), tolerance = 1e-12)
  g <- glance(f)
  expect_equal(nrow(g), 1L)
  expect_identical(g$t_rd, "3, 4, 5")
  out <- capture.output(print(f))
  expect_match(out[1], "averaged over RD periods 3, 4, 5", fixed = TRUE)
  expect_match(out[2], "aggregated weights", fixed = TRUE)
  expect_identical(glance(.fit(d, t_rd = 3, comparisons = .comps, h = 0.4))$t_rd, 3)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(f), "ggplot")
})

test_that("min_variance with one RD period scales exactly with rd_weights (pinned and not pinned)", {
  r <- function(...) suppressMessages(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id",
                                            t_rd = 3, h = 0.3, weighting = "min_variance", ...))
  for (comps in list(2, c(1, 2))) {                  # one comparison: pinned; two: not pinned
    a1 <- r(comparisons = comps); a2 <- r(comparisons = comps, rd_weights = 2)
    expect_equal(a2$weights_detail$pinned, length(comps) == 1L)
    expect_equal(a2$coef, 2 * a1$coef, tolerance = 1e-12)
    expect_equal(a2$estimates[, "est"], 2 * a1$estimates[, "est"], tolerance = 1e-12)
    expect_equal(a2$weights_by_target, a1$weights_by_target, tolerance = 1e-12)
    expect_equal(a2$weights_detail$pilot, a1$weights_detail$pilot, tolerance = 1e-12)
  }
})

test_that("joint rule, PC, equal pilots: Veff / h0 = c' Sigma c with RD-RD covariances built by hand", {
  d <- .agg_data("pc")
  phi <- c(0.2, 0.3, 0.5)
  coef <- c(stats::setNames(phi, .rd), stats::setNames(-rep(sum(phi) / 2, 2), .comps))
  per <- names(coef); h0 <- 0.4; b0 <- 0.55
  plist <- stats::setNames(lapply(per, function(k) { s <- d[d$year == as.numeric(k), ]
    data.frame(y = s$Y, x = s$R, id = s$id) }), per)
  pb <- stats::setNames(rep(list(c(h = h0, b = b0)), length(per)), per)
  jb <- .bw_joint(plist, coef, scheme = "pc", pilot_bws = pb)
  fits <- lapply(per, function(k) rd_period(plist[[k]]$y, plist[[k]]$x, h = h0, b = b0, id = plist[[k]]$id))
  same <- function(ft, fs) sum(vapply(c("+", "-"), function(sd) {
    a <- ft$sides[[sd]]; b <- fs$sides[[sd]]; m <- match(a$id, b$id); ok <- !is.na(m)
    sum(a$g[ok] * b$g[m[ok]]) }, numeric(1)))
  K <- length(per); S <- matrix(0, K, K)
  for (i in seq_len(K)) for (j in seq_len(K)) S[i, j] <- if (i == j) fits[[i]]$V_D else same(fits[[i]], fits[[j]])
  expect_equal(jb$Veff / h0, drop(t(coef) %*% S %*% coef), tolerance = 1e-10)
  expect_true(abs(S[1, 2]) > 0)                                     # an RD-RD pair is live
})

test_that("character and factor period labels give the same aggregate as numeric labels", {
  d <- .agg_data("pc")
  f  <- .fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(0.2, 0.3, 0.5), h = 0.4)
  dc <- transform(d, year = paste0("y", year))
  fc <- .fit(dc, t_rd = c("y3", "y4", "y5"), comparisons = c("y1", "y2"),
             rd_weights = c(y5 = 0.5, y3 = 0.2, y4 = 0.3), h = 0.4)
  expect_equal(fc$estimates, f$estimates, tolerance = 1e-12)
  ff <- .fit(dc, t_rd = factor(c("y3", "y4", "y5")), comparisons = c("y1", "y2"),
             rd_weights = c(0.2, 0.3, 0.5), h = 0.4)
  expect_equal(ff$estimates, f$estimates, tolerance = 1e-12)
  expect_error(.fit(d, t_rd = .rd, comparisons = .comps, rd_weights = c(0, 0, 0)), "all zero")
})

test_that("linear-trend weights are exact with calendar years and invariant to shifting the labels", {
  w <- function(comps, t) unname(.rddid_weights("linear", comps, t))
  expect_equal(w(c(1999, 2000), 2004), c(-4, 5), tolerance = 1e-14)
  expect_equal(w(c(1999, 2000), 2001), c(-1, 2), tolerance = 1e-14)
  # three comparison periods: the OLS line, same weights whatever the origin of the time axis
  for (shift in c(0, 1000, 1990, 1e5))
    expect_equal(w(c(1, 2, 4) + shift, 6 + shift), w(c(1, 2, 4), 6), tolerance = 1e-12)
  # weights sum to one and reproduce t_rd (Corollary 1's constraints)
  ww <- w(c(1996, 1998, 1999, 2000), 2003)
  expect_equal(sum(ww), 1, tolerance = 1e-13)
  expect_equal(sum(ww * c(1996, 1998, 1999, 2000)), 2003, tolerance = 1e-9)
})
