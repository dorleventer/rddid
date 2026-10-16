# Behaviour fixed after the 2026-10-16 code and UX reviews: argument checks, error headers,
# period handling, and what the printouts and tidiers report.

d  <- rddid_sim
dp <- rddid_sim_pv
fit <- function(...) suppressMessages(rddid(d, y = "Y", x = "R", time = "year", id = "id", h = 0.5, ...))
caught <- function(expr) tryCatch(expr, error = function(e) e)

test_that("errors name the user's call, also when raised inside a helper", {
  for (e in list(caught(fit(t_rd = 3, trend = "quadratic")),
                 caught(fit(t_rd = 3, trend = c(0.5, 0.3, 0.2))),
                 caught(fit(t_rd = 3, kernel = "gauss")),
                 caught(rddid(d, y = "Y", x = "r", time = "year", id = "id", t_rd = 3)))) {
    expect_identical(as.character(conditionCall(e)[[1L]]), "rddid")
  }
  e <- caught(fit(t_rd = 3, trend = c(0.5, 0.3, 0.2)))
  expect_match(conditionMessage(e), "one entry per comparison period (2); got 3", fixed = TRUE)
  e <- caught(do.call(rddid, list(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                                  kernel = "gauss")))
  expect_identical(conditionCall(e), quote(rddid()))     # named, without the data
})

test_that("a missing column names its argument and lists the columns", {
  e <- caught(rddid(d, y = "Y", x = "r", time = "yr", id = "id", t_rd = 3))
  expect_match(conditionMessage(e), "not a column of `data`: `x` = \"r\", `time` = \"yr\"", fixed = TRUE)
  expect_match(conditionMessage(e), "The columns are id, year, R", fixed = TRUE)
})

test_that("a factor running variable gets a conversion hint", {
  df <- dp; df$R <- factor(round(df$R, 1))
  expect_error(rddid(df, y = "Y", x = "R", time = "year", id = "id", t_rd = 3),
               "it is a factor (as.numeric(as.character(data$R))", fixed = TRUE)
})

test_that("rd_weights: none with a single RD period; names on every entry or none", {
  expect_error(fit(t_rd = 3, rd_weights = -2), "leave it out")
  s4 <- rbind(d, transform(d[d$year == 3, ], year = 4))
  f4 <- function(...) suppressMessages(rddid(s4, y = "Y", x = "R", time = "year", id = "id",
                                             h = 0.5, t_rd = c(4, 3), comparisons = 1:2, ...))
  expect_error(f4(rd_weights = c(`3` = 0.2, 0.8)), "named on some entries only")
  expect_equal(f4(rd_weights = c(`3` = 0.2, `4` = 0.8))$rd_weights, c(`4` = 0.8, `3` = 0.2))
  expect_error(fit(t_rd = 3, comparisons = 1:2, trend = c(`2` = 0.9, 0.1)), "named on some entries only")
})

test_that("periods given as text or factor are read in the period column's type", {
  a <- fit(t_rd = 3, trend = "linear")
  expect_equal(fit(t_rd = "3", trend = "linear")$estimates, a$estimates, tolerance = 0)
  expect_equal(fit(t_rd = factor(3), comparisons = factor(1:2), trend = "linear")$estimates,
               a$estimates, tolerance = 0)
  dc <- d; dc$year <- as.character(dc$year)
  expect_error(suppressMessages(rddid(dc, y = "Y", x = "R", time = "year", id = "id", t_rd = "3",
                                      trend = "linear")),
               "the period column `year` must be numeric (it is character)", fixed = TRUE)
})

test_that("rd_homog / rd_trendcell name a comparison period with no observed outcome", {
  na <- dp; na$Y[na$year == 1] <- NA
  for (f in list(rd_homog, rd_trendcell))
    expect_error(suppressMessages(f(na, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)),
                 "no row with an observed `Y` in comparison period(s) 1", fixed = TRUE)
})

test_that("the tests ask for `id` and for `t_rd` with the same message everywhere", {
  expect_error(rd_typecont(dp, x = "R", time = "year", t_rd = 3), "`id` (the unit identifier)",
               fixed = TRUE)
  for (f in list(function(...) rddid(dp, y = "Y", ...), function(...) rd_compstable(dp, ...),
                 function(...) rd_homog(dp, y = "Y", ...), function(...) rd_trendcell(dp, y = "Y", ...)))
    expect_error(f(x = "R", time = "year", id = "id"), "`t_rd` (the RD period) is required.",
                 fixed = TRUE)
})

test_that("rd_checks(): level is a confidence level, as in rddid()", {
  ch <- suppressMessages(rd_checks(dp, y = "Y", x = "R", time = "year", id = "id", t_rd = 3))
  expect_identical(ch$level, 0.95)
  out <- capture.output(print(ch))
  expect_true(any(grepl("type continuity +not rejected", out)))      # p = 0.868
  expect_error(rd_checks(dp, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, level = 95),
               "strictly between 0 and 1")
})

test_that("summary(): time order for text periods; the robust coefficients under min_variance", {
  dc <- d; dc$per <- c("early", "mid", "late")[dc$year]
  dc$per <- factor(dc$per, levels = c("early", "mid", "late"))
  g <- suppressMessages(rddid(dc, y = "Y", x = "R", time = "per", id = "id", t_rd = "late", h = 0.5))
  expect_identical(summary(g)$per_period$period, c("early", "mid", "late"))
  mv <- fit(t_rd = 3, weighting = "min_variance")
  tab <- summary(mv)$per_period
  expect_equal(sum(tab$coef_bc * tab$jump_bc), mv$estimates["Robust", "est"], tolerance = 1e-12)
  expect_equal(sum(tab$coef * tab$jump), mv$estimates["Conventional", "est"], tolerance = 1e-12)
  expect_null(summary(fit(t_rd = 3))$per_period$coef_bc)
})

test_that("tidy() of a test: df only for Wald statistics; test and assumption columns", {
  tc <- rd_typecont(dp, x = "R", time = "year", id = "id", t_rd = 3)
  td <- tidy(tc)
  expect_identical(td$test, "rd_typecont")
  expect_identical(td$assumption, "type continuity")
  expect_true(is.na(td$df))                                      # sup-t: no chi-squared df
  hg <- suppressMessages(rd_homog(dp, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, h = 0.5))
  expect_identical(tidy(hg)$df, hg$df)
  s4 <- rbind(d, transform(d[d$year == 3, ], year = 4))
  g <- suppressMessages(rddid(s4, y = "Y", x = "R", time = "year", id = "id", t_rd = c(3, 4), h = 0.5))
  expect_identical(glance(g)$rd_weights, "0.5, 0.5")
})

test_that("rd_trendcell(trend = \"linear\") is titled as a linear-trend test", {
  tr <- suppressMessages(rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id",
                                      t_rd = 3, h = 0.5))
  expect_output(print(tr), "Test of constant within-type confounding")
  expect_identical(tidy(tr)$assumption, "constant within-type confounding")
})
