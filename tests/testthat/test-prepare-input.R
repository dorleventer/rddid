# The shared input step (.prepare_input(), input_checks.R): the same checks and messages in
# rddid() and the four tests, and the same handling of missing values.

d <- rddid_sim_pv
no_call <- function(obj) { obj <- unclass(obj); obj[setdiff(names(obj), "call")] }

calls <- list(
  rddid         = function(dat, ...) rddid(dat, y = "Y", x = "R", time = "year", id = "id", h = 0.5, ...),
  rd_typecont   = function(dat, ...) rd_typecont(dat, x = "R", time = "year", id = "id", h = 0.5, ...),
  rd_compstable = function(dat, ...) rd_compstable(dat, x = "R", time = "year", id = "id", h = 0.5, ...),
  rd_homog      = function(dat, ...) rd_homog(dat, y = "Y", x = "R", time = "year", id = "id", h = 0.5, ...),
  rd_trendcell  = function(dat, ...) rd_trendcell(dat, y = "Y", x = "R", time = "year", id = "id", h = 0.5, ...)
)

test_that("every entry point gives the same message for the same input error, under its own name", {
  dup <- rbind(d, d[d$year == 1, ][1:2, ])
  chr <- d; chr$R <- as.character(chr$R)
  # the error's call is the exported function's, so the header reads "Error in <fn>(...)"
  caught <- function(expr) tryCatch(expr, error = function(e) e)
  for (fn in names(calls)) {
    f <- calls[[fn]]
    e <- caught(f(d, t_rd = 7))
    expect_identical(conditionMessage(e), "`t_rd` = 7 is not a period in `year`.")
    expect_identical(as.character(conditionCall(e)[[1L]]), fn)
    expect_error(f(d, t_rd = 3, comparisons = c(1, 9)), "^comparison period\\(s\\) not in `year`: 9$")
    expect_error(f(d, t_rd = 3, comparisons = c(1, 1)), "contains a period more than once")
    expect_error(f(dup, t_rd = 3, comparisons = 1:2), "^2 row\\(s\\) repeat a \\(id, year\\) pair")
    expect_error(f(chr, t_rd = 3), "running variable column 'R' must be numeric")
    expect_error(f(as.list(d), t_rd = 3), "`data` must be a data frame")
  }
})

test_that("the four tests say plainly that they take one RD period", {
  for (fn in setdiff(names(calls), "rddid"))
    expect_error(calls[[fn]](d, t_rd = c(2, 3), comparisons = 1), "must be a single period")
  expect_error(rddid(d, y = "Y", x = "R", time = "year", id = "id"), "`t_rd` (the RD period) is required",
               fixed = TRUE)
})

test_that("a period with no complete row is named", {
  na_x <- d; na_x$R[na_x$year == 1] <- NA
  expect_message(expect_error(calls$rddid(na_x, t_rd = 3), "no complete row in period(s) 1", fixed = TRUE),
                 "3000|1000")
})

test_that("rddid(): rows with a missing value are dropped, with a message", {
  na_y <- d; na_y$Y[c(5, 1500)] <- NA
  expect_message(f_na <- calls$rddid(na_y, t_rd = 3), "2 row(s) with a missing value", fixed = TRUE)
  f_drop <- calls$rddid(d[-c(5, 1500), ], t_rd = 3)
  expect_equal(f_na$estimates, f_drop$estimates, tolerance = 0)
})

test_that("rddid(): duplicates in a period that is not used do not matter", {
  dup <- rbind(d, d[d$year == 1, ][1:2, ])
  expect_equal(calls$rddid(dup, t_rd = 3, comparisons = 2)$estimates,
               calls$rddid(d, t_rd = 3, comparisons = 2)$estimates, tolerance = 0)
})

test_that("rd_homog / rd_trendcell: a missing outcome drops the row from the fits, not the type", {
  # the types come from the running variable: a missing outcome in t_rd (not used in the fits)
  # changes nothing; one in comparison period 1 removes the row from period 1's cells only
  na_rd <- d; na_rd$Y[na_rd$year == 3] <- NA
  idx   <- which(d$year == 1)[c(3, 40, 700)]
  na_c  <- d; na_c$Y[idx] <- NA
  tab   <- c(rd_homog = "period_type_jumps", rd_trendcell = "cell_period_jumps")
  for (fn in c("rd_homog", "rd_trendcell")) {
    base <- calls[[fn]](d, t_rd = 3)
    expect_message(r_rd <- calls[[fn]](na_rd, t_rd = 3), "1000 row(s) with a missing `Y`", fixed = TRUE)
    expect_equal(no_call(r_rd), no_call(base), tolerance = 0)
    expect_message(r_c <- calls[[fn]](na_c, t_rd = 3), "3 row(s) with a missing `Y`", fixed = TRUE)
    j_c <- r_c[[tab[[fn]]]]; j_b <- base[[tab[[fn]]]]
    expect_equal(j_c[j_c$period == 2, ], j_b[j_b$period == 2, ], tolerance = 0)
    expect_equal(sum(j_c$n[j_c$period == 1]), sum(j_b$n[j_b$period == 1]) - 3L)
  }
})

test_that("the kernel name is not case sensitive in any entry point", {
  for (fn in names(calls))
    expect_equal(no_call(suppressMessages(calls[[fn]](d, t_rd = 3, kernel = "Triangular"))),
                 no_call(suppressMessages(calls[[fn]](d, t_rd = 3))), tolerance = 0)
})

test_that("a period outside t_rd and comparisons does not change any test", {
  # a year 4 that only half the units have: before, the tests built types from every period in
  # the data, so the units missing from year 4 were dropped (rd_homog: 1000 -> 500 units)
  d4 <- rbind(d, transform(d[d$year == 3 & d$id <= 500, ], year = 4))
  for (fn in c("rd_compstable", "rd_homog", "rd_trendcell")) {
    set.seed(1); r  <- suppressMessages(calls[[fn]](d,  t_rd = 3, comparisons = 1:2))
    set.seed(1); r4 <- suppressMessages(calls[[fn]](d4, t_rd = 3, comparisons = 1:2))
    expect_equal(no_call(r4), no_call(r), tolerance = 0)
  }
})
