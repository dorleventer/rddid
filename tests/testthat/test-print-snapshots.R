# Snapshot tests of the console text (the numerical gate compares objects, not printouts).
# After an intended change to a printout, run testthat::snapshot_review() and accept.

local_width <- function() testthat::local_reproducible_output(width = 100, unicode = FALSE)

test_that("print and summary of a fit", {
  local_width()
  fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_snapshot(print(fit))
  expect_snapshot(print(summary(fit)))
  expect_snapshot(print(coef(fit)))
  expect_snapshot(print(confint(fit)))
})

test_that("print of the four tests", {
  local_width()
  expect_snapshot(print(rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)))
  expect_snapshot(print(rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)))
  expect_snapshot(print(rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)))
  expect_snapshot(print(rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)))
})

test_that("the main error and message texts", {
  local_width()
  expect_snapshot(error = TRUE, rddid(rddid_sim, y = "Yy", x = "R", time = "year", id = "id", t_rd = 3))
  expect_snapshot(error = TRUE, rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                                      comparisons = c(1, 7)))
  expect_snapshot(error = TRUE, rd_typecont(rddid_sim, x = "R", time = "year", id = "id"))
  expect_snapshot(rddid(rddid_sim, y = "Y", x = "R", time = "year", t_rd = 3, bwselect = "cct")$scheme)
})
