# Snapshot tests of the console text (the numerical gate compares objects, not printouts).
# After an intended change to a printout, run testthat::snapshot_review() and accept.

local_width <- function() testthat::local_reproducible_output(width = 100, unicode = FALSE)
# The sup-t and sWald p-values are simulated: the draws go through an eigendecomposition whose eigenvector signs depend
# on the platform's LAPACK, so their third decimal differs across platforms. Mask it.
mask_sim_p <- function(x) {
  sim <- grepl("critical values", x)
  x[sim] <- gsub("p = ([0-9]\\.[0-9]{2})[0-9]", "p = \\1*", x[sim])
  x[sim] <- sub("^(\\s+\\S.*\\s[0-9.]+\\s+)([0-9]\\.[0-9]{2})[0-9](\\s+(sup-t|sWald),)", "\\1\\2*\\3", x[sim])
  x
}

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
  expect_snapshot(print(rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)),
                  transform = mask_sim_p)
  expect_snapshot(print(rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)),
                  transform = mask_sim_p)
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

test_that("print of rd_checks()", {
  skip_on_cran()
  expect_snapshot(print(rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id",
                                  t_rd = 3)), transform = mask_sim_p)
})
