# The plot methods: each returns a ggplot object that builds without error, and the objects
# carry the fields the methods read.

test_that("fits and tests store what the plot methods need", {
  fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  expect_named(fit$data, c("3", "1", "2"))
  expect_equal(nrow(fit$data[["3"]]), 1000L)
  tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id")
  expect_equal(dim(tc$fits), c(length(tc$meta$type_values), length(tc$meta$periods)))
  expect_equal(tc$meta$c, 0)
  expect_true(all(c("id", "R_1", "side_1", "R_3", "side_3") %in% names(tc$sides)))
  cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
  pr <- cs$pairs[[1]]
  expect_named(pr$sample, c("x_trd", "type_trd", "x_t0", "type_t0"))
  expect_equal(length(pr$fits), length(pr$type_values))
})

test_that("every plot method returns a ggplot that builds", {
  skip_if_not_installed("ggplot2")
  built <- function(p) { expect_s3_class(p, "ggplot"); expect_silent(ggplot2::ggplot_build(p)) }
  fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  built(plot(fit))
  tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id")
  built(plot(tc))
  built(plot(tc, t_rd = 3, comparison = 2))
  expect_error(plot(tc, t_rd = 3, comparison = 3), "two different periods")
  cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
  built(plot(cs))
  built(plot(cs, pair = 2))
  built(plot(cs, pair = "3::1"))
  hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  built(plot(hg))
  tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  built(plot(tr))
  built(plot_switchers(rddid_sim_pv, x = "R", time = "year", id = "id", periods = c(1, 3)))
  built(plot_switchers(rddid_sim_pv, x = "R", time = "year", id = "id"))
  expect_error(plot_switchers(rddid_sim_pv, x = "R", time = "year", id = "id", periods = 1),
               "two values")
})

test_that("the binned-means helper bins over the requested range", {
  b <- rddid:::.bin_means(seq(-1, 1, length.out = 201), rep(1, 201), bins = 4, range = c(-1, 1))
  expect_equal(nrow(b), 4L)
  expect_equal(b$y, rep(1, 4))
  expect_equal(sum(b$n), 201L)
})
