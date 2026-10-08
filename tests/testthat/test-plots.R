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

test_that("the binned-means helper bins each side of the cutoff separately", {
  x <- seq(-1, 1, length.out = 201)
  b <- rddid:::.bin_sides(x, rep(1, 201), cutoff = 0, bins = 4, lo = -1, hi = 1)
  expect_equal(nrow(b), 8L)                      # 4 bins per side
  expect_true(all(b$x[1:4] < 0) && all(b$x[5:8] >= 0))   # no bin straddles the cutoff
  expect_equal(b$y, rep(1, 8))
  expect_equal(sum(b$n), 201L)
})

test_that("degenerate panels (no unit changes side) give a clear error from every test", {
  for (f in list(function() rd_typecont(rddid_sim, x = "R", time = "year", id = "id"),
                 function() rd_compstable(rddid_sim, x = "R", time = "year", id = "id", t_rd = 3),
                 function() rd_homog(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3),
                 function() rd_trendcell(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)))
    expect_error(f(), "no unit changes side")
})

test_that("glance() is one row under every bandwidth rule", {
  for (bw in c("joint", "cct", "iter")) {
    g <- generics::glance(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                                bwselect = bw))
    expect_equal(nrow(g), 1L)
  }
  expect_error(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, level = 95),
               "between 0 and 1")
  expect_error(rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                     comparisons = c(1, 7)), "not in `year`")
})

test_that("summary() per-period table is aligned with its period labels (time order)", {
  fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
  tab <- summary(fit)$per_period
  expect_equal(tab$period, c("1", "2", "3"))
  for (i in seq_len(nrow(tab))) {
    f <- fit$fits[[tab$period[i]]]
    expect_equal(tab$jump[i], f$D)
    expect_equal(tab$h[i], unname(f$h))
    expect_equal(tab$n[i], as.integer(f$n))
    expect_equal(tab$coef[i], unname(fit$coef[tab$period[i]]))
  }
  expect_equal(sum(tab$coef * tab$jump), unname(coef(fit)[["Conventional"]]))
})
