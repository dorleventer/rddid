# App. B.4 "Bandwidth estimation": the joint rule's pilots d and b for the aggregate, and the
# common h at b. The reference is rdrobust's own selector (rdbwselect, bwselect = "mserd",
# vce = "hc1"): with one period the joint rule is that selector, and in a balanced panel with a
# running variable fixed over time it is that selector on the combined outcome sum_t coef_t Y_t.

.jp_panel <- function(n = 2000, seed = 21, xdist = c("unif", "t")) {
  xdist <- match.arg(xdist)
  set.seed(seed)
  # a heavy-tailed running variable makes IQR / 1.349 < sd, so the rule of thumb takes the IQR
  # branch (rdrobust: type-2 quantiles). No ties: rdrobust's mass-point adjustment would count
  # distinct values instead of rows
  R <- if (xdist == "unif") runif(n, -1, 1) else 0.4 * rt(n, df = 3)
  ue <- rnorm(n, 0, 0.4)
  rows <- lapply(1:3, function(t) {
    m <- 0.3 * t + 0.5 * R + (0.4 + 0.5 * t) * R^2 - 0.3 * R^3 + 0.4 * (R >= 0) + (t == 3) * (R >= 0)
    data.frame(id = seq_len(n), year = t, R = R, Y = m + 0.3 * rnorm(n) + ue)
  })
  do.call(rbind, rows)
}
.jp_plist <- function(d) stats::setNames(lapply(c(3, 1, 2), function(tv) {
  r <- d$year == tv; data.frame(y = d$Y[r], x = d$R[r], id = d$id[r]) }), c("3", "1", "2"))
.jp_coef <- c(`3` = 1, `1` = -0.5, `2` = -0.5)

test_that("B.4 joint rule, one period — h and b equal rdbwselect (mserd, hc1), p = 1, 2, both branches of the rule of thumb", {
  for (xd in c("unif", "t")) for (p in 1:2) {
    q <- p + 1L
    d <- .jp_panel(xdist = xd); d3 <- d[d$year == 3, ]
    pl <- list(`3` = data.frame(y = d3$Y, x = d3$R, id = d3$id))
    jb <- .bw_joint(pl, c(`3` = 1), "3", scheme = "cs", p = p, q = q)
    rb <- rdrobust::rdbwselect(d3$Y, d3$R, p = p, q = q, vce = "hc1", bwselect = "mserd")$bws
    expect_equal(jb$h, rb[1, 1], tolerance = 1e-8, info = paste(xd, p))
    expect_equal(jb$b, rb[1, 3], tolerance = 1e-8, info = paste(xd, p))
  }
})

test_that("B.4 joint rule, balanced PC panel — d, b and h are rdbwselect on the combined outcome at the same h_rot", {
  d <- .jp_panel(); pl <- .jp_plist(d)
  yc <- 0
  for (k in names(.jp_coef)) yc <- yc + .jp_coef[[k]] * pl[[k]]$y[match(pl[["1"]]$id, pl[[k]]$id)]
  xc <- pl[["1"]]$x
  # rdbwselect's own rule of thumb on the single running variable (the pooled one differs by the
  # n - 1 of sd over 3n rows, ~1e-4); pass it to both
  iqr <- unname(quantile(xc, 0.75, type = 2) - quantile(xc, 0.25, type = 2))
  hr  <- 2.576 * min(sd(xc), iqr / 1.349) * length(xc)^(-1 / 5)
  pil <- .bw_pilots(pl, .jp_coef, scheme = "pc", h_rot = hr)
  jb  <- .bw_joint(pl, .jp_coef, "3", scheme = "pc",
                   constants = .bw_constants(pl, .jp_coef, scheme = "pc", b = pil$b, h_rot = hr))
  rb  <- rdrobust::rdbwselect(yc, xc, vce = "hc1", bwselect = "mserd")$bws
  expect_equal(pil$b, rb[1, 3], tolerance = 1e-8)
  expect_equal(jb$h, rb[1, 1], tolerance = 1e-8)
  # d enters only through b; rdbwselect does not return it, so check it through the b stage:
  # at a d moved by 10% the b stage moves
  st <- .bw_stage(pl, .jp_coef, "pc", o = 2L, nu = 2L, oB = 3L, h_V = hr, h_B = 1.1 * pil$d,
                  reg_const = 3)
  expect_false(isTRUE(all.equal(min(st$bw, max(abs(xc))), pil$b, tolerance = 1e-6)))
})

test_that("B.4 joint rule — the h stage of .bw_stage at b is .bw_joint's h (two implementations of eq:common_h_opt)", {
  for (sc in c("cs", "pc", "pv")) {
    d <- .jp_panel(); pl <- .jp_plist(d)
    if (sc == "cs") for (k in names(pl)) pl[[k]]$id <- paste(k, pl[[k]]$id)
    if (sc == "pv") for (k in names(pl)) pl[[k]]$x <- pl[[k]]$x + rnorm(nrow(pl[[k]]), 0, 0.05)
    jb <- .bw_joint(pl, .jp_coef, "3", scheme = sc)
    st <- .bw_stage(pl, .jp_coef, sc, o = 1L, nu = 0L, oB = 2L, h_V = jb$h_rot, h_B = jb$b,
                    reg_const = 3)
    expect_equal(st$bw, jb$h, tolerance = 1e-8, info = sc)
  }
})

test_that("B.4 joint rule — regularization moves b and h but not d; reg_const scales it", {
  d <- .jp_panel(); pl <- .jp_plist(d)
  on  <- .bw_pilots(pl, .jp_coef, scheme = "pc")
  off <- .bw_pilots(pl, .jp_coef, scheme = "pc", regularize = FALSE)
  expect_identical(on$d, off$d)
  expect_equal(off$stage_b$R, 0)
  expect_true(on$stage_b$R > 0)
  expect_true(on$b < off$b)
  r6 <- .bw_pilots(pl, .jp_coef, scheme = "pc", reg_const = 6)
  expect_equal(r6$stage_b$R, 2 * on$stage_b$R, tolerance = 1e-12)
})

test_that("B.4 joint rule through rddid() — one h and one b in every period, the bias correction at that b", {
  d <- .jp_panel()
  f <- rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, comparisons = c(1, 2),
             bwselect = "joint", scheme = "pc")
  bw <- f$bandwidth
  expect_true(all(bw$h_by_period == bw$h))
  expect_true(all(bw$b_by_period == bw$b))
  expect_true(is.finite(bw$d) && bw$d > 0)
  for (k in names(f$fits)) expect_identical(unname(f$fits[[k]]$b), bw$b)
})
