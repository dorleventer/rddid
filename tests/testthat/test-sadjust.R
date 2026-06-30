# Tests for rd_sadjust(): within-period composition term S_t (Prop sdisc).
# DGP mirrors test-adjust.R's .adj_draw (P=2, binary partner-side types): the
# unit's type is its side in the partner period, exactly the structure S_t needs.

.sadj_draw <- function(n, delta, a1_lo, a1_hi, a2_lo, a2_hi,
                       rho = 0.6, s = 0.5, sigma = 1, tau = 5) {
  eta <- rnorm(n); e1 <- rnorm(n, 0, sqrt(1 - rho^2)); e2 <- rnorm(n, 0, sqrt(1 - rho^2))
  R1 <- rho * eta + e1 + delta; R2 <- rho * eta + e2
  v1 <- as.numeric(R1 >= 0); v2 <- as.numeric(R2 >= 0)
  a1 <- ifelse(v2 == 1, a1_hi, a1_lo); a2 <- ifelse(v1 == 1, a2_hi, a2_lo)
  Y1 <- s * R1 + a1 * v1 + rnorm(n, 0, sigma)
  Y2 <- s * R2 + a2 * v2 + tau * v2 + rnorm(n, 0, sigma)
  id <- seq_len(n)
  rbind(data.frame(id, time = 1L, x = R1, y = Y1),
        data.frame(id, time = 2L, x = R2, y = Y2))
}

test_that("rd_sadjust returns the documented object structure", {
  set.seed(201)
  d <- .sadj_draw(3000, delta = 1, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_sadjust(d, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 1, t_rd = 2, h = 0.6, se = "none")
  expect_s3_class(r, "rd_sadjust")
  expect_true(all(c("comparison", "period", "role", "S", "S_bc", "S_se",
                    "S_bc_se", "D", "D_bc", "adjusted", "adjusted_bc",
                    "n", "n_control_switchers") %in% names(r$table)))
  # one comparison -> one pair -> a comparison row and an rd row
  expect_equal(nrow(r$table), 2L)
  expect_setequal(r$table$role, c("comparison", "rd"))
})

test_that("S_t is produced at the RD period (control side observed there)", {
  set.seed(202)
  d <- .sadj_draw(3000, delta = 1, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_sadjust(d, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 1, t_rd = 2, h = 0.6, se = "none")
  rd_row <- r$table[r$table$role == "rd", ]
  expect_equal(nrow(rd_row), 1L)
  expect_true(is.finite(rd_row$S))
  expect_true(is.finite(rd_row$adjusted))
})

test_that("binary-type identity S = dpi * (mu_above - mu_below) holds", {
  set.seed(203)
  d  <- .sadj_draw(3000, delta = 1, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  h  <- 0.6
  r  <- rd_sadjust(d, y = "y", x = "x", time = "time", id = "id",
                   comparisons = 1, t_rd = 2, h = h, se = "none")
  # reconstruct S for the comparison row (focal = period 1, partner = period 2)
  d1 <- d[d$time == 1, ]
  side2 <- stats::setNames(as.integer(d[d$time == 2, "x"] >= 0),
                           as.character(d[d$time == 2, "id"]))
  d1$type <- side2[as.character(d1$id)]
  d1 <- d1[stats::complete.cases(d1$y, d1$x, d1$type), ]
  ab <- d1$type == 1; be <- d1$type == 0
  dpi   <- rd_period(d1$type, d1$x, h = h, c = 0, kernel = "triangular")$D
  mu_ab <- rd_period(d1$y[ab], d1$x[ab], h = h, c = 0,
                     kernel = "triangular")$sides[["-"]]$beta0
  mu_be <- rd_period(d1$y[be], d1$x[be], h = h, c = 0,
                     kernel = "triangular")$sides[["-"]]$beta0
  S_manual <- dpi * (mu_ab - mu_be)
  S_pkg <- r$table$S[r$table$role == "comparison"]
  expect_equal(S_pkg, S_manual, tolerance = 1e-8)
})

test_that("cluster bootstrap returns finite, positive SEs", {
  set.seed(204)
  d <- .sadj_draw(2000, delta = 1, a1_lo = 1, a1_hi = 3, a2_lo = 1, a2_hi = 3)
  r <- rd_sadjust(d, y = "y", x = "x", time = "time", id = "id",
                  comparisons = 1, t_rd = 2, h = 0.7, se = "bootstrap", B = 80)
  expect_true(all(is.finite(r$table$S_se)))
  expect_true(all(r$table$S_se > 0))
  expect_true(is.numeric(r$n_boot_ok) && r$n_boot_ok > 0)
})
