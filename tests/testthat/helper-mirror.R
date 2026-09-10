# Helpers for the mirror-invariance tests (dev/atu_estimand_plan.md, S1, S4
# items 1-3). Mirroring the running variable (x -> -x, cutoff 0) is applied by
# hand here, independent of any package `estimand` argument (there is none in
# this worktree).

mirror_x <- function(d, x = "R") {
  d[[x]] <- -d[[x]]
  d
}

# ---- dgp_a: verbatim from vignettes/rddid-estimation.Rmd, chunk `dgp` -------
# 3 periods, time-invariant running variable R, PC (panel, constant-side) scheme.
m_a <- function(r, theta) r + (theta / 2) * r^2 * (r >= 0)
dgp_a <- function(n = 2000, alpha = c(1, 1, 1), theta = c(2, 2, 2), tau = 1, seed = 1) {
  set.seed(seed)
  R <- runif(n, -1, 1)                  # time-invariant running variable, cutoff at 0
  u <- rnorm(n, 0, 0.5)                 # unit effect
  do.call(rbind, lapply(1:3, function(t) {
    V <- as.integer(R >= 0)             # confounding treatment: sharp RD every period
    W <- V * (t == 3)                    # treatment of interest: sharp RD in the RD period only
    data.frame(id = seq_len(n), t = t, R = R,
               Y = m_a(R, theta[t]) + alpha[t] * V + tau * W + u + rnorm(n, 0, 0.5))
  }))
}

# ---- dgp_b: verbatim from vignettes/rddid-validation-tests.Rmd, chunk `dgp` -
# 2 or 3 periods, time-varying running variable R, PV (panel, side-switching)
# scheme.
m_b <- function(r) r + r^2 * (r >= 0)
dgp_b <- function(n = 4000, d = 0, alpha = c(1, 1), kappa = 0, sort_p = 0, sort_w = 0.5,
                  tau = 1, periods = 2, gamma = 0, seed = 1) {
  set.seed(seed)
  eta <- rnorm(n)                                          # latent level of the running variable
  R <- sapply(seq_len(periods), function(t) eta + d * (t - 1) + rnorm(n, 0, 0.5))
  if (sort_p > 0) {                                        # sorting: some units just below the RD-period
    V1 <- R[, 1] >= 0; r2 <- R[, periods]                  # cutoff that were above in period 1 move above
    mover <- r2 > -sort_w & r2 < 0 & V1 & runif(n) < sort_p
    R[mover, periods] <- -R[mover, periods]
  }
  V <- R >= 0
  do.call(rbind, lapply(seq_len(periods), function(t) {
    other <- if (periods == 2) V[, 3 - t] else V[, periods]  # the unit's type: side in the other period
    a <- alpha[other + 1] + gamma * t                        # within-type confounding discontinuity
    data.frame(id = seq_len(n), t = t, R = R[, t],
               Y = m_b(R[, t]) + a * V[, t] + tau * V[, t] * (t == periods) + kappa * eta + rnorm(n, 0, 0.5))
  }))
}

S1   <- dgp_b(d = 0.5)          # drift: composition changes across periods (PV)
S0_3 <- dgp_b(periods = 3)      # three periods, everything holds (PV)
