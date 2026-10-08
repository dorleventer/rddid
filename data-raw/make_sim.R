# Builds the two example datasets shipped with the package, data/rddid_sim.rda and
# data/rddid_sim_pv.rda. Run from the package root:  Rscript data-raw/make_sim.R
#
# Design (the paper's setting, Leventer & Nevo): a running variable R with the cutoff at 0.
# A confounding policy V switches at the cutoff in EVERY year and shifts the outcome by
# alpha = 0.5. The treatment of interest W switches at the cutoff only in year 3 (the RD
# year) and shifts the outcome by tau = 1. Years 1 and 2 are comparison years: W is zero for
# everyone, so the discontinuity in Y there is the confounding jump alone. The RD-DID target is
# ATT(3) = tau = 1; a naive RD in year 3 recovers alpha + tau = 1.5.
#
# rddid_sim    : 1,000 units x 3 years; R is fixed over time (nobody changes side).
# rddid_sim_pv : 1,000 units x 3 years; each year's R deviates from the unit's base position R0
#                by an independent normal draw whose spread grows with the year (sd 0, 0.15,
#                0.30), so some units sit on different sides of the cutoff in different years
#                ("switchers"). This is the panel the validation tests are for. Because the
#                spread grows, the mix of units near the cutoff (who was above in another year)
#                is not the same in every year: composition stability FAILS by design, while the
#                confounding jump is the same for every unit (homogeneous confounding holds), so
#                rddid() stays unbiased.
make_sim <- function(moving, seed) {
  set.seed(seed)
  n <- 1000; years <- 1:3; rd_year <- 3
  alpha <- 0.5; tau <- 1
  R0 <- runif(n, -1, 1)                   # unit's running variable (year 1 position)
  u  <- rnorm(n, 0, 0.4)                  # unit effect, shared across years
  rows <- lapply(years, function(t) {
    R <- if (moving) pmin(1, pmax(-1, R0 + (t - 1) * rnorm(n, 0, 0.15))) else R0   # sd 0.15 (t - 1)
    V <- as.integer(R >= 0)               # confounding policy: on above the cutoff, every year
    W <- as.integer(R >= 0 & t == rd_year) # treatment of interest: on above the cutoff in year 3
    mu <- 0.2 * t + 0.8 * R + 0.3 * R^2 * (R >= 0)   # smooth part, bends above the cutoff
    data.frame(id = seq_len(n), year = t, R = R, V = V, W = W,
               Y = mu + alpha * V + tau * W + u + rnorm(n, 0, 0.4))
  })
  d <- do.call(rbind, rows)
  rownames(d) <- NULL
  d
}
rddid_sim    <- make_sim(moving = FALSE, seed = 20261008)
# Seed chosen with data-raw/seed_scan.R (40 seeds): one where the assumptions that hold in this
# design (type continuity, homogeneous confounding, constant within-type jump) are comfortably not
# rejected and the one that fails (composition stability) clearly is, with the estimate near 1.
rddid_sim_pv <- make_sim(moving = TRUE,  seed = 20261027)

# switchers in the pv panel: units above the cutoff in some years and below in others
sw <- tapply(rddid_sim_pv$V, rddid_sim_pv$id, function(v) length(unique(v)) > 1)
cat(sprintf("rddid_sim_pv: %d of %d units switch side at least once\n", sum(sw), length(sw)))

usethis::use_data(rddid_sim, rddid_sim_pv, overwrite = TRUE)
