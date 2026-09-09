# rddid

Estimation and inference for regression-discontinuity
difference-in-discontinuities (RD-DID) designs, where a confounding
policy switches at the same cutoff as the treatment of interest.
Implements the framework of Leventer and Nevo.

Paper: <https://arxiv.org/abs/2408.05847>

## Installation

``` r

# install.packages("devtools")
devtools::install_github("dorleventer/rddid")
```

## Quick start

The data are long, one row per unit-period, with an outcome, a running
variable, a period and a unit id. Below, periods 1 and 2 are comparison
periods (the confounding discontinuity is 1 in both) and period 3 is the
RD period, where the treatment of interest adds 1 at the cutoff.

``` r

library(rddid)

set.seed(1)
n <- 2000
R <- runif(n, -1, 1)                 # time-invariant running variable, cutoff at 0
u <- rnorm(n, 0, 0.5)                # unit effect
dat <- do.call(rbind, lapply(1:3, function(t) {
  V <- as.integer(R >= 0)            # confounding treatment: sharp RD in every period
  W <- V * (t == 3)                  # treatment of interest: sharp RD in period 3 only
  data.frame(id = seq_len(n), t = t, R = R,
             Y = R + R^2 * (R >= 0) + 1 * V + 1 * W + u + rnorm(n, 0, 0.5))
}))

rddid(dat, y = "Y", x = "R", time = "t", id = "id",
      t_rd = 3, comparisons = c(1, 2), weights = "constant")
#> RD-DID estimate of ATT(t_RD)
#>   RD period: 3   comparison periods: 1, 2
#>   weights: constant [0.5, 0.5]
#>   bandwidth: period-specific joint AMSE (coord. descent, 5 iters)
#>   sampling scheme: PC (auto-detected)
#> 
#>                    Estimate   Std.Err.   95% CI
#>   Conventional      1.09413    0.09746   [  0.90312,   1.28514]
#>   Robust            1.12176    0.11715   [  0.89215,   1.35138]
#> 
#>   SEs by scheme (Robust): CS=0.18902  PC=0.11715  PV=0.11715
```

The vignettes walk through the estimator, its bandwidth rules and
sampling schemes, and the composition validation tests:
<https://dorleventer.github.io/rddid/>.
