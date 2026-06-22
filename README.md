# rddid

<!-- badges: start -->
<!-- badges: end -->

**Treatment-effect estimation in regression-discontinuity
difference-in-discontinuities (RD-DID) designs.**

Standard RD requires the mean potential outcome to be continuous at the cutoff.
That fails when a *confounding* policy switches at the same threshold as the
treatment of interest, biasing the discontinuity. `rddid` implements the RD-DID
framework of **Leventer and Nevo**: comparison periods — where the confounding
is present at the cutoff but the treatment of interest is uniform — identify and
net out the bias, recovering the causal effect

<!-- math -->
> ATT(t_RD) = D_{t_RD} − Σ_t w_t D_t,

where each `D_t` is a standard local-linear RD discontinuity and `{w_t}` are
trend weights.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("dorleventer/rddid")
```

> **Note.** This is a from-scratch rewrite (v0.1) targeting the current
> estimation framework. The original work-in-progress package is preserved at
> the git tag [`archive/v0`](https://github.com/dorleventer/rddid/releases/tag/archive/v0).

## Quick start

``` r
library(rddid)

# long data: one row per unit-period, with outcome y, running variable x,
# period `time`, and unit `id`.
fit <- rddid(data, y = "y", x = "x", time = "time", id = "id",
             t_rd = 3,                 # period with the RD treatment
             weights = "constant",     # or "linear", or a numeric vector
             bwselect = "joint")       # or "cct"
fit
```

## What it does

**Estimation.** The period-`t_rd` discontinuity is estimated by standard
local-linear RD, conventional and robust bias-corrected (Calonico, Cattaneo and
Titiunik 2014), and aggregated across periods. The per-period engine is
validated against [`rdrobust`](https://github.com/rdpackages/rdrobust) to
machine precision.

**Weights.** `"constant"` (equal weights; a constant confounding trend),
`"linear"` (the line through the comparison discontinuities extrapolated to
`t_rd`), or a custom numeric vector.

**Two bandwidth rules.**

| `bwselect` | Rule |
|---|---|
| `"cct"` | Within-period MSE-optimal bandwidth (Imbens–Kalyanaraman / CCT) applied to each `D_t` on its own. Never degenerates; leaves the cross-period bias cancellation unexploited. |
| `"joint"` | One common bandwidth `h* = (V / B²)^(1/5) n^(-1/5)` minimising the asymptotic MSE of the **aggregate** estimator. The bias constant `B = b_{t_RD} − Σ w_t b_t` is signed, so comparison-period biases can offset the RD-period bias. |

**Three sampling schemes.** Standard errors are reported under repeated
cross-section (`CS`), panel with a time-constant running variable (`PC`), and
panel with a time-varying running variable (`PV`). The applicable scheme is
auto-detected from the id / side structure; all three are always returned.

## Status

v0.1 — estimation and inference complete and verified. Simulation study,
empirical vignette, and inference with data-driven weights are in progress.

## Reference

Leventer, D. and Nevo, D. *Regression discontinuity with a confounded cutoff.*
