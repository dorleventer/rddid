
<!-- README.md is generated from README.Rmd. Edit README.Rmd, then run devtools::build_readme(). -->

# rddid

<!-- badges: start -->

[![R-CMD-check](https://github.com/dorleventer/rddid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dorleventer/rddid/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

A treatment of interest switches on at a cutoff of a running variable in
one period, the **RD period**. A **confounding policy** switches at the
same cutoff, in every period, so the jump in the outcome at the cutoff
in the RD period mixes the treatment effect with the **confounding
jump**. In the **comparison periods** the treatment of interest is
uniform at the cutoff (nobody treated, or everybody treated), so,
provided the treatment of interest has no anticipation or carry-over
effects there (which the paper assumes), the jump there *is* the
confounding jump. `rddid()` estimates the jump in every period by
local-linear RD and subtracts a weighted average of the
comparison-period jumps from the RD-period jump. How the weights are set
is the **confounding-trend assumption**: constant (equal weights) or
linear in time.

Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data.*
[arXiv:2408.05847](https://arxiv.org/abs/2408.05847).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("dorleventer/rddid")
```

## Quick start

`rddid_sim` is a simulated panel: 1,000 units in three years, a
confounding jump of 0.5 in every year, and a treatment effect of 1 in
year 3, the RD period.

``` r
library(rddid)
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, running variable fixed over time: no unit changes side of the cutoff (detected from the data)
#>   Bandwidth: common h = 0.2672 (rule "joint", AMSE-optimal for the aggregate)
#>   Pilot bandwidth b (period = value): 3 = 0.3951, 1 = 0.4102, 2 = 0.3868
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0927     0.1264    8.64   <0.001   [0.8450, 1.3405]
#>   Robust (bias-corrected)      1.1414     0.1494    7.64   <0.001   [0.8486, 1.4342]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

Years 1 and 2 are the comparison periods, with equal weights.
`Conventional` is the local-linear estimate with its conventional
standard error; `Robust (bias-corrected)` is the bias-corrected estimate
with its robust standard error. A plain RD in year 3 would target 1.5,
the effect plus the confounding jump.

When the running variable moves over time, so that units can change side
of the cutoff between periods, four tests check the assumptions this
adds, for example type continuity:

``` r
rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id")
#> Test of a continuous type distribution  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, running variable varies over time: some units change side (detected from the data)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Joint Wald chi-squared(9) = 4.782,  p = 0.853
#>     Period 1: chi-squared(3) = 1.531,  p = 0.675
#>     Period 2: chi-squared(3) = 0.382,  p = 0.944
#>     Period 3: chi-squared(3) = 2.202,  p = 0.532
```

The other three are `rd_compstable()`, `rd_homog()` and
`rd_trendcell()`.

## Learn more

- [Get
  started](https://dorleventer.github.io/rddid/articles/rddid-estimation.html):
  the estimate, its printout, and the four tests.
- [Reference](https://dorleventer.github.io/rddid/reference/index.html):
  every function, with examples.
