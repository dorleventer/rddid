
<!-- README.md is generated from README.Rmd. Edit README.Rmd, then run devtools::build_readme(). -->

# didisc

<!-- badges: start -->

[![R-CMD-check](https://github.com/dorleventer/didisc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dorleventer/didisc/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

A treatment of interest switches on at a cutoff of a running variable in
one period, the **RD period**. A **confounding policy** switches at the
same cutoff, in every period, so the jump in the outcome at the cutoff
in the RD period mixes the treatment effect with the **confounding
jump**. In the **comparison periods** the treatment of interest is
uniform at the cutoff (nobody treated, or everybody treated), so,
provided the treatment of interest has no anticipation or carry-over
effects there (which the paper assumes), the jump there *is* the
confounding jump. `didisc()` estimates the jump in every period by
local-linear RD and subtracts a weighted average of the
comparison-period jumps from the RD-period jump. How the weights are set
is the **confounding-trend assumption**: constant (equal weights) or
linear in time.

Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data.*
[arXiv:2408.05847](https://arxiv.org/abs/2408.05847).
`citation("didisc")` returns this reference, with a BibTeX entry.

Compared with a single `rdrobust` fit, `didisc()` runs one local-linear
RD per period with the same conventions (`rd_period()` reproduces
`rdrobust` at a given bandwidth) and then nets the comparison-period
jumps out of the RD-period jump, with standard errors that account for
the same units appearing in several periods and bandwidth rules chosen
for the combined estimate.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("dorleventer/didisc")
```

## Quick start

`didisc_sim` is a simulated panel: 1,000 units in three years, a
confounding jump of 0.5 in every year, and a treatment effect of 1 in
year 3, the RD period.

``` r
library(didisc)
fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2672 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4102, 2 = 0.3868, 3 = 0.3951
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0927     0.1264    8.64   <0.001   [0.8450, 1.3405]
#>   Robust (bias-corrected)      1.1414     0.1494    7.64   <0.001   [0.8486, 1.4342]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

The ATT is the effect of the treatment on the units just above the
cutoff, the ones treated in year 3. Years 1 and 2 are the comparison
periods, with equal weights. `h` is the main bandwidth, one for every
year, chosen for the RD-DID estimate; `b` is the pilot bandwidth of the
bias correction in each year. `Conventional` is the local-linear
estimate with its conventional standard error; `Robust (bias-corrected)`
is the bias-corrected estimate with its robust standard error. A plain
RD in year 3 would target 1.5, the effect plus the confounding jump.
[Get
started](https://dorleventer.github.io/didisc/articles/didisc-estimation.html)
reads every line.

When the running variable moves over time, so that units can change side
of the cutoff between periods, four tests check the assumptions this
adds, for example type continuity:

``` r
rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
#> Test of type continuity  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Joint sup-t (t critical values) max |t| = 1.232 over 9 jumps,  p = 0.868  (largest: period:type 1:+-)
#>     Period 1: sup-t (t critical values) max |t| = 1.232 over 3 jumps,  p = 0.500
#>     Period 2: sup-t (t critical values) max |t| = 0.582 over 3 jumps,  p = 0.901
#>     Period 3: sup-t (t critical values) max |t| = 1.048 over 3 jumps,  p = 0.631
#>   All joint statistics: x$joint_tests
```

The other three are `rd_compstable()`, `rd_homog()` and
`rd_trendcell()`. `plot()` draws any fit or test (one RD plot per
period; the binned type shares with their fitted lines; the within-type
jumps), and `plot_switchers()` shows who changes side of the cutoff
between two periods.

## Learn more

- [Get
  started](https://dorleventer.github.io/didisc/articles/didisc-estimation.html):
  the estimate, its printout, and the four tests.
- [Checking the identification
  assumptions](https://dorleventer.github.io/didisc/articles/didisc-validation-tests.html):
  the four tests, one at a time.
- [Plots of the estimate and the
  checks](https://dorleventer.github.io/didisc/articles/didisc-plots.html):
  `plot()` on a fit and on each test.
- [Reference](https://dorleventer.github.io/didisc/reference/index.html):
  every function, with examples.
