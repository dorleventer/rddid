# Bandwidth rules and sampling schemes

This article covers the
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
arguments that set the bandwidths, the main bandwidth `h` (for the point
estimate) and the pilot bandwidth `b` (for the bias correction), the
standard error under the three sampling schemes, and the choice of the
comparison-period weights. It uses `rddid_sim`, a simulated panel with
RD period 3, comparison periods 1 and 2 and a true effect of 1; [Get
started](https://dorleventer.github.io/rddid/articles/rddid-estimation.md)
introduces the data and the estimator.

``` r

library(rddid)
```

## Bandwidth rules

`bwselect` picks the rule. Whatever the rule,
`fit$bandwidth$h_by_period` and `fit$bandwidth$b_by_period` hold the
bandwidths used in each period, named by period (stored with the RD
period first; the printout and
[`summary()`](https://rdrr.io/r/base/summary.html) list them in time
order).

### The common bandwidth: `bwselect = "joint"` (the default)

``` r

fit_joint <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_joint$bandwidth$h_by_period
#>         3         1         2 
#> 0.2672096 0.2672096 0.2672096
fit_joint$bandwidth$b_by_period
#>         3         1         2 
#> 0.3950840 0.4101520 0.3868401
```

One main bandwidth for every period, chosen to minimize the asymptotic
mean squared error of the RD-DID estimate, not of each period’s jump.
The bias and variance of each period’s jump are estimated at that
period’s own CCT pilot bandwidths, so the rule does not depend on which
period is labeled `t_rd`. The pilot bandwidth of each period keeps that
period’s CCT ratio of pilot to main bandwidth. For year 3, from its own
CCT pair:

``` r

y3  <- subset(rddid_sim, year == 3)
bw3 <- rd_bw_cct(y3$Y, y3$R)
fit_joint$bandwidth$h * bw3[["b"]] / bw3[["h"]]
#> [1] 0.395084
```

which is the year-3 entry of `b_by_period` above.

### Per-period CCT bandwidths: `bwselect = "cct"`

``` r

fit_cct <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                 bwselect = "cct")
fit_cct$bandwidth$h_by_period
#>         3         1         2 
#> 0.4140584 0.2637172 0.3255516
```

Each period gets its own MSE-optimal bandwidth from
[`rdrobust::rdbwselect`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
(Calonico, Cattaneo and Titiunik, 2014), computed by
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md),
ignoring the other periods. Each bandwidth is optimal for its own jump;
the common bandwidth is optimal for the difference of jumps that
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
reports.

### Period-specific bandwidths for the aggregate: `bwselect = "iter"`

``` r

fit_iter <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                  bwselect = "iter")
fit_iter$bandwidth$h_by_period
#>         3         1         2 
#> 0.2913714 0.2348810 0.2916165
```

One bandwidth per period, found by coordinate descent on the asymptotic
mean squared error of the RD-DID estimate. The paper does not use this
rule; it is kept for simulations. `start` sets where the descent starts:
`"hstar"` (the default) starts every period at the common bandwidth,
`"cct"` starts each period at its own CCT bandwidth, and a named vector
gives one starting value per period.

``` r

fit_iter_cct <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                      bwselect = "iter", start = "cct")
fit_iter_cct$bandwidth$h_by_period
#>         3         1         2 
#> 0.2913710 0.2348758 0.2916176
```

### Fixed bandwidths

``` r

fit_fixed <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                   h = 0.3)
fit_fixed$bandwidth$method
#> [1] "fixed"
fit_fixed$bandwidth$h_by_period
#>   3   1   2 
#> 0.3 0.3 0.3
```

A given `h` is used in every period and bypasses bandwidth selection;
`b` defaults to `h`.

`regularize` (default `TRUE`) adds an `rdrobust`-style regularization
term to the squared bias in the `"joint"` and `"iter"` rules, so that a
near-zero estimated curvature does not blow the bandwidth up;
`reg_const` (default 3) is its constant.

### The rules side by side

The four fits above, with the main bandwidth in each year and the two
estimates (standard errors in parentheses):

| rule          | h, year 1 | h, year 2 | h, year 3 | Conventional  | Robust        |
|:--------------|----------:|----------:|----------:|:--------------|:--------------|
| joint         |     0.267 |     0.267 |     0.267 | 1.093 (0.126) | 1.141 (0.149) |
| cct           |     0.264 |     0.326 |     0.414 | 1.040 (0.117) | 1.096 (0.141) |
| iter          |     0.235 |     0.292 |     0.291 | 1.100 (0.125) | 1.141 (0.147) |
| fixed h = 0.3 |     0.300 |     0.300 |     0.300 | 1.068 (0.120) | 1.118 (0.161) |

The rules give different bandwidths, from 0.235 to 0.414, and, on these
data, conventional estimates from 1.040 to 1.100, with standard errors
of about 0.12; the truth is 1.

## Sampling schemes

The sampling scheme describes how the data were collected:

- `"cs"`, repeated cross-section: different units in each period;
- `"pc"`, panel, running variable fixed over time: the same units, each
  on the same side of the cutoff in every period;
- `"pv"`, panel, running variable varies over time: some units change
  side between periods (switchers).

With `scheme = "auto"` (the default)
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
reads the scheme off the data: no `id`, or no unit observed in more than
one period, gives `"cs"`; units observed in more than one period, each
on the same side of the cutoff in every period, give `"pc"`; any unit on
different sides in different periods gives `"pv"`. Without `id`,
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
treats every row as a different unit and says so; a panel needs `id`.

Three data sets, one per scheme: a repeated cross-section made from
`rddid_sim` by keeping each unit in one year only, `rddid_sim` itself,
and `rddid_sim_pv`, whose running variable drifts between years.

``` r

sim_cs <- rddid_sim[rddid_sim$id %% 3 == rddid_sim$year %% 3, ]  # each unit in one year
fit_cs <- rddid(sim_cs, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_pc <- fit_joint
fit_pv <- rddid(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
c(cs = fit_cs$scheme, pc = fit_pc$scheme, pv = fit_pv$scheme)
#>   cs   pc   pv 
#> "cs" "pc" "pv"
```

[`summary()`](https://rdrr.io/r/base/summary.html) ends with the robust
standard error under each scheme:

``` r

summary(fit_pv)
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Bandwidth: common h = 0.2367 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.3889, 2 = 0.3847, 3 = 0.3933
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9654     0.1933    4.99   <0.001   [0.5866, 1.3442]
#>   Robust (bias-corrected)      0.9823     0.2292    4.29   <0.001   [0.5332, 1.4315]
#> 
#>   Per-period local-linear fits, in time order (estimate = sum of coef x jump):
#>   period   role          coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison    -0.5   1000   0.2367   0.3889     0.7519    0.1593     0.7893    0.1893
#>   2        comparison    -0.5   1000   0.2367   0.3847     0.7327    0.1824     0.7974    0.2231
#>   3        RD               1   1000   0.2367   0.3933     1.7077    0.1559     1.7757    0.1818
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pv"; the others are for comparison):
#>     cross-section 0.2333   panel, no unit changes side 0.2254   panel, some change side 0.2292
```

The same line for the three data sets:

``` r

se_tab <- rbind(fit_cs$estimates["Robust", c("se_cs", "se_pc", "se_pv")],
                fit_pc$estimates["Robust", c("se_cs", "se_pc", "se_pv")],
                fit_pv$estimates["Robust", c("se_cs", "se_pc", "se_pv")])
dimnames(se_tab) <- list(c("cross-section data", "panel data, R fixed", "panel data, R varies"),
                         c("cross-section s.e.", "panel, R fixed s.e.", "panel, R varies s.e."))
knitr::kable(se_tab, digits = 3)
```

|  | cross-section s.e. | panel, R fixed s.e. | panel, R varies s.e. |
|:---|---:|---:|---:|
| cross-section data | 0.388 | 0.388 | 0.388 |
| panel data, R fixed | 0.239 | 0.149 | 0.149 |
| panel data, R varies | 0.233 | 0.225 | 0.229 |

- **Cross-section data.** No unit appears in two periods, so the jumps
  of different periods are independent and the three formulas coincide
  (0.388).
- **Panel data, running variable fixed** (`rddid_sim`). Every unit
  enters every year’s fit, on the same side of the cutoff. Its unit
  effect makes the yearly jumps positively correlated, and because the
  comparison jumps are subtracted, that correlation lowers the variance.
  The panel standard error is 0.149; the cross-section one, which
  ignores the correlation, is 0.239. With nobody changing side, the
  formula for a moving running variable reduces to the one for a fixed
  running variable.
- **Panel data, running variable varies** (`rddid_sim_pv`). Units move,
  so a unit can be above the cutoff in one year and below it in another.
  The cross-period covariances are then small (in large samples they
  vanish relative to the variances); here the three standard errors
  range from 0.225 to 0.233.

The scheme describes how the data were sampled; the printout reports the
column for that scheme, and the other columns show the other formulas on
the same data. Setting `scheme` replaces the detection. The scheme sets
which standard error is reported and, under the `"joint"` (default) and
`"iter"` bandwidth rules, it also enters the bandwidth: these rules
weigh bias against the asymptotic variance of the estimate, and that
variance has a cross-period covariance term under `"pc"` only (the term
is zero under `"cs"` and negligible in large samples under `"pv"`). On
`rddid_sim`:

``` r

fit_as_cs <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                   scheme = "cs")
#> rddid(): scheme = "cs" as requested; the data look like "pc" (panel, no unit changes side of the cutoff).
rbind(pc = c(h = fit_pc$bandwidth$h, coef(fit_pc)),
      cs = c(h = fit_as_cs$bandwidth$h, coef(fit_as_cs)))
#>            h Conventional   Robust
#> pc 0.2672096     1.092709 1.141380
#> cs 0.3228792     1.047013 1.082019
```

Treating the panel as a cross-section changes the common bandwidth and
so the estimate. Only with a fixed `h` or `bwselect = "cct"` do the
bandwidths not depend on the scheme; then the scheme leaves the estimate
untouched and changes only the reported standard error.

## Choosing the comparison-period weights

Under `trend = "constant"` any comparison-period weights that sum to one
remove the confounding jump, and under `"linear"` any weights that also
extrapolate the comparison periods to the RD period do. All of them
estimate the same effect; they differ in precision. By default
(`weighting = "ols"`)
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) uses
equal weights under `"constant"` and the least-squares line under
`"linear"`. `weighting = "min_variance"` uses the minimum-variance
weights instead: the allowed weights that make the variance of the
conventional estimate smallest, computed from the estimated variances of
the yearly jumps and, in a panel, their covariances (the same quantities
behind the conventional standard error).

``` r

fit_ols <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_mv  <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                 weighting = "min_variance")
fit_mv$weights
#>         1         2 
#> 0.4523365 0.5476635
c(ols = fit_ols$estimates["Conventional", "se"],
  min_variance = fit_mv$estimates["Conventional", "se"])
#>          ols min_variance 
#>    0.1264043    0.1255210
```

The weights on years 1 and 2 are 0.452 and 0.548, and the standard error
barely moves. In `rddid_sim` the two comparison years are equally noisy
(their jumps have the same variance), so the weights stay near one half.
The tilt toward year 2 comes from the panel: the same units appear every
year, and the estimated year-2 jump covaries more with the year-3 jump
than the year-1 jump does, so weighting it more cancels more of the
year-3 noise. The gain is larger when the comparison periods differ in
precision, for example in sample size near the cutoff; in a repeated
cross-section with a constant trend the minimum-variance weights are
inverse-variance weights.

Under the default bandwidth rule the weights and the common bandwidth
are chosen in three steps: the common bandwidth with the `"ols"`
weights, the minimum-variance weights at that bandwidth, and the
bandwidths (common `h` and pilot `b`) again with those weights (here
0.267 and then 0.270). With `bwselect = "cct"` or a fixed `h` the
bandwidths do not depend on the weights, so the weights are computed
once. The standard errors treat the weights as known; in large samples,
estimating them does not change the distribution of the estimate,
because every allowed set of weights cancels the confounding jump
exactly, so an error in the weights multiplies only the estimation error
of the jumps. `fit_mv$weights_detail` keeps the `"ols"` weights, the
estimated variances and covariances used, and the bandwidths at which
they were estimated. With one comparison period under `"constant"`, or
two under `"linear"`, there is a single set of allowed weights, and
`"min_variance"` returns it.

## References

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295–2326.

Leventer, D. and D. Nevo (2024). Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data. arXiv:2408.05847.
<https://arxiv.org/abs/2408.05847>
