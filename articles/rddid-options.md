# Bandwidth rules and sampling schemes

This article covers the
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
arguments that set the bandwidths, the main bandwidth `h` (for the point
estimate) and the pilot bandwidth `b` (for the bias correction), and the
standard error under the three sampling schemes. It uses `rddid_sim`, a
simulated panel with RD period 3, comparison periods 1 and 2 and a true
effect of 1; [Get
started](https://dorleventer.github.io/rddid/articles/rddid-estimation.md)
introduces the data and the estimator.

``` r

library(rddid)
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
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

``` r

estse <- function(f, row) paste0(fmt(f$estimates[row, "est"]), " (", fmt(f$estimates[row, "se"]), ")")
fits <- list(joint = fit_joint, cct = fit_cct, iter = fit_iter, "fixed h = 0.3" = fit_fixed)
tab <- data.frame(
  rule = names(fits),
  t(sapply(fits, function(f) f$bandwidth$h_by_period[c("1", "2", "3")])),
  Conventional = sapply(fits, estse, row = "Conventional"),
  Robust = sapply(fits, estse, row = "Robust"),
  check.names = FALSE, row.names = NULL
)
names(tab)[2:4] <- c("h, year 1", "h, year 2", "h, year 3")
knitr::kable(tab, digits = 3)
```

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
treats every row as a different unit and says so; for a panel, pass
`id`.

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
#>   Sampling scheme: panel, running variable varies over time: some units change side (detected from the data)
#>   Bandwidth: common h = 0.2367 (rule "joint", AMSE-optimal for the aggregate)
#>   Pilot bandwidth b (period = value): 1 = 0.3889, 2 = 0.3847, 3 = 0.3933
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9654     0.1933    4.99   <0.001   [0.5866, 1.3442]
#>   Robust (bias-corrected)      0.9823     0.2292    4.29   <0.001   [0.5332, 1.4315]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
#> 
#>   Per-period local-linear fits, in time order (estimate = sum of coef x jump):
#>   period   role          coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison    -0.5   1000   0.2367   0.3933     1.7077    0.1559     1.7757    0.1818
#>   2        comparison    -0.5   1000   0.2367   0.3889     0.7519    0.1593     0.7893    0.1893
#>   3        RD               1   1000   0.2367   0.3847     0.7327    0.1824     0.7974    0.2231
#> 
#>   Robust s.e. under each sampling scheme:  cross-section 0.2333   panel, fixed R 0.2254   panel, varying R 0.2292
#>   (the printed s.e. is the one for scheme "pv"; the others are shown for comparison)
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

The scheme describes how the data were sampled; it is not a choice among
the three columns. Setting `scheme` replaces the detection. The scheme
sets which standard error is reported and, under the `"joint"` (default)
and `"iter"` bandwidth rules, it also enters the bandwidth: these rules
weigh bias against the asymptotic variance of the estimate, and that
variance has a cross-period covariance term under `"pc"` only (the term
is zero under `"cs"` and negligible in large samples under `"pv"`). On
`rddid_sim`:

``` r

fit_as_cs <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                   scheme = "cs")
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

## References

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295–2326.

Leventer, D. and D. Nevo (2024). Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data. arXiv:2408.05847.
<https://arxiv.org/abs/2408.05847>
