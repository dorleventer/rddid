# Bandwidth rules and sampling schemes

This article covers the
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
arguments that set the bandwidths, the main bandwidth `h` (for the point
estimate) and the bias-correction bandwidth `b`, the standard error
under the three sampling schemes, and the choice of the
comparison-period weights. It uses `didisc_sim`, a simulated panel with
RD period 3, comparison periods 1 and 2 and a true effect of 1; [Get
started](https://dorleventer.github.io/didisc/articles/didisc-estimation.md)
introduces the data and the estimator.

``` r

library(didisc)
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

fit_joint <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_joint$bandwidth$h_by_period
#>         3         1         2 
#> 0.3985921 0.3985921 0.3985921
fit_joint$bandwidth$b_by_period
#>       3       1       2 
#> 0.59356 0.59356 0.59356
```

One main bandwidth for every period, chosen to minimize the asymptotic
mean squared error of the diff-in-disc estimate, not of each period’s
jump. As in rdrobust, everything that depends only on the running
variable and the residuals is estimated at one rule-of-thumb bandwidth
shared by every period (`fit_joint$bandwidth$h_rot`): each period’s
variance, the covariances between periods in a panel, and the design
part of each period’s bias. Only the curvature of the conditional mean,
the other part of the bias, needs its own fit, at the bias-correction
bandwidth `b`. `b` is chosen the same way, for the same difference of
jumps one derivative up, with its own curvature from fits at a bandwidth
`d` (`fit_joint$bandwidth$d`), as rdrobust chooses its `b` for one
period. So one `b` serves every period:

``` r

c(h = fit_joint$bandwidth$h, b = fit_joint$bandwidth$b, d = fit_joint$bandwidth$d)
#>         h         b         d 
#> 0.3985921 0.5935600 0.9750591
```

These inputs to the rule are the same whichever period is the RD period.
With a single period the rule is
[`rdrobust::rdbwselect()`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
(`bwselect = "mserd"`, `vce = "hc1"`); in a panel whose running variable
does not move, it is that selector applied to the combined outcome, the
RD-period outcome minus the weighted comparison outcomes.

### Per-period CCT bandwidths: `bwselect = "cct"`

``` r

fit_cct <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                 bwselect = "cct")
fit_cct$bandwidth$h_by_period
#>         3         1         2 
#> 0.4140584 0.2637172 0.3255516
```

Each period gets its own MSE-optimal bandwidth from
[`rdrobust::rdbwselect`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
(Calonico, Cattaneo and Titiunik, 2014), computed by
[`rd_bw_cct()`](https://dorleventer.github.io/didisc/reference/rd_bw_cct.md),
ignoring the other periods. Each bandwidth is optimal for its own jump;
the common bandwidth is optimal for the difference of jumps that
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
reports.

### Fixed bandwidths

``` r

fit_fixed <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
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
term to the squared bias in the `"joint"` rule, so that a near-zero
estimated curvature does not blow the bandwidth up; `reg_const` (default
3) is its constant.

### The rules side by side

The three fits above, with the main bandwidth in each year and the two
estimates (standard errors in parentheses):

| rule          | h, year 1 | h, year 2 | h, year 3 | Conventional  | Robust        |
|:--------------|----------:|----------:|----------:|:--------------|:--------------|
| joint         |     0.399 |     0.399 |     0.399 | 1.008 (0.106) | 1.025 (0.126) |
| cct           |     0.264 |     0.326 |     0.414 | 1.040 (0.117) | 1.096 (0.141) |
| fixed h = 0.3 |     0.300 |     0.300 |     0.300 | 1.068 (0.120) | 1.118 (0.161) |

The rules give different bandwidths, from 0.264 to 0.414, and, on these
data, conventional estimates from 1.008 to 1.068, with standard errors
of about 0.11; the truth is 1.

## Sampling schemes

The sampling scheme describes how the data were collected:

- `"cs"`, repeated cross-section: different units in each period;
- `"pc"`, panel, running variable fixed over time: the same units, each
  on the same side of the cutoff in every period;
- `"pv"`, panel, running variable varies over time: some units change
  side between periods (switchers).

With `scheme = "auto"` (the default)
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
reads the scheme off the data: no `id`, or no unit observed in more than
one period, gives `"cs"`; units observed in more than one period, each
on the same side of the cutoff in every period, give `"pc"`; any unit on
different sides in different periods gives `"pv"`. Without `id`,
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
treats every row as a different unit and says so; a panel needs `id`.

Three data sets, one per scheme: a repeated cross-section made from
`didisc_sim` by keeping each unit in one year only, `didisc_sim` itself,
and `didisc_sim_pv`, whose running variable drifts between years.

``` r

sim_cs <- didisc_sim[didisc_sim$id %% 3 == didisc_sim$year %% 3, ]  # each unit in one year
fit_cs <- didisc(sim_cs, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_pc <- fit_joint
fit_pv <- didisc(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
c(cs = fit_cs$scheme, pc = fit_pc$scheme, pv = fit_pv$scheme)
#>   cs   pc   pv 
#> "cs" "pc" "pv"
```

[`summary()`](https://rdrr.io/r/base/summary.html) ends with the robust
standard error under each scheme:

``` r

summary(fit_pv)
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Bandwidth: common h = 0.2417 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Bias-correction bandwidth b = 0.3683 (also chosen for the diff-in-disc estimate)
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9687     0.1910    5.07   <0.001   [0.5943, 1.3431]
#>   Robust (bias-corrected)      0.9735     0.2331    4.18   <0.001   [0.5166, 1.4305]
#> 
#>   Per-period local-linear fits:
#>                                                        Conventional         Robust (bias-corrected)
#>   Period   Role             n        h        b          Jump   Std. err.          Jump   Std. err.
#>   1        comparison    1000   0.2417   0.3683        0.7501      0.1574        0.7989      0.1921
#>   2        comparison    1000   0.2417   0.3683        0.7263      0.1797        0.8056      0.2263
#>   3        RD            1000   0.2417   0.3683        1.7069      0.1546        1.7758      0.1849
#> 
#>   Std. err. of the estimate under each sampling scheme (the estimate table at the top uses PV):
#>                                                    Conventional   Robust (bias-corrected)
#>   RCS  repeated cross-section                            0.1953                    0.2371
#>   PC   panel, no unit changes side of the cutoff         0.1882                    0.2291
#>   PV   panel, some units change side of the cutoff       0.1910                    0.2331
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
| cross-section data | 0.409 | 0.409 | 0.409 |
| panel data, R fixed | 0.200 | 0.126 | 0.126 |
| panel data, R varies | 0.237 | 0.229 | 0.233 |

- **Cross-section data.** No unit appears in two periods, so the jumps
  of different periods are independent and the three formulas coincide
  (0.409).
- **Panel data, running variable fixed** (`didisc_sim`). Every unit
  enters every year’s fit, on the same side of the cutoff. Its unit
  effect makes the yearly jumps positively correlated, and because the
  comparison jumps are subtracted, that correlation lowers the variance.
  The panel standard error is 0.126; the cross-section one, which
  ignores the correlation, is 0.200. With nobody changing side, the
  formula for a moving running variable reduces to the one for a fixed
  running variable.
- **Panel data, running variable varies** (`didisc_sim_pv`). Units move,
  so a unit can be above the cutoff in one year and below it in another.
  The cross-period covariances are then small (in large samples they
  vanish relative to the variances); here the three standard errors
  range from 0.229 to 0.237.

The scheme describes how the data were sampled; the printout reports the
column for that scheme, and the other columns show the other formulas on
the same data. Setting `scheme` replaces the detection. The scheme sets
which standard error is reported and, under the `"joint"` (default)
bandwidth rule, it also enters the bandwidth: this rule weighs bias
against the asymptotic variance of the estimate, and that variance has a
cross-period covariance term under `"pc"` only (the term is zero under
`"cs"` and negligible in large samples under `"pv"`). On `didisc_sim`:

``` r

fit_as_cs <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                   scheme = "cs")
#> didisc(): scheme = "cs" as requested; the data look like "pc" (panel, no unit changes side of the cutoff).
rbind(pc = c(h = fit_pc$bandwidth$h, coef(fit_pc)),
      cs = c(h = fit_as_cs$bandwidth$h, coef(fit_as_cs)))
#>            h Conventional   Robust
#> pc 0.3985921     1.007537 1.024654
#> cs 0.4196584     1.004214 1.021247
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
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
uses equal weights under `"constant"` and the least-squares line under
`"linear"`. `weighting = "min_variance"` uses the minimum-variance
weights instead: the allowed weights that make the variance of the
estimate smallest, computed from the estimated variances of the yearly
jumps and, in a panel, their covariances. The conventional estimate gets
the weights that minimize its variance; the bias-corrected (robust)
estimate gets its own, from the variances and covariances of the
bias-corrected jumps.

``` r

fit_ols <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_mv  <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                 weighting = "min_variance")
fit_mv$weights                      # the conventional estimate's weights
#>         1         2 
#> 0.4794084 0.5205916
fit_mv$weights_detail$weights_bc    # the bias-corrected estimate's weights
#>         1         2 
#> 0.4652808 0.5347192
se <- cbind(ols = fit_ols$estimates$se, min_variance = fit_mv$estimates$se)
rownames(se) <- rownames(fit_ols$estimates)
se                                  # standard errors under the two weightings
#>                    ols min_variance
#> Conventional 0.1061519    0.1047821
#> Robust       0.1256024    0.1241868
```

The conventional weights on years 1 and 2 are 0.479 and 0.521, the
bias-corrected ones 0.465 and 0.535; the standard errors barely move. In
`didisc_sim` the two comparison years are equally noisy (their jumps
have the same variance), so the weights stay near one half. The tilt
toward year 2 comes from the panel: the same units appear every year,
and the estimated year-2 jump covaries more with the year-3 jump than
the year-1 jump does, so weighting it more cancels more of the year-3
noise. The gain is larger when the comparison periods differ in
precision, for example in sample size near the cutoff; in a repeated
cross-section with a constant trend the minimum-variance weights are
inverse-variance weights.

Under the default bandwidth rule the weights and the common bandwidth
are chosen in three steps: the common bandwidth with the `"ols"`
weights, the minimum-variance weights at that bandwidth, and the
bandwidths (common `h` and bias-correction `b`) again with those weights
(here 0.399 and then 0.409). With `bwselect = "cct"` or a fixed `h` the
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
