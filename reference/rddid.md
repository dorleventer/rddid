# Estimate the effect of a treatment at a cutoff shared with a confounding policy

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

## Usage

``` r
rddid(
  data,
  y,
  x,
  time,
  id = NULL,
  t_rd,
  comparisons = NULL,
  trend = "constant",
  rd_weights = NULL,
  weighting = c("ols", "min_variance"),
  estimand = c("att", "atu"),
  bwselect = c("joint", "cct"),
  h = NULL,
  b = NULL,
  scheme = c("auto", "cs", "pc", "pv"),
  c = 0,
  p = 1L,
  q = 2L,
  kernel = "triangular",
  level = 0.95,
  regularize = TRUE,
  reg_const = 3
)
```

## Arguments

- data:

  a data frame in long format, one row per unit and period: a repeated
  cross-section (different units in each period) or a panel (the same
  units in several periods; it need not be balanced). Only the rows of
  `t_rd` and `comparisons` are used; among them, rows with a missing
  outcome, running variable, period or id are dropped, with a message.

- y:

  name of the outcome column (a string).

- x:

  name of the running-variable column (a string).

- time:

  name of the period column (a string). The column is usually numeric (a
  year); character or factor labels work as well, except under
  `trend = "linear"` (in `rddid()` or
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)),
  where the line is fitted on the period values: these then need to be
  numeric and are the time scale of the line (so 2015, 2017, 2018 are
  unequally spaced);
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) converts
  character labels such as `"2019"`. In `rddid()` the default comparison
  periods are every period other than `t_rd`, in sorted order of the
  values (alphabetical for character or factor labels), the order a
  numeric `trend` follows.

- id:

  name of the unit-identifier column (a string), needed for panel
  standard errors. With `NULL` (default) every row is treated as a
  different unit, which gives repeated cross-section standard errors
  (with a message).

- t_rd:

  the RD period: the value of `time` in which the treatment of interest
  switches on at the cutoff. A vector of several RD periods estimates
  the weighted average of their effects (see "Several RD periods").

- comparisons:

  the comparison periods: values of `time` in which the treatment of
  interest is uniform at the cutoff (nobody treated, or, with
  `estimand = "atu"`, everybody treated). `NULL` (default) uses every
  period other than `t_rd`; pass the periods explicitly when the data
  contain periods that belong in neither.

- trend:

  the confounding-trend assumption, which restricts the
  comparison-period weights (`weighting` picks them within it):
  `"constant"` (default; the confounding jump is the same in every
  period; equal weights under the default `weighting`) or `"linear"`
  (the confounding jump moves linearly in time; needs at least two
  comparison periods). A numeric vector gives the weights directly, one
  per entry of `comparisons` in that order. Weights that sum to one
  cancel a constant confounding jump; `rddid()` warns when they do not.
  With several RD periods a numeric vector applies to every RD period,
  and a numeric matrix gives one row of weights per RD period (rows in
  `t_rd` order or named by period; columns in `comparisons` order or
  named by period).

- rd_weights:

  the weights \\\phi_t\\ of the RD periods in the average, one per entry
  of `t_rd` (in that order, or named by period). `NULL` (default) gives
  equal weights. Any finite weights are allowed; weights that sum to one
  make the estimand a weighted average. Only for several RD periods (a
  single RD period has weight 1).

- weighting:

  how the comparison-period weights are chosen within the `trend`
  assumption: `"ols"` (default) uses the weights `trend` implies (equal
  weights for `"constant"`, the least-squares line for `"linear"`);
  `"min_variance"` uses the minimum-variance weights, the
  comparison-period weights the confounding-trend assumption allows that
  make the variance of the estimate smallest (see Details). For a single
  RD period and a character `trend`.

- estimand:

  `"att"` (default) when the comparison periods are uniformly untreated:
  the ATT, the effect of the treatment on the units just above the
  cutoff that are treated in the RD period. `"atu"` when the comparison
  periods are uniformly treated: the ATU, the effect for the units just
  below the cutoff, which are untreated in the RD period. The numbers
  are the same either way; see "Targeting the ATU".

- bwselect:

  the bandwidth rule, used when `h` is not given: `"joint"` (default;
  the common bandwidth, one `h` for every period, chosen to minimize the
  asymptotic mean squared error of the RD-DID estimate), `"cct"`
  (per-period CCT bandwidths: each period's own MSE-optimal bandwidth of
  Calonico, Cattaneo and Titiunik, from
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)).
  See "Bandwidth rules".

- h:

  the main bandwidth (point estimate). If given, it is used in every
  period and `bwselect` is ignored.

- b:

  the pilot bandwidth (bias correction), used with a given `h`; defaults
  to `h`.

- scheme:

  the sampling scheme, which sets the standard error: `"cs"` (repeated
  cross-section: different units in each period), `"pc"` (panel, running
  variable fixed over time: the same units, each on the same side of the
  cutoff in every period), `"pv"` (panel, running variable varies over
  time: some units change side between periods), or `"auto"` (default),
  which reads it off the data: no repeated `id` gives `"cs"`, repeated
  units that never change side give `"pc"`, any unit that changes side
  gives `"pv"`.

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- p:

  order of the local polynomial for the point estimate (default 1, local
  linear).

- q:

  order of the local polynomial for the bias correction (default 2);
  must exceed `p`.

- kernel:

  the kernel: `"triangular"` (default), `"epanechnikov"` or `"uniform"`.

- level:

  confidence level of the reported intervals (default 0.95).

- regularize:

  logical. If `TRUE` (default), the `"joint"` rule adds a regularization
  term to the estimated squared bias, as rdrobust does, so that a
  near-zero estimated bias cannot make the bandwidth very large. Not
  used with a given `h` or `bwselect = "cct"`.

- reg_const:

  the regularization constant: the multiple of the estimated variance of
  the bias constants added to the squared bias (default 3).

## Value

An object of class `"rddid"`, a list with:

- `estimates`:

  a data frame with two rows, `Conventional` (the local-linear estimate
  with its conventional standard error) and `Robust` (the bias-corrected
  estimate with its robust standard error, printed as "Robust
  (bias-corrected)"), and columns `est`, `se`, `ci_l`, `ci_u`, `z`, `p`
  (estimate, standard error, confidence limits, z statistic and p-value
  under `scheme`) and `se_cs`, `se_pc`, `se_pv` (that row's standard
  error under each sampling scheme, at the same bandwidths).

- `coef`:

  named numeric vector, the coefficient of each period's jump in the
  estimate: \\\phi_t\\ for each RD period (1 for a single RD period),
  minus its (aggregated) weight for each comparison period.
  [`coef()`](https://rdrr.io/r/stats/coef.html) returns the estimate
  itself.

- `coef_bc`:

  the same for the bias-corrected estimate; it differs from `coef` only
  under `weighting = "min_variance"`.

- `weights`:

  named numeric vector of the comparison-period weights; with several RD
  periods, the aggregated weights \\w^\phi_s\\; under
  `weighting = "min_variance"`, those of the conventional estimate (the
  bias-corrected estimate's are in `weights_detail$weights_bc`).

- `rd_weights`:

  named numeric vector of the RD-period weights \\\phi_t\\ (1 for a
  single RD period).

- `weights_by_target`:

  numeric matrix of the comparison-period weights for each RD period,
  one row per RD period.

- `weights_type`:

  `"constant"`, `"linear"`, or `"custom"` for numeric weights.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `t_rd`, `comparisons`:

  the RD period(s) and the comparison periods used.

- `scheme`:

  the sampling scheme behind `se`, `ci_l`, `ci_u`, `z` and `p`;
  `scheme_detected` is the scheme read off the data and
  `scheme_requested` the argument as passed.

- `bandwidth`:

  a list: `method` (the `bwselect` value, or `"fixed"` when `h` is
  given); `h_by_period` and `b_by_period` (named numeric vectors, the
  main and pilot bandwidth used in each period, whatever the rule); for
  `"fixed"` and `"joint"`, the common main bandwidth `h` and the pilot
  `b` (one number for `"fixed"`, one per period for `"joint"`); and
  intermediate quantities of the rule (`bws` for `"cct"`; `B`, `Veff`,
  `reg`, `pilot_bws` for `"joint"`).

- `fits`:

  named list of the per-period
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  fits, the RD period(s) first (index by name).

- `n_by_period`:

  the number of observations in each period.

- `level`, `c`, `p`, `q`, `kernel`:

  as passed.

- `data`:

  the per-period data used (a named list of data frames with columns
  `y`, `x`, `id`), for
  [`plot.rddid()`](https://dorleventer.github.io/rddid/reference/plot.rddid.md).

- `call`:

  the matched call.

- `weighting`:

  `"ols"` or `"min_variance"`, as passed.

- `weights_detail`:

  `NULL` for `"ols"`; for `"min_variance"` a list with `pilot` (the
  `"ols"` weights), `Xi` (the estimated covariance matrix of the
  comparison-period jumps), `xi` (their estimated covariances with the
  RD-period jump), `h` (the main bandwidths at which these were
  estimated; under `"joint"`, the first-step bandwidths, not the final
  ones in `bandwidth`) and `pinned` (`TRUE` when the assumption allows
  only one set of weights), and `weights_bc`, `Xi_bc`, `xi_bc`, the same
  for the bias-corrected estimate.

## Details

### The estimate

\$\$\mathrm{ATT}(t\_{\mathrm{RD}}) = D\_{t\_{\mathrm{RD}}} - \sum_t w_t
D_t,\$\$ where \\D_t\\ is the jump in the outcome at the cutoff in
period \\t\\, estimated by a local-linear RD on that period's
observations; \\t\_{\mathrm{RD}}\\ is the RD period (`t_rd`); the sum
runs over the comparison periods; and \\w_t\\ is the weight of
comparison period \\t\\.
[`summary()`](https://rdrr.io/r/base/summary.html) lists every \\D_t\\
with its coefficient in the sum.

### Several RD periods

With a vector `t_rd`, `rddid()` estimates the weighted average of the RD
periods' effects, \\\mathrm{ATT}^\phi = \sum_t \phi_t \mathrm{ATT}(t)\\,
with weights \\\phi_t\\ from `rd_weights` (equal by default; any finite
weights are allowed, as they define the estimand). Each period's jump is
estimated once: each RD period \\t\\ enters with \\\phi_t\\, and each
comparison period \\s\\ with the aggregated weight \\w^\phi_s = \sum_t
\phi_t w_s(t)\\, where \\w_s(t)\\ is the weight that RD period \\t\\
alone gives comparison period \\s\\ under the trend assumption. The
standard errors include the covariances between RD periods, and
`"joint"` chooses one bandwidth for the average, not for any single RD
period. `weighting = "min_variance"` is available for a single RD period
only.

### The confounding-trend assumption

With `trend = "constant"` the confounding jump is the same in every
period, so the comparison periods get equal weights, \\w_t = 1/m\\ with
\\m\\ comparison periods. With `trend = "linear"` the confounding jump
moves linearly in time; the weights then extrapolate the least-squares
line through the comparison-period jumps to the RD period (with two
comparison periods, the line through them). They sum to one and can be
negative: with comparison periods 1 and 2 and RD period 3 they are -1
and 2. A numeric `trend` supplies the weights directly. These are the
weights under the default `weighting = "ols"`;
`weighting = "min_variance"` chooses others the same assumption allows
(see below).

### Minimum-variance weights

All weights the confounding-trend assumption allows estimate the same
effect; they differ in precision. `weighting = "min_variance"` picks,
for each row of `estimates`, the ones with the smallest variance of that
estimate. For the conventional estimate it uses the estimated variances
of the per-period jumps and, in a panel (`scheme` `"pc"` or `"pv"`),
their covariances (the quantities behind the conventional standard
error). In a repeated cross-section with `trend = "constant"` these are
inverse-variance weights: a comparison period whose jump is estimated
more precisely gets more weight. Under `bwselect = "joint"` the weights
and the bandwidth are chosen in three steps: the common bandwidth with
the `"ols"` weights, the minimum-variance weights at that bandwidth, and
the bandwidths (common `h` and pilot `b`) again with those weights. With
`"cct"` or a fixed `h` the bandwidths do not depend on the weights, so
the weights are computed once. The bias-corrected estimate gets its own
weights, computed in the same way from the variances and covariances of
the bias-corrected jumps (those behind the robust standard error), at
the same bandwidths; they are in `weights_detail$weights_bc`, and the
common bandwidth is chosen with the conventional ones. The standard
errors treat the weights as known; in large samples, estimating them
does not change the distribution of the estimate, because every allowed
set of weights cancels the confounding jump exactly, so an error in the
weights multiplies only the estimation error of the jumps. With one
comparison period under `"constant"`, or two under `"linear"`, the
assumption allows only one set of weights; `rddid()` returns it and sets
`weights_detail$pinned` to `TRUE`.

### Sampling schemes

The scheme sets the standard error and, under `"joint"`, also the
bandwidth, and with it the estimate, because those rules balance bias
against the variance under the scheme in use; with a fixed `h` or
`"cct"` it changes only the standard error. In a repeated cross-section
(`"cs"`) the periods' samples are independent and the variance is the
weighted sum of the period variances. In a panel the same units appear
in several periods, so the period jumps are correlated and the variance
adds their covariances, computed by matching units on `id`: under `"pc"`
from units on the same side of the cutoff in both periods; under `"pv"`
also from the units that change side, which enter with the opposite
sign. The standard errors under all three schemes, at the bandwidths
actually used, are in `estimates` and in
[`summary()`](https://rdrr.io/r/base/summary.html).

### Bandwidth rules

- `"joint"` (default), the common bandwidth: one `h` in every period,
  chosen to minimize the asymptotic mean squared error of the RD-DID
  estimate, not of each period's jump, so the biases of the period jumps
  can partly cancel.

- `"cct"`, per-period CCT bandwidths: each period gets its own
  MSE-optimal bandwidth from
  [`rdrobust::rdbwselect()`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html),
  as if it were a stand-alone RD (Calonico, Cattaneo and Titiunik,
  2014).

- A numeric `h` is used in every period, with pilot bandwidth `b`
  (default `h`).

The bandwidths used in each period are in `bandwidth$h_by_period` and
`bandwidth$b_by_period`, and in
[`summary()`](https://rdrr.io/r/base/summary.html).

### Targeting the ATU

When the treatment of interest is uniformly present in the comparison
periods (everybody at the cutoff is treated), the same difference of
jumps identifies the ATU: the effect for the units just below the
cutoff, which are untreated in the RD period. The paper shows that this
design is the ATT design with the two sides of the cutoff exchanged (the
running variable mirrored around the cutoff). Mirroring changes neither
the jump estimates nor their standard errors nor the bandwidth rules, so
`estimand = "atu"` returns the same numbers as `"att"` and only labels
the output. The four tests of the assumptions take the same argument;
only
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
computes differently under `"atu"`.

### Technical details

The `"joint"` rule is AMSE-optimal: it minimizes the asymptotic mean
squared error (AMSE) of the RD-DID estimate over a common `h`, in the
way rdrobust chooses a bandwidth for one period. Each period's bias
constant is estimated from a fit at that period's own CCT bandwidths.
Each period's variance constant, and under `"pc"` the covariances
between periods, are estimated at one rule-of-thumb bandwidth shared by
every period (the preliminary bandwidth of rdrobust, computed from the
running variables of all periods pooled; `bandwidth$h_rot` in the
result). These are combined, with the coefficients of the sum above,
into the AMSE, which is then minimized in closed form. Each period's
pilot bandwidth keeps that period's CCT ratio \\b/h\\. Nothing in the
rule depends on which period is the RD period, so an aggregate of
several RD periods gets the same bandwidth whatever the order of `t_rd`.
With `regularize = TRUE` the rule adds `reg_const` times the estimated
variance of the aggregate bias constant, including its covariances
between periods, to the squared bias, as rdrobust does, so a near-zero
estimated bias cannot make the bandwidth very large. Standard errors use
the HC1 convention of rdrobust. The derivations are in the paper's
appendix on estimation and bandwidth choice.

## References

Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data.*
arXiv:2408.05847. <https://arxiv.org/abs/2408.05847>

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295-2326.

## See also

[`summary.rddid()`](https://dorleventer.github.io/rddid/reference/summary.rddid.md)
and
[rddid-methods](https://dorleventer.github.io/rddid/reference/rddid-methods.md)
([`coef()`](https://rdrr.io/r/stats/coef.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html)) for the fitted object;
the tests of the assumptions
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
and
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md);
the example data
[rddid_sim](https://dorleventer.github.io/rddid/reference/rddid_sim.md)
and
[rddid_sim_pv](https://dorleventer.github.io/rddid/reference/rddid_sim_pv.md).

Other RD-DID estimation:
[`plot.rddid()`](https://dorleventer.github.io/rddid/reference/plot.rddid.md),
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md),
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)

## Examples

``` r
# rddid_sim: confounding jump 0.5 in every year, treatment effect 1 in year 3
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit            # the estimate, its standard error and confidence interval
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2734 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4197, 2 = 0.3958, 3 = 0.4043
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0903     0.1251    8.72   <0.001   [0.8451, 1.3355]
#>   Robust (bias-corrected)      1.1387     0.1478    7.70   <0.001   [0.8489, 1.4284]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
summary(fit)   # the jump in every period, and the s.e. under each scheme
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2734 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4197, 2 = 0.3958, 3 = 0.4043
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0903     0.1251    8.72   <0.001   [0.8451, 1.3355]
#>   Robust (bias-corrected)      1.1387     0.1478    7.70   <0.001   [0.8489, 1.4284]
#> 
#>   Per-period local-linear fits (Conventional = sum of coef x jump; Robust = sum of coef x jump (bc)):
#>   period   role           coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison   -0.500   1000   0.2734   0.4197     0.5152    0.1672     0.4653    0.1932
#>   2        comparison   -0.500   1000   0.2734   0.3958     0.6308    0.1667     0.6269    0.1991
#>   3        RD            1.000   1000   0.2734   0.4043     1.6633    0.1624     1.6848    0.1912
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
#>     cross-section 0.2362   panel, no unit changes side 0.1478   panel, some change side 0.1478
coef(fit)
#> Conventional       Robust 
#>     1.090323     1.138689 
confint(fit, "Robust")
#>            2.5 %   97.5 %
#> Robust 0.8489493 1.428429
# the confounding jump moves linearly in time
rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
      trend = "linear")
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (linear confounding trend)
#>   Comparison weights: 1: -1.000, 2: 2.000
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.308 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4728, 2 = 0.4459, 3 = 0.4554
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9196     0.2416    3.81   <0.001   [0.4462, 1.3931]
#>   Robust (bias-corrected)      0.8810     0.2874    3.06    0.002   [0.3176, 1.4443]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
# minimum-variance comparison-period weights
rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
      weighting = "min_variance")
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights (minimum variance): 1: 0.453, 2: 0.547; Robust row: 1: 0.462, 2: 0.538
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.277 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4252, 2 = 0.401, 3 = 0.4096
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0828     0.1242    8.72   <0.001   [0.8394, 1.3262]
#>   Robust (bias-corrected)      1.1290     0.1467    7.70   <0.001   [0.8415, 1.4165]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
# several RD periods: the average of their effects (year 3 copied as a
# second RD year, 4, for illustration)
sim4 <- rbind(rddid_sim, transform(rddid_sim[rddid_sim$year == 3, ], year = 4))
rddid(sim4, y = "Y", x = "R", time = "year", id = "id", t_rd = c(3, 4))
#> RD-DID estimate of the ATT averaged over RD periods 3, 4
#>   RD-period weights: 3: 0.500, 4: 0.500
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights (summed over the RD periods): 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2734 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4197, 2 = 0.3958, 3 = 0.4043, 4 = 0.4043
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0903     0.1251    8.72   <0.001   [0.8451, 1.3355]
#>   Robust (bias-corrected)      1.1387     0.1478    7.70   <0.001   [0.8489, 1.4284]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
# comparison periods uniformly treated: add estimand = "atu" (same numbers)
```
