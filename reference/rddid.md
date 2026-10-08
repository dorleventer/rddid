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
  estimand = c("att", "atu"),
  bwselect = c("joint", "iter", "cct"),
  h = NULL,
  b = NULL,
  scheme = c("auto", "cs", "pc", "pv"),
  c = 0,
  p = 1L,
  q = 2L,
  kernel = "triangular",
  level = 0.95,
  start = "hstar",
  regularize = TRUE,
  reg_const = 3,
  weights = NULL
)
```

## Arguments

- data:

  a data frame in long format, one row per unit and period: a repeated
  cross-section (different units in each period) or a panel (the same
  units in several periods; it need not be balanced).

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
  switches on at the cutoff.

- comparisons:

  the comparison periods: values of `time` in which the treatment of
  interest is uniform at the cutoff (nobody treated, or, with
  `estimand = "atu"`, everybody treated). `NULL` (default) uses every
  period other than `t_rd`; pass the periods explicitly when the data
  contain periods that are neither (a second RD period, say).

- trend:

  the confounding-trend assumption, which sets the comparison-period
  weights: `"constant"` (default; the confounding jump is the same in
  every period, equal weights) or `"linear"` (the confounding jump moves
  linearly in time; needs at least two comparison periods). A numeric
  vector gives the weights directly, one per entry of `comparisons` in
  that order. Weights that sum to one cancel a constant confounding
  jump; `rddid()` warns when they do not.

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
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)),
  or `"iter"` (a separate bandwidth in each period, chosen together for
  the RD-DID estimate; not used in the paper, kept for simulations). See
  "Bandwidth rules".

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

- start:

  where the `"iter"` rule starts: `"hstar"` (default; the common
  bandwidth in every period), `"cct"` (each period's CCT bandwidth), or
  a named numeric vector or list with one starting bandwidth per period.
  Used only with `bwselect = "iter"`.

- regularize:

  logical. If `TRUE` (default), the `"joint"` and `"iter"` rules add a
  regularization term to the estimated squared bias, as rdrobust does,
  so that a near-zero estimated bias cannot make the bandwidth very
  large. Not used with a given `h` or `bwselect = "cct"`.

- reg_const:

  the regularization constant: the multiple of the estimated variance of
  the bias constants added to the squared bias (default 3).

- weights:

  the old name of `trend`, still accepted with a message. It is not a
  vector of observation weights; `rddid()` has none.

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
  estimate: 1 for the RD period, minus its weight for each comparison
  period. [`coef()`](https://rdrr.io/r/stats/coef.html) returns the
  estimate itself.

- `weights`:

  named numeric vector of the comparison-period weights.

- `weights_type`:

  `"constant"`, `"linear"`, or `"custom"` for numeric weights.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `t_rd`, `comparisons`:

  the RD period and the comparison periods used.

- `scheme`:

  the sampling scheme behind `se`, `ci_l`, `ci_u`, `z` and `p`;
  `scheme_detected` is the scheme read off the data and
  `scheme_requested` the argument as passed.

- `bandwidth`:

  a list: `method` (the `bwselect` value, or `"fixed"` when `h` is
  given); `h_by_period` and `b_by_period` (named numeric vectors, the
  main and pilot bandwidth used in each period, whatever the rule); for
  `"fixed"` and `"joint"`, the common main bandwidth `h` and the pilot
  `b` (one number for `"fixed"`, one per period for `"joint"`); for
  `"iter"`, the number of iterations `niter`; and intermediate
  quantities of the rule (`bws` for `"cct"` and `"iter"`; `B`, `Veff`,
  `reg`, `pilot_bws` for `"joint"`).

- `fits`:

  named list of the per-period
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  fits, the RD period first (index by name).

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

### The confounding-trend assumption

With `trend = "constant"` the confounding jump is the same in every
period, so the comparison periods get equal weights, \\w_t = 1/m\\ with
\\m\\ comparison periods. With `trend = "linear"` the confounding jump
moves linearly in time; the weights then extrapolate the least-squares
line through the comparison-period jumps to the RD period (with two
comparison periods, the line through them). They sum to one and can be
negative: with comparison periods 1 and 2 and RD period 3 they are -1
and 2. A numeric `trend` supplies the weights directly.

### Sampling schemes

The scheme sets the standard error and, under `"joint"` and `"iter"`,
also the bandwidth, and with it the estimate, because those rules
balance bias against the variance under the scheme in use; with a fixed
`h` or `"cct"` it changes only the standard error. In a repeated
cross-section (`"cs"`) the periods' samples are independent and the
variance is the weighted sum of the period variances. In a panel the
same units appear in several periods, so the period jumps are correlated
and the variance adds their covariances, computed by matching units on
`id`: under `"pc"` from units on the same side of the cutoff in both
periods; under `"pv"` also from the units that change side, which enter
with the opposite sign. The standard errors under all three schemes, at
the bandwidths actually used, are in `estimates` and in
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

- `"iter"`: a separate bandwidth in each period, chosen together to
  minimize the asymptotic mean squared error of the RD-DID estimate; not
  used in the paper, kept for simulations.

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
squared error (AMSE) of the RD-DID estimate over a common `h`. Each
period is first fitted at its own CCT bandwidths to estimate that
period's bias constant and variance constant; these are combined, with
the coefficients of the sum above (and, under `"pc"`, the covariances
between periods), into the AMSE, which is then minimized in closed form.
Each period's pilot bandwidth keeps that period's CCT ratio \\b/h\\.
Each period's bias and variance constants are estimated at that period's
own CCT pilot bandwidths, not at the RD period's, so an aggregate of
several RD periods gets the same bandwidth whichever of them carries the
`t_rd` label. The `"iter"` rule minimizes the same objective over one
bandwidth per period by coordinate descent, starting from `start`. With
`regularize = TRUE` both rules add `reg_const` times the estimated
variance of the bias constants to the squared bias, as rdrobust does, so
a near-zero estimated bias cannot make the bandwidth very large.
Standard errors use the HC1 convention of rdrobust. The derivations are
in the paper's appendix on estimation and bandwidth choice.

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
summary(fit)   # the jump in every period, and the s.e. under each scheme
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
#>   Per-period local-linear fits, in time order (estimate = sum of coef x jump):
#>   period   role          coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison    -0.5   1000   0.2672   0.4102     0.5119    0.1688     0.4639    0.1948
#>   2        comparison    -0.5   1000   0.2672   0.3868     0.6306    0.1686     0.6240    0.2012
#>   3        RD               1   1000   0.2672   0.3951     1.6640    0.1641     1.6854    0.1931
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
#>     cross-section 0.2385   panel, no unit changes side 0.1494   panel, some change side 0.1494
coef(fit)
#> Conventional       Robust 
#>     1.092709     1.141380 
confint(fit, "Robust")
#>            2.5 %   97.5 %
#> Robust 0.8485735 1.434187
# the confounding jump moves linearly in time
rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
      trend = "linear")
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (linear confounding trend; weights -1, 2)
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2672 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4102, 2 = 0.3869, 3 = 0.3951
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9147     0.2582    3.54   <0.001   [0.4087, 1.4207]
#>   Robust (bias-corrected)      0.9012     0.3068    2.94    0.003   [0.3000, 1.5025]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
# comparison periods uniformly treated: add estimand = "atu" (same numbers)
```
