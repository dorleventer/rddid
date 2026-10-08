# Test of constant within-type confounding

When the running variable moves over time, some units are above the
cutoff in one period and below it in another. A unit's **type** is the
side of the cutoff it is on in the other period(s); by default here, its
side in the RD period, which stays fixed across the comparison periods.
The confounding-trend assumption of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) must
then hold within each type. It concerns the RD period, where the
confounding jump is not observed separately, so, like a pre-trends check
in difference-in-differences, `rd_trendcell()` tests it across the
comparison periods: the null is that **within each type, the confounding
jump is the same in every comparison period** (with `trend = "linear"`:
moves linearly in time). A rejection means the comparison periods do not
support the trend assumption, and the estimate of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
under that assumption can be biased; if the jumps move linearly,
consider `rddid(trend = "linear")`. With a running variable fixed over
time (as in
[rddid_sim](https://dorleventer.github.io/rddid/reference/rddid_sim.md))
the types are degenerate and the test is not informative.

## Usage

``` r
rd_trendcell(
  data,
  y,
  x,
  time,
  id,
  t_rd = NULL,
  comparisons = NULL,
  estimand = c("att", "atu"),
  trend = c("constant", "linear"),
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  min_n = 10L,
  bc = TRUE,
  type_by = c("rd_side", "pattern"),
  p = 1L,
  q = 2L
)
```

## Arguments

- data:

  a data frame in long format, one row per unit and period, from a panel
  (it need not be balanced). A unit's type is read from the periods in
  which it is observed; a unit missing from a period that its type needs
  is left out of the cells that use that type.

- y:

  name of the outcome column (a string).

- x:

  name of the running-variable column (a string).

- time:

  name of the period column (a string).

- id:

  name of the unit-identifier column (a string). Required: types are
  read across periods.

- t_rd:

  the RD period. Required with `type_by = "rd_side"` (the default),
  where a unit's side of the cutoff in the RD period is its type. The
  test does not use the RD period's outcomes.

- comparisons:

  the comparison periods in which the test runs, taken in time order.
  `NULL` (default) uses every period other than `t_rd` (every period if
  `t_rd` is `NULL`).

- estimand:

  `"att"` (default) or `"atu"`, as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).
  Label only: the test is the same either way.

- trend:

  the trend assumption tested within each type: `"constant"` (default;
  the confounding jump is the same in every comparison period) or
  `"linear"` (it moves linearly in time; needs at least three comparison
  periods). Use the `trend` of the
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  call being checked. With `"linear"` the second differences are taken
  in time (the period values), so unequally spaced comparison periods
  are handled.

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- h:

  a bandwidth to use, as both main and pilot bandwidth, in every cell.
  If given, `bwselect` is ignored.

- bwselect:

  the bandwidth rule when `h` is not given: `"cct"` (default; each
  cell's own CCT bandwidths from
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md))
  or `"rot"` (the rule of thumb `0.2` times the range of the running
  variable within the cell).

- kernel:

  the kernel: `"triangular"` (default), `"epanechnikov"` or `"uniform"`.

- scheme:

  the sampling scheme, which sets the covariance across comparison
  periods in the test: `"cs"`, `"pc"` or `"pv"` (as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)),
  or `"auto"` (default), which reads it off the comparison periods by
  the rule of
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).
  See Details.

- min_n:

  the minimum number of observations on each side of the cutoff for a
  (type, period) cell to enter the test (default 10). Smaller cells are
  dropped, with a message listing them.

- bc:

  logical. `TRUE` (default): test the bias-corrected jumps with their
  robust variance, as in the `Robust` row of
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  `FALSE`: the conventional jumps and variances.

- type_by:

  how types are defined: `"rd_side"` (default; the unit's side of the
  cutoff in the RD period, which needs `t_rd`) or `"pattern"` (its sides
  in all periods other than the RD period). Either way the type is fixed
  across the comparison periods.

- p, q:

  orders of the local polynomials in every cell, for the point estimate
  and the bias correction (defaults 1 and 2; `q` must exceed `p`). The
  per-cell CCT bandwidths are chosen for the order-`p` fit.

## Value

An object of class `"rd_trendcell"`, a list with:

- `statistic`, `df`, `p_value`:

  the Wald statistic, its degrees of freedom (the number of positive
  directions of the covariance used) and its chi-squared p-value; `NA`,
  `0`, `NA` when `trend = "linear"` and no type has jumps in three
  comparison periods.

- `scheme`:

  the sampling scheme used.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `call`:

  the matched call.

- `cell_period_jumps`:

  a data frame with one row per (type, comparison period) cell that was
  fitted: `cell` (the type: the unit's side(s), `"+"` above and `"-"`
  below the cutoff), `period`, `jump` (the cell's confounding jump,
  bias-corrected when `bc = TRUE`), `se`, `n` (observations in the cell)
  and `reference` (`TRUE` for each type's first comparison period under
  `trend = "constant"`).

- `contrasts`:

  named numeric vector of the tested differences, stacked across types.

- `cov_matrix`:

  the estimated covariance matrix of `contrasts` (block-diagonal by
  type).

- `bc`, `trend`:

  as passed.

- `comparisons`:

  the comparison periods, in time order.

## Details

### What is estimated

Each unit gets one type, fixed across the comparison periods. In each
comparison period and for each type, a local-linear RD of the outcome on
the running variable gives that type's confounding jump. Under
`trend = "constant"` each type's jump in every comparison period is
compared with its jump in the first one (one difference fewer than the
number of comparison periods, per type). Under `trend = "linear"` the
test uses the second differences of the time-ordered jumps (two fewer
per type), so it needs at least three comparison periods; if no type has
jumps in three, the function returns `statistic = NA` and `df = 0`, with
a message. All differences are tested jointly by a Wald statistic. Their
covariance is estimated and can be numerically indefinite, so the
statistic uses only its positive directions, and `df` counts them.

### Shared units and the sampling scheme

A unit's type is fixed, so different types are different units and their
jumps are independent. Within a type the same units can appear in
several comparison periods, and the covariance follows `scheme`,
matching units on `id`: none under `"cs"`; from the units on the same
side of the cutoff in both periods under `"pc"`; under `"pv"` also from
the units that change side, with the opposite sign. `"auto"` reads the
scheme off the comparison periods, by the rule of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).

### Options

`type_by = "rd_side"` (default) types each unit by its side of the
cutoff in the RD period, so `t_rd` is required; units not observed in
the RD period are left out. `type_by = "pattern"` types units by their
sides in all periods in `data` other than the RD period (other than the
first comparison period if `t_rd` is `NULL`). A (type, period) cell with
fewer than `min_n` observations on either side of the cutoff is dropped,
and a message lists the dropped cells. `bc = TRUE` (default) tests the
bias-corrected jumps with their robust variance, as in the `Robust` row
of [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
`bc = FALSE` uses the conventional jumps and variances. With
`bwselect = "cct"` (default) each cell gets its own CCT bandwidths from
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md),
computed on that cell's outcome and running variable. With
`bwselect = "rot"` the rule of thumb is `h = b = 0.2` times the range of
the running variable within each cell. A numeric `h` is used as both
bandwidths in every cell.

### ATU designs

When the comparison periods are uniformly treated, the comparison-period
jumps are the confounding jumps among treated units. The computation is
the same, so `estimand` only labels the output.

## References

Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data.*
arXiv:2408.05847. <https://arxiv.org/abs/2408.05847>

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295-2326.

## See also

[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) for
the estimate;
[rddid_sim_pv](https://dorleventer.github.io/rddid/reference/rddid_sim_pv.md)
for example data; `tidy()` in
[rddid-tidiers](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
for a one-row summary.

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md),
[`plot.rd_homog()`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md),
[`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
# rddid_sim_pv: the running variable moves, so some units change side
tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
tr          # the Wald test, then each type's jump in each comparison period
#> Test of a constant within-type confounding discontinuity  [rd_trendcell()]
#>   H0: within each type, the confounding jump is the same in every comparison period
#>   Sampling scheme: panel, running variable varies over time: some units change side (detected from the data)
#>   Comparison periods: 1, 2   Trend: constant
#> 
#>   Wald chi-squared(2) = 0.078,  p = 0.962
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Type       Period             jump       s.e.       n
#>     +          1                0.7374     0.3582     490  (reference)
#>     +          2                0.7356     0.3588     490
#>     -          1                0.8985     0.2213     510  (reference)
#>     -          2                1.0282     0.4007     510
# trend = "linear" needs three comparison periods; with two it is not testable
rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
             trend = "linear")
#> rd_trendcell: linear trend is not testable -- no cell has 3 or more comparison periods (degrees of freedom = 0). Returning an object with df = 0, statistic = NA, p_value = NA.
#> Test of a constant within-type confounding discontinuity  [rd_trendcell()]
#>   H0: within each type, the confounding jump moves linearly across the comparison periods
#>   Sampling scheme: panel, running variable varies over time: some units change side (detected from the data)
#>   Comparison periods: 1, 2   Trend: linear
#> 
#>   Wald: not testable (df = 0)
#>     (a linear within-type trend needs at least 3 comparison periods to be testable)
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Type       Period             jump       s.e.       n
#>     +          1                0.7374     0.3582     490
#>     +          2                0.7356     0.3588     490
#>     -          1                0.8985     0.2213     510
#>     -          2                1.0282     0.4007     510
```
