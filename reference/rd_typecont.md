# Test of type continuity

When the running variable moves over time, some units are above the
cutoff in one period and below it in another. A unit's **type** is the
side of the cutoff it is on in the other period(s); with two periods,
"above in the other period" or "below in the other period".
`rd_typecont()` tests the null that **the share of each type jumps by
zero at the cutoff, in every period**. A rejection means that units sort
across the cutoff by type, so the jump in the outcome at the cutoff can
reflect who is on each side as well as the policies, and the estimate of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) can
be biased. With a running variable fixed over time (as in
[rddid_sim](https://dorleventer.github.io/rddid/reference/rddid_sim.md))
every unit's type is its own side, the shares jump from 0 to 1 by
construction, and the test is not informative.

## Usage

``` r
rd_typecont(
  data,
  x,
  time,
  id,
  t_rd = NULL,
  comparisons = NULL,
  estimand = c("att", "atu"),
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  bc = TRUE
)
```

## Arguments

- data:

  a data frame in long format, one row per unit and period, from a panel
  (it need not be balanced). A unit's type is read from the periods in
  which it is observed; a unit missing from a period that its type needs
  is left out of the cells that use that type.

- x:

  name of the running-variable column (a string).

- time:

  name of the period column (a string). The column is usually numeric (a
  year); character or factor labels work as well, except under
  `trend = "linear"` (in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) or
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)),
  where the line is fitted on the period values: these then need to be
  numeric and are the time scale of the line (so 2015, 2017, 2018 are
  unequally spaced);
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) converts
  character labels such as `"2019"`. In
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  the default comparison periods are every period other than `t_rd`, in
  sorted order of the values (alphabetical for character or factor
  labels), the order a numeric `trend` follows.

- id:

  name of the unit-identifier column (a string). Required: types are
  read across periods.

- t_rd, comparisons:

  optional: the RD period and the comparison periods, to restrict the
  test to these periods. With `comparisons = NULL` (default) every
  period in `data` enters, whatever `t_rd`; since the test treats all
  periods alike, `t_rd` alone changes nothing and only lets you write
  the same call as for
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).

- estimand:

  `"att"` (default) or `"atu"`, as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).
  Label only: the test is the same either way.

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- h:

  a bandwidth to use, as both main and pilot bandwidth, in every
  regression of the test. If given, `bwselect` is ignored.

- bwselect:

  the bandwidth rule when `h` is not given: `"cct"` (default; each
  regression's own CCT bandwidths from
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md))
  or `"rot"` (the rule of thumb `0.5 * IQR(x)`, the same in every
  regression).

- kernel:

  the kernel: `"triangular"` (default), `"epanechnikov"` or `"uniform"`.

- scheme:

  the sampling scheme, which sets the covariance across periods in the
  test: `"cs"`, `"pc"` or `"pv"` (as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)),
  or `"auto"` (default), which reads it off the data as
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  does. See Details.

- bc:

  logical. `TRUE` (default): test the bias-corrected jumps with their
  robust variance, as in the `Robust` row of
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  `FALSE`: the conventional jumps and variances.

## Value

An object of class `"rd_typecont"`, a list with:

- `statistic`, `df`, `p_value`:

  the joint Wald statistic over all periods, its degrees of freedom and
  its chi-squared p-value.

- `scheme`:

  the sampling scheme used; `scheme_requested` is the argument as
  passed.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `fits`:

  the per-cell
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  fits, a list-matrix indexed by type and period (`NULL` where a cell
  could not be fitted).

- `data`:

  the typed data by period: for each period a data frame with `id`, `R`
  (the running variable) and `type`.

- `sides`:

  one row per unit with its running variable (`R_<period>`) and side of
  the cutoff (`side_<period>`, `"+"`/`"-"`) in every period;
  [`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md)
  reads it.

- `call`:

  the matched call.

- `per_period`:

  a list by period; each element holds `ll_wald`, that period's own Wald
  test (`stat`, `df`, `p`).

- `ll_wald`:

  the joint test again, as a list (`stat`, `df`, `p`).

- `meta`:

  a list with `periods`, `type_values` (the types, written as the sides
  in the other periods in time order, e.g. `"+-"`), `h` (the common
  bandwidth, `NA` with `bwselect = "cct"`), `bwselect`, `scheme`, `bc`
  and `estimand`.

## Details

### What is estimated

In each period and for each type, a local-linear RD of the indicator
"the unit is of this type" on the running variable estimates the jump in
that type's share at the cutoff. The shares sum to one within a period,
so one reference type per period is dropped (the type below the cutoff
in every other period), and the remaining jumps are tested jointly by a
Wald statistic, chi-squared with as many degrees of freedom as
independent jumps tested. Each period's own Wald test is reported as
well. A unit's type in a period needs its side in every other period, so
a unit missing from some period is left out of the regressions of the
other periods. The test treats all periods alike: `t_rd` and
`comparisons` only select which periods enter.

### Shared units and the sampling scheme

Within a period the type-indicator regressions use the same units, so
their covariance always enters. Across periods it follows `scheme`,
matching units on `id`: none under `"cs"`; from the units on the same
side of the cutoff in both periods under `"pc"`; under `"pv"` also from
the units that change side, with the opposite sign. `"auto"` reads the
scheme off the data as
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
does.

### Options

`bc = TRUE` (default) tests the bias-corrected jumps with their robust
variance, as in the `Robust` row of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
`bc = FALSE` uses the conventional jumps and variances. With
`bwselect = "cct"` (default) each (period, type) regression gets its own
CCT bandwidths from
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md).
With `bwselect = "rot"` the rule of thumb is `h = b = 0.5 * IQR(x)`, the
interquartile range of the running variable over all periods used
(`sd(x)` if that is zero), the same in every regression. A numeric `h`
is used as both bandwidths in every regression.

### ATU designs

The null treats the two sides of the cutoff alike, so the test is the
same under `estimand = "att"` and `"atu"`; `estimand` only labels the
output.

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
for example data;
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) in
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
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)

## Examples

``` r
# rddid_sim_pv: the running variable moves, so some units change side
tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
tc          # the null, the joint Wald test, then each period's own test
#> Test of a continuous type distribution  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Joint Wald chi-squared(9) = 4.782,  p = 0.853
#>     Period 1: chi-squared(3) = 1.531,  p = 0.675
#>     Period 2: chi-squared(3) = 0.382,  p = 0.944
#>     Period 3: chi-squared(3) = 2.202,  p = 0.532
tc$p_value
#> [1] 0.8529133
```
