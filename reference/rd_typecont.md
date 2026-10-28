# Test of type continuity

When the running variable moves over time, some units are above the
cutoff in one period and below it in another. A unit's **type** is the
side of the cutoff it is on in the other period(s); with two periods,
"above in the other period" or "below in the other period".
`rd_typecont()` tests the null that **the share of each type jumps by
zero at the cutoff, in every period**. A rejection means that units sort
across the cutoff by type, so the jump in the outcome at the cutoff can
reflect who is on each side as well as the policies, and the estimate of
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
can be biased. With a running variable fixed over time (as in
[didisc_sim](https://dorleventer.github.io/didisc/reference/didisc_sim.md))
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
  bc = TRUE,
  joint = c("supt", "swald", "wald"),
  studentize = TRUE
)
```

## Arguments

- data:

  a data frame in long format, one row per unit and period, from a panel
  (it need not be balanced). Only the rows of the periods used (see
  `t_rd`, `comparisons`) are read, and a unit's type comes from those
  periods; a unit missing from a period that its type needs is left out
  of the cells that use that type. Rows with a missing running variable,
  period or id are dropped, with a message.

- x:

  name of the running-variable column (a string).

- time:

  name of the period column (a string). The column is usually numeric (a
  year); character or factor labels work as well, except under
  `trend = "linear"` (in
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
  or
  [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)),
  where the line is fitted on the period values: these then need to be
  numeric and are the time scale of the line (so 2015, 2017, 2018 are
  unequally spaced);
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) converts
  character labels such as `"2019"`. In
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
  the default comparison periods are every period other than `t_rd`, in
  sorted order of the values (alphabetical for character or factor
  labels), the order a numeric `trend` follows.

- id:

  name of the unit-identifier column (a string). Required: types are
  read across periods.

- t_rd, comparisons:

  the periods the test uses: `t_rd` and the `comparisons`;
  `comparisons = NULL` (default) uses every period in `data` (other than
  `t_rd`, which is then used too). The test treats all periods alike, so
  `t_rd` is optional: it lets you write the same call as for
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md).

- estimand:

  `"att"` (default) or `"atu"`, as in
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md).
  Label only: the test is the same either way.

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- h:

  a bandwidth to use, as both main and pilot bandwidth, in every
  regression of the test. If given, `bwselect` is ignored.

- bwselect:

  the bandwidth rule when `h` is not given: `"cct"` (default; each
  regression's own CCT bandwidths from
  [`rd_bw_cct()`](https://dorleventer.github.io/didisc/reference/rd_bw_cct.md))
  or `"rot"` (the rule of thumb `0.5 * IQR(x)`, the same in every
  regression).

- kernel:

  the kernel: `"triangular"` (default), `"epanechnikov"` or `"uniform"`.

- scheme:

  the sampling scheme, which sets the covariance across periods in the
  test: `"cs"`, `"pc"` or `"pv"` (as in
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)),
  or `"auto"` (default), which reads it off the data as
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
  does. See Details.

- bc:

  logical. `TRUE` (default): test the bias-corrected jumps with their
  robust variance, as in the `Robust` row of
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md);
  `FALSE`: the conventional jumps and variances.

- joint:

  the joint statistic: `"supt"` (default; the largest \|t\| against a
  sup-t critical value), `"swald"` (the sum of squared t's) or `"wald"`
  (the chi-squared Wald statistic). See "The joint statistic" in
  Details.

- studentize:

  logical, for `joint = "supt"` or `"swald"`. `TRUE` (default): t-based
  critical values that account for the noise in each jump's estimated
  standard error; `FALSE`: normal critical values.

## Value

An object of class `"rd_typecont"`, a list with:

- `statistic`, `df`, `p_value`:

  the joint test over all periods under `joint` and `studentize`: the
  statistic (the largest \|t\| for `"supt"`, the sum of squared t's for
  `"swald"`, the chi-squared statistic for `"wald"`), the number of
  jumps tested (`"supt"`, `"swald"`) or the degrees of freedom
  (`"wald"`), and the p-value.

- `joint_method`, `studentize`:

  as used; `joint_tests` has every statistic (a data frame with
  `method`, `reference` (`"t"`, `"normal"` or `"chi-squared"`),
  `statistic`, `df`, `p_value`).

- `t_stats`, `nu`:

  the t-statistics of the jumps tested and the degrees of freedom of
  their variances, named `"<period>:<type>"`.

- `scheme`:

  the sampling scheme used; `scheme_requested` is the argument as
  passed.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `fits`:

  the per-cell
  [`rd_period()`](https://dorleventer.github.io/didisc/reference/rd_period.md)
  fits, a list-matrix indexed by type and period (`NULL` where a cell
  could not be fitted).

- `data`:

  the typed data by period: for each period a data frame with `id`, `R`
  (the running variable) and `type`.

- `sides`:

  one row per unit with its running variable (`R_<period>`) and side of
  the cutoff (`side_<period>`, `"+"`/`"-"`) in every period;
  [`plot.rd_typecont()`](https://dorleventer.github.io/didisc/reference/plot.rd_typecont.md)
  reads it.

- `call`:

  the matched call.

- `per_period`:

  a list by period; each element holds `test`, that period's own test
  under `joint` (`stat`, `df`, `p`), and `ll_wald`, its Wald test.

- `ll_wald`:

  the joint Wald test, as a list (`stat`, `df`, `p`). A type share that
  is constant near the cutoff on both sides with no jump carries no
  information and is left out (`dropped` counts them); one that jumps
  deterministically makes the test undefined (`NA`).

- `meta`:

  a list with `periods`, `type_values` (the types, written as the sides
  in the other periods in time order, e.g. `"+-"`), `h` (the common
  bandwidth, `NA` with `bwselect = "cct"`), `bwselect`, `scheme`, `bc`,
  `estimand`, `joint` and `studentize`.

## Details

### What is estimated

In each period and for each type, a local-linear RD of the indicator
"the unit is of this type" on the running variable estimates the jump in
that type's share at the cutoff. The shares sum to one within a period,
so one reference type per period is dropped (the type below the cutoff
in every other period), and the remaining jumps are tested jointly (see
"The joint statistic"). Each period's own test is reported as well. A
unit's type in a period needs its side in every other period, so a unit
missing from some period is left out of the regressions of the other
periods. The test treats all periods alike: `t_rd` and `comparisons`
only select which periods enter.

### The joint statistic

Two arguments choose it. `joint` picks the statistic: `"supt"` (default)
rejects when the largest \|t\| of the jumps exceeds a sup-t critical
value, and the printout names the (period, type) cell with the largest
\|t\|; `"swald"` uses the sum of the squared t's (the sWald test of
Fusejima, Ishihara and Sawada, 2025); `"wald"` is the chi-squared Wald
statistic. `studentize` picks the reference distribution of sup-t and
sWald: with `TRUE` (default) each t is compared as a t-statistic, its
standard error being an estimate with the degrees of freedom of its
robust variance (computed from the jump's regression weights); with
`FALSE` as a standard normal. The critical values come from the
estimated correlation of the jumps.

Why these defaults: sup-t and sWald use the covariance only through its
correlations, so they avoid inverting it, an error that grows with the
number of jumps; t critical values account for the noise in each jump's
estimated standard error. The simulation evidence (size and power of
each statistic) is in the article "Checking the identification
assumptions".

Every statistic is reported in `joint_tests`, so any of them can be read
off one call. The simulated critical values use 100,000 draws under a
fixed seed, so a result is reproducible, and the caller's random-number
state is left as it was. With one jump, sup-t and sWald are the
two-sided t (or normal) test.

### Shared units and the sampling scheme

Within a period the type-indicator regressions use the same units, so
their covariance always enters. Across periods it follows `scheme`,
matching units on `id`: none under `"cs"`; from the units on the same
side of the cutoff in both periods under `"pc"`; under `"pv"` also from
the units that change side, with the opposite sign. `"auto"` reads the
scheme off the data as
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
does.

### Options

`bc = TRUE` (default) tests the bias-corrected jumps with their robust
variance, as in the `Robust` row of
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md);
`bc = FALSE` uses the conventional jumps and variances. With
`bwselect = "cct"` (default) each (period, type) regression gets its own
CCT bandwidths from
[`rd_bw_cct()`](https://dorleventer.github.io/didisc/reference/rd_bw_cct.md).
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

[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
for the estimate;
[didisc_sim_pv](https://dorleventer.github.io/didisc/reference/didisc_sim_pv.md)
for example data;
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) in
[didisc-tidiers](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
for a one-row summary.

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/didisc/reference/plot.rd_compstable.md),
[`plot.rd_homog()`](https://dorleventer.github.io/didisc/reference/plot.rd_homog.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/didisc/reference/plot.rd_trendcell.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/didisc/reference/plot.rd_typecont.md),
[`plot_switchers()`](https://dorleventer.github.io/didisc/reference/plot_switchers.md),
[`rd_checks()`](https://dorleventer.github.io/didisc/reference/rd_checks.md),
[`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)

## Examples

``` r
# didisc_sim_pv: the running variable moves, so some units change side
tc <- rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
tc          # the null, the joint sup-t test, then each period's own test
#> Test of type continuity  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#>   Type in a period: a unit's sides of the cutoff in the other periods, in time order
#>     ("+" at or above, "-" below).
#> 
#>   Joint sup-t (t critical values) max |t| = 1.232 over 9 jumps,  p = 0.869  (largest: period:type 1:+-)
#>     Period 1: sup-t (t critical values) max |t| = 1.232 over 3 jumps,  p = 0.499
#>     Period 2: sup-t (t critical values) max |t| = 0.582 over 3 jumps,  p = 0.902
#>     Period 3: sup-t (t critical values) max |t| = 1.048 over 3 jumps,  p = 0.631
#>   All joint statistics: x$joint_tests
tc$p_value
#> [1] 0.8688813
tc$joint_tests   # every statistic side by side
#>   method   reference statistic df   p_value
#> 1   supt           t  1.231775  9 0.8688813
#> 2   supt      normal  1.231775  9 0.8609214
#> 3  swald           t  5.110313  9 0.7953520
#> 4  swald      normal  5.110313  9 0.7861821
#> 5   wald chi-squared  4.781633  9 0.8529133
# the chi-squared Wald, or sup-t with normal critical values:
rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, joint = "wald")
#> Test of type continuity  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#>   Type in a period: a unit's sides of the cutoff in the other periods, in time order
#>     ("+" at or above, "-" below).
#> 
#>   Joint Wald chi-squared(9) = 4.782,  p = 0.853
#>     Period 1: chi-squared(3) = 1.531,  p = 0.675
#>     Period 2: chi-squared(3) = 0.382,  p = 0.944
#>     Period 3: chi-squared(3) = 2.202,  p = 0.532
#>   All joint statistics: x$joint_tests
rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, studentize = FALSE)
#> Test of type continuity  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#>   Type in a period: a unit's sides of the cutoff in the other periods, in time order
#>     ("+" at or above, "-" below).
#> 
#>   Joint sup-t (normal critical values) max |t| = 1.232 over 9 jumps,  p = 0.861  (largest: period:type 1:+-)
#>     Period 1: sup-t (normal critical values) max |t| = 1.232 over 3 jumps,  p = 0.490
#>     Period 2: sup-t (normal critical values) max |t| = 0.582 over 3 jumps,  p = 0.900
#>     Period 3: sup-t (normal critical values) max |t| = 1.048 over 3 jumps,  p = 0.628
#>   All joint statistics: x$joint_tests
```
