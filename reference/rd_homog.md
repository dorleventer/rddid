# Test of homogeneous confounding

When the running variable moves over time, some units are above the
cutoff in one period and below it in another. A unit's **type** is the
side of the cutoff it is on in the other period(s); by default here, its
side in the RD period. `rd_homog()` tests the null that **in each
comparison period the confounding jump is the same for every type**,
using the comparison periods only, where the jump in the outcome is the
confounding jump. Homogeneous confounding and composition stability
([`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md))
are alternatives: the estimate of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
needs one of the two (together with type continuity,
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)).
A rejection here alone therefore does not invalidate the estimate; if
composition stability is rejected as well, the comparison periods mix
the types differently from the RD period and the estimate can be biased.
With a running variable fixed over time (as in
[rddid_sim](https://dorleventer.github.io/rddid/reference/rddid_sim.md))
the types are degenerate and the test is not informative.

## Usage

``` r
rd_homog(
  data,
  y,
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

  the comparison periods in which the test runs. `NULL` (default) uses
  every period other than `t_rd` (every period if `t_rd` is `NULL`).

- estimand:

  `"att"` (default) or `"atu"`, as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).
  Label only: the test is the same either way.

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
  (period, type) cell to enter the test (default 10). Smaller cells are
  dropped, with a message listing them.

- bc:

  logical. `TRUE` (default): test the bias-corrected jumps with their
  robust variance, as in the `Robust` row of
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  `FALSE`: the conventional jumps and variances.

- type_by:

  how types are defined: `"rd_side"` (default; the unit's side of the
  cutoff in the RD period, which needs `t_rd`) or `"pattern"` (its sides
  in all other periods in `data`).

- p, q:

  orders of the local polynomials in every cell, for the point estimate
  and the bias correction (defaults 1 and 2; `q` must exceed `p`). The
  CCT bandwidths are always chosen for a local-linear fit.

## Value

An object of class `"rd_homog"`, a list with:

- `statistic`, `df`, `p_value`:

  the Wald statistic, its degrees of freedom (the number of positive
  directions of the covariance used) and its chi-squared p-value.

- `scheme`:

  the sampling scheme used.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `call`:

  the matched call.

- `period_type_jumps`:

  a data frame with one row per (comparison period, type) cell that was
  fitted: `period`, `type` (the unit's side(s), `"+"` above and `"-"`
  below the cutoff), `jump` (the cell's confounding jump, bias-corrected
  when `bc = TRUE`), `se`, `n` (observations in the cell) and
  `reference` (`TRUE` for the reference type).

- `contrasts`:

  named numeric vector of the tested differences (type minus reference),
  stacked across periods.

- `cov_matrix`:

  the estimated covariance matrix of `contrasts`.

- `bc`:

  as passed.

- `comparisons`:

  the comparison periods that contributed a difference.

## Details

### What is estimated

In each comparison period the units are split by type, and a
local-linear RD of the outcome on the running variable within each type
gives that type's confounding jump. Within each period every type's jump
is compared with the jump of a reference type (the type below the cutoff
in the other period(s)), and these differences are tested jointly across
the comparison periods by a Wald statistic. The null is equality across
types within each period, not equality across periods (that is
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)).
The covariance of the differences is estimated and can be numerically
indefinite, so the Wald statistic uses only its positive directions, and
`df` counts them.

### Shared units and the sampling scheme

Within a period the types are different units, so their jumps are
independent. Across comparison periods the same units can appear in
both, and the covariance follows `scheme`, matching units on `id`: none
under `"cs"`; from the units on the same side of the cutoff in both
periods under `"pc"`; under `"pv"` also from the units that change side,
with the opposite sign. `"auto"` reads the scheme off the comparison
periods, by the rule of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).

### Options

`type_by = "rd_side"` (default) types each unit by its side of the
cutoff in the RD period, so `t_rd` is required and each comparison
period gives one difference; units not observed in the RD period are
left out. `type_by = "pattern"` types units by their sides in all other
periods in `data`. A (period, type) cell with fewer than `min_n`
observations on either side of the cutoff is dropped, and a message
lists the dropped cells; a period with fewer than two types is skipped.
`bc = TRUE` (default) tests the bias-corrected jumps with their robust
variance, as in the `Robust` row of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
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
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
# rddid_sim_pv: the running variable moves, so some units change side
hc <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
hc          # the Wald test, then each type's jump in each comparison period
#> Test of homogeneous confounding  [rd_homog()]
#>   H0: in each comparison period the confounding jump is the same for every type
#>   Sampling scheme: panel, running variable varies over time: some units change side (detected from the data)
#>   Comparison periods: 1, 2
#> 
#>   Wald chi-squared(2) = 0.465,  p = 0.793
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Period     Type               jump       s.e.       n
#>     1          -                0.8985     0.2213     510  (reference)
#>     1          +                0.7374     0.3582     490
#>     2          -                1.0282     0.4007     510  (reference)
#>     2          +                0.7356     0.3588     490
hc$period_type_jumps
#>   period type      jump        se   n reference
#> 1      1    - 0.8985337 0.2213343 510      TRUE
#> 2      1    + 0.7373828 0.3582032 490     FALSE
#> 3      2    - 1.0281886 0.4007468 510      TRUE
#> 4      2    + 0.7355842 0.3588090 490     FALSE
```
