# Test of composition stability

When the running variable moves over time, some units are above the
cutoff in one period and below it in another. A unit's **type** is the
side of the cutoff it is on in the other period(s); with two periods,
"above in the other period" or "below in the other period".
`rd_compstable()` tests the null that **the share of each type among the
units just above the cutoff is the same in the RD period and in each
comparison period**. Composition stability and homogeneous confounding
([`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md))
are alternatives: the estimate of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
needs one of the two (together with type continuity,
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)).
A rejection here alone therefore does not invalidate the estimate; if
homogeneous confounding is rejected as well, the comparison periods mix
the types differently from the RD period and the estimate can be biased.
With a running variable fixed over time (as in
[rddid_sim](https://dorleventer.github.io/rddid/reference/rddid_sim.md))
the types are degenerate and the test is not informative.

## Usage

``` r
rd_compstable(
  data,
  x,
  time,
  id,
  t_rd,
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

- t_rd:

  the RD period: the value of `time` in which the treatment of interest
  switches on at the cutoff.

- comparisons:

  the comparison periods: values of `time` in which the treatment of
  interest is uniform at the cutoff (nobody treated, or, with
  `estimand = "atu"`, everybody treated). `NULL` (default) uses every
  period other than `t_rd`; pass the periods explicitly when the data
  contain periods that are neither (a second RD period, say).

- estimand:

  `"att"` (default) or `"atu"`, as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).
  Under `"atu"` the test is on the shares among the units below the
  cutoff; see "ATU designs" in Details.

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

  the sampling scheme for the covariance between the two groups of a
  pair: `"auto"` (default; `"pv"` for a pair in which some unit is above
  the cutoff in both periods, `"cs"` otherwise), `"cs"`, `"pc"` or
  `"pv"` (as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)).
  See Details.

- bc:

  logical. `TRUE` (default): test the bias-corrected jumps with their
  robust variance, as in the `Robust` row of
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  `FALSE`: the conventional jumps and variances.

## Value

An object of class `"rd_compstable"`, a list with:

- `statistic`, `df`, `p_value`:

  the joint test over pairs: the sum of the pair Wald statistics, the
  sum of their degrees of freedom, and the chi-squared p-value
  (approximate; see Details).

- `scheme`:

  the sampling scheme used: one value if every pair used the same,
  `"mixed"` otherwise; `scheme_requested` is the argument as passed.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `call`:

  the matched call.

- `t_rd`, `comparisons`:

  the RD period and the comparison periods used.

- `pairs`:

  a list with one element per pair, named
  `"<RD period>::<comparison period>"`, each holding `ll_wald` (the
  pair's Wald test: `stat`, `df`, `p`); `jumps` and `jump_se` (the
  tested share jumps and their standard errors, named by type: the
  unit's sides, `1` above and `0` below the cutoff, with its side in the
  other period of the pair last); `type_values` (the types present);
  `scheme` (the scheme used for the pair); and `n_trd`, `n_t0`, `n_both`
  (the number of units above the cutoff in the RD period, in the
  comparison period, and in both; below the cutoff under
  `estimand = "atu"`); `fits` (each type's
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  fit on the reflected sample, `NULL` where it failed) and `sample` (the
  reflected sample: `x_trd`, `type_trd`, `x_t0`, `type_t0`), which feed
  [`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md).

- `joint`:

  the joint test again, as `ll_wald` (`stat`, `df`, `p`).

- `meta`:

  a list with `t_rd`, `comparisons`, `h` (the common bandwidth, `NA`
  with `bwselect = "cct"`), `bwselect`, `c` (the cutoff as passed,
  before any mirroring), `bc` and `estimand`.

## Details

### What is estimated

The test runs on each pair of the RD period and one comparison period.
Take the units above the cutoff in each of the two periods and stack
them into one artificial sample: the RD-period units at their distance
above the cutoff, the comparison-period units reflected to the same
distance below an artificial cutoff at zero. In a local-linear RD of a
type indicator on this reflected running variable, the jump at zero is
the type's share among the RD-period units just above the cutoff minus
its share among the comparison-period units just above the cutoff. With
two periods the type is the unit's side in the other period of the pair,
and there is one jump per pair; with more periods the type also records
the unit's sides in the remaining periods. The shares sum to one, so one
reference type is dropped and the remaining jumps are tested jointly by
a Wald statistic, one test per pair. The joint test over pairs adds up
the pair statistics and degrees of freedom, which treats the pairs as
independent although they share the RD-period units, so its p-value is
approximate (the paper's test is per pair). A pair with fewer than three
units in either group is skipped with a warning.

### Shared units and the sampling scheme

A unit above the cutoff in both periods of a pair appears on both sides
of the artificial cutoff. Under `"pv"` the covariance between its two
appearances, matched on `id`, is subtracted from the variance of the
jump; `"cs"` and `"pc"` treat the two groups as independent and give the
same test. `"auto"` (default) uses `"pv"` for a pair in which some unit
is above the cutoff in both periods and `"cs"` otherwise; the scheme is
reported per pair.

### Options

`bc = TRUE` (default) tests the bias-corrected jumps with their robust
variance, as in the `Robust` row of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
`bc = FALSE` uses the conventional jumps and variances. With
`bwselect = "cct"` (default) each type-indicator regression in the
reflected sample gets its own CCT bandwidths from
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md).
With `bwselect = "rot"` the rule of thumb is `h = b = 0.5 * IQR(x)`, the
interquartile range of the running variable in `data` (`sd(x)` if that
is zero), the same in every regression. A numeric `h` is used as both
bandwidths in every regression.

### ATU designs

With `estimand = "atu"` (comparison periods uniformly treated) the
running variable is mirrored around the cutoff before the construction
above, so the test is on the units *below* the cutoff: the null becomes
that the share of each type among the units below the cutoff is the same
in the RD period and in each comparison period. The type keeps its
meaning, so the reported jump is the change in the share of below-cutoff
units that are above the cutoff in the other period. This is the only
one of the four tests whose computation changes with `estimand`. Units
exactly at the cutoff count as above it in the original design and
cannot be placed in the mirrored one, so `"atu"` stops with an error if
any `x == c`; put the cutoff between support points (e.g. `c = 4999.5`
for integer populations).

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
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
# rddid_sim_pv: the running variable moves, so some units change side
cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
cs          # one Wald test per (RD period, comparison period) pair, then their sum
#> Test of composition stability  [rd_compstable()]
#>   H0: the share of each type among the units just above the cutoff is the same in the RD period and in each comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Pair 3::1: chi-squared(3) = 21.982,  p = <0.001
#>     n above the cutoff: 490 (RD period), 508 (comparison), 430 in both
#>   Pair 3::2: chi-squared(3) = 6.137,  p = 0.105
#>     n above the cutoff: 490 (RD period), 497 (comparison), 418 in both
#> 
#>   Joint over pairs (sum of chi-squared): chi-squared(6) = 28.119,  p = <0.001
cs$pairs[["3::1"]]$jumps
#>         01         10         11 
#> -0.2519944 -0.2544709  0.3613727 
# comparison periods uniformly treated: the shares below the cutoff
rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3,
              estimand = "atu")
#> Test of composition stability  [rd_compstable()]
#>   H0: the share of each type among the units just below the cutoff is the same in the RD period and in each comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Estimand: ATU (the units below the cutoff are the ones untreated in the RD period, so the test is on their shares (mirrored design))
#>   RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Pair 3::1: chi-squared(3) = 6.206,  p = 0.102
#>     n below the cutoff: 510 (RD period), 492 (comparison), 432 in both
#>   Pair 3::2: chi-squared(3) = 3.164,  p = 0.367
#>     n below the cutoff: 510 (RD period), 503 (comparison), 431 in both
#> 
#>   Joint over pairs (sum of chi-squared): chi-squared(6) = 9.370,  p = 0.154
```
