# Run the four tests of the identification assumptions

When the running variable moves over time, the estimate of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) is
identified under type continuity, constant within-type confounding, and
either composition stability or homogeneous confounding (Leventer and
Nevo, 2024). `rd_checks()` runs the four tests,
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
and
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
with one set of options, and prints one table. Each test's own result is
kept, so every detail stays available.

## Usage

``` r
rd_checks(
  data,
  y,
  x,
  time,
  id,
  t_rd,
  comparisons = NULL,
  estimand = c("att", "atu"),
  c = 0,
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  bc = TRUE,
  joint = c("supt", "swald", "wald"),
  studentize = TRUE,
  outcome_bwselect = c("pooled", "cct", "rot"),
  typecont = list(),
  compstable = list(),
  homog = list(),
  trendcell = list(),
  level = 0.95
)
```

## Arguments

- data, y, x, time, id:

  as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  `id` is required.

- t_rd:

  the RD period (required).

- comparisons, estimand, c, kernel, scheme, bc:

  shared by the four tests; defaults as in each test.

- joint:

  the joint statistic of
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  and
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md):
  `"supt"` (default), `"swald"` or `"wald"` (the chi-squared Wald).

- studentize:

  for `joint = "supt"` or `"swald"`: `TRUE` (default; t critical values)
  or `FALSE` (normal critical values).

- outcome_bwselect:

  the bandwidth rule of
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  and
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md):
  `"pooled"` (default; one bandwidth for every cell, the CCT bandwidth
  of the outcome on the pooled comparison periods), `"cct"` (one CCT
  bandwidth per cell) or `"rot"` (a rule of thumb per cell). The share
  tests' own `bwselect` (per-fit CCT by default) can be changed in
  `typecont` / `compstable`; their rule of thumb is a different one
  (`0.5 * IQR(x)` in every regression).

- typecont, compstable, homog, trendcell:

  lists of further arguments for one test, which override the shared
  ones.

- level:

  the confidence level of the reading in the printout (default 0.95, as
  in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)):
  a test rejects when its p-value is below `1 - level`.

## Value

An object of class `"rd_checks"`: a list with `table` (one row per test:
`test`, `assumption`, `statistic`, `p_value`, `method`, `bandwidth`,
`error`), the four test results `typecont`, `compstable`, `homog`,
`trendcell` (`NULL` where a test failed), `level` and `call`.

## Details

The arguments shared by the four tests (`t_rd`, `comparisons`,
`estimand`, `c`, `kernel`, `scheme`, `bc`) are set once. `joint` and
`studentize` apply to the two tests on type shares, `outcome_bwselect`
to the two tests on the outcome; see the help of each test for what they
do. Anything else, or a different value for one test only, goes in that
test's list, for example `typecont = list(joint = "wald")` or
`homog = list(bwselect = "cct", min_n = 20)`. A test that cannot run on
the data (for example, no unit changes side of the cutoff) gives a row
with its error message and `NA` results.

The printout ends with a reading at the confidence level `level`, as in
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md): a
test rejects its assumption when its p-value is below `1 - level`. Not
rejecting an assumption does not establish it; constant within-type
confounding concerns the RD period and is checked in the comparison
periods only (a pre-trends check).

## See also

[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md);
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) in
[rddid-tidiers](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md).

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md),
[`plot.rd_homog()`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md),
[`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
ch <- rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
ch
#> Tests of the identification assumptions  [rd_checks()]
#> 
#>   Assumption                          Statistic        p   Statistic type             Bandwidth
#>   type continuity                         1.232    0.869   sup-t, t critical values   CCT per cell
#>   composition stability                   3.086    0.016   sup-t, t critical values   CCT per cell
#>   homogeneous confounding                 0.247    0.884   Wald chi-squared(2)        pooled CCT, h = 0.21
#>   constant within-type confounding        1.445    0.485   Wald chi-squared(2)        pooled CCT, h = 0.21
#> 
#>   Reading at level 0.95, rejecting when p < 0.05 (not rejecting does not establish an assumption):
#>     type continuity                    not rejected
#>     composition stability              rejected     } identification needs
#>     homogeneous confounding            not rejected } only one of these two
#>     constant within-type confounding   not rejected (checked in the comparison periods)
#>   Each test's details: x$typecont, x$compstable, x$homog, x$trendcell.
# switch options: the chi-squared Wald and per-cell bandwidths everywhere
rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
          joint = "wald", outcome_bwselect = "cct")
#> Tests of the identification assumptions  [rd_checks()]
#> 
#>   Assumption                          Statistic        p   Statistic type             Bandwidth
#>   type continuity                         4.782    0.853   Wald chi-squared(9)        CCT per cell
#>   composition stability                  23.647   <0.001   Wald chi-squared(6)        CCT per cell
#>   homogeneous confounding                 0.465    0.793   Wald chi-squared(2)        CCT per cell
#>   constant within-type confounding        0.078    0.962   Wald chi-squared(2)        CCT per cell
#> 
#>   Reading at level 0.95, rejecting when p < 0.05 (not rejecting does not establish an assumption):
#>     type continuity                    not rejected
#>     composition stability              rejected     } identification needs
#>     homogeneous confounding            not rejected } only one of these two
#>     constant within-type confounding   not rejected (checked in the comparison periods)
#>   Each test's details: x$typecont, x$compstable, x$homog, x$trendcell.
```
