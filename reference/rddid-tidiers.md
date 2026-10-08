# Tidy output for RD-DID fits and validation tests

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods
following the `broom` conventions, so that fits and tests can be passed
to table makers such as `modelsummary`. The generics come from the
`generics` package; the methods are registered when that package (or
`broom`) is loaded.

## Usage

``` r
# S3 method for class 'rddid'
tidy(x, conf.int = TRUE, conf.level = NULL, ...)

# S3 method for class 'rddid'
glance(x, ...)

# S3 method for class 'rd_typecont'
tidy(x, ...)

# S3 method for class 'rd_compstable'
tidy(x, ...)

# S3 method for class 'rd_homog'
tidy(x, ...)

# S3 method for class 'rd_trendcell'
tidy(x, ...)
```

## Arguments

- x:

  an object returned by
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md),
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  or
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md).

- conf.int, conf.level:

  for [`tidy()`](https://generics.r-lib.org/reference/tidy.html) of a
  fit: include the confidence interval (default `TRUE`) and at which
  level (default the fit's `level`; another value recomputes the
  interval from the estimate and its standard error), as `modelsummary`
  and other table makers pass them.

- ...:

  unused.

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
data frame with one row per estimate (`term`, `estimate`, `std.error`,
`statistic`, `p.value`, `conf.low`, `conf.high`) for a fit, and one row
per test (`test`, `statistic`, `df`, `p.value`) for a validation test;
for
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
that row is the joint test over pairs, whose p-value is approximate (see
its Details).
[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
one-row data frame describing the fit (`nobs`, `t_rd`, `comparisons`,
`trend`, `weights`, `bwselect`, `h` (the common bandwidth under
`"joint"`/fixed `h`, `NA` otherwise), `scheme`, `level`).

## Examples

``` r
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
tidy(fit)
#>           term estimate std.error statistic      p.value  conf.low conf.high
#> 1 Conventional 1.092709 0.1264043  8.644556 5.401507e-18 0.8449613  1.340457
#> 2       Robust 1.141380 0.1493939  7.640071 2.171020e-14 0.8485735  1.434187
glance(fit)
#>   nobs t_rd comparisons    trend  weights bwselect         h scheme level
#> 1 3000    3        1, 2 constant 0.5, 0.5    joint 0.2672096     pc  0.95
tidy(rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3))
#>              test statistic df   p.value
#> 1 type continuity  4.781633  9 0.8529133
```
