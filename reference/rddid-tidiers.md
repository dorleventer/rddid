# Tidy output for RD-DID fits and validation tests

`tidy()` and `glance()` methods following the `broom` conventions, so
that fits and tests can be passed to table makers such as
`modelsummary`. The generics come from the `generics` package; the
methods are registered when that package (or `broom`) is loaded.

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

  for `tidy()` of a fit: include the confidence interval (default
  `TRUE`) and at which level (default the fit's `level`; another value
  recomputes the interval from the estimate and its standard error), as
  `modelsummary` and other table makers pass them.

- ...:

  unused.

## Value

`tidy()` returns a data frame with one row per estimate (`term`,
`estimate`, `std.error`, `statistic`, `p.value`, `conf.low`,
`conf.high`) for a fit, and one row per test (`test`, `statistic`, `df`,
`p.value`) for a validation test; for
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
that row is the joint test over pairs, whose p-value is approximate (see
its Details). `glance()` returns a one-row data frame describing the fit
(`nobs`, `t_rd`, `comparisons`, `trend`, `weights`, `bwselect`, `h` (the
common bandwidth under `"joint"`/fixed `h`, `NA` otherwise), `scheme`,
`level`).

## Examples

``` r
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
if (requireNamespace("generics", quietly = TRUE)) {
  generics::tidy(fit)
  generics::glance(fit)
  generics::tidy(rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3))
}
#>              test statistic df   p.value
#> 1 type continuity  4.781633  9 0.8529133
```
