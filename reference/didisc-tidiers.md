# Tidy output for diff-in-disc fits and validation tests

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods
following the `broom` conventions, so that fits and tests can be passed
to table makers such as `modelsummary`. The generics come from the
`generics` package; the methods are registered when that package (or
`broom`) is loaded.

## Usage

``` r
# S3 method for class 'didisc'
tidy(x, conf.int = TRUE, conf.level = NULL, ...)

# S3 method for class 'didisc'
glance(x, ...)

# S3 method for class 'rd_typecont'
tidy(x, ...)

# S3 method for class 'rd_compstable'
tidy(x, ...)

# S3 method for class 'rd_homog'
tidy(x, ...)

# S3 method for class 'rd_trendcell'
tidy(x, ...)

# S3 method for class 'rd_checks'
tidy(x, ...)
```

## Arguments

- x:

  an object returned by
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md),
  [`rd_typecont()`](https://dorleventer.github.io/didisc/reference/rd_typecont.md),
  [`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md),
  [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md)
  or
  [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md).

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
per test for a validation test: `test` (the function, as in
[`rd_checks()`](https://dorleventer.github.io/didisc/reference/rd_checks.md)),
`assumption`, `statistic` (the largest \|t\| under `joint = "supt"` for
[`rd_typecont()`](https://dorleventer.github.io/didisc/reference/rd_typecont.md)
and
[`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md),
whose row is the joint test over pairs; the chi-squared statistic for
[`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md)
and
[`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)),
`df` (the chi-squared degrees of freedom of a Wald statistic; `NA` for
sup-t and sWald), `p.value` and `method`.
[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
one-row data frame describing the fit (`nobs`, `t_rd`, `rd_weights`,
`comparisons`, `trend`, `weights`, `bwselect`, `h` (the common bandwidth
under `"joint"`/fixed `h`, `NA` otherwise), `scheme`, `level`,
`weighting`).

## Examples

``` r
fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
tidy(fit)
#>           term estimate std.error statistic      p.value  conf.low conf.high
#> 1 Conventional 1.092242 0.1259386  8.672816 4.215617e-18 0.8454070  1.339077
#> 2       Robust 1.141077 0.1488256  7.667206 1.757833e-14 0.8493838  1.432769
glance(fit)
#>   nobs t_rd comparisons rd_weights    trend  weights bwselect         h scheme
#> 1 3000    3        1, 2          1 constant 0.5, 0.5    joint 0.2694485     pc
#>   level weighting
#> 1  0.95       ols
tidy(rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3))
#>          test      assumption statistic df   p.value                   method
#> 1 rd_typecont type continuity  1.231775 NA 0.8688813 sup-t, t critical values
```
