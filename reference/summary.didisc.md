# Summary of a diff-in-disc fit

[`summary()`](https://rdrr.io/r/base/summary.html) of a
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
object adds, to what [`print()`](https://rdrr.io/r/base/print.html)
shows, a per-period table (the local-linear jump in every period with
its standard error, conventional and bias-corrected, its bandwidths and
sample size) and a table of the standard error of the estimate,
conventional and robust, under each of the three sampling schemes (RCS,
PC, PV).

## Usage

``` r
# S3 method for class 'didisc'
summary(object, ...)

# S3 method for class 'summary.didisc'
print(x, digits = 4, ...)
```

## Arguments

- object:

  an object of class `"didisc"`.

- ...:

  unused.

- x:

  a `"summary.didisc"` object.

- digits:

  number of decimals in the printed tables.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.didisc"`: a list with `fit` (the object) and
`per_period`, a data frame with one row per period in time order and
columns `period`, `role` (`"RD"` or `"comparison"`), `weight` (the
period's signed weight in the conventional estimate: the RD period's
weight in `rd_weights`, 1 for a single RD period, and minus its weight
for a comparison period), `n`, `h`, `b`, `jump` and `se` (the
conventional local-linear jump and its standard error), `jump_bc`,
`se_rb` (bias-corrected jump, robust standard error) and, under
`weighting = "min_variance"`, `weight_bc` (the signed weights of the
bias-corrected estimate, which has its own weights).

## Examples

``` r
fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
summary(fit)
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.3986 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b = 0.5936 (also chosen for the diff-in-disc estimate)
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0075     0.1062    9.49   <0.001   [0.7995, 1.2156]
#>   Robust (bias-corrected)      1.0247     0.1256    8.16   <0.001   [0.7785, 1.2708]
#> 
#>   Per-period local-linear fits:
#>                                                        Conventional         Robust (bias-corrected)
#>   Period   Role             n        h        b          Jump   Std. err.          Jump   Std. err.
#>   1        comparison    1000   0.3986   0.5936        0.5768      0.1393        0.5716      0.1672
#>   2        comparison    1000   0.3986   0.5936        0.6308      0.1387        0.6558      0.1655
#>   3        RD            1000   0.3986   0.5936        1.6113      0.1364        1.6383      0.1624
#> 
#>   Std. err. of the estimate under each sampling scheme (the estimate table at the top uses PC):
#>                                                    Conventional   Robust (bias-corrected)
#>   RCS  repeated cross-section                            0.1681                    0.2005
#>   PC   panel, no unit changes side of the cutoff         0.1062                    0.1256
#>   PV   panel, some units change side of the cutoff       0.1062                    0.1256
summary(fit)$per_period
#>   period       role weight    n         h       b      jump        se   jump_bc
#> 1      1 comparison   -0.5 1000 0.3985921 0.59356 0.5767858 0.1393244 0.5715896
#> 2      2 comparison   -0.5 1000 0.3985921 0.59356 0.6308135 0.1387137 0.6557607
#> 3      3         RD    1.0 1000 0.3985921 0.59356 1.6113369 0.1363900 1.6383291
#>       se_rb
#> 1 0.1671519
#> 2 0.1654705
#> 3 0.1623620
```
