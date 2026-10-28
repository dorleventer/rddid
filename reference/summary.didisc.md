# Summary of a diff-in-disc fit

[`summary()`](https://rdrr.io/r/base/summary.html) of a
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
object adds, to what [`print()`](https://rdrr.io/r/base/print.html)
shows, a per-period table (the local-linear jump in every period with
its bandwidths and sample size, and its signed weight in the estimate)
and the standard error of the estimate under each of the three sampling
schemes.

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
#>   Bandwidth: common h = 0.2694 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4136, 2 = 0.3901, 3 = 0.3984
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0922     0.1259    8.67   <0.001   [0.8454, 1.3391]
#>   Robust (bias-corrected)      1.1411     0.1488    7.67   <0.001   [0.8494, 1.4328]
#> 
#>   Per-period local-linear fits (Conventional = sum of weight x jump; Robust = sum of weight x jump (bc)):
#>   period   role         weight      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison   -0.500   1000   0.2694   0.4136     0.5129    0.1683     0.4642    0.1942
#>   2        comparison   -0.500   1000   0.2694   0.3901     0.6306    0.1679     0.6249    0.2004
#>   3        RD            1.000   1000   0.2694   0.3984     1.6640    0.1635     1.6856    0.1924
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
#>     cross-section 0.2377   panel, no unit changes side 0.1488   panel, some change side 0.1488
summary(fit)$per_period
#>   period       role weight    n         h         b      jump        se
#> 1      1 comparison   -0.5 1000 0.2694485 0.4135886 0.5128736 0.1682545
#> 2      2 comparison   -0.5 1000 0.2694485 0.3900813 0.6305524 0.1678831
#> 3      3         RD    1.0 1000 0.2694485 0.3983944 1.6639551 0.1635322
#>     jump_bc     se_rb
#> 1 0.4641595 0.1942384
#> 2 0.6248868 0.2004214
#> 3 1.6855998 0.1924041
```
