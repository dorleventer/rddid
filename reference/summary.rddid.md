# Summary of an RD-DID fit

[`summary()`](https://rdrr.io/r/base/summary.html) of an
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
object adds, to what [`print()`](https://rdrr.io/r/base/print.html)
shows, a per-period table (the local-linear jump in every period with
its bandwidths and sample size, and its coefficient in the estimate) and
the standard error of the estimate under each of the three sampling
schemes.

## Usage

``` r
# S3 method for class 'rddid'
summary(object, ...)

# S3 method for class 'summary.rddid'
print(x, digits = 4, ...)
```

## Arguments

- object:

  an object of class `"rddid"`.

- ...:

  unused.

- x:

  a `"summary.rddid"` object.

- digits:

  number of decimals in the printed tables.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.rddid"`: a list with `fit` (the object) and
`per_period`, a data frame with one row per period in time order and
columns `period`, `role` (`"RD"` or `"comparison"`), `coef` (the
period's coefficient in the conventional estimate: the RD period's
weight in `rd_weights`, 1 for a single RD period, and minus its weight
for a comparison period), `n`, `h`, `b`, `jump` and `se` (the
conventional local-linear jump and its standard error), `jump_bc`,
`se_rb` (bias-corrected jump, robust standard error) and, under
`weighting = "min_variance"`, `coef_bc` (the coefficients of the
bias-corrected estimate, which has its own weights).

## Examples

``` r
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
summary(fit)
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2734 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4197, 2 = 0.3958, 3 = 0.4043
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0903     0.1251    8.72   <0.001   [0.8451, 1.3355]
#>   Robust (bias-corrected)      1.1387     0.1478    7.70   <0.001   [0.8489, 1.4284]
#> 
#>   Per-period local-linear fits (Conventional = sum of coef x jump; Robust = sum of coef x jump (bc)):
#>   period   role           coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison   -0.500   1000   0.2734   0.4197     0.5152    0.1672     0.4653    0.1932
#>   2        comparison   -0.500   1000   0.2734   0.3958     0.6308    0.1667     0.6269    0.1991
#>   3        RD            1.000   1000   0.2734   0.4043     1.6633    0.1624     1.6848    0.1912
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
#>     cross-section 0.2362   panel, no unit changes side 0.1478   panel, some change side 0.1478
summary(fit)$per_period
#>   period       role coef    n         h         b      jump        se   jump_bc
#> 1      1 comparison -0.5 1000 0.2734096 0.4196687 0.5152062 0.1672386 0.4652979
#> 2      2 comparison -0.5 1000 0.2734096 0.3958158 0.6307688 0.1666649 0.6268708
#> 3      3         RD  1.0 1000 0.2734096 0.4042511 1.6633103 0.1624416 1.6847733
#>       se_rb
#> 1 0.1932246
#> 2 0.1990750
#> 3 0.1911622
```
