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
`per_period`, a data frame with one row per period and columns `period`,
`role` (`"RD"` or `"comparison"`), `coef` (its coefficient in the
estimate, 1 for the RD period and minus its weight for a comparison
period), `n`, `h`, `b`, `jump` and `se` (the conventional local-linear
jump and its standard error) and `jump_bc`, `se_rb` (bias-corrected
jump, robust standard error).

## Examples

``` r
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
summary(fit)
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, running variable fixed over time (detected from the data)
#>   Bandwidth: common h = 0.2672 (rule "joint", AMSE-optimal for the aggregate)
#>   Pilot bandwidth b by period: 0.3951, 0.4102, 0.3868
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0927     0.1264    8.64   <0.001   [0.8450, 1.3405]
#>   Robust (bias-corrected)      1.1414     0.1494    7.64   <0.001   [0.8486, 1.4342]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
#> 
#>   Per-period local-linear fits (estimate = sum of coef x jump):
#>   period   role          coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   3        RD               1   1000   0.2672   0.3951     1.6640    0.1641     1.6854    0.1931
#>   1        comparison    -0.5   1000   0.2672   0.4102     0.5119    0.1688     0.4639    0.1948
#>   2        comparison    -0.5   1000   0.2672   0.3868     0.6306    0.1686     0.6240    0.2012
#> 
#>   Robust s.e. under each sampling scheme:  cross-section 0.2385   panel, fixed R 0.1494   panel, varying R 0.1494
#>   (the printed s.e. uses "pc"; set scheme= to choose another)
summary(fit)$per_period
#>   period       role coef    n         h         b      jump        se   jump_bc
#> 1      3         RD  1.0 1000 0.2672096 0.3950840 1.6639868 0.1641410 1.6853502
#> 2      1 comparison -0.5 1000 0.2672096 0.4101520 0.5119443 0.1688159 0.4639253
#> 3      2 comparison -0.5 1000 0.2672096 0.3868401 0.6306109 0.1685612 0.6240145
#>       se_rb
#> 1 0.1931055
#> 2 0.1948448
#> 3 0.2011668
```
