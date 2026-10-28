# Coefficients, confidence intervals and sample size of an RD-DID fit

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the two estimates
of an
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
fit, `Conventional` (local-linear) and `Robust` (bias-corrected);
[`confint()`](https://rdrr.io/r/stats/confint.html) their confidence
intervals under the fit's sampling scheme;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of observations
used.

## Usage

``` r
# S3 method for class 'didisc'
coef(object, ...)

# S3 method for class 'didisc'
confint(object, parm = c("Conventional", "Robust"), level = NULL, ...)

# S3 method for class 'didisc'
nobs(object, ...)
```

## Arguments

- object:

  an object of class `"didisc"`.

- ...:

  unused.

- parm:

  which rows of the estimate table: `"Conventional"` (local-linear
  estimate, conventional standard error), `"Robust"` (bias-corrected
  estimate, robust standard error), or both (default).

- level:

  confidence level; `NULL` (default) returns the interval stored in the
  object (at the `level` given to
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named numeric vector,
the `Conventional` and `Robust` estimates;
[`confint()`](https://rdrr.io/r/stats/confint.html) a matrix with one
row per `parm` and columns giving the lower and upper limits;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of unit-period
rows used across all periods.

## Examples

``` r
fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
coef(fit)
#> Conventional       Robust 
#>     1.092242     1.141077 
confint(fit)
#>                  2.5 %   97.5 %
#> Conventional 0.8454070 1.339077
#> Robust       0.8493838 1.432769
confint(fit, "Robust", level = 0.9)
#>              5 %     95 %
#> Robust 0.8962803 1.385873
nobs(fit)
#> [1] 3000
```
