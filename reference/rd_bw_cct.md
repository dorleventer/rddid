# CCT bandwidths for one period (building block)

The building block behind `rddid(bwselect = "cct")`, the pilot fits of
the `"joint"` and `"iter"` rules, and the default bandwidths of the four
tests of the assumptions: the MSE-optimal main and pilot bandwidths of
Calonico, Cattaneo and Titiunik (2014) for a local-linear RD in one
period, from
[`rdrobust::rdbwselect()`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
with `bwselect = "mserd"`. Like
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
it takes vectors, not a data frame and column names.

## Usage

``` r
rd_bw_cct(y, x, c = 0, p = 1L, kernel = "triangular")
```

## Arguments

- y:

  the outcome (a numeric vector).

- x:

  the running variable (a numeric vector, same length as `y`).

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- p:

  order of the local polynomial (default 1, local linear).

- kernel:

  the kernel: `"triangular"` (default), `"epanechnikov"` or `"uniform"`.

## Value

A named numeric vector `c(h = , b = )`: the main bandwidth (point
estimate) and the pilot bandwidth (bias correction).

## Details

If `rdbwselect()` fails, or returns a main bandwidth that is not a
positive number, `rd_bw_cct()` falls back to `h = b = 0.5 * IQR(x)`
(`sd(x)` if the interquartile range is zero) and says so in a message.

## References

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295-2326.

## See also

Other RD-DID estimation:
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
# each year of rddid_sim gets its own bandwidths (what bwselect = "cct" uses)
sapply(split(rddid_sim, rddid_sim$year), function(d) rd_bw_cct(y = d$Y, x = d$R))
#>           1         2         3
#> h 0.2637172 0.3255516 0.4140584
#> b 0.4047913 0.4713019 0.6122080
```
