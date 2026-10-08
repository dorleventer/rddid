# Plot the per-period RD fits behind an RD-DID estimate

One panel per period, the RD period first: the outcome averaged within
`bins` equal-width bins of the running variable, and the two
local-linear fits of that period drawn over their bandwidth on each side
of the cutoff. The jump between the two lines at the cutoff is the
period's discontinuity \\D_t\\; the RD-DID estimate is the RD-period
jump minus the weighted comparison-period jumps (see
[`summary()`](https://rdrr.io/r/base/summary.html)).

## Usage

``` r
# S3 method for class 'rddid'
plot(x, bins = 20L, ...)
```

## Arguments

- x:

  an object returned by
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).

- bins:

  number of equal-width bins for the binned means (default 20, over the
  running variable's range in each period).

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed).

## See also

Other RD-DID estimation:
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md),
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
```
