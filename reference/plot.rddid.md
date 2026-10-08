# Plot the per-period RD fits behind an RD-DID estimate

One panel per period, in time order: the outcome averaged within `bins`
equal-width bins of the running variable on each side of the cutoff, and
the two local-linear fits of that period drawn over their bandwidth on
each side. The jump between the two lines at the cutoff is the period's
discontinuity \\D_t\\;
[`summary()`](https://rdrr.io/r/base/summary.html) lists every \\D_t\\
and the RD-DID estimate they combine into. Green panels are comparison
periods, pink the RD period.

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

  number of equal-width bins on each side of the cutoff (default 20,
  over the running variable's range in each period).

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed), which can be changed with
`+` and saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

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
