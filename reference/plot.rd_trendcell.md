# Plot a constant-within-type-confounding test: each type's confounding jump over time

One panel per type: its local-linear jump in each comparison period with
its 95% interval, and a dashed reference line, the type's average jump
(`trend = "constant"`) or the least-squares line through its jumps
(`trend = "linear"`). Under the null the points sit on the dashed line,
up to sampling error.

## Usage

``` r
# S3 method for class 'rd_trendcell'
plot(x, ...)
```

## Arguments

- x:

  an object returned by
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md).

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed).

## See also

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md),
[`plot.rd_homog()`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md),
[`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(tr)
```
