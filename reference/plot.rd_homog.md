# Plot a homogeneous-confounding test: the confounding jump of each type, by comparison period

The paper's homogeneous-confounding figure: for each comparison period,
the local-linear jump in the outcome at the cutoff within each type
(point) with its 95% interval, types side by side. Under the null the
types' jumps coincide within each period.

## Usage

``` r
# S3 method for class 'rd_homog'
plot(x, ...)
```

## Arguments

- x:

  an object returned by
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md).

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed).

## See also

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md),
[`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(hg)
```
