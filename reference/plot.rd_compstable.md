# Plot a composition-stability test: the reflected sample of one pair of periods

The paper's composition-stability figure. The units above the cutoff in
the comparison period are placed to the left of an artificial cutoff at
their mirrored distance \\-(R\_{t_0} - c)\\, the units above the cutoff
in the RD period to the right at \\R\_{t\_{RD}} - c\\; the outcome is
whether the unit is above the cutoff in the other period of the pair.
The binned share is drawn with the local-linear fit on each side, at the
bandwidth rule the test used. Under the null the two lines meet at the
artificial cutoff: the units just above the cutoff are the same mix in
both periods. With `estimand = "atu"` the same picture is drawn for the
units below the cutoff.

## Usage

``` r
# S3 method for class 'rd_compstable'
plot(x, pair = 1L, bins = 20L, ...)
```

## Arguments

- x:

  an object returned by
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md).

- pair:

  which pair to draw: an index or a name of `x$pairs` (default the
  first).

- bins:

  number of equal-width bins for the binned shares (default 20, over the
  bandwidth on each side).

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed).

## Details

With more than two periods the test's types also record the unit's sides
in the remaining periods; the picture shows the pairwise version.

## See also

Other tests of the assumptions:
[`plot.rd_homog()`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md),
[`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(cs)
```
