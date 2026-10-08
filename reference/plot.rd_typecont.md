# Plot a type-continuity test: one comparison period against one RD period

The paper's type-continuity figure. Two panels: in the comparison
period, the share of units that are above the cutoff in the RD period,
against the comparison-period running variable; in the RD period, the
share that are above the cutoff in the comparison period, against the
RD-period running variable. In each panel the share is averaged within
`bins` equal-width bins on each side of the cutoff inside the bandwidth,
and the local-linear fit is drawn on each side, at the bandwidth rule
the test used. Under the null the two lines of a panel meet at the
cutoff: who a unit is in the other period does not jump there.

## Usage

``` r
# S3 method for class 'rd_typecont'
plot(x, t_rd = NULL, comparison = NULL, bins = 20L, ...)
```

## Arguments

- x:

  an object returned by
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md).

- t_rd, comparison:

  the RD period and the comparison period to draw (values of the period
  variable). Defaults: the `t_rd` given to
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  (else the last period), and the first other period.

- bins:

  number of equal-width bins on each side of the cutoff (default 20).

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed).

## Details

With more than two periods the test itself uses the full sign pattern of
the other periods as the type; the picture shows the pairwise version,
one pair of periods at a time.

## See also

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md),
[`plot.rd_homog()`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md),
[`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(tc, comparison = 1)
```
