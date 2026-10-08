# Plot the switchers: the running variable in one period against another

Each unit observed in both periods is a point; the dashed lines are the
cutoff. Units in the off-diagonal quadrants changed side of the cutoff
between the two periods (the "switchers" that the composition tests are
about); the title gives their share.

## Usage

``` r
plot_switchers(data, x, time, id, periods = NULL, c = 0, ...)
```

## Arguments

- data:

  a data frame in long format, one row per unit and period: a repeated
  cross-section (different units in each period) or a panel (the same
  units in several periods; it need not be balanced).

- x:

  name of the running-variable column (a string).

- time:

  name of the period column (a string).

- id:

  name of the unit-identifier column (a string), needed for panel
  standard errors. With `NULL` (default) every row is treated as a
  different unit, which gives repeated cross-section standard errors
  (with a message).

- periods:

  the two periods to compare (values of `time`); default the first two
  in `data`.

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- ...:

  unused.

## Value

A ggplot object (ggplot2 must be installed).

## See also

Other tests of the assumptions:
[`plot.rd_compstable()`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md),
[`plot.rd_homog()`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md),
[`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md),
[`plot.rd_typecont()`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md),
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE))
  plot_switchers(rddid_sim_pv, x = "R", time = "year", id = "id", periods = c(1, 3))
```
