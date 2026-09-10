# Test of a constant within-type confounding discontinuity

Pre-trends check for the constant within-type confounding assumption
(Section 4.4 of Leventer and Nevo), in the difference-in-differences
sense: the assumption concerns the RD period, where the confounding is
not separately observed, so the test asks whether the per-cell outcome
RD discontinuity is constant (or linear) **across the comparison
periods**. In a comparison period \\t_0\\ the discontinuity equals the
confounding, \\D\_{t_0}(k) = \alpha\_{t_0,0}(k)\\, where \\k\\ is the
unit's cell (its side of the cutoff in \\t\_{\mathrm{RD}}\\ under the
default `type_by = "rd_side"`). The function estimates the jump per cell
via
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
then forms a joint Wald test that the within-cell jumps conform to the
hypothesised trend \\g_0\\ across comparison periods.

## Usage

``` r
rd_trendcell(
  data,
  y,
  x,
  time,
  id,
  comparisons = NULL,
  t_rd = NULL,
  estimand = c("att", "atu"),
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  min_n = 10L,
  bc = TRUE,
  type_by = c("rd_side", "pattern"),
  trend = c("constant", "linear"),
  ...
)
```

## Arguments

- data:

  a long data frame, one row per unit-period. A unit's type in period
  \\t\\ is read from its running variable in the other period(s); units
  unobserved there are dropped from period \\t\\, so the panel need not
  be balanced.

- y, x, time:

  Column names (character strings) for the outcome, running variable,
  and period indicator.

- id:

  Column name for the unit identifier. Required (the test needs a panel
  to define a fixed cell per unit).

- comparisons:

  Values of `time` to use as comparison periods. Defaults to all periods
  except `t_rd` (if `t_rd` is supplied) or all periods (if
  `t_rd = NULL`).

- t_rd:

  Value of `time` for the RD period. Required for `type_by = "rd_side"`.
  Under `type_by = "pattern"`, if supplied, the cell is the sign pattern
  of all comparison periods from \\t\_{\mathrm{RD}}\\'s perspective; if
  `NULL`, the pattern is taken from the first comparison period's
  perspective.

- estimand:

  `"att"` (default) or `"atu"`. Label only: under `"atu"` the
  within-type comparison-period discontinuities are the confounding
  discontinuities among TREATED units, \\\alpha\_{t_0,1}(v)\\; the
  estimates and test are numerically identical to the `"att"` call.

- c:

  Cutoff value for the running variable (default 0).

- h:

  Main bandwidth. If `NULL` (default), bandwidth is chosen according to
  `bwselect`. An explicit numeric value overrides `bwselect` and is used
  directly for every (cell, period) combination.

- bwselect:

  Bandwidth selection when `h = NULL`: `"cct"` (default) computes a
  per-(cell, period) CCT MSE-optimal bandwidth via
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  on that cell's outcome and running variable; `"rot"` uses the 0.2 ×
  range rule of thumb (the original behavior). Ignored when `h` is
  supplied.

- kernel:

  Kernel for the local-linear RD: `"triangular"` (default),
  `"epanechnikov"`, or `"uniform"`.

- scheme:

  Sampling scheme for the cross-period covariance: `"auto"` (detect from
  the id/side structure, default), `"cs"` (repeated cross-section, no
  cross-period terms), `"pc"` (panel, time-constant running variable),
  or `"pv"` (panel, time-varying running variable).

- min_n:

  Minimum number of observations per cell-side before a (cell, period)
  pair is included. Pairs with fewer than `min_n` obs on either side are
  silently dropped (default 10).

- bc:

  Use robust bias-corrected per-cell jumps and variances (Calonico,
  Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test with the
  bias-corrected
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  estimator; `FALSE` uses the conventional local-linear jumps and
  variances.

- type_by:

  How a unit's type is defined. `"rd_side"` (default) = the unit's side
  of the cutoff in the RD period, the partition of the paper's Section
  4.4. `"pattern"` = the sign pattern of the other periods' running
  variables. The cell is **fixed** across comparison periods for both
  choices.

- trend:

  Trend form for the null hypothesis. `"constant"` (default) tests equal
  per-cell jumps across all comparison periods; `"linear"` tests that
  the per-cell jumps lie on a line in time, using second-difference
  contrasts (requires \\\geq 3\\ comparison periods per cell).

- ...:

  Further arguments passed to
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  (e.g. `p`, `q`, `b`).

## Value

An object of class `"rd_trendcell"`, a list with:

- `statistic`:

  Joint Wald chi-squared statistic (`NA` when `df = 0`).

- `df`:

  Degrees of freedom (number of positive eigenvalues used). `0` when
  `trend = "linear"` and no cell has \\\geq 3\\ comparison periods.

- `p_value`:

  p-value from the chi-squared distribution (`NA` when `df = 0`).

- `cell_period_jumps`:

  Data frame with one row per (cell, comparison period): cell, period,
  jump estimate, SE, number of observations, and a flag indicating
  whether this is the reference period for that cell (under
  `trend = "constant"` only).

- `contrasts`:

  Named numeric vector of jump contrasts stacked across cells.

- `cov_matrix`:

  Estimated covariance matrix of the contrasts; block-diagonal by cell.

- `scheme`:

  Sampling scheme used.

- `bc`:

  Whether bias-corrected jumps were used.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `trend`:

  Trend form used (`"constant"` or `"linear"`).

- `comparisons`:

  Comparison periods actually used (character).

- `call`:

  The matched call.

**`trend = "linear"`** requires at least 3 comparison periods per cell
to be informative. With only 2 comparison periods the within-cell linear
trend is just-identified (any two points define a line), so no
second-difference contrast exists and the test returns `df = 0`.

## Details

### Null hypothesis

\$\$H_0 : D\_{t_0}(k) \text{ is } g_0 \text{ in } t_0, \text{ jointly
for all cells } k,\$\$ where \\g_0\\ is `trend = "constant"` (equal
jumps across all comparison periods) or `trend = "linear"` (jumps lie on
a line in \\t_0\\).

Under `trend = "constant"`, the contrasts for cell \\k\\ with
\\\|T_0\|\\ comparison periods are \\D\_{t_0}(k) - D\_{t_1}(k)\\ for
\\t_0 \neq t_1\\ (reference = first comparison period), yielding
\\\|T_0\| - 1\\ contrasts per cell.

Under `trend = "linear"`, the testable contrasts are the **second
differences** of the time-ordered per-cell jumps: \\D\_{t\_{j+1}}(k) - 2
D\_{t_j}(k) + D\_{t\_{j-1}}(k)\\, yielding \\\|T_0\| - 2\\ contrasts per
cell. **This requires at least 3 comparison periods per cell**; if no
cell reaches this threshold the function returns an object with
`df = 0`, `statistic = NA`, and a message.

## References

Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
Designs." Working paper.

Calonico, S., Cattaneo, M. D., and Titiunik, R. (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica*, 82(6), 2295-2326.

## See also

[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Three-period panel: period 3 is RD, periods 1 and 2 are comparison.
# Constant within-cell confounding across periods 1 and 2 (H0 holds).
set.seed(1)
n <- 600
r_rd   <- runif(n, -1, 1)
r_comp <- runif(n, -1, 1)   # same running variable both comparison periods
cell   <- ifelse(r_rd >= 0, "+", "-")
# constant confounding: same jump in both comparison periods within each cell
alpha  <- ifelse(cell == "+", 0.5, 0.3)
y1 <- 0.3 * r_comp + alpha * (r_comp >= 0) + rnorm(n, 0, 0.3)
y2 <- 0.3 * r_comp + alpha * (r_comp >= 0) + rnorm(n, 0, 0.3)
y3 <- 0.3 * r_rd + 1.0 * (r_rd >= 0) + rnorm(n, 0, 0.3)  # RD period
dat <- data.frame(
  id   = rep(seq_len(n), 3),
  time = rep(1:3, each = n),
  x    = c(r_comp, r_comp, r_rd),
  y    = c(y1, y2, y3)
)
rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
             comparisons = 1:2, t_rd = 3, h = 0.3)
} # }
```
