# Test of homogeneous confounding

Wald test of the homogeneous-confounding assumption (Section 4.4 of
Leventer and Nevo): in the comparison periods, the outcome RD
discontinuity is the same across types. In a comparison period \\t_0\\
the discontinuity equals the confounding, \\D\_{t_0}(v) =
\alpha\_{t_0,0}(v)\\, with \\v\\ the unit's type (by default its side of
the cutoff in the RD period). The function estimates the jump by type
via
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
then forms a joint Wald test that the jumps are equal across types and
comparison periods.

## Usage

``` r
rd_homog(
  data,
  y,
  x,
  time,
  id,
  comparisons = NULL,
  t_rd = NULL,
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  min_n = 10L,
  bc = TRUE,
  type_by = c("rd_side", "pattern"),
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

  Column name for the unit identifier. Required (the test uses the other
  periods' running variables to assign types, which needs a panel).

- comparisons:

  Values of `time` to use as comparison periods. Defaults to all periods
  except `t_rd` (if `t_rd` is supplied) or all periods (if
  `t_rd = NULL`).

- t_rd:

  Value of `time` for the RD period. Used only to exclude it from
  comparison periods when `comparisons = NULL`; the RD period is **not**
  used in the test itself.

- c:

  Cutoff value for the running variable (default 0).

- h:

  Main bandwidth. If `NULL` (default), bandwidth is chosen according to
  `bwselect`. An explicit numeric value overrides `bwselect` and is used
  directly for every cell.

- bwselect:

  Bandwidth selection when `h = NULL`: `"cct"` (default) computes a
  per-cell CCT MSE-optimal bandwidth via
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

  Minimum number of observations per type-side before that type is
  included. Types with fewer than `min_n` obs on either side in a given
  period are silently dropped (default 10).

- bc:

  Use robust bias-corrected per-type jumps and variances (Calonico,
  Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test with the
  bias-corrected
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  estimator; `FALSE` uses the conventional local-linear jumps and
  variances.

- type_by:

  How a unit's type is defined. `"rd_side"` (default) = the unit's side
  of the cutoff in the RD period, the partition of the paper's Section
  4.4. `"pattern"` = the sign pattern of the other periods' running
  variables.

- ...:

  Further arguments passed to
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  (e.g. `p`, `q`, `b`).

## Value

An object of class `"rd_homog"`, a list with:

- `statistic`:

  Joint Wald chi-squared statistic.

- `df`:

  Degrees of freedom (total contrasts used).

- `p_value`:

  p-value from the chi-squared distribution.

- `period_type_jumps`:

  Data frame with one row per (comparison period, type) cell: period,
  type, jump estimate, SE, number of observations, and whether it was
  the reference type.

- `contrasts`:

  Named numeric vector of jump contrasts (non-reference minus
  reference), stacked across periods.

- `cov_matrix`:

  Estimated covariance matrix of the contrasts.

- `scheme`:

  Sampling scheme used.

- `comparisons`:

  Comparison periods actually used.

- `call`:

  The matched call.

## Details

### Null hypothesis

\\H_0 :\\ the outcome RD jump is equal across all type cells
\\\mathbf{v}\_{-t_0}\\, jointly across all comparison periods:
\\D\_{t_0}(\mathbf{v}) = D\_{t_0}(\mathbf{v}') \\ \forall \\ \mathbf{v},
\mathbf{v}', t_0\\.

The Wald statistic is formed from the per-type jump contrasts \\\Delta =
D\_{t_0,v} - D\_{t_0,\mathrm{ref}}\\, stacked across comparison periods,
with covariance estimated from the
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
influence vectors. Within a period, cross-type covariance is zero (types
partition the sample). Across periods, id-matched covariance is used
under the detected/requested sampling scheme (see `scheme`).

## References

Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
Designs." Working paper.

Imbens, G. W. and Lemieux, T. (2008). Regression discontinuity designs:
A guide to practice. *Journal of Econometrics*, 142(2), 615-635.

Calonico, S., Cattaneo, M. D., and Titiunik, R. (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica*, 82(6), 2295-2326.

## See also

[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-period panel: period 1 is RD, period 2 is comparison.
# Type in period 2 = sign of R_{i,1} (the running variable of the other period).
set.seed(1)
n <- 500
r1 <- runif(n, -1, 1)
r2 <- runif(n, -1, 1)
# homogeneous confounding in comparison period (period 2)
y2 <- 0.3 * r2 + 0.4 * (r2 >= 0) + rnorm(n, 0, 0.3)
dat <- data.frame(
  id   = rep(seq_len(n), 2),
  time = rep(1:2, each = n),
  x    = c(r1, r2),
  y    = c(rnorm(n), y2)
)
rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
         comparisons = 2, t_rd = 1, h = 0.3)
} # }
```
