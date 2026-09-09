# Test the continuity of the type distribution

Wald test of the continuous-type-distribution assumption (Section 4.4 of
Leventer and Nevo). The "type" of unit \\i\\ in period \\t\\ is the sign
pattern of its running variables in the OTHER periods,
\\\mathbf{V}\_{i,-t} = (1\\R\_{i,s} \ge c\\)\_{s \ne t}\\.

## Usage

``` r
rd_typecont(
  data,
  x,
  time,
  id,
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  bc = TRUE,
  ...
)
```

## Arguments

- data:

  a long data frame, one row per unit-period. A unit's type in period
  \\t\\ is read from its running variable in the other period(s); units
  unobserved there are dropped from period \\t\\, so the panel need not
  be balanced.

- x:

  column name (string) for the running variable.

- time:

  column name (string) for the period.

- id:

  column name (string) for the unit identifier.

- c:

  cutoff (default 0).

- h:

  bandwidth. If `NULL`, the bandwidth is determined by `bwselect`; an
  explicit numeric value overrides `bwselect` and is used directly.

- bwselect:

  bandwidth selection rule when `h = NULL`: `"cct"` (default) computes a
  per-cell CCT MSE-optimal bandwidth via
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  for each (period, type) RD; `"rot"` uses the `0.5 * IQR(x)` rule of
  thumb applied to the full sample (the previous default behaviour).
  Ignored when `h` is supplied explicitly.

- kernel:

  kernel for the local-linear RD: `"triangular"` (default),
  `"epanechnikov"`, or `"uniform"`.

- scheme:

  covariance scheme for the joint Wald: `"auto"` detects from the data
  (same logic as
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)),
  or one of `"cs"`, `"pc"`, `"pv"`.

- bc:

  use robust bias-corrected jumps and variances in the LL-Wald
  (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the
  test with the bias-corrected
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  estimator; `FALSE` uses the conventional local-linear jumps and
  variances.

- ...:

  currently unused.

## Value

An object of class `"rd_typecont"`, a named list with:

- ll_wald:

  list with `stat` (chi-square), `df`, `p`.

- per_period:

  named list (by period) of the per-period components the joint test
  aggregates; each entry has `ll_wald` (that period's own LL-Wald
  `stat`/`df`/`p`, restricted to its kept contrasts).

- meta:

  list with `periods`, `type_values`, `h` (NA when `bwselect = "cct"`),
  `bwselect`, `scheme`, `bc`.

## Details

For each period \\t\\ and each type value \\v\\, estimate the
local-linear RD jump of the type indicator \\1\\\mathbf{V}\_{i,-t} =
v\\\\ on the running variable \\R\_{i,t}\\. The jump
\\\hat\pi\_{t,(+)}(v) - \hat\pi\_{t,(-)}(v)\\ is the output of
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
with a binary outcome. The joint Wald statistic across all (period,
type) pairs drops one reference type per period, because within each
period the type indicators sum to 1 (the full block is singular); df =
number of kept contrasts. The covariance is built from the per-unit
influence vectors returned by
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
using the same within-period and cross-period id-matching as the main
estimator, scheme-aware.
