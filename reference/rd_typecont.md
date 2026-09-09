# Test the continuity of the type distribution (assumption `ass:type-cont`)

Runs four tests for the continuous-type-distribution assumption
(`ass:type-cont`, Section 4.4 of Leventer and Nevo). Section 4.4 states
the LL-Wald; the McCrary variants are run in the paper's simulation
appendix and the permutation test in its Section 6 validation table (see
`dev/tests_map.md`). The "type" of unit \\i\\ in period \\t\\ is the
sign pattern of its running variables in the OTHER periods,
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
  q = NULL,
  S = 499L,
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  bc = TRUE,
  ...
)
```

## Arguments

- data:

  a long data frame, one row per unit × period (balanced panel).

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

- q:

  number of observations nearest the cutoff on each side for the
  Canay-Kamat permutation test. `NULL` (default) selects `q` per period
  by the Canay & Kamat (2018) rule of thumb; this is the recommended
  choice, since a fixed `q` over-rejects in finite samples when the type
  distribution varies steeply in the running variable at the cutoff.
  Supply an integer to force a fixed `q` on every period. The per-period
  `q` actually used is returned in `meta$q_used`.

- S:

  number of permutation replications for the Canay-Kamat test (default
  499).

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
  variances. Does not affect the Canay-Kamat permutation test (which
  operates on raw near-cutoff indicator means).

- ...:

  currently unused.

## Value

An object of class `"rd_typecont"`, a named list with:

- ll_wald:

  list with `stat` (chi-square), `df`, `p`.

- per_period:

  named list (by period) of the per-period components the joint test
  aggregates; each entry has `ll_wald` (that period's own LL-Wald
  `stat`/`df`/`p`, restricted to its kept contrasts) and `ck_p`
  (per-period Canay-Kamat p-value, read off the same permutation draws
  as the joint statistic, or `NA` for a period with no active type
  cell).

- ck_perm:

  list with `stat` (observed sum of \|mean diffs\|), `p`.

- mccrary_within:

  data frame with columns `period`, `type`, `p_raw` (per-type McCrary
  p-value); plus `p_bonf` (Bonferroni-adjusted p-value for the period,
  minimum × number of types).

- mccrary_pooled:

  data frame with columns `period`, `p` (per-period McCrary p-value on
  pooled sample).

- meta:

  list with `periods`, `type_values`, `h` (NA when `bwselect = "cct"`),
  `bwselect`, `q` (`"rot"` when the Canay-Kamat rule of thumb is used,
  otherwise the integer supplied), `q_used` (per-period integer `q`
  actually used), `S`, `scheme`, `bc`.

## Details

**Four tests, in order of interpretive strength:**

1.  **LL-Wald** (necessary AND sufficient): for each period \\t\\ and
    each type value \\v\\, estimate the local-linear RD jump of the type
    indicator \\1\\\mathbf{V}\_{i,-t} = v\\\\ on the running variable
    \\R\_{i,t}\\. The jump \\\hat\pi\_{t,(+)}(v) - \hat\pi\_{t,(-)}(v)\\
    is the output of
    [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
    with a binary outcome. The joint Wald statistic across all (period,
    type) pairs drops one reference type per period, because within each
    period the type indicators sum to 1 (the full block is singular); df
    = number of kept contrasts. The covariance is built from the
    per-unit influence vectors returned by
    [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
    using the same within-period and cross-period id-matching as the
    main estimator, scheme-aware.

2.  **Canay-Kamat permutation** (necessary AND sufficient): uses the
    \\q\\ observations nearest the cutoff on each side to test whether
    the type distribution is the same on both sides, via approximate
    sign randomisation (Canay and Kamat 2018). Run jointly over all
    types and periods: the test statistic is the sum of absolute mean
    differences across (period, type) pairs.

3.  **McCrary within-type** (sufficient, NOT necessary): McCrary (2008)
    density test run separately for each type value within each period,
    then combined across types via Bonferroni (minimum p-value times
    number of tests). If every within-type density is continuous at
    \\c\\ then the type shares are continuous (sufficiency), but shares
    can remain continuous even when all within-type densities jump by a
    common proportional factor (not necessary).

4.  **McCrary pooled** (NEITHER sufficient NOR necessary):
    McCrary (2008) density test on the pooled sample, run per period.
    Neither sufficient (type shares can jump while pooled density is
    smooth) nor necessary (pooled density can jump while shares stay
    continuous).
