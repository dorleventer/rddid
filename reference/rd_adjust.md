# Composition-adjusted RD-DID estimator (Theorem 3)

Estimates the composition-adjusted ATT for a time-varying running
variable, the route of Section 4.4 / Theorem `thm:adjust` to take when
composition stability (`ass:comp-stable`,
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md))
is rejected but the within-type confounding trend (`ass:trend-cell`,
[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md))
is credible. Reports the per-type RD-DID effects
`ATT(t_rd | v_comp = a)`, the RD-period composition shares, and their
share-weighted aggregate, alongside the unadjusted estimator. Standard
errors are a unit-level cluster bootstrap.

## Usage

``` r
rd_adjust(
  data,
  y,
  x,
  time,
  id,
  t_rd,
  comparisons = NULL,
  c = 0,
  h,
  b = h,
  weights = "constant",
  kernel = "triangular",
  min_n = 10L,
  p = 1L,
  q = 2L,
  se = c("bootstrap", "none"),
  B = 500L
)
```

## Arguments

- data:

  long data frame, one row per unit-period.

- y, x, time, id:

  column names (strings).

- t_rd:

  RD-period value of `time`.

- comparisons:

  comparison-period values of `time`; `NULL` uses all others.

- c:

  cutoff (default 0).

- h:

  common bandwidth for every block (jumps and shares); cells are thin,
  so a bandwidth wider than the aggregate is appropriate.

- b:

  pilot bandwidth for bias correction (defaults to `h`).

- weights:

  `"constant"`, `"linear"`, or a numeric vector over `comparisons` (the
  trend `g0`, as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)).

- kernel:

  `"triangular"` (default), `"epanechnikov"`, `"uniform"`.

- min_n:

  minimum observations per side of any block (default 10).

- p, q:

  point / bias-correction polynomial orders (default 1, 2).

- se:

  `"bootstrap"` (default) or `"none"`.

- B:

  bootstrap replications (default 500).

## Value

object of class `"rd_adjust"`.
