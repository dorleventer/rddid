# RD-DID estimation and inference

Estimates a treatment effect in a regression-discontinuity
difference-in-discontinuities design. The period-\\t\_{\mathrm{RD}}\\
discontinuity is contaminated by a confounding policy that switches at
the same cutoff; comparison periods, where the confounding is present
but the treatment of interest is uniform at the cutoff, identify and net
out that confounding. The estimator is \\\widehat{\att} = \widehat
D\_{t\_{\mathrm{RD}}} - \sum_t w_t \widehat D_t\\, with each \\\widehat
D_t\\ a standard local-linear RD.

## Usage

``` r
rddid(
  data,
  y,
  x,
  time,
  id = NULL,
  t_rd,
  comparisons = NULL,
  weights = "constant",
  bwselect = c("joint", "cct", "iter"),
  h = NULL,
  b = NULL,
  start = "hstar",
  scheme = c("auto", "cs", "pc", "pv"),
  regularize = TRUE,
  reg_const = 3,
  c = 0,
  p = 1L,
  q = 2L,
  kernel = "triangular",
  level = 0.95
)
```

## Arguments

- data:

  a long data frame, one row per unit-period.

- y, x, time:

  column names (strings) for the outcome, running variable, and period.

- id:

  column name for the unit id; `NULL` (default) treats every row as a
  distinct unit (repeated cross-section).

- t_rd:

  the value of `time` identifying the RD (treated-at-cutoff) period.

- comparisons:

  values of `time` to use as comparison periods; `NULL` (default) uses
  every other period present.

- weights:

  `"constant"` (equal weights; constant confounding trend), `"linear"`
  (line through the comparison discontinuities extrapolated to `t_rd`),
  or a numeric vector over `comparisons`.

- bwselect:

  `"joint"` (default; a single common AMSE-optimal bandwidth for the
  aggregate estimator), `"cct"` (per-period MSE-optimal bandwidths via
  `rdrobust`), or `"iter"` (period-specific bandwidths chosen jointly by
  coordinate descent on the aggregate AMSE). Ignored if `h` is supplied.

- h, b:

  optional common point / pilot bandwidths; if `h` is given it is used
  for every period (with `b` defaulting to `h`).

- start:

  seed for the iterative (`bwselect = "iter"`) coordinate descent.
  `"hstar"` (default) starts all periods at the common joint-optimal
  h\*; `"cct"` starts each period at its own CCT/IK pilot h; or supply a
  named numeric vector/list with one entry per period. Ignored unless
  `bwselect = "iter"`.

- scheme:

  `"auto"` (detect from the id/side structure) or one of `"cs"`, `"pc"`,
  `"pv"`; selects which sampling-scheme variance is reported as the
  headline standard error. All three are always returned.

- regularize:

  logical; if `TRUE` (default) the joint AMSE-optimal bandwidth adds an
  `rdrobust`-style regularization term to the squared bias so a
  near-zero estimated curvature cannot blow the bandwidth up. Ignored if
  `h` is supplied.

- reg_const:

  regularization constant for `regularize` (default 3, matching the CCT
  convention).

- c:

  cutoff (default 0).

- p, q:

  point / bias-correction polynomial orders (default 1, 2).

- kernel:

  `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.

- level:

  confidence level (default 0.95).

## Value

An object of class `"rddid"` with the conventional and robust
bias-corrected estimates, standard errors under all three sampling
schemes, confidence intervals at the recommended scheme, the per-period
fits, the weights, and the bandwidth(s) used.
