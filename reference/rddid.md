# RD-DID estimation and inference

Estimates a treatment effect in a regression-discontinuity
difference-in-discontinuities design. The period-\\t\_{\mathrm{RD}}\\
discontinuity is contaminated by a confounding policy that switches at
the same cutoff; comparison periods, where the confounding is present
but the treatment of interest is uniform at the cutoff, identify and net
out that confounding. The estimator is \\\widehat{\mathrm{ATT}} =
\widehat D\_{t\_{\mathrm{RD}}} - \sum_t w_t \widehat D_t\\, with each
\\\widehat D_t\\ a standard local-linear RD. This estimates the ATT, or
the ATU when the comparison periods are uniformly treated
(`estimand = "atu"`).

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
  estimand = c("att", "atu"),
  bwselect = c("iter", "joint", "cct"),
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

  a long data frame, one row per unit-period (repeated cross-section or
  panel; a panel need not be balanced).

- y, x, time:

  column names (strings) for the outcome, running variable, and period.

- id:

  column name for the unit id; `NULL` (default) treats every row as a
  distinct unit (repeated cross-section).

- t_rd:

  the value of `time` identifying the RD (treated-at-cutoff) period.

- comparisons:

  values of `time` to use as comparison periods: periods in which the
  treatment of interest does not switch at the cutoff. `NULL` (default)
  uses every other period present, so with more than one RD period pass
  `comparisons` explicitly.

- weights:

  `"constant"` (equal weights; constant confounding trend), `"linear"`
  (line through the comparison discontinuities extrapolated to `t_rd`),
  or a numeric vector over `comparisons`.

- estimand:

  `"att"` (default) when the treatment of interest is uniformly ZERO in
  the comparison periods (targets the ATT), `"atu"` when it is uniformly
  ONE (targets the ATU). See "Targeting the ATU" below.

- bwselect:

  `"iter"` (default; period-specific bandwidths chosen jointly by
  coordinate descent on the aggregate AMSE, started at the common
  joint-optimal bandwidth), `"joint"` (a single common AMSE-optimal
  bandwidth for the aggregate estimator), or `"cct"` (per-period CCT
  MSE-optimal bandwidths,
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)).
  Ignored if `h` is supplied. Both joint rules estimate each period's
  bias and variance constants at that period's own CCT pilot, so neither
  depends on which period is labelled `t_rd`; under `"joint"` the pilot
  `b_t` keeps each period's CCT ratio `b_t^CCT / h_t^CCT` (Appendix B.4
  of the paper).

- h, b:

  optional common point / pilot bandwidths; if `h` is given it is used
  for every period (with `b` defaulting to `h`).

- start:

  starting point of the iterative (`bwselect = "iter"`) coordinate
  descent. `"hstar"` (default) starts all periods at the common
  joint-optimal h\*; `"cct"` starts each period at its own CCT h; or
  supply a named numeric vector/list with one entry per period. Ignored
  unless `bwselect = "iter"`.

- scheme:

  sampling scheme: `"cs"` (repeated cross-section), `"pc"` (panel,
  time-constant running variable), `"pv"` (panel, time-varying running
  variable), or `"auto"` (default), which reads it off the data: no `id`
  gives `"cs"`; units observed in more than one period, each on the same
  side of the cutoff in every period, give `"pc"`; any unit on different
  sides in different periods gives `"pv"`. The scheme selects which
  standard error and CI are printed; all three are always returned.

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

An object of class `"rddid"`, a list with:

- `estimates`:

  matrix with rows `Conventional` and `Robust` (the bias-corrected
  estimate with its robust variance) and columns `est`, `se`, `ci_l`,
  `ci_u` (at `scheme`), and `se_cs`, `se_pc`, `se_pv`.

- `scheme`:

  the sampling scheme used; `scheme_requested` is the argument as
  passed.

- `weights`, `weights_type`:

  the comparison-period weights and their kind.

- `estimand`:

  `"att"` or `"atu"`, as passed.

- `bandwidth`:

  list with `method` (the `bwselect` value, or `"fixed"`), `h`, `b`, and
  `niter` for `"iter"`.

- `fits`:

  named list of
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  objects by period, the RD period first (index by name).

- `t_rd`, `comparisons`, `level`, `call`:

  as passed.

## Details

### Targeting the ATU

When the treatment of interest is uniformly present (equal to one) in
the comparison periods rather than uniformly absent, the same difference
of discontinuities identifies the ATU: Leventer and Nevo, Section 6,
show that the ATU design is the ATT design with the sides of the cutoff
exchanged (mirror \\\tilde R = c - R\\ and apply the ATT procedure
unchanged). The point estimate, standard errors, and bandwidth rules are
numerically the same either way, so `estimand` only labels the output
here; the argument is also passed through to the validation tests
([`?rd_typecont`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
[`?rd_homog`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
[`?rd_trendcell`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md),
[`?rd_compstable`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)),
where only
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
computes differently.
