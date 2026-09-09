# Within-period composition-adjusted RD estimator (Proposition `sdisc`)

Estimates the within-period composition term \\S_t\\ that contaminates
the observed period-\\t\\ discontinuity \\D_t\\ when units sort across
the cutoff between regimes, and reports the composition-free
discontinuity \\D_t - S_t\\. With binary types (the unit's side of the
cutoff in the partner period) the term collapses to \\S_t =
\Delta\pi\_{above}\\(\mu\_{above}^{(-)} - \mu\_{below}^{(-)})\\, where
\\\Delta\pi\_{above}\\ is the cutoff jump in the above-type share and
\\\mu_a^{(-)}\\ is the below-cutoff (control-side) outcome limit among
type-\\a\\ units. Because \\\mu_a^{(-)}\\ uses the control side it is
observed in every period, so \\S_t\\ is estimable even at the RD period.

## Usage

``` r
rd_sadjust(
  data,
  y,
  x,
  time,
  id,
  comparisons,
  t_rd,
  c = 0,
  h = NULL,
  bwselect = "cct",
  kernel = "triangular",
  se = c("bootstrap", "none"),
  B = 499L
)
```

## Arguments

- data:

  long data frame, one row per unit-period.

- y, x, time, id:

  column names (strings) for outcome, running variable, period, and unit
  id.

- comparisons:

  comparison-period values of `time`; the pair `{t0, t_rd}` is formed
  for each.

- t_rd:

  RD-period value of `time`.

- c:

  cutoff (default 0).

- h:

  optional fixed numeric bandwidth used for every fit; if `NULL`
  (default) a per-(period, type, fit) CCT bandwidth is selected.

- bwselect:

  bandwidth selector label when `h` is `NULL` (default `"cct"`);
  retained for forward compatibility.

- kernel:

  `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.

- se:

  `"bootstrap"` (default) or `"none"`.

- B:

  bootstrap replications (default 499).

## Value

An object of class `"rd_sadjust"`: a list with a tidy per-period `table`
(columns `comparison`, `period`, `role`, `S`, `S_bc`, `S_se`, `S_bc_se`,
`D`, `D_bc`, `adjusted`, `adjusted_bc`, `n`, `n_control_switchers`) and
meta fields.

## Details

Uses the pairwise (P = 2) design: for each comparison period in
`comparisons`, the pair `{t0, t_rd}` yields \\S\_{t0}\\ (type = the
unit's `t_rd` side) and \\S\_{t\\rd}\\ (type = the unit's `t0` side).
Standard errors are a unit-level cluster bootstrap at bandwidths held
fixed at the point-estimate values.
