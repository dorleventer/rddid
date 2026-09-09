# Correction-form composition-adjusted RD-DID ATT family (Theorem 3)

Assembles the family of composition-adjusted
average-treatment-on-the-treated (ATT) estimators in their CORRECTION
form from the three per-block primitives of the package – the period
discontinuity \\D_t\\
([`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)),
the within-period composition term \\S_t\\ (Proposition `sdisc`,
[`rd_sadjust()`](https://dorleventer.github.io/rddid/reference/rd_sadjust.md)),
and the cross-period composition term \\C\\ (Proposition `cdisc`,
[`rd_c()`](https://dorleventer.github.io/rddid/reference/rd_c.md)) – and
exposes the full \\D\\/\\S\\/\\C\\ decomposition. The headline
quantities are \$\$\mathrm{ATT}\_{\mathrm{unadj}} = D\_{t\_{RD}} -
g_0(\\D\_{t_0}\\),\$\$ \$\$\mathrm{ATT}\_{s} =
(D\_{t\_{RD}}-S\_{t\_{RD}}) - g_0(\\D\_{t_0}-S\_{t_0}\\),\$\$
\$\$\mathrm{ATT}\_{c} = D\_{t\_{RD}} - g_0(\\D\_{t_0}\\) -
C\_{\mathrm{agg}},\$\$ \$\$\mathrm{ATT}\_{sc} =
(D\_{t\_{RD}}-S\_{t\_{RD}}) - g_0(\\D\_{t_0}-S\_{t_0}\\) -
C\_{\mathrm{agg}},\$\$ where \\g_0(\\x\_{t_0}\\) = \sum\_{t_0} w\_{t_0}
x\_{t_0}\\ is the trend extrapolation with the package trend weights and
\\C\_{\mathrm{agg}} = \sum\_{t_0} w\_{t_0} C\_{t_0,t\_{RD}}\\.

## Usage

``` r
rd_att(
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
  trend = "constant",
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

  optional fixed numeric bandwidth used for every block; if `NULL`
  (default) a per-block CCT bandwidth is selected once on the original
  sample.

- bwselect:

  bandwidth selector label when `h` is `NULL` (default `"cct"`);
  retained for forward compatibility.

- kernel:

  `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.

- trend:

  trend extrapolation weights: `"constant"` (default), `"linear"`, or a
  numeric vector with one entry per comparison period (as in
  [`rd_adjust()`](https://dorleventer.github.io/rddid/reference/rd_adjust.md)).

- se:

  `"bootstrap"` (default) or `"none"`.

- B:

  bootstrap replications (default 499).

## Value

An object of class `"rd_att"`: a list with

- `att` – the ATT table (one row per `estimator` in `{unadj, s, c, sc}`,
  columns `conv`, `bc`, `se`, `bc_se`);

- `per_period` – the per-period decomposition table (`period`, `role`,
  `D`, `D_bc`, `D_se`, `D_bc_se`, `S`, `S_bc`, `S_se`, `S_bc_se`,
  `adj_s_D`, `adj_s_D_bc`, `adj_s_D_se`, `adj_s_D_bc_se`);

- `cterm` – the cross-period composition table (`comparison`, `C`,
  `C_bc`, `C_se`, `C_bc_se`);

- `cagg` – the aggregate `Cagg` (conv + bc + SE);

- meta fields (`t_rd`, `comparisons`, `weights`, `trend`, `c`, `kernel`,
  `bwselect`, `se`, `B`, `n_boot_ok`, `call`).

## Details

The estimator uses the pairwise (\\P = 2\\) binary-type design: each
comparison period `t0` forms the pair `{t0, t_rd}`. \\S\_{t_0}\\ uses
partner `t_rd`; \\S\_{t\_{RD}}\\ is the trend-weighted average over
partners, so the single-comparison case (\\w = 1\\) reduces exactly to
\\\mathrm{adj}\text{-}c = D\_{t\_{RD}} - D\_{t_0} - C\_{t_0,t\_{RD}}\\
and \\\mathrm{adj}\text{-}sc = (D\_{t\_{RD}}-S\_{t\_{RD}}) -
(D\_{t_0}-S\_{t_0}) - C\_{t_0,t\_{RD}}\\.

All standard errors come from ONE shared unit-level cluster bootstrap:
the per-block CCT bandwidths are chosen once on the original sample and
held fixed; each replication resamples units with replacement, rebuilds
the panel with fresh sequential cluster ids, and recomputes every
primitive at the fixed bandwidths before re-forming every ATT. Because
all quantities are recomputed on the same draws, the
SD-over-replications standard errors carry the cross-primitive
covariances, so the combination SEs are valid.

This is the correction form and does not reimplement or depend on the
within-type (reweighting) `att_adj` of
[`rd_adjust()`](https://dorleventer.github.io/rddid/reference/rd_adjust.md),
which is retained as-is.
