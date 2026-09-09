# Cross-period composition term C (Proposition `cdisc`)

Estimates the cross-period composition term \\C\_{t,s}\\ between two
time periods as a first-class primitive: the Kitagawa (between-group)
part of the change in the confounding discontinuity from period `t` to
period `s`. With binary sides \\a\in\\0,1\\\\ (the unit's side of the
cutoff) it is
\$\$C\_{t,s}=\sum_a\bigl(\pi\_{s,(+)}(v_t{=}a)-\pi\_{t,(+)}(v_s{=}a)\bigr)\\
D_t(v_s{=}a),\$\$ the linear-\\g_0\\ estimator of Appendix
`app:est-adjust`, where \\D_t(v_s{=}a)\\ is the period-`t` discontinuity
among units with side `a` in the partner period `s`, and
\\\pi\_{\cdot,(+)}\\ are above-cutoff cell shares (above-cutoff
intercepts of a local-linear RD of the side indicator).

## Usage

``` r
rd_c(
  data,
  y,
  x,
  time,
  id,
  t,
  s,
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

- t:

  first (comparison-role) period value of `time`.

- s:

  second (RD-role) period value of `time`.

- c:

  cutoff (default 0).

- h:

  optional fixed numeric bandwidth used for every block; if `NULL`
  (default) a per-block CCT bandwidth is selected.

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

An object of class `"rd_c"`: a list with the conventional and
bias-corrected composition term (`C`, `C_bc`), their cluster-bootstrap
standard errors (`C_se`, `C_bc_se`), the ordered periods (`t`, `s`), a
per-cell `cells` data frame (`label`, `dpi`, `dpi_bc`, `D`, `D_bc`,
`pi_s_plus`, `pi_t_plus`, `contrib`, `contrib_bc`), and meta fields.

## Details

In the composition-adjusted ATT (Theorem `thm:adjust`) the ordered pair
is `(t, s) = (t0, t_rd)`: `t` is the comparison period, `s` the RD
period, and the correction-form ATT subtracts \\C\_{t_0,t\_{RD}}\\ from
the unadjusted aggregate. \\C\_{t,s}\\ is the pure CROSS-period term and
is distinct from the within-period term \\S_t\\ of
[`rd_sadjust()`](https://dorleventer.github.io/rddid/reference/rd_sadjust.md);
the two are combined only at the ATT level. \\C\_{t,s}\\ is not
symmetric in `(t, s)`: the comparison jump \\D\\ is taken at `t`, the
share weight \\\pi\_{s,(+)}\\ at `s`.

Bandwidths default to a per-block CCT selection, made once on the
original sample and held fixed across the unit-level cluster bootstrap
that supplies the standard error (recomputing CCT on the thin cells is
unstable).
