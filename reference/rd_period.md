# Single-period local-linear RD discontinuity

Estimates the period-\\t\\ outcome discontinuity \\D_t =
\beta^{(0)}\_{(+)} - \beta^{(0)}\_{(-)}\\ by a standard local-linear
regression discontinuity, conventional and robust-bias-corrected
(Calonico, Cattaneo and Titiunik 2014). At a given bandwidth pair (`h`,
`b`) it reproduces the Conventional and Bias-Corrected estimates and the
Conventional and Robust standard errors of rdrobust to machine
precision. This is the per-period engine that
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
aggregates across periods.

## Usage

``` r
rd_period(
  y,
  x,
  h,
  b = h,
  id = NULL,
  c = 0,
  p = 1L,
  q = 2L,
  kernel = "triangular"
)
```

## Arguments

- y:

  outcome vector.

- x:

  running variable.

- h:

  main (point-estimate) bandwidth.

- b:

  pilot (bias-correction) bandwidth; defaults to `h`.

- id:

  optional unit identifiers, needed only when this period will be
  combined with others under panel sampling. If `NULL`, sequential ids
  are assigned and no cross-period matching is possible.

- c:

  cutoff (default 0).

- p:

  point-estimate polynomial order (default 1, local linear).

- q:

  bias-correction polynomial order (default 2); must exceed `p`.

- kernel:

  `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.

## Value

An object of class `"rd_period"`: a list with the conventional and
bias-corrected discontinuity (`D`, `D_bc`) and their variances (`V_D`,
`V_D_bc`), the per-period plug-in bias and variance constants used for
joint bandwidth selection (`b_const = (p+1)! (D - D_bc) / h^{p+1}` and
`v_const = n h V_D`, the plug-ins for the constants `B_t` and `V_t` of
`eq:per-period-orders` in Appendix B.3), the effective sample size `n`,
the bandwidths, and a per-side list `sides` holding, for each side `"+"`
and `"-"`, the active units' `id`, their conventional / bias-corrected
`g` vectors, the conventional intercept at the cutoff (`beta0`), and the
conventional local-linear `slope`. The fitted line on a side is
`beta0 + slope * (x - c)`.

## Details

The function returns, for each side of the cutoff, the per-unit
influence weight on the intercept times its local-linear residual, `g`.
These `g` vectors are the single primitive the rest of the package
reuses:

- the discontinuity variance is \\V(\hat D_t) = \sum g\_{(+)}^2 + \sum
  g\_{(-)}^2\\;

- any cross-period covariance is a merge of two periods' `g` vectors on
  shared unit `id` (see the internal covariance helpers);

- the bandwidth constants are read off `V_D` and the conventional /
  bias-corrected gap.

Variances use the HC1 finite-sample convention of rdrobust
(`vce = "hc1"`): residuals are scaled by `sqrt(n_s / (n_s - k))` with
`n_s` the side's active sample and `k` the number of fitted
coefficients; the BC variance uses the residuals of the order-`q` pilot
fit at `b`, as in rdrobust. The code follows the matrix form of Appendix
B of Leventer and Nevo; the object-by-object map is `dev/appB_map.md` in
the source repository.
