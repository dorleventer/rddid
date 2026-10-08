# Local-linear RD in one period (building block)

The building block that
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) runs
in every period: the jump in the outcome at the cutoff in one period,
estimated by local-linear RD, both conventional and bias-corrected, with
standard errors. Unlike
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) it
takes vectors (`y`, `x`, `id`), not a data frame and column names, and
the bandwidths must be given, for example from
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md).
At a given pair (`h`, `b`) it reproduces the conventional and
bias-corrected estimates and the conventional and robust standard errors
of
[`rdrobust::rdrobust()`](https://rdrr.io/pkg/rdrobust/man/rdrobust.html)
with `vce = "hc1"`.

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

  the outcome (a numeric vector).

- x:

  the running variable (a numeric vector, same length as `y`).

- h:

  the main bandwidth (point estimate).

- b:

  the pilot bandwidth (bias correction); defaults to `h`.

- id:

  optional unit identifiers (a vector, same length as `y`), needed only
  to combine this period with others in a panel. With `NULL` (default)
  the observations are numbered `1, 2, ...`, so they cannot be matched
  across periods.

- c:

  the cutoff (default 0). A unit with `x >= c` is above the cutoff.

- p:

  order of the local polynomial for the point estimate (default 1, local
  linear).

- q:

  order of the local polynomial for the bias correction (default 2);
  must exceed `p`.

- kernel:

  the kernel: `"triangular"` (default), `"epanechnikov"` or `"uniform"`.

## Value

An object of class `"rd_period"`, a list with:

- `D`, `V_D`:

  the conventional jump and its variance.

- `D_bc`, `V_D_bc`:

  the bias-corrected jump and its robust variance.

- `b_const`, `v_const`:

  plug-in constants used by the bandwidth rules (an estimate of the bias
  constant, from the gap between the conventional and bias-corrected
  jumps, and `n * h * V_D`).

- `n`:

  the number of observations supplied.

- `h`, `b`, `c`, `p`, `q`, `kernel`:

  as passed.

- `sides`:

  a list with one element per side of the cutoff, `"+"` (above) and
  `"-"` (below), each holding the `id` of the observations used, their
  conventional and bias-corrected `g` vectors (`g`, `g_bc`) and `g_diff`
  (for the variance of the estimated bias), the conventional intercept
  at the cutoff `beta0`, its bias-corrected version `beta0_bc`, and the
  conventional `slope`. The fitted line on a side is
  `beta0 + slope * (x - c)`.

## Details

On each side of the cutoff the function keeps, for every observation
within the main or pilot bandwidth, its influence on the intercept times
its residual (`g`). These vectors are what the rest of the package
reuses: the variance of the jump is the sum of the squared `g` on both
sides, and in a panel the covariance between two periods' jumps sums the
products of `g` over the units present in both periods, matched on `id`.
Variances use the HC1 convention of rdrobust: residuals are scaled by
\\\sqrt{n_s / (n_s - k)}\\, with \\n_s\\ the observations used on that
side and \\k\\ the number of fitted coefficients; the bias-corrected
variance uses the residuals of the order-`q` pilot fit at `b`.

## References

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295-2326.

## See also

Other RD-DID estimation:
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
# the RD period (year 3) of rddid_sim, where the jump is 0.5 + 1 = 1.5
d3 <- rddid_sim[rddid_sim$year == 3, ]
bw <- rd_bw_cct(y = d3$Y, x = d3$R)
fit <- rd_period(y = d3$Y, x = d3$R, h = bw[["h"]], b = bw[["b"]], id = d3$id)
fit
#> Single-period RD (p=1, h=0.4141, b=0.6122, kernel=triangular, n=1000)
#>   D (conventional)   = +1.6108  (se 0.1337)
#>   D (bias-corrected) = +1.6397  (se 0.1599)
c(jump = fit$D_bc, se = sqrt(fit$V_D_bc))
#>      jump        se 
#> 1.6397286 0.1599211 
```
