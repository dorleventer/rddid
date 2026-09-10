# Bandwidth rules and sampling schemes

## Setup

This vignette covers the
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
arguments that control the bandwidths ($`h_t`$, the period-$`t`$
bandwidth, and $`b_t`$, the pilot bandwidth for the bias correction) and
the standard errors reported under the three sampling schemes. For the
estimator, the weights, and the per-period output, see
[`vignette("rddid-estimation")`](https://dorleventer.github.io/rddid/articles/rddid-estimation.md).

Both parts use the same one-sided-curvature DGP: periods 1 and 2 are
comparison periods, period 3 is the RD period, $`\alpha_t = 1`$ in every
period and $`\tau = 1`$, so $`D_1 = D_2 = 1`$, $`D_3 = 2`$,
$`\mathrm{ATT} = 1`$. Part 1 varies the curvature $`\theta_t`$ across
periods so that the per-period bandwidths differ; Part 2 varies the
sampling scheme.

``` r

library(rddid)
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
# per-period bandwidths h_t of a fitted rddid object, ordered by period
# (r$fits lists the RD period first)
bw_by_period <- function(r) {
  h <- vapply(r$fits, function(f) unname(f$h), numeric(1))
  round(h[order(as.numeric(names(h)))], 3)
}
# bandwidth method, h_t by period, and the estimates
peek <- function(r) {
  cat("method:", r$bandwidth$method, "\n"); print(bw_by_period(r))
  print(round(r$estimates[, c("est", "se")], 3))
}

m_a <- function(r, theta) r + (theta / 2) * r^2 * (r >= 0)
dgp_a <- function(n = 2000, alpha = c(1, 1, 1), theta = c(2, 2, 2), tau = 1,
                  scheme = c("pc", "cs", "pv"), seed = 1) {
  scheme <- match.arg(scheme)
  set.seed(seed)
  R0 <- runif(n, -1, 1)                 # running variable at baseline, cutoff at 0
  u0 <- rnorm(n, 0, 0.5)                # unit effect
  do.call(rbind, lapply(1:3, function(t) {
    if (scheme == "cs") {               # new units every period
      R <- runif(n, -1, 1); u <- rnorm(n, 0, 0.5); id <- (t - 1) * n + seq_len(n)
    } else {                            # same units every period
      R <- if (scheme == "pv") R0 + rnorm(n, 0, 0.1) else R0
      u <- u0; id <- seq_len(n)
    }
    V <- as.integer(R >= 0)             # confounding treatment: sharp RD every period
    W <- V * (t == 3)                    # treatment of interest: RD period only
    data.frame(id = id, t = t, R = R,
               Y = m_a(R, theta[t]) + alpha[t] * V + tau * W + u + rnorm(n, 0, 0.5))
  }))
}
```

## Bandwidth rules

Data for this section have a curvature $`\theta_t`$ that differs across
periods, under the `pc` sampling scheme:

``` r

dat1 <- dgp_a(theta = c(1, 6, 3))
```

### `bwselect = "cct"`

Per-period CCT bandwidths: each $`h_t`$ (and pilot $`b_t`$) minimizes
the MSE of its own $`\widehat D_t`$, ignoring the other periods.

``` r

r_cct <- rddid(dat1, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, bwselect = "cct")
peek(r_cct)
#> method: cct 
#>     1     2     3 
#> 0.390 0.328 0.398 
#>                est    se
#> Conventional 1.073 0.089
#> Robust       1.107 0.106
```

The per-period $`h_t`$ sit in `r_cct$fits[[k]]$h`: $`h_1 =
0.390`$, $`h_2 = 0.328`$, $`h_3 = 0.398`$; the narrowest is period 2,
the most curved.

### `bwselect = "joint"`

A single common $`h^\ast`$ minimizes the asymptotic MSE of the aggregate
$`\widehat{\mathrm{ATT}}(t_{\mathrm{RD}})`$, using one bandwidth for
every period.

``` r

r_joint <- rddid(dat1, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, bwselect = "joint")
peek(r_joint)
#> method: joint 
#>     1     2     3 
#> 0.305 0.305 0.305 
#>               est   se
#> Conventional 1.09 0.10
#> Robust       1.12 0.12
```

`r_joint$bandwidth$h` is the common $`h^\ast = 0.305`$, applied to every
period’s `fits[[k]]$h`.

### `bwselect = "iter"` (the default)

Period-specific bandwidths targeting the same aggregate MSE as
`"joint"`, but with one $`h_t`$ per period, found by coordinate descent
started at $`h^\ast`$.

``` r

r_iter <- rddid(dat1, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, bwselect = "iter")
peek(r_iter)
#> method: iter 
#>     1     2     3 
#> 0.335 0.263 0.317 
#>                est    se
#> Conventional 1.094 0.098
#> Robust       1.124 0.117
```

`start` sets the starting point of the descent: `start = "hstar"`
(default) starts every period at the common $`h^\ast`$; `start = "cct"`
starts each period at its own CCT $`h_t`$; or a named vector supplies a
per-period starting point.

``` r

r_iter_cct <- rddid(dat1, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                    bwselect = "iter", start = "cct")
bw_by_period(r_iter_cct)
#>     1     2     3 
#> 0.334 0.263 0.317
```

`regularize` (default `TRUE`) adds an `rdrobust`-style regularization
term to the squared-bias part of the `"joint"` and `"iter"` objectives,
so a near-zero estimated curvature does not blow the bandwidth up;
`reg_const` (default 3) is its constant. $`h`$ and $`b`$ can also be
fixed by hand, bypassing bandwidth selection entirely:

``` r

r_fixed <- rddid(dat1, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, h = 0.3)
r_fixed$bandwidth$method
#> [1] "fixed"
```

``` r

estse <- function(r, row) paste0(fmt(r$estimates[row, "est"]), " (", fmt(r$estimates[row, "se"]), ")")
tab <- data.frame(
  rule = c("cct", "joint", "iter"),
  h_1 = c(r_cct$fits[["1"]]$h, r_joint$bandwidth$h, r_iter$fits[["1"]]$h),
  h_2 = c(r_cct$fits[["2"]]$h, r_joint$bandwidth$h, r_iter$fits[["2"]]$h),
  h_3 = c(r_cct$fits[["3"]]$h, r_joint$bandwidth$h, r_iter$fits[["3"]]$h),
  Conventional = sapply(list(r_cct, r_joint, r_iter), estse, row = "Conventional"),
  Robust = sapply(list(r_cct, r_joint, r_iter), estse, row = "Robust"),
  truth = 1
)
knitr::kable(tab, digits = 3)
```

| rule  |   h_1 |   h_2 |   h_3 | Conventional  | Robust        | truth |
|:------|------:|------:|------:|:--------------|:--------------|------:|
| cct   | 0.390 | 0.328 | 0.398 | 1.073 (0.089) | 1.107 (0.106) |     1 |
| joint | 0.305 | 0.305 | 0.305 | 1.090 (0.100) | 1.120 (0.120) |     1 |
| iter  | 0.335 | 0.263 | 0.317 | 1.094 (0.098) | 1.124 (0.117) |     1 |

The three rules give different per-period bandwidths — from
$`h_3 = 0.398`$ under `"cct"` down to $`h_3 = 0.317`$ under `"iter"` —
and, on these data, Conventional estimates within one standard error of
each other (1.073 to 1.094, each with SE around 0.089).

## Sampling schemes

The three schemes, `cs` (repeated cross-section), `pc` (panel,
time-constant running variable) and `pv` (panel, time-varying running
variable), are defined and their detection rule stated in the Get
started vignette. They differ in the standard error: under `pc` the same
units enter every period’s local regression, so the per-period
discontinuity estimates are correlated through the unit effect; under
`pv` units move across the window, so the shared-unit covariance is of
smaller order than the variance; under `cs` there is no cross-period
covariance.

For the `cs` data, `id = NULL` treats every row as a distinct unit:

``` r

dat_cs <- dgp_a(scheme = "cs")
dat_pc <- dgp_a(scheme = "pc")
dat_pv <- dgp_a(scheme = "pv")
r_cs <- rddid(dat_cs, y = "Y", x = "R", time = "t", id = NULL, t_rd = 3, bwselect = "cct")
r_pc <- rddid(dat_pc, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, bwselect = "cct")
r_pv <- rddid(dat_pv, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, bwselect = "cct")
c(cs = r_cs$scheme, pc = r_pc$scheme, pv = r_pv$scheme)
#>   cs   pc   pv 
#> "cs" "pc" "pv"
```

``` r

se_tab <- data.frame(
  data = c("cs", "pc", "pv"),
  rbind(r_cs$estimates["Conventional", c("se_cs", "se_pc", "se_pv")],
        r_pc$estimates["Conventional", c("se_cs", "se_pc", "se_pv")],
        r_pv$estimates["Conventional", c("se_cs", "se_pc", "se_pv")]))
knitr::kable(se_tab, digits = 3, row.names = FALSE)
```

| data | se_cs | se_pc | se_pv |
|:-----|------:|------:|------:|
| cs   | 0.147 | 0.147 | 0.147 |
| pc   | 0.143 | 0.089 | 0.089 |
| pv   | 0.139 | 0.129 | 0.140 |

On the `pc` data, `se_pc` (0.089) is smaller than `se_cs` (0.143); on
the `pv` data, `se_pv` (0.140) is close to `se_cs` (0.139) while `se_pc`
(0.129) is smaller; on the `cs` data all three coincide, at 0.147. The
scheme is a description of how the data were sampled, not a choice among
the three columns.

Passing `scheme` explicitly picks which of the three the print method’s
SE and CI correspond to, rather than detecting it from the data:

``` r

r_pv_explicit <- rddid(dat_pv, y = "Y", x = "R", time = "t", id = "id", t_rd = 3,
                       bwselect = "cct", scheme = "pv")
r_pv_explicit$scheme
#> [1] "pv"
r_pv_explicit$estimates["Conventional", c("se", "ci_l", "ci_u")]
#>                     se      ci_l     ci_u
#> Conventional 0.1398743 0.4860585 1.034356
```
