# CCT (MSE-optimal) bandwidth for a single local-linear RD

Returns the Calonico–Cattaneo–Titiunik MSE-optimal bandwidths `h` and
`b` for a single local-linear RD using rdrobust. Falls back gracefully
when rdrobust is unavailable, the call fails, or the returned bandwidth
is non-positive/non-finite.

## Usage

``` r
rd_bw_cct(y, x, c = 0, p = 1L, kernel = "triangular")
```

## Arguments

- y:

  Outcome vector.

- x:

  Running variable vector.

- c:

  Cutoff (default 0).

- p:

  Polynomial order (default 1L, local linear).

- kernel:

  Kernel type: `"triangular"` (default), `"epanechnikov"`, or
  `"uniform"`.

## Value

A named numeric vector `c(h = ..., b = ...)` with the main and pilot
bandwidths. When the CCT computation is unavailable, both equal
`0.5 * IQR(x)` (or `sd(x)` if IQR is zero), and a message is emitted
naming the reason.
