#' Kernel weights for local-linear RD
#'
#' Triangular (default), Epanechnikov, or uniform kernel weights, matching the
#' conventions of \pkg{rdrobust}. Returned weights are *not* divided by the
#' bandwidth, since only relative weights enter the weighted least squares.
#'
#' @param u scaled distance `(x - c) / h`.
#' @param kernel one of `"triangular"`, `"epanechnikov"`, `"uniform"`.
#' @return numeric vector of weights, zero outside `[-1, 1]`.
#' @keywords internal
#' @noRd
.rd_kweight <- function(u, kernel = "triangular") {
  kernel <- match.arg(tolower(kernel),
                       c("triangular", "epanechnikov", "uniform"))
  inwin <- abs(u) <= 1
  switch(kernel,
    triangular   = (1 - abs(u)) * inwin,
    epanechnikov = 0.75 * (1 - u^2) * inwin,
    uniform      = 0.5 * inwin)
}

#' Inverse of a weighted Gram matrix via Cholesky, given the square-root design
#'
#' The columns of `x` are powers of the centred running variable, so
#' `crossprod(x)` has a condition number of order 1e4–1e6 at orders 2–3 for
#' typical bandwidths. The inverse is computed after scaling each column of `x`
#' to unit norm and undoing the scaling afterwards (exact algebra, `G = D G* D`),
#' which keeps the result stable to ~1e-13 across BLAS implementations instead
#' of ~1e-10.
#' @keywords internal
#' @noRd
.qrXXinv <- function(x) {
  s  <- sqrt(colSums(x^2))
  Gi <- chol2inv(chol(crossprod(x / rep(s, each = nrow(x)))))
  Gi / outer(s, s)
}
