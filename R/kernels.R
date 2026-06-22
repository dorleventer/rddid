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
#' @keywords internal
#' @noRd
.qrXXinv <- function(x) {
  chol2inv(chol(crossprod(x)))
}
