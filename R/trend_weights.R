# trend_weights.R -- comparison-period weights implied by the confounding-trend assumption
# (constant: equal weights; linear: the least-squares line through the comparison jumps extrapolated to the RD
# period; numeric: as given). Used by rddid().

#' Period coefficients from a trend / weighting scheme
#' @keywords internal
#' @noRd
.rddid_weights <- function(weights, comps, t_rd) {
  m <- length(comps)
  if (is.numeric(weights)) {
    if (length(weights) != m)
      stop("numeric `trend` weights must have one entry per comparison period.")
    if (!is.null(names(weights)) && all(nzchar(names(weights)))) {
      if (!setequal(names(weights), as.character(comps)))
        stop("the names of the numeric `trend` weights must be the comparison periods: ",
             paste(comps, collapse = ", "))
      weights <- weights[as.character(comps)]     # match by name, not by position
    }
    w <- unname(weights)
    if (abs(sum(w) - 1) > 1e-8)
      warning("comparison weights do not sum to 1 (constant-trend admissibility).")
  } else {
    if (is.character(weights) && length(weights) == 1L && weights %in% c("ols", "min_variance"))
      stop("\"", weights, "\" is a value of `weighting`, not of `trend` (\"constant\" or ",
           "\"linear\"); use weighting = \"", weights, "\".", call. = FALSE)
    weights <- match.arg(weights, c("constant", "linear"))
    if (weights == "constant") {
      w <- rep(1 / m, m)                 # equal weights: constant confounding trend
    } else {
      if (m < 2) stop("linear weights need at least 2 comparison periods.")
      if (!is.numeric(comps) || !is.numeric(t_rd))
        stop("trend = \"linear\" needs numeric period values (the line is fitted on the period ",
             "values); the `time` column is not numeric.")
      X <- cbind(1, comps)               # extrapolate the line through the D_t to t_rd
      w <- as.numeric(c(1, t_rd) %*% solve(crossprod(X), t(X)))
    }
  }
  stats::setNames(w, as.character(comps))
}

#' Minimum-variance comparison weights from the per-period fits
#'
#' Among the weights the confounding-trend assumption admits (they sum to one; under
#' `"linear"` they also extrapolate the comparison periods to the RD period), the ones that
#' minimize the variance of the RD-DID estimate. The variances and covariances are those of
#' the conventional per-period jumps under `scheme`, the same objects the standard errors
#' use. Closed form: w = Xi^{-1} [xi + L'(L Xi^{-1} L')^{-1}(l - L Xi^{-1} xi)], with Xi the
#' covariance matrix of the comparison jumps, xi their covariances with the RD-period jump and
#' L w = l the constraints (dates centred at the RD period, which leaves the weights unchanged
#' and keeps L Xi^{-1} L' well conditioned with calendar years). When the constraints pin the
#' weights (as many comparison periods as constraint rows) the `w_pilot` weights are returned.
#' @return list: `w` (named weights), `Xi`, `xi`, `pinned`.
#' @keywords internal
#' @noRd
.mv_weights <- function(fits, comps, t_rd, trend, scheme, w_pilot) {
  keys <- as.character(comps)
  rd <- as.character(t_rd)
  m <- length(keys)
  Xi <- matrix(0, m, m, dimnames = list(keys, keys))
  xi <- stats::setNames(numeric(m), keys)
  for (i in seq_len(m)) {
    Xi[i, i] <- fits[[keys[i]]]$V_D
    xi[i] <- .cov_scheme(fits[[keys[i]]], fits[[rd]], scheme)
    if (i < m) for (j in (i + 1L):m)
      Xi[i, j] <- Xi[j, i] <- .cov_scheme(fits[[keys[i]]], fits[[keys[j]]], scheme)
  }
  k <- if (trend == "constant") 1L else 2L          # number of constraint rows
  if (m == k)                                        # a single admissible weight vector
    return(list(w = w_pilot, Xi = Xi, xi = xi, pinned = TRUE))
  if (inherits(tryCatch(chol(Xi), error = function(e) e), "error"))
    stop("the estimated covariance matrix of the comparison-period jumps is not positive ",
         "definite, so the minimum-variance weights are not defined (are two comparison periods ",
         "identical?); drop one of them or use `weighting = \"ols\"`.", call. = FALSE)
  L <- if (k == 1L) matrix(1, 1L, m) else rbind(1, as.numeric(comps) - as.numeric(t_rd))
  l <- if (k == 1L) 1 else c(1, 0)
  list(w = stats::setNames(.mv_solve(Xi, xi, L, l), keys), Xi = Xi, xi = xi, pinned = FALSE)
}

#' The minimum-variance weights given the covariance blocks and the constraints L w = l
#' @keywords internal
#' @noRd
.mv_solve <- function(Xi, xi, L, l) {
  Xi_xi <- solve(Xi, xi)
  Xi_Lt <- solve(Xi, t(L))
  zeta <- solve(L %*% Xi_Lt, l - L %*% Xi_xi)
  as.numeric(Xi_xi + Xi_Lt %*% zeta)
}
