# trend_weights.R -- comparison-period weights implied by the confounding-trend assumption
# (constant: equal weights; linear: the least-squares line through the comparison jumps extrapolated to the RD
# period; numeric: as given), the RD-period weights of an aggregate over several RD periods,
# and the minimum-variance weights. Used by didisc().

#' Comparison-period weights from a trend / weighting scheme
#' @keywords internal
#' @noRd
.didisc_weights <- function(weights, comps, t_rd) {
  m <- length(comps)
  if (is.numeric(weights)) {
    if (length(weights) != m)
      .abort("numeric `trend` weights need one entry per comparison period (", m, "); got ",
             length(weights), ".")
    weights <- .match_names(weights, comps, "the numeric `trend` weights", "comparison periods")
    w <- unname(weights)
    if (abs(sum(w) - 1) > 1e-8)
      .warn("comparison weights do not sum to 1 (constant-trend admissibility).")
  } else {
    if (weights == "constant") {
      w <- rep(1 / m, m)                 # equal weights: constant confounding trend
    } else {
      if (m < 2) .abort("trend = \"linear\" needs at least 2 comparison periods.")
      # extrapolate the line through the D_t to t_rd, with time centred at t_rd: shifting the
      # time axis changes neither the fitted line nor its value at t_rd, but the normal
      # equations in raw calendar years (~2000) lose ~1e-9 to cancellation
      X <- cbind(1, comps - t_rd)
      w <- as.numeric(c(1, 0) %*% solve(crossprod(X), t(X)))
    }
  }
  stats::setNames(w, as.character(comps))
}

#' RD-period weights phi of the aggregate over several RD periods
#'
#' `NULL` gives equal weights 1/J. Any finite weights are allowed (they are the
#' researcher's choice of estimand); named weights are matched to `t_rd` by name.
#' @keywords internal
#' @noRd
.didisc_rd_weights <- function(rd_weights, t_rd) {
  J <- length(t_rd)
  if (is.null(rd_weights)) return(stats::setNames(rep(1 / J, J), as.character(t_rd)))
  if (J == 1L)
    .abort("`rd_weights` weights several RD periods; with the single RD period `t_rd` = ", t_rd,
           " leave it out.")
  if (!is.numeric(rd_weights) || length(rd_weights) != J || any(!is.finite(rd_weights)))
    .abort("`rd_weights` must be finite numbers, one per RD period in `t_rd` (", J, "); got ",
           length(rd_weights), " value(s).")
  if (all(rd_weights == 0))
    .abort("`rd_weights` are all zero, so the estimand is identically zero.")
  rd_weights <- .match_names(rd_weights, t_rd, "`rd_weights`", "RD periods")
  stats::setNames(unname(rd_weights), as.character(t_rd))
}

#' Comparison weights for each RD period (one row per target, one column per comparison)
#'
#' Row t holds the weights w_{t0}(t) of Corollary cor:trend for target t: `.didisc_weights()`
#' applied to each RD period. A numeric vector `trend` gives the same weights to every target;
#' a numeric matrix gives one row per RD period (rows in `t_rd` order or named by period,
#' columns in `comps` order or named by period).
#' @return numeric matrix, J x m, with dimnames (RD periods, comparison periods).
#' @keywords internal
#' @noRd
.didisc_weights_by_target <- function(trend, comps, t_rd) {
  J <- length(t_rd)
  if (is.matrix(trend)) {
    if (!is.numeric(trend) || nrow(trend) != J || ncol(trend) != length(comps))
      .abort("a matrix `trend` must be numeric with one row per RD period (", J, ") and one ",
             "column per comparison period (", length(comps), ").")
    if (!is.null(rownames(trend)))
      trend <- trend[.match_names(stats::setNames(seq_len(J), rownames(trend)), t_rd,
                                  "the rows of a matrix `trend`", "RD periods"), , drop = FALSE]
    rows <- lapply(seq_len(J), function(j)
      .didisc_weights(stats::setNames(trend[j, ], colnames(trend)), comps, t_rd[j]))
  } else if (is.numeric(trend)) {
    rows <- rep(list(.didisc_weights(trend, comps, t_rd[1])), J)   # same weights for every target
  } else {
    rows <- lapply(seq_len(J), function(j) .didisc_weights(trend, comps, t_rd[j]))
  }
  W <- do.call(rbind, rows)
  dimnames(W) <- list(as.character(t_rd), as.character(comps))
  W
}

#' Minimum-variance comparison weights from the per-period fits
#'
#' Among the weights the confounding-trend assumption admits (they sum to one; under
#' `"linear"` they also extrapolate the comparison periods to the RD period), the ones that
#' minimize the variance of the diff-in-disc estimate. The variances and covariances are those of
#' the conventional per-period jumps under `scheme`, the same objects the standard errors
#' use. Closed form: w = Xi^{-1} [xi + L'(L Xi^{-1} L')^{-1}(l - L Xi^{-1} xi)], with Xi the
#' covariance matrix of the comparison jumps, xi their covariances with the RD-period jump and
#' L w = l the constraints (dates centred at the RD period, which leaves the weights unchanged
#' and keeps L Xi^{-1} L' well conditioned with calendar years). When the constraints pin the
#' weights (as many comparison periods as constraint rows) the `w_pilot` weights are returned.
#' With `bc = TRUE` the same objects for the bias-corrected jumps (the quantities behind the
#' robust standard error), which gives the minimum-variance weights of the bias-corrected estimate.
#' @return list: `w` (named weights), `Xi`, `xi`, `pinned`.
#' @keywords internal
#' @noRd
.mv_weights <- function(fits, comps, t_rd, trend, scheme, w_pilot, bc = FALSE) {
  keys <- as.character(comps)
  rd <- as.character(t_rd)
  m <- length(keys)
  Xi <- matrix(0, m, m, dimnames = list(keys, keys))
  xi <- stats::setNames(numeric(m), keys)
  for (i in seq_len(m)) {
    Xi[i, i] <- if (bc) fits[[keys[i]]]$V_D_bc else fits[[keys[i]]]$V_D
    xi[i] <- .cov_scheme(fits[[keys[i]]], fits[[rd]], scheme, bc = bc)
    if (i < m) for (j in (i + 1L):m)
      Xi[i, j] <- Xi[j, i] <- .cov_scheme(fits[[keys[i]]], fits[[keys[j]]], scheme, bc = bc)
  }
  k <- if (trend == "constant") 1L else 2L          # number of constraint rows
  if (m == k)                                        # a single admissible weight vector
    return(list(w = w_pilot, Xi = Xi, xi = xi, pinned = TRUE))
  if (inherits(tryCatch(chol(Xi), error = function(e) e), "error"))
    .abort("the estimated covariance matrix of the comparison-period jumps is not positive ",
           "definite, so the minimum-variance weights are not defined (are two comparison periods ",
           "identical?); drop one of them or use `weighting = \"ols\"`.")
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

#' Order `v` by `keys` when it is named (every entry named, the names being exactly `keys`);
#' an unnamed `v` stays in the order given. Names on some entries only are an error.
#' @keywords internal
#' @noRd
.match_names <- function(v, keys, what, key_label) {
  nm <- names(v)
  if (is.null(nm) || all(!nzchar(nm))) return(v)
  if (any(!nzchar(nm)))
    .abort(what, " are named on some entries only: name every entry by its period, or none.")
  if (!setequal(nm, as.character(keys)) || anyDuplicated(nm))
    .abort("the names of ", what, " must be the ", key_label, ": ", paste(keys, collapse = ", "), ".")
  v[as.character(keys)]
}
