# Internal helpers shared by rd_typecont() and rd_compstable().
# All functions are prefixed `.` and are not exported.
#
# Design contract:
#   .build_types(data, x, time, id, c)       → list($wide, $period_types, $periods)
#   .joint_wald(thetas, Sigma)               → list(stat, df, p)
#   .wald_eigen(Delta, Sigma)               → list(stat, df, p)
#   .cell_bandwidth(y, x, c, kernel, ...)   → c(h, b)

# ---------------------------------------------------------------------------
# .build_types
# ---------------------------------------------------------------------------
#' Build the per-period type vectors from a long panel
#'
#' For each period t the "type" of unit i is the sign pattern of the OTHER
#' periods' running variables: V_{i,-t} = (1{R_{i,s} >= c})_{s != t}, recorded
#' as a "+"/"-" string (e.g. "+-").  Units exactly at the cutoff are treated as
#' above it (V_i = 1{R_i >= c}).  Units not observed in every period are
#' dropped from each period's type frame (their sign pattern is undefined).
#'
#' This is the single canonical implementation, shared by [rd_typecont()]
#' and [rd_homog()].
#'
#' @param data long data frame with one row per unit × period.
#' @param x,time,id column name strings for running variable, period, period id.
#' @param c cutoff.
#' @return A list:
#'   \item{wide}{data frame with columns id, one `R_<period>` per period, one
#'     `side_<period>` ("+"/"-"/NA) per period.}
#'   \item{period_types}{named list (one element per period t): a data frame
#'     with columns `id`, `R` (running variable in period t), and `type`
#'     (sign-pattern string of the other periods, e.g. "+-").}
#'   \item{periods}{character vector of period labels.}
#' @keywords internal
#' @noRd
.build_types <- function(data, x, time, id, c = 0) {
  periods <- sort(unique(data[[time]]))
  plab    <- as.character(periods)
  n_per   <- length(periods)

  # Pivot to wide: one row per unit, R and side per period.  Side is "+" if
  # R >= c (treated at the cutoff), "-" if R < c, NA if the unit is unobserved.
  wide <- data.frame(id = unique(data[[id]]))
  for (k in seq_along(periods)) {
    sub <- data[data[[time]] == periods[k], , drop = FALSE]
    Rk  <- sub[[x]][match(wide$id, sub[[id]])]
    wide[[paste0("R_",    plab[k])]] <- Rk
    wide[[paste0("side_", plab[k])]] <-
      ifelse(is.na(Rk), NA_character_, ifelse(Rk >= c, "+", "-"))
  }

  # For each period t, type = sign pattern of the OTHER periods' sides.  Units
  # with any unobserved other period (NA side) are dropped from that period.
  period_types <- stats::setNames(vector("list", n_per), plab)
  for (k in seq_along(periods)) {
    other_idx <- setdiff(seq_along(periods), k)
    Rk        <- wide[[paste0("R_", plab[k])]]
    if (length(other_idx) == 0L) {
      type <- rep("all", nrow(wide))
    } else {
      side_mat <- as.matrix(wide[, paste0("side_", plab[other_idx]), drop = FALSE])
      type <- apply(side_mat, 1L, function(r)
        if (anyNA(r)) NA_character_ else paste(r, collapse = ""))
    }
    keep <- !is.na(type) & !is.na(Rk)
    period_types[[k]] <- data.frame(
      id   = wide$id[keep],
      R    = Rk[keep],
      type = type[keep],
      stringsAsFactors = FALSE
    )
  }

  list(wide = wide, period_types = period_types, periods = plab)
}


# ---------------------------------------------------------------------------
# .joint_wald
# ---------------------------------------------------------------------------
#' Joint Wald test via Moore-Penrose inverse
#'
#' Given a vector of estimates theta and their covariance Sigma (possibly
#' singular), computes the Wald statistic theta' Sigma^+ theta where Sigma^+
#' is the Moore-Penrose pseudo-inverse, with degrees of freedom = rank(Sigma).
#'
#' @param theta numeric vector of jump estimates.
#' @param Sigma numeric square covariance matrix (same length as theta).
#' @return list with elements `stat` (chi-square statistic), `df` (rank of
#'   Sigma), `p` (p-value from chi-square distribution).
#' @keywords internal
#' @noRd
.joint_wald <- function(theta, Sigma) {
  # Moore-Penrose pseudo-inverse via SVD, truncating near-zero singular values.
  # Sigma is structurally rank-deficient here (the per-period type indicators
  # sum to 1, so each period contributes one exact-zero direction).  Use the
  # MASS::ginv relative tolerance sqrt(eps)*max(sv): a tighter tolerance leaves
  # a structural-zero singular value just above the cut on some LAPACK builds,
  # and its 1/sv blows the Wald statistic up (platform-dependent false rejects).
  sv <- svd(Sigma)
  tol <- sqrt(.Machine$double.eps) * max(sv$d)
  keep <- sv$d > tol
  df <- sum(keep)
  if (df == 0L) return(list(stat = 0, df = 0L, p = 1))
  d_inv <- ifelse(keep, 1 / sv$d, 0)
  Sigma_pinv <- sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
  stat <- as.numeric(t(theta) %*% Sigma_pinv %*% theta)
  p    <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  list(stat = stat, df = df, p = p)
}


# ---------------------------------------------------------------------------
# .mccrary
#
# ---------------------------------------------------------------------------
# .cell_bandwidth
# ---------------------------------------------------------------------------
#' Per-cell bandwidth for type-indicator and outcome RDs
#'
#' Single helper encapsulating the three bandwidth cases used in the four test
#' functions.  Priority: (1) explicit `h` overrides everything; (2) `bwselect
#' = "rot"`: use `rot_val` when supplied (full-sample pre-computed, as in
#' [rd_typecont()] and [rd_compstable()]) or compute 0.2 × range of `x` per
#' cell (as in [rd_homog()] and [rd_trendcell()]); (3) `bwselect = "cct"`:
#' call [rd_bw_cct()] and return its `h` and `b`.
#'
#' @param y outcome vector for the cell.
#' @param x running variable vector for the cell.
#' @param c cutoff.
#' @param kernel kernel name (passed to [rd_bw_cct()]).
#' @param h explicit bandwidth override; if non-`NULL`, overrides `bwselect`.
#' @param bwselect `"cct"` (default) or `"rot"`.
#' @param rot_val precomputed rule-of-thumb bandwidth for use when
#'   `bwselect = "rot"`.  If `NULL`, the 0.2 × range rule is applied to `x`.
#' @return Named numeric vector of length 2: `h` (main bandwidth) and `b`
#'   (pilot bandwidth for bias correction).
#' @keywords internal
#' @noRd
.cell_bandwidth <- function(y, x, c, kernel,
                            h = NULL, bwselect = "cct",
                            rot_val = NULL) {
  if (!is.null(h)) return(c(h = h, b = h))
  if (bwselect == "rot") {
    hw <- if (!is.null(rot_val)) {
      rot_val
    } else {
      rng <- diff(range(x[is.finite(x)], na.rm = TRUE))
      0.2 * rng
    }
    return(c(h = hw, b = hw))
  }
  # bwselect = "cct"
  bw <- rd_bw_cct(y, x, c = c, kernel = kernel)
  c(h = bw[["h"]], b = bw[["b"]])
}


# ---------------------------------------------------------------------------
# .wald_eigen
# ---------------------------------------------------------------------------
#' Conservative eigen-decomposition Wald test
#'
#' Given a contrast vector `Delta` and its (possibly indefinite) covariance
#' matrix `Sigma`, computes the Wald statistic \eqn{\Delta' \Sigma^+ \Delta}
#' using an eigen pseudo-inverse that **drops non-positive eigenvalue
#' directions**.  This is the conservative choice for contrast covariances
#' that are linear combinations of estimated covariances and can be
#' numerically indefinite.
#'
#' This is deliberately DISTINCT from [.joint_wald()], which uses an
#' SVD-based Moore-Penrose pseudo-inverse and is designed for structurally
#' rank-deficient but positive-semidefinite matrices.
#'
#' @param Delta numeric contrast vector.
#' @param Sigma numeric square matrix (covariance of `Delta`).
#' @return list with elements `stat` (Wald statistic), `df` (number of
#'   positive eigenvalue directions retained), and `p` (chi-square p-value).
#'   Returns `list(stat = NA_real_, df = 0L, p = NA_real_)` when `df = 0`.
#' @keywords internal
#' @noRd
.wald_eigen <- function(Delta, Sigma) {
  K   <- length(Delta)
  ev  <- eigen(Sigma, symmetric = TRUE)
  tol <- max(abs(ev$values)) * K * .Machine$double.eps^0.5
  pos <- ev$values > tol
  df  <- sum(pos)
  if (df == 0L) return(list(stat = NA_real_, df = 0L, p = NA_real_))
  Sigma_inv <- ev$vectors[, pos, drop = FALSE] %*%
               diag(1 / ev$values[pos], nrow = sum(pos)) %*%
               t(ev$vectors[, pos, drop = FALSE])
  W  <- as.numeric(t(Delta) %*% Sigma_inv %*% Delta)
  pv <- stats::pchisq(W, df = df, lower.tail = FALSE)
  list(stat = W, df = df, p = pv)
}
