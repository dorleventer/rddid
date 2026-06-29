# Internal helpers shared by rd_typecont() (A7) and rd_compstable() (A8).
# All functions are prefixed `.` and are not exported.
#
# Design contract:
#   .build_types(data, x, time, id, c)       → list($wide, $period_types, $periods)
#   .joint_wald(thetas, Sigma)               → list(stat, df, p)
#   .wald_eigen(Delta, Sigma)               → list(stat, df, p)
#   .cell_bandwidth(y, x, c, kernel, ...)   → c(h, b)
#   .q_rot(z, type, c)                       → integer q (Canay-Kamat rule of thumb)
#   .mccrary(x, h)                           → scalar p-value

# ---------------------------------------------------------------------------
# .q_rot
# ---------------------------------------------------------------------------
#' Canay-Kamat (2018) rule-of-thumb number of nearest observations per side
#'
#' The approximate permutation test uses the `q` observations nearest the cutoff
#' on each side and is valid with `q` FIXED as `n` grows.  With a fixed `q` it
#' over-rejects in finite samples when the covariate's conditional distribution
#' varies steeply in the running variable at the cutoff (Canay & Kamat 2018,
#' Model 7).  Their feasible rule of thumb (eq. 15 and Appendix D.1) adapts `q`:
#' \deqn{q = \lceil \hat f(0)\, \hat\sigma_z \sqrt{1-\hat\rho^2}\; n^{0.9}/\log n
#'        \rceil,} bounded to \eqn{[10,\, n^{0.9}/\log n]}.  \eqn{\hat f(0)} is a
#' Gaussian-kernel density of the running variable at the cutoff and
#' \eqn{\hat\rho} measures the running-variable / covariate association — a
#' steeper relationship (large \eqn{\rho}) shrinks `q`.  Here the covariate is
#' the categorical type, so we use the correlation ratio
#' \eqn{\eta=\sqrt{SS_{\text{between}}/SS_{\text{total}}}} of the running
#' variable grouped by type (the generalisation of \eqn{|\rho|} to a categorical
#' covariate; it equals \eqn{|\rho|} when the type is binary).
#'
#' @param z numeric running variable.
#' @param type covariate (the per-unit type label) the same length as `z`.
#' @param c cutoff.
#' @return integer `q`.
#' @keywords internal
#' @noRd
.q_rot <- function(z, type, c = 0) {
  n  <- length(z)
  ub <- n^0.9 / log(n)
  sd_z <- stats::sd(z)
  if (n < 20L || !is.finite(sd_z) || sd_z == 0) return(max(10L, 1L))
  dens <- stats::density(z)
  f0   <- stats::approx(dens$x, dens$y, xout = c, rule = 2)$y
  # correlation ratio eta of z by type (categorical generalisation of |rho|)
  grand  <- mean(z)
  ss_tot <- sum((z - grand)^2)
  ss_bet <- sum(vapply(split(z, type), function(g)
    length(g) * (mean(g) - grand)^2, numeric(1)))
  rho2 <- if (ss_tot > 0) min(1, ss_bet / ss_tot) else 0
  qhat <- f0 * sd_z * sqrt(max(0, 1 - rho2)) * ub
  as.integer(max(10, min(ceiling(qhat), floor(ub))))
}

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
#' This is the single canonical implementation, shared by [rd_typecont()] (A7)
#' and [rd_homog()] (A9).
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
# Note: there is no shared single-cell Canay-Kamat helper.  rd_typecont() and
# rd_compstable() each run a *joint* Canay-Kamat permutation that draws one
# shared per-period unit-level shuffle across all type columns, so the test is
# inlined in those functions rather than factored out here.
# ---------------------------------------------------------------------------
#' McCrary (2008) density-discontinuity test
#'
#' Bins the running variable into a histogram (bin width = 2*sd(x)*n^(-1/2)),
#' fits a local-linear density estimator on each side of the cutoff (triangular
#' kernel, bandwidth h), and tests log f_(+)/f_(-) = 0.  Sandwich standard
#' error uses the Poisson bin-count variance approximation from McCrary (2008).
#'
#' @param x numeric running variable, centred at cutoff 0.
#' @param h bandwidth for the local-linear density fit.
#' @return scalar p-value (NA if fewer than 30 observations or numerical
#'   failure).
#' @keywords internal
#' @noRd
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


.mccrary <- function(x, h) {
  n <- length(x)
  if (n < 30L) return(NA_real_)
  binw <- 2 * stats::sd(x) * n^(-1 / 2)
  lo   <- floor(min(x) / binw) * binw
  hi   <- ceiling(max(x) / binw) * binw
  brks <- seq(lo, hi, by = binw)
  mids <- brks[-length(brks)] + binw / 2
  cnt  <- tabulate(findInterval(x, brks, rightmost.closed = TRUE),
                   nbins = length(mids))
  Y    <- cnt / (n * binw)

  fit_side <- function(keep) {
    g <- mids[keep]; y <- Y[keep]
    w <- pmax(0, 1 - abs(g) / h)
    use <- w > 0
    g <- g[use]; y <- y[use]; w <- w[use]
    if (length(g) < 3L) return(NULL)
    Z    <- cbind(1, g)
    Zw   <- Z * w
    XtWX <- crossprod(Z, Zw)
    XtWXi <- tryCatch(solve(XtWX), error = function(e) NULL)
    if (is.null(XtWXi)) return(NULL)
    beta <- XtWXi %*% crossprod(Zw, y)
    vj   <- pmax(y, 1e-8) / (n * binw)        # Poisson bin variance
    Vb   <- XtWXi %*% crossprod(Zw, vj * Zw) %*% XtWXi
    c(f0 = beta[1L], var = Vb[1L, 1L])
  }

  Rr <- fit_side(mids > 0)
  Ll <- fit_side(mids < 0)
  if (is.null(Rr) || is.null(Ll)) return(NA_real_)
  if (Rr["f0"] <= 0 || Ll["f0"] <= 0)  return(NA_real_)

  theta <- log(Rr["f0"]) - log(Ll["f0"])
  se    <- sqrt(Rr["var"] / Rr["f0"]^2 + Ll["var"] / Ll["f0"]^2)
  unname(2 * stats::pnorm(-abs(theta / se)))
}
