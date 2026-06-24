# Internal helpers shared by rd_typecont() (A7) and rd_compstable() (A8).
# All functions are prefixed `.` and are not exported.
#
# Design contract:
#   .build_types(data, x, time, id, c)   → list with $wide, $types per period
#   .joint_wald(thetas, Sigma, df)        → list(stat, df, p)
#   .ck_perm(x, g, q, S)                 → scalar p-value
#   .mccrary(x, h)                       → scalar p-value

# ---------------------------------------------------------------------------
# .build_types
# ---------------------------------------------------------------------------
#' Build the per-period type vectors from a long panel
#'
#' For each period t the "type" of unit i is the sign pattern of the OTHER
#' periods' running variables: V_{i,-t} = (1{R_{i,s} >= c})_{s != t}.
#' With P periods this produces 2^(P-1) possible types per period, encoded
#' as an integer (binary, with the first other-period as the least significant
#' bit).
#'
#' @param data long data frame with one row per unit × period.
#' @param x,time,id column name strings for running variable, period, period id.
#' @param c cutoff.
#' @return A list:
#'   \item{wide}{data frame with columns id, one `R_<period>` per period, one
#'     `side_<period>` (0/1) per period.}
#'   \item{period_types}{named list (one element per period t): a data frame
#'     with columns `id`, `R` (running variable in period t), `type` (integer
#'     encoding V_{i,-t}), and `type_str` (e.g. "01").}
#'   \item{periods}{character vector of period labels.}
#' @keywords internal
#' @noRd
.build_types <- function(data, x, time, id, c = 0) {
  periods <- sort(unique(data[[time]]))
  plab    <- as.character(periods)
  n_per   <- length(periods)

  # Pivot to wide: one row per unit, one column per period for R and for side
  wide <- data.frame(id = unique(data[[id]]))
  for (k in seq_along(periods)) {
    sub <- data[data[[time]] == periods[k], , drop = FALSE]
    m   <- match(wide$id, sub[[id]])
    wide[[paste0("R_",    plab[k])]] <- sub[[x]][m]
    wide[[paste0("side_", plab[k])]] <- as.integer(sub[[x]][m] >= c)
  }

  # For each period t, build the type = binary encoding of the OTHER periods' sides
  period_types <- stats::setNames(vector("list", n_per), plab)
  for (k in seq_along(periods)) {
    other_idx <- setdiff(seq_along(periods), k)
    # encode: sum side_{other[j]} * 2^(j-1)
    side_mat <- as.matrix(wide[, paste0("side_", plab[other_idx]), drop = FALSE])
    type_int <- as.integer(side_mat %*% (2L ^ (seq_along(other_idx) - 1L)))

    # type_str for labeling
    type_str <- apply(side_mat, 1L, paste, collapse = "")

    period_types[[k]] <- data.frame(
      id       = wide$id,
      R        = wide[[paste0("R_", plab[k])]],
      type     = type_int,
      type_str = type_str
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
  # Moore-Penrose pseudo-inverse via SVD, truncating near-zero singular values
  sv <- svd(Sigma)
  tol <- max(dim(Sigma)) * .Machine$double.eps * max(sv$d)
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
# .ck_perm
# ---------------------------------------------------------------------------
#' Canay-Kamat (2018) approximate sign randomisation test
#'
#' Takes the q observations nearest the cutoff on each side, and tests the
#' null that the covariate distribution is the same on both sides by permuting
#' the side label. The test statistic is the absolute difference in covariate
#' means across sides. Works on a binary covariate (the type indicator for a
#' single type value).
#'
#' For the joint test over multiple (period, type) pairs, call once with the
#' concatenated covariate and x vectors (or call per-cell and combine — the
#' caller decides).
#'
#' @param x numeric running variable, centred at the cutoff (cutoff = 0).
#' @param g numeric covariate (e.g. type indicator).
#' @param q number of observations to use from each side.
#' @param S number of permutation replications.
#' @return scalar p-value.
#' @keywords internal
#' @noRd
.ck_perm <- function(x, g, q, S = 499L) {
  # Nearest q on each side
  pos  <- x >= 0
  neg  <- x <  0
  rt   <- order(x[pos])[seq_len(min(q, sum(pos)))]
  lt   <- order(-x[neg])[seq_len(min(q, sum(neg)))]
  gr   <- g[pos][rt]
  gl   <- g[neg][lt]
  pool <- c(gr, gl)
  nr   <- length(gr)
  if (nr < 1L || (length(pool) - nr) < 1L) return(NA_real_)
  obs  <- abs(mean(gr) - mean(gl))
  perm <- replicate(S, {
    idx <- sample.int(length(pool), nr)
    abs(mean(pool[idx]) - mean(pool[-idx]))
  })
  (1 + sum(perm >= obs)) / (S + 1)
}


# ---------------------------------------------------------------------------
# .mccrary
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
