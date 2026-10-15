# assumption_tests_helpers.R -- internals shared by the four assumption tests: .build_types()
# (a unit's type = its sides of the cutoff in the other periods), .cell_bandwidth() (CCT or rule of
# thumb per cell), .joint_wald() and .wald_eigen() (two Wald statistics; see the notes on each),
# .gauss_joint() / .joint_tests() (sup-t and sWald, which do not invert the covariance) and
# .common_bandwidth() (one pooled CCT bandwidth for the outcome tests).
# The code-to-paper map for the tests is dev/tests_map.md (developer material, not shipped).

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
  if (anyDuplicated(data[, c(id, time)]))
    stop("some (", id, ", ", time, ") pairs appear more than once: each unit may appear once per ",
         "period.", call. = FALSE)
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
#'   Sigma), `p` (p-value from chi-square distribution); `dropped` (the number of
#'   shares left out for having no variation near the cutoff) when it is positive.
#' @keywords internal
#' @noRd
.joint_wald <- function(theta, Sigma, scale = 1) {
  if (length(theta) == 0L) return(list(stat = NA_real_, df = 0L, p = NA_real_))
  # A share whose jump has a standard error that is zero to working precision is constant on
  # each side of the cutoff within the window (the type indicators are 0/1, so the scale is 1).
  # - Its jump is not zero: the share jumps deterministically (a type defined by the side of
  #   the cutoff in the tested period, one switcher at the edge of the window). The statistic
  #   would be rounding noise blown up, so the test is undefined.
  # - Its jump is zero too: the type is absent (or universal) near the cutoff on both sides.
  #   The share carries no information and its row of Sigma is zero, so it is left out exactly
  #   and the remaining shares are tested (2026-10-09: this case used to make the whole test
  #   NA as well).
  tol0    <- sqrt(.Machine$double.eps) * scale
  zero_se <- sqrt(pmax(diag(Sigma), 0)) <= tol0
  if (any(zero_se & abs(theta) > tol0))
    return(list(stat = NA_real_, df = 0L, p = NA_real_, degenerate = TRUE))
  n_dropped <- sum(zero_se)
  if (n_dropped > 0L) {
    theta <- theta[!zero_se]
    Sigma <- Sigma[!zero_se, !zero_se, drop = FALSE]
    if (length(theta) == 0L)
      return(list(stat = NA_real_, df = 0L, p = NA_real_, dropped = n_dropped))
  }
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
  if (df == 0L) return(list(stat = NA_real_, df = 0L, p = NA_real_))   # nothing testable
  d_inv <- ifelse(keep, 1 / sv$d, 0)
  Sigma_pinv <- sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
  stat <- as.numeric(t(theta) %*% Sigma_pinv %*% theta)
  p    <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  out  <- list(stat = stat, df = df, p = p)
  if (n_dropped > 0L) out$dropped <- n_dropped
  out
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
                            h = NULL, bwselect = "cct", p = 1L) {
  if (!is.null(h)) return(c(h = h, b = h))
  if (bwselect == "rot") {
    rng <- diff(range(x[is.finite(x)], na.rm = TRUE))
    hw  <- 0.2 * rng
    return(c(h = hw, b = hw))
  }
  # bwselect = "cct": MSE-optimal for the order-p fit the cell will run (2026-10-08 fix:
  # the bandwidth used to be chosen at p = 1 whatever p the fit used)
  bw <- rd_bw_cct(y, x, c = c, p = p, kernel = kernel)
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
.wald_eigen <- function(Delta, Sigma, scale = NULL) {
  K   <- length(Delta)
  # A covariance that is zero to working precision (no residual variation in the cells, e.g. a
  # constant outcome) has eigenvalues made of rounding noise; the relative tolerance below
  # would keep them and turn noise into a chi-squared statistic. Compare the standard errors
  # with the outcome's scale first (2026-10-08 fix).
  if (!is.null(scale) && is.finite(scale) &&
      max(sqrt(pmax(diag(Sigma), 0))) <= sqrt(.Machine$double.eps) * scale)
    return(list(stat = NA_real_, df = 0L, p = NA_real_, degenerate = TRUE))
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

# ---------------------------------------------------------------------------
# .gauss_joint, .joint_tests: joint tests that do not invert the estimated covariance
# ---------------------------------------------------------------------------
# Why (T1652 size study, rd-did/docs/test_size_morning_briefing.md): with many contrasts the
# chi-squared Wald over-rejects because it inverts a noisy estimated covariance (7% at 9 df,
# 15% at 28 df, n = 1,000). The max-|z| test with a sup-t critical value and the sum of squared
# t-statistics (sWald, Fusejima, Ishihara & Sawada 2025, J. Econometrics) use the covariance only
# through the correlations of the standardized jumps, never its inverse, and hold size.

#' Run `expr` under a fixed seed and RNG kind, leaving the caller's RNG state untouched
#' @keywords internal
#' @noRd
.with_seed <- function(seed, expr) {
  genv     <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = genv, inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    # .Random.seed encodes the RNG kinds, so restoring it restores them; without a seed the kinds
    # are reset explicitly (quietly: an old "Rounding" sampler would warn on every call)
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else {
      suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
      if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
    }
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}

#' Degrees of freedom of the robust variance of a jump (studentized joint tests)
#'
#' The variance of D = (+) intercept - (-) intercept is a sum of squared weighted residuals,
#' sum_i a_i^2 e_i^2. Under a homoskedastic working model its degrees of freedom are
#' (sum_i a_i^2)^2 / sum_i a_i^4 (Satterthwaite; the Bell and McCaffrey 2002 small-sample df
#' without the residual-maker correction, which changed nu by under 1% in a check). `a` are the
#' intercept weights of rd_period() on both sides (bias-corrected weights when `bc`).
#' @keywords internal
#' @noRd
.jump_df <- function(fit, bc = TRUE) {
  a  <- c(if (bc) fit$sides$`+`$a_bc else fit$sides$`+`$a,
          if (bc) fit$sides$`-`$a_bc else fit$sides$`-`$a)
  a2 <- a^2
  sum(a2)^2 / sum(a2^2)
}

#' sup-t and sWald tests of theta = 0 from the standardized jumps, normal and studentized
#'
#' t_k = theta_k / se_k with se = sqrt(diag(Sigma)). sup-t: max_k |t_k|; sWald: sum_k t_k^2. Each
#' is compared with its distribution when Z ~ N(0, R), R the correlation matrix of Sigma (negative
#' eigenvalues of an indefinite estimate set to zero), under two references:
#'   normal       T_k = Z_k;
#'   studentized  T_k = Z_k / sqrt(S_k), S_k ~ chisq(nu_k) / nu_k independent: the estimated
#'                standard error of each jump is itself noisy (T1682/T1687: this noise, with
#'                nu ~ 50-120, is what leaves the normal sup-t at 6-7% instead of 5%).
#' With one jump both are exact: the two-sided normal test, and the two-sided t(nu) test. Otherwise
#' the references are simulated with `n_sim` draws under a fixed seed (reproducible; the caller's
#' RNG state is restored), p = (1 + #{draws >= observed}) / (n_sim + 1). Zero-SE jumps are handled
#' as in .joint_wald().
#' @param nu degrees of freedom of each jump's variance (.jump_df()); NULL = no studentized version.
#' @return list(supt_t, supt_normal, swald_t, swald_normal), each list(stat, df, p) and, for the
#'   sup-t ones, `cell` (the index in `theta` of the largest |t|); `df` is the number of jumps.
#' @keywords internal
#' @noRd
.gauss_joint <- function(theta, Sigma, nu = NULL, scale = 1, n_sim = 100000L, seed = 20261015L) {
  empty <- function(extra = list()) c(list(stat = NA_real_, df = 0L, p = NA_real_), extra)
  na_out <- function(extra = list()) {
    s <- empty(c(list(cell = NA_integer_), extra)); w <- empty(extra)
    list(supt_t = s, supt_normal = s, swald_t = w, swald_normal = w)
  }
  if (length(theta) == 0L) return(na_out())
  if (is.null(nu)) nu <- rep(Inf, length(theta))
  tol0    <- sqrt(.Machine$double.eps) * scale
  se      <- sqrt(pmax(diag(Sigma), 0))
  zero_se <- se <= tol0
  if (any(zero_se & abs(theta) > tol0)) return(na_out(list(degenerate = TRUE)))
  idx <- which(!zero_se)
  if (length(idx) == 0L) return(na_out(list(dropped = sum(zero_se))))
  tt <- theta[idx] / se[idx]
  nu <- nu[idx]
  K  <- length(tt)
  m_obs <- max(abs(tt))
  s_obs <- sum(tt^2)
  cell  <- idx[which.max(abs(tt))]
  if (K == 1L) {
    p_n <- 2 * stats::pnorm(-abs(tt))
    p_t <- if (is.finite(nu)) 2 * stats::pt(-abs(tt), df = nu) else p_n
    out <- list(supt_t = list(stat = m_obs, df = 1L, p = p_t, cell = cell),
                supt_normal = list(stat = m_obs, df = 1L, p = p_n, cell = cell),
                swald_t = list(stat = s_obs, df = 1L, p = p_t),
                swald_normal = list(stat = s_obs, df = 1L, p = p_n))
  } else {
    R <- Sigma[idx, idx, drop = FALSE] / tcrossprod(se[idx])
    e <- eigen((R + t(R)) / 2, symmetric = TRUE)
    L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), K)
    n_max_n <- n_max_t <- n_sum_n <- n_sum_t <- 0
    .with_seed(seed, {
      left <- n_sim
      while (left > 0) {   # in chunks, so K x n_sim never sits in memory at once
        b  <- min(left, 20000L)
        Z  <- L %*% matrix(stats::rnorm(K * b), K, b)
        fin <- is.finite(nu)          # nu = Inf: the standard error is taken as known (S = 1)
        S  <- matrix(1, K, b)
        if (any(fin))
          S[fin, ] <- matrix(stats::rchisq(sum(fin) * b, df = rep(nu[fin], b)), sum(fin), b) / nu[fin]
        Zt <- Z / sqrt(S)
        mn <- abs(Z[1L, ]);  for (k in seq_len(K)[-1L]) mn <- pmax(mn, abs(Z[k, ]))
        mt <- abs(Zt[1L, ]); for (k in seq_len(K)[-1L]) mt <- pmax(mt, abs(Zt[k, ]))
        n_max_n <- n_max_n + sum(mn >= m_obs)
        n_max_t <- n_max_t + sum(mt >= m_obs)
        n_sum_n <- n_sum_n + sum(colSums(Z^2) >= s_obs)
        n_sum_t <- n_sum_t + sum(colSums(Zt^2) >= s_obs)
        left <- left - b
      }
    })
    pp <- function(n) (1 + n) / (n_sim + 1)
    out <- list(supt_t = list(stat = m_obs, df = K, p = pp(n_max_t), cell = cell),
                supt_normal = list(stat = m_obs, df = K, p = pp(n_max_n), cell = cell),
                swald_t = list(stat = s_obs, df = K, p = pp(n_sum_t)),
                swald_normal = list(stat = s_obs, df = K, p = pp(n_sum_n)))
  }
  if (any(zero_se)) out <- lapply(out, function(z) c(z, list(dropped = sum(zero_se))))
  out
}

#' The joint tests of theta = 0: the chi-squared Wald (.joint_wald()) and sup-t / sWald with normal
#' and studentized references (.gauss_joint()); returns a named list of list(stat, df, p, ...).
#' @keywords internal
#' @noRd
.joint_tests <- function(theta, Sigma, nu = NULL, scale = 1) {
  if (length(theta) == 0L) {
    empty <- list(stat = NA_real_, df = 0L, p = NA_real_)
    return(list(wald = empty, supt_t = c(empty, list(cell = NA_integer_)),
                supt_normal = c(empty, list(cell = NA_integer_)), swald_t = empty,
                swald_normal = empty))
  }
  c(list(wald = .joint_wald(theta, Sigma, scale = scale)),
    .gauss_joint(theta, Sigma, nu = nu, scale = scale))
}

#' Stop unless `studentize` is a single TRUE or FALSE
#' @keywords internal
#' @noRd
.check_studentize <- function(studentize, fn) {
  if (!is.logical(studentize) || length(studentize) != 1L || is.na(studentize))
    stop(fn, ": `studentize` must be TRUE or FALSE.", call. = FALSE)
  invisible(TRUE)
}

#' The element of .joint_tests() for a choice of `joint` and `studentize`
#' @keywords internal
#' @noRd
.pick_joint <- function(tests, joint, studentize) {
  if (joint == "wald") return(tests$wald)
  tests[[paste0(joint, if (studentize) "_t" else "_normal")]]
}

#' The joint tests as a data frame (method, reference, statistic, df, p_value), for the results
#' @keywords internal
#' @noRd
.joint_tests_df <- function(tests) {
  keys <- c("supt_t", "supt_normal", "swald_t", "swald_normal", "wald")
  data.frame(method    = c("supt", "supt", "swald", "swald", "wald"),
             reference = c("t", "normal", "t", "normal", "chi-squared"),
             statistic = vapply(tests[keys], function(z) z$stat, numeric(1)),
             df        = vapply(tests[keys], function(z) as.integer(z$df), integer(1)),
             p_value   = vapply(tests[keys], function(z) z$p, numeric(1)),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Common bandwidth of rd_homog() and rd_trendcell() under bwselect = "common": the CCT
#' MSE-optimal h of the outcome on the running variable, pooling every comparison period (all
#' units, all types), used as main and pilot bandwidth in every cell.
#' @keywords internal
#' @noRd
.common_bandwidth <- function(data, y, x, time, comparisons, cutoff, kernel, p) {
  # Why (T1652): a bandwidth chosen on each thin cell's own data makes every cell variance
  # 10-15% too small (the variance formula treats h as fixed), so the tests reject 7-8% at
  # n = 1,000; one h chosen on the pooled comparison periods (n(T - 1) rows instead of ~n/2) is
  # far less noisy and gives nominal size (5.6% / 5.9%), at about 8 points of power.
  pooled <- data[data[[time]] %in% comparisons, , drop = FALSE]
  pooled <- pooled[is.finite(pooled[[y]]) & is.finite(pooled[[x]]), , drop = FALSE]
  unname(rd_bw_cct(pooled[[y]], pooled[[x]], c = cutoff, p = p, kernel = kernel)[["h"]])
}

#' Rule-of-thumb bandwidth of rd_typecont() and rd_compstable(): half the IQR of the running
#' variable (its sd if the IQR is zero), used as both h and b in every cell
#' @keywords internal
#' @noRd
.rot_bandwidth_iqr <- function(all_x) {
  h <- 0.5 * stats::IQR(all_x)
  if (h <= 0) h <- stats::sd(all_x)
  h
}

#' Stop when no unit changes side of the cutoff between the periods used: every unit then has
#' the same type in every period, and the tests of the assumptions for a time-varying running
#' variable are not defined (2026-10-08: they used to return degenerate statistics)
#' @keywords internal
#' @noRd
.stop_if_no_switchers <- function(wide, periods, fn) {
  side_cols <- paste0("side_", periods)
  side_cols <- side_cols[side_cols %in% names(wide)]
  if (length(side_cols) < 2L) return(invisible(FALSE))
  sides <- wide[, side_cols, drop = FALSE]
  observed <- rowSums(!is.na(sides))
  if (!any(observed >= 2L))
    stop(fn, ": no unit is observed in two of the periods used (the ids never repeat, as in a ",
         "repeated cross-section); the tests of the assumptions need a panel.", call. = FALSE)
  switches <- apply(sides, 1L, function(s) { s <- s[!is.na(s)]; length(unique(s)) > 1L })
  if (!any(switches))
    stop(fn, ": no unit changes side of the cutoff between the periods, so every unit has the ",
         "same type in every period and the test is not defined. The tests are for a running ",
         "variable that varies over time (see ?", fn, ").", call. = FALSE)
  invisible(TRUE)
}
