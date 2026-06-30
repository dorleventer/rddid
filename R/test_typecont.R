#' Test the continuity of the type distribution (Assumption A7)
#'
#' Runs the four tests for Assumption A7 ("continuity of the type distribution")
#' from Leventer and Nevo.  The "type" of unit \eqn{i} in period \eqn{t} is the
#' sign pattern of its running variables in the OTHER periods,
#' \eqn{\mathbf{V}_{i,-t} = (1\{R_{i,s} \ge c\})_{s \ne t}}.
#'
#' **Four tests, in order of interpretive strength:**
#' \enumerate{
#'   \item **LL-Wald** (necessary AND sufficient): for each period \eqn{t} and
#'     each type value \eqn{v}, estimate the local-linear RD jump of the type
#'     indicator \eqn{1\{\mathbf{V}_{i,-t} = v\}} on the running variable
#'     \eqn{R_{i,t}}.  The jump \eqn{\hat\pi_{t,(+)}(v) - \hat\pi_{t,(-)}(v)}
#'     is the output of [rd_period()] with a binary outcome.  The joint Wald
#'     statistic across all (period, type) pairs uses a Moore-Penrose
#'     pseudo-inverse because within each period the type indicators sum to 1
#'     (the block is singular); df = rank of the covariance matrix.  The
#'     covariance is built from the per-unit influence vectors returned by
#'     [rd_period()], using the same within-period and cross-period id-matching
#'     as the main estimator, scheme-aware.
#'   \item **Canay-Kamat permutation** (necessary AND sufficient): uses the
#'     \eqn{q} observations nearest the cutoff on each side to test whether the
#'     type distribution is the same on both sides, via approximate sign
#'     randomisation (Canay and Kamat 2018).  Run jointly over all types and
#'     periods: the test statistic is the sum of absolute mean differences
#'     across (period, type) pairs.
#'   \item **McCrary within-type** (sufficient, NOT necessary): McCrary (2008)
#'     density test run separately for each type value within each period, then
#'     combined across types via Bonferroni (minimum p-value times number of
#'     tests).  If every within-type density is continuous at \eqn{c} then the
#'     type shares are continuous (sufficiency), but shares can remain continuous
#'     even when all within-type densities jump by a common proportional factor
#'     (not necessary).
#'   \item **McCrary pooled** (NEITHER sufficient NOR necessary): McCrary (2008)
#'     density test on the pooled sample, run per period.  Neither sufficient
#'     (type shares can jump while pooled density is smooth) nor necessary
#'     (pooled density can jump while shares stay continuous).
#' }
#'
#' @param data a long data frame, one row per unit × period (balanced panel).
#' @param x column name (string) for the running variable.
#' @param time column name (string) for the period.
#' @param id column name (string) for the unit identifier.
#' @param c cutoff (default 0).
#' @param h bandwidth.  If `NULL`, the bandwidth is determined by `bwselect`;
#'   an explicit numeric value overrides `bwselect` and is used directly.
#' @param bwselect bandwidth selection rule when `h = NULL`: `"cct"` (default)
#'   computes a per-cell CCT MSE-optimal bandwidth via [rd_bw_cct()] for each
#'   (period, type) RD; `"rot"` uses the `0.5 * IQR(x)` rule of thumb applied
#'   to the full sample (the previous default behaviour).  Ignored when `h` is
#'   supplied explicitly.
#' @param q number of observations nearest the cutoff on each side for the
#'   Canay-Kamat permutation test. `NULL` (default) selects `q` per period by
#'   the Canay & Kamat (2018) rule of thumb; this is the
#'   recommended choice, since a fixed `q` over-rejects in finite samples when
#'   the type distribution varies steeply in the running variable at the cutoff.
#'   Supply an integer to force a fixed `q` on every period. The per-period `q`
#'   actually used is returned in `meta$q_used`.
#' @param S number of permutation replications for the Canay-Kamat test
#'   (default 499).
#' @param kernel kernel for the local-linear RD: `"triangular"` (default),
#'   `"epanechnikov"`, or `"uniform"`.
#' @param scheme covariance scheme for the joint Wald: `"auto"` detects from
#'   the data (same logic as [rddid()]), or one of `"cs"`, `"pc"`, `"pv"`.
#' @param bc use robust bias-corrected jumps and variances in the LL-Wald
#'   (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test
#'   with the bias-corrected [rddid()] estimator; `FALSE` uses the conventional
#'   local-linear jumps and variances. Does not affect the Canay-Kamat
#'   permutation test (which operates on raw near-cutoff indicator means).
#' @param ... currently unused.
#'
#' @return An object of class `"rd_typecont"`, a named list with:
#'   \item{ll_wald}{list with `stat` (chi-square), `df`, `p`.}
#'   \item{per_period}{named list (by period) of the per-period components the
#'     joint test aggregates; each entry has `ll_wald` (that period's own
#'     LL-Wald `stat`/`df`/`p`, restricted to its kept contrasts) and `ck_p`
#'     (per-period Canay-Kamat p-value, read off the same permutation draws as
#'     the joint statistic, or `NA` for a period with no active type cell).}
#'   \item{ck_perm}{list with `stat` (observed sum of |mean diffs|), `p`.}
#'   \item{mccrary_within}{data frame with columns `period`, `type`,
#'     `p_raw` (per-type McCrary p-value); plus `p_bonf` (Bonferroni-adjusted
#'     p-value for the period, minimum × number of types).}
#'   \item{mccrary_pooled}{data frame with columns `period`, `p` (per-period
#'     McCrary p-value on pooled sample).}
#'   \item{meta}{list with `periods`, `type_values`, `h` (NA when
#'     `bwselect = "cct"`), `bwselect`, `q` (`"rot"` when the Canay-Kamat
#'     rule of thumb is used, otherwise the integer supplied), `q_used`
#'     (per-period integer `q` actually used), `S`, `scheme`, `bc`.}
#' @export
rd_typecont <- function(data, x, time, id,
                        c = 0,
                        h = NULL,
                        bwselect = c("cct", "rot"),
                        q = NULL,
                        S = 499L,
                        kernel = "triangular",
                        scheme = c("auto", "cs", "pc", "pv"),
                        bc = TRUE,
                        ...) {
  scheme   <- match.arg(scheme)
  bwselect <- match.arg(bwselect)

  # ----- input checks -------------------------------------------------------
  for (nm in c(x, time, id)) {
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")
  }
  data  <- data[stats::complete.cases(data[, c(x, time, id)]), , drop = FALSE]
  periods <- sort(unique(data[[time]]))
  plab    <- as.character(periods)
  P       <- length(periods)
  if (P < 2L) stop("need at least 2 periods to define a type.")

  # ----- default bandwidth --------------------------------------------------
  # h_rot: IQR-based pilot, computed whenever h is NULL (regardless of bwselect)
  # so that scheme detection and McCrary always have a finite window.
  # h_aux: the bandwidth actually used for scheme detection + McCrary tests:
  #   equals h when h is explicit; equals h_rot otherwise.
  # When bwselect = "rot" and h = NULL, h is set to h_rot (current behaviour).
  # When bwselect = "cct" and h = NULL, h stays NULL; per-cell CCT is computed
  # inside the LL-Wald loop below.
  h_rot <- NULL
  if (is.null(h)) {
    all_x <- data[[x]]
    h_rot <- 0.5 * stats::IQR(all_x)
    if (h_rot <= 0) h_rot <- stats::sd(all_x)
    if (bwselect == "rot") h <- h_rot
  }
  h_aux <- if (!is.null(h)) h else h_rot

  # ----- build types --------------------------------------------------------
  # Per-period frames (id, R, type) from the shared canonical builder in
  # R/test_helpers.R.  `type` is the "+"/"-" sign-pattern string of the other
  # periods; units at the cutoff are treated as above it.
  pt <- .build_types(data, x, time, id, c = c)$period_types

  # All type values that appear anywhere across all periods
  all_type_values <- sort(unique(unlist(lapply(pt, `[[`, "type"))))
  n_types <- length(all_type_values)

  # ----- detect scheme -----------------------------------------------------
  # Classify from the in-window units of each period via the shared primitive
  # (side = 1{R >= c}, treated at the cutoff).
  if (scheme == "auto") {
    long <- do.call(rbind, lapply(plab, function(k) {
      df_k  <- pt[[k]]
      inwin <- abs(df_k$R - c) <= h_aux
      data.frame(period = k,
                 id     = df_k$id[inwin],
                 side   = as.integer(df_k$R[inwin] >= c))
    }))
    use_scheme <- .scheme_from_long(long)
  } else {
    use_scheme <- scheme
  }

  # ----- (1) LL-Wald -------------------------------------------------------
  # For each (period t, type value v): run rd_period on the type indicator
  # index: (t-1)*n_types + v_rank
  type_rank <- stats::setNames(seq_along(all_type_values), all_type_values)
  fits_by_pt <- vector("list", P * n_types)   # row-major: period varies fast
  dim(fits_by_pt) <- c(n_types, P)
  dimnames(fits_by_pt) <- list(type  = as.character(all_type_values),
                                period = plab)

  theta <- numeric(P * n_types)   # jump estimates (stacked)
  idx   <- 0L

  for (ki in seq_along(plab)) {
    df_k   <- pt[[plab[ki]]]
    for (vi in seq_along(all_type_values)) {
      idx <- idx + 1L
      v   <- all_type_values[vi]
      y_v <- as.numeric(df_k$type == v)
      # Per-cell bandwidth: CCT when h = NULL and bwselect = "cct"; otherwise
      # h is non-NULL (explicit or pre-set from h_rot for bwselect = "rot").
      bw     <- .cell_bandwidth(y_v, df_k$R, c, kernel, h, bwselect)
      h_cell <- bw[["h"]]
      b_cell <- bw[["b"]]
      fit <- tryCatch(
        rd_period(y = y_v, x = df_k$R, h = h_cell, b = b_cell, id = df_k$id,
                  c = c, p = 1L, q = 2L, kernel = kernel),
        error = function(e) NULL
      )
      fits_by_pt[vi, ki] <- list(fit)  # use [ to allow NULL without error
      theta[idx] <- if (is.null(fit)) NA_real_ else if (bc) fit$D_bc else fit$D
    }
  }

  # Build the covariance matrix Sigma (P*n_types x P*n_types)
  N <- P * n_types
  Sigma <- matrix(0, N, N)

  for (a_t in seq_along(plab)) {
    for (a_v in seq_along(all_type_values)) {
      row_a <- (a_t - 1L) * n_types + a_v
      fit_a <- fits_by_pt[[a_v, a_t]]
      if (is.null(fit_a)) next

      for (b_t in seq_along(plab)) {
        for (b_v in seq_along(all_type_values)) {
          row_b <- (b_t - 1L) * n_types + b_v
          fit_b <- fits_by_pt[[b_v, b_t]]
          if (is.null(fit_b)) next

          if (a_t == b_t) {
            # Within-period covariance: units share the same running variable,
            # so the same unit contributes to both type-indicator RDs on the
            # same side — the same-side ("pc") component, scheme-independent.
            Sigma[row_a, row_b] <- .cross_cov(fit_a, fit_b, bc = bc)$pc
          } else {
            # Cross-period covariance: scheme-dependent.
            Sigma[row_a, row_b] <- .cov_scheme(fit_a, fit_b, use_scheme, bc = bc)
          }
        }
      }
    }
  }

  # The n_types type-indicator jumps sum to zero within each period (the
  # indicators sum to 1), so the full P*n_types covariance is rank-deficient by
  # construction. Testing all of them through a pseudo-inverse is numerically
  # fragile — a structural-zero singular value is kept on some LAPACK builds and
  # its 1/sv inflates the statistic (platform-dependent p-values). Instead drop
  # one (reference) type per period: with k present types we keep k-1, an
  # equivalent full-rank test (the dropped jump is minus the sum of the rest).
  keep <- logical(length(theta))
  for (ki in seq_along(plab)) {
    rows_k  <- (ki - 1L) * n_types + seq_len(n_types)
    present <- rows_k[!is.na(theta[rows_k])]
    if (length(present) >= 2L) keep[present[-length(present)]] <- TRUE
  }
  ok_idx <- which(keep)
  if (length(ok_idx) == 0L) {
    ll_result <- list(stat = 0, df = 0L, p = 1)
  } else {
    ll_result <- .joint_wald(theta[ok_idx], Sigma[ok_idx, ok_idx, drop = FALSE])
  }

  # Per-period LL-Wald: restrict the kept (full-rank) contrasts to each period's
  # own rows. With binary types this is the single type-share jump in that period
  # => chi^2(1); these are the components the joint test aggregates (the joint is
  # not their sum, since it also carries the cross-period covariance).
  per_period_wald <- stats::setNames(vector("list", length(plab)), plab)
  for (ki in seq_along(plab)) {
    rows_k <- (ki - 1L) * n_types + seq_len(n_types)
    keep_k <- rows_k[keep[rows_k]]
    per_period_wald[[ki]] <- if (length(keep_k) == 0L)
      list(stat = 0, df = 0L, p = 1)
    else
      .joint_wald(theta[keep_k], Sigma[keep_k, keep_k, drop = FALSE])
  }

  # ----- (2) Canay-Kamat permutation ----------------------------------------
  # Joint test over all (period, type) pairs: observed statistic = sum of
  # absolute mean differences across all cells.
  # The joint permutation: for each period independently, permute the side
  # label among the 2q nearest observations (the canonical CK construction).
  # We sum the observed stats and compare to the sum under permutation.

  ck_obs_total <- 0
  # ck_period_parts: one entry per period, grouping all type cells together.
  # Each entry stores the pooled indicator matrix (n_near × n_active_types) and
  # the number of right-side units (nr), so the permutation loop can draw ONE
  # unit-level index per period and apply it identically to every type column.
  ck_period_parts <- vector("list", length(plab))
  ck_obs_by_period <- stats::setNames(rep(NA_real_, length(plab)), plab)  # per-period observed CK stat
  q_used <- stats::setNames(integer(length(plab)), plab)  # per-period q actually used

  for (ki in seq_along(plab)) {
    df_k <- pt[[plab[ki]]]
    # Number of nearest observations per side: the Canay-Kamat rule of thumb
    # (per period) when q is NULL, otherwise the user-supplied fixed q.
    q_k <- if (is.null(q)) .q_rot(df_k$R, df_k$type, c) else as.integer(q)
    q_used[ki] <- q_k
    # Work with the q_k nearest on each side (pool all types, then split)
    pos <- df_k$R >= c
    neg <- df_k$R <  c
    rt  <- order(df_k$R[pos])[seq_len(min(q_k, sum(pos)))]
    lt  <- order(-(df_k$R[neg]))[seq_len(min(q_k, sum(neg)))]

    R_near  <- c(df_k$R[pos][rt],  df_k$R[neg][lt])
    id_near <- c(df_k$id[pos][rt], df_k$id[neg][lt])
    side_near <- c(rep(1L, length(rt)), rep(0L, length(lt)))
    type_near <- df_k$type[match(id_near, df_k$id)]

    nr   <- sum(side_near == 1L)
    npool <- length(side_near)

    # Build indicator matrix: one column per active type value
    g_mat_list <- list()
    obs_vec    <- numeric(0)
    obs_period <- 0
    for (vi in seq_along(all_type_values)) {
      v      <- all_type_values[vi]
      g_near <- as.numeric(type_near == v)
      gr     <- g_near[side_near == 1L]
      gl     <- g_near[side_near == 0L]
      if (length(gr) == 0L || length(gl) == 0L) next
      obs_cell <- abs(mean(gr) - mean(gl))
      ck_obs_total <- ck_obs_total + obs_cell
      obs_period   <- obs_period + obs_cell
      g_mat_list[[length(g_mat_list) + 1L]] <- g_near
      obs_vec <- c(obs_vec, obs_cell)
    }

    if (length(g_mat_list) > 0L) {
      # Matrix: npool rows × n_active_types columns
      g_mat <- matrix(unlist(g_mat_list), nrow = npool, ncol = length(g_mat_list))
      ck_period_parts[[ki]] <- list(g_mat = g_mat, nr = nr, npool = npool, period = plab[ki])
      ck_obs_by_period[ki]  <- obs_period
    }
  }
  # Drop NULL entries (periods with no active type cells)
  ck_period_parts <- Filter(Negate(is.null), ck_period_parts)

  # Generate the null distribution: for each replication, draw ONE unit-level
  # permutation per period (over the 2q nearest units) and apply it identically
  # to all type columns for that period.  This is the canonical CK construction:
  # within each period the type cells share the SAME permuted side assignment.
  # NB: do NOT reseed here. The permutation must inherit whatever RNG state the
  # caller set (e.g. set.seed() before the call) so the CK p-value is
  # reproducible; an internal set.seed(NULL) would re-init from entropy and make
  # every call irreproducible.
  # Draw per-period contributions in the SAME order as before; the joint null is
  # their colSums, so the joint CK p-value is byte-identical to the original
  # single-vector loop, while per-period nulls are read off the same draws.
  ck_p_by_period <- stats::setNames(rep(NA_real_, length(plab)), plab)
  if (length(ck_period_parts) == 0L) {
    ck_perm_dist <- rep(0, S)
    ck_p <- 1
  } else {
    ck_perm_mat <- replicate(S, {
      vapply(ck_period_parts, function(part) {
        g_mat <- part$g_mat; nr <- part$nr; npool <- part$npool
        # One shared draw for all type columns in this period
        idx <- sample.int(npool, nr)
        # Contribution = sum over columns of |mean(pool[idx,]) - mean(pool[-idx,])|
        perm_right <- colMeans(g_mat[idx,  , drop = FALSE])
        perm_left  <- colMeans(g_mat[-idx, , drop = FALSE])
        sum(abs(perm_right - perm_left))
      }, numeric(1))
    })
    if (is.null(dim(ck_perm_mat)))
      ck_perm_mat <- matrix(ck_perm_mat, nrow = length(ck_period_parts))
    ck_perm_dist <- colSums(ck_perm_mat)                    # joint null distribution
    ck_p <- (1 + sum(ck_perm_dist >= ck_obs_total)) / (S + 1)
    for (j in seq_along(ck_period_parts)) {
      pj <- ck_period_parts[[j]]$period
      ck_p_by_period[pj] <- (1 + sum(ck_perm_mat[j, ] >= ck_obs_by_period[pj])) / (S + 1)
    }
  }

  # ----- (3) McCrary within-type -------------------------------------------
  mcc_within_rows <- list()
  for (ki in seq_along(plab)) {
    df_k <- pt[[plab[ki]]]
    p_vec <- numeric(n_types)
    for (vi in seq_along(all_type_values)) {
      v    <- all_type_values[vi]
      x_v  <- df_k$R[df_k$type == v]
      # Centre x_v at cutoff (already: cutoff = c, so subtract c)
      p_vec[vi] <- .mccrary(x_v - c, h_aux)
    }
    for (vi in seq_along(all_type_values)) {
      mcc_within_rows[[length(mcc_within_rows) + 1L]] <- data.frame(
        period = plab[ki],
        type   = as.character(all_type_values[vi]),
        p_raw  = p_vec[vi],
        stringsAsFactors = FALSE
      )
    }
    # Attach bonferroni to last row of this period — will be added post-hoc
    # Actually store p_bonf per period separately below
  }
  mcc_within <- do.call(rbind, mcc_within_rows)
  rownames(mcc_within) <- NULL

  # Add per-period Bonferroni column
  mcc_within$p_bonf <- NA_real_
  for (ki in plab) {
    rows_k <- mcc_within$period == ki
    p_k    <- mcc_within$p_raw[rows_k]
    p_min_k <- if (all(is.na(p_k))) NA_real_ else min(p_k, na.rm = TRUE)
    bonf_k  <- if (is.na(p_min_k)) NA_real_ else pmin(1, p_min_k * n_types)
    mcc_within$p_bonf[rows_k] <- bonf_k
  }

  # ----- (4) McCrary pooled ------------------------------------------------
  mcc_pooled <- data.frame(
    period = plab,
    p      = vapply(seq_along(plab), function(ki) {
      df_k <- pt[[plab[ki]]]
      .mccrary(df_k$R - c, h_aux)
    }, numeric(1L)),
    stringsAsFactors = FALSE
  )

  # Per-period components (the building blocks behind the joint tests): each
  # period's own LL-Wald (chi^2 with its kept contrasts) and Canay-Kamat p.
  per_period <- stats::setNames(lapply(seq_along(plab), function(ki) {
    list(ll_wald = per_period_wald[[ki]], ck_p = unname(ck_p_by_period[plab[ki]]))
  }), plab)

  # ----- assemble output ---------------------------------------------------
  structure(
    list(
      ll_wald        = ll_result,
      per_period     = per_period,
      ck_perm        = list(stat = ck_obs_total, p = ck_p),
      mccrary_within = mcc_within,
      mccrary_pooled = mcc_pooled,
      meta = list(
        periods      = plab,
        type_values  = all_type_values,
        h            = if (!is.null(h)) h else NA_real_,
        bwselect     = bwselect,
        q            = if (is.null(q)) "rot" else as.integer(q),
        q_used       = q_used,
        S            = S,
        scheme       = use_scheme,
        bc           = bc
      )
    ),
    class = "rd_typecont"
  )
}


#' @export
print.rd_typecont <- function(x, ...) {
  cat("Type-continuity tests (Assumption A7)\n")
  h_str <- if (is.na(x$meta$h)) paste0("per-cell ", toupper(x$meta$bwselect)) else sprintf("%.4g", x$meta$h)
  cat(sprintf("  Periods: %s   Types: %s   h=%s   bwselect=%s   scheme=%s\n\n",
              paste(x$meta$periods, collapse = ", "),
              paste(x$meta$type_values, collapse = ", "),
              h_str, x$meta$bwselect, toupper(x$meta$scheme)))

  cat("(1) LL-Wald [nec & suff]:\n")
  cat(sprintf("    chi2(%.0f) = %.4f   p = %.4f\n\n",
              x$ll_wald$df, x$ll_wald$stat, x$ll_wald$p))

  cat("(2) Canay-Kamat permutation [nec & suff]:\n")
  cat(sprintf("    stat = %.4f   p = %.4f\n\n",
              x$ck_perm$stat, x$ck_perm$p))

  cat("(3) McCrary within-type [suff, NOT nec]  (Bonferroni per period):\n")
  for (k in unique(x$mccrary_within$period)) {
    sub <- x$mccrary_within[x$mccrary_within$period == k, ]
    cat(sprintf("    period %s  p_bonf=%.4f  (by type: %s)\n",
                k, sub$p_bonf[1L],
                paste(sprintf("t%s=%.3f", sub$type, sub$p_raw), collapse = ", ")))
  }
  cat("\n")

  cat("(4) McCrary pooled [NEITHER suff NOR nec]:\n")
  for (i in seq_len(nrow(x$mccrary_pooled))) {
    cat(sprintf("    period %s  p = %.4f\n",
                x$mccrary_pooled$period[i], x$mccrary_pooled$p[i]))
  }
  invisible(x)
}
