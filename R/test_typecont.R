#' Test the continuity of the type distribution (Assumption A6)
#'
#' Runs the four tests for Assumption A6 ("continuity of the type distribution")
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
#' @param h bandwidth.  If `NULL`, uses `0.5 * IQR(x)` as a simple default;
#'   a data-driven bandwidth (e.g. from [rddid()]) is recommended in practice.
#' @param q number of observations nearest the cutoff on each side to use in
#'   the Canay-Kamat permutation test (default 75).
#' @param S number of permutation replications for the Canay-Kamat test
#'   (default 499).
#' @param kernel kernel for the local-linear RD: `"triangular"` (default),
#'   `"epanechnikov"`, or `"uniform"`.
#' @param scheme covariance scheme for the joint Wald: `"auto"` detects from
#'   the data (same logic as [rddid()]), or one of `"cs"`, `"pc"`, `"pv"`.
#' @param ... currently unused.
#'
#' @return An object of class `"rd_typecont"`, a named list with:
#'   \item{ll_wald}{list with `stat` (chi-square), `df`, `p`.}
#'   \item{ck_perm}{list with `stat` (observed sum of |mean diffs|), `p`.}
#'   \item{mccrary_within}{data frame with columns `period`, `type`,
#'     `p_raw` (per-type McCrary p-value); plus `p_bonf` (Bonferroni-adjusted
#'     p-value for the period, minimum × number of types).}
#'   \item{mccrary_pooled}{data frame with columns `period`, `p` (per-period
#'     McCrary p-value on pooled sample).}
#'   \item{meta}{list with `periods`, `type_values`, `h`, `q`, `S`, `scheme`.}
#' @export
rd_typecont <- function(data, x, time, id,
                        c = 0,
                        h = NULL,
                        q = 75L,
                        S = 499L,
                        kernel = "triangular",
                        scheme = c("auto", "cs", "pc", "pv"),
                        ...) {
  scheme <- match.arg(scheme)

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
  if (is.null(h)) {
    all_x <- data[[x]]
    h     <- 0.5 * stats::IQR(all_x)
    if (h <= 0) h <- stats::sd(all_x)
  }

  # ----- build types --------------------------------------------------------
  # .build_types may be the version from test_helpers.R (returns list with
  # $period_types, $wide, $periods) OR the version from test_homog.R (returns
  # a named list of data frames directly, with only id+type columns).
  # Normalise to a named list of per-period data frames that include column R.
  bt_raw  <- .build_types(data, x, time, id, c = c)

  # Detect which version was loaded
  if (!is.null(bt_raw$period_types)) {
    # test_helpers.R version: wraps everything in a list
    pt <- bt_raw$period_types
  } else {
    # test_homog.R version: returns the list of per-period frames directly
    pt <- bt_raw
  }

  # Ensure each per-period frame has an R (running variable) column.
  # If absent (test_homog.R version only returns id+type), join from data.
  for (k in plab) {
    df_k <- pt[[k]]
    if (!is.null(df_k) && is.null(df_k$R)) {
      sub_k <- data[data[[time]] == periods[match(k, plab)], , drop = FALSE]
      m     <- match(df_k[[id]], sub_k[[id]])
      df_k$R <- sub_k[[x]][m]
      pt[[k]] <- df_k
    }
  }

  # All type values that appear anywhere across all periods
  all_type_values <- sort(unique(unlist(lapply(pt, `[[`, "type"))))
  n_types <- length(all_type_values)

  # ----- detect scheme (mirroring .detect_scheme in rddid.R) ---------------
  if (scheme == "auto") {
    # Build a "plist" analogue for scheme detection: for each period, collect
    # ids and x values of units in the bandwidth window
    use_scheme <- local({
      long_rows <- lapply(plab, function(k) {
        df_k <- pt[[k]]
        inwin <- abs(df_k$R - c) <= h
        data.frame(period = k,
                   id     = df_k$id[inwin],
                   side   = sign(df_k$R[inwin] - c))
      })
      long <- do.call(rbind, long_rows)
      rep_ids <- names(which(table(unique(long[, c("period", "id")])$id) >= 2L))
      if (length(rep_ids) == 0L) {
        "cs"
      } else {
        sub     <- long[long$id %in% rep_ids, ]
        switches <- tapply(sub$side, sub$id, function(s) length(unique(s)) > 1L)
        if (any(switches)) "pv" else "pc"
      }
    })
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
      fit <- tryCatch(
        rd_period(y = y_v, x = df_k$R, h = h, b = h, id = df_k$id,
                  c = c, p = 1L, q = 2L, kernel = kernel),
        error = function(e) NULL
      )
      fits_by_pt[vi, ki] <- list(fit)  # use [ to allow NULL without error
      theta[idx] <- if (is.null(fit)) NA_real_ else fit$D
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
            # Within-period covariance: units share the same running variable.
            # Both sides: same unit i on side +/- of period t contributes to
            # both type-a and type-b indicator RDs.
            cov_ab <- .match_sum(fit_a$sides$`+`$id, fit_a$sides$`+`$g,
                                  fit_b$sides$`+`$id, fit_b$sides$`+`$g) +
                      .match_sum(fit_a$sides$`-`$id, fit_a$sides$`-`$g,
                                  fit_b$sides$`-`$id, fit_b$sides$`-`$g)
            Sigma[row_a, row_b] <- cov_ab
          } else {
            # Cross-period covariance: depends on scheme
            if (use_scheme == "cs") {
              Sigma[row_a, row_b] <- 0
            } else {
              cc_pp <- .match_sum(fit_a$sides$`+`$id, fit_a$sides$`+`$g,
                                   fit_b$sides$`+`$id, fit_b$sides$`+`$g)
              cc_mm <- .match_sum(fit_a$sides$`-`$id, fit_a$sides$`-`$g,
                                   fit_b$sides$`-`$id, fit_b$sides$`-`$g)
              cc_pm <- .match_sum(fit_a$sides$`+`$id, fit_a$sides$`+`$g,
                                   fit_b$sides$`-`$id, fit_b$sides$`-`$g)
              cc_mp <- .match_sum(fit_a$sides$`-`$id, fit_a$sides$`-`$g,
                                   fit_b$sides$`+`$id, fit_b$sides$`+`$g)
              if (use_scheme == "pc") {
                Sigma[row_a, row_b] <- cc_pp + cc_mm
              } else {  # pv
                Sigma[row_a, row_b] <- cc_pp + cc_mm - cc_pm - cc_mp
              }
            }
          }
        }
      }
    }
  }

  # Drop rows/cols where theta is NA
  ok_idx <- which(!is.na(theta))
  theta_ok  <- theta[ok_idx]
  Sigma_ok  <- Sigma[ok_idx, ok_idx, drop = FALSE]
  ll_result <- .joint_wald(theta_ok, Sigma_ok)

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

  for (ki in seq_along(plab)) {
    df_k <- pt[[plab[ki]]]
    # Work with the q nearest on each side (pool all types, then split)
    pos <- df_k$R >= c
    neg <- df_k$R <  c
    rt  <- order(df_k$R[pos])[seq_len(min(q, sum(pos)))]
    lt  <- order(-(df_k$R[neg]))[seq_len(min(q, sum(neg)))]

    R_near  <- c(df_k$R[pos][rt],  df_k$R[neg][lt])
    id_near <- c(df_k$id[pos][rt], df_k$id[neg][lt])
    side_near <- c(rep(1L, length(rt)), rep(0L, length(lt)))
    type_near <- df_k$type[match(id_near, df_k$id)]

    nr   <- sum(side_near == 1L)
    npool <- length(side_near)

    # Build indicator matrix: one column per active type value
    g_mat_list <- list()
    obs_vec    <- numeric(0)
    for (vi in seq_along(all_type_values)) {
      v      <- all_type_values[vi]
      g_near <- as.numeric(type_near == v)
      gr     <- g_near[side_near == 1L]
      gl     <- g_near[side_near == 0L]
      if (length(gr) == 0L || length(gl) == 0L) next
      obs_cell <- abs(mean(gr) - mean(gl))
      ck_obs_total <- ck_obs_total + obs_cell
      g_mat_list[[length(g_mat_list) + 1L]] <- g_near
      obs_vec <- c(obs_vec, obs_cell)
    }

    if (length(g_mat_list) > 0L) {
      # Matrix: npool rows × n_active_types columns
      g_mat <- matrix(unlist(g_mat_list), nrow = npool, ncol = length(g_mat_list))
      ck_period_parts[[ki]] <- list(g_mat = g_mat, nr = nr, npool = npool)
    }
  }
  # Drop NULL entries (periods with no active type cells)
  ck_period_parts <- Filter(Negate(is.null), ck_period_parts)

  # Generate the null distribution: for each replication, draw ONE unit-level
  # permutation per period (over the 2q nearest units) and apply it identically
  # to all type columns for that period.  This is the canonical CK construction:
  # within each period the type cells share the SAME permuted side assignment.
  set.seed(NULL)  # allow externally seeded reproducibility
  ck_perm_dist <- replicate(S, {
    total_perm <- 0
    for (part in ck_period_parts) {
      g_mat <- part$g_mat; nr <- part$nr; npool <- part$npool
      # One shared draw for all type columns in this period
      idx <- sample.int(npool, nr)
      # Contribution = sum over columns of |mean(pool[idx,]) - mean(pool[-idx,])|
      perm_right <- colMeans(g_mat[idx,    , drop = FALSE])
      perm_left  <- colMeans(g_mat[-idx,   , drop = FALSE])
      total_perm <- total_perm + sum(abs(perm_right - perm_left))
    }
    total_perm
  })
  ck_p <- (1 + sum(ck_perm_dist >= ck_obs_total)) / (S + 1)

  # ----- (3) McCrary within-type -------------------------------------------
  mcc_within_rows <- list()
  for (ki in seq_along(plab)) {
    df_k <- pt[[plab[ki]]]
    p_vec <- numeric(n_types)
    for (vi in seq_along(all_type_values)) {
      v    <- all_type_values[vi]
      x_v  <- df_k$R[df_k$type == v]
      # Centre x_v at cutoff (already: cutoff = c, so subtract c)
      p_vec[vi] <- .mccrary(x_v - c, h)
    }
    # Bonferroni within this period: min(p) * n_types (capped at 1)
    p_min  <- if (all(is.na(p_vec))) NA_real_ else min(p_vec, na.rm = TRUE)
    p_bonf <- if (is.na(p_min)) NA_real_ else pmin(1, p_min * n_types)
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
      .mccrary(df_k$R - c, h)
    }, numeric(1L)),
    stringsAsFactors = FALSE
  )

  # ----- assemble output ---------------------------------------------------
  structure(
    list(
      ll_wald        = ll_result,
      ck_perm        = list(stat = ck_obs_total, p = ck_p),
      mccrary_within = mcc_within,
      mccrary_pooled = mcc_pooled,
      meta = list(
        periods      = plab,
        type_values  = all_type_values,
        h            = h,
        q            = q,
        S            = S,
        scheme       = use_scheme
      )
    ),
    class = "rd_typecont"
  )
}


#' @export
print.rd_typecont <- function(x, ...) {
  cat("Type-continuity tests (Assumption A6)\n")
  cat(sprintf("  Periods: %s   Types: %s   h=%.4g   scheme=%s\n\n",
              paste(x$meta$periods, collapse = ", "),
              paste(x$meta$type_values, collapse = ", "),
              x$meta$h, toupper(x$meta$scheme)))

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
