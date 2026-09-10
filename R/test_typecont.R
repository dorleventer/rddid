#' Test the continuity of the type distribution
#'
#' Wald test of the continuous-type-distribution assumption (Section 4.4 of
#' Leventer and Nevo). The "type" of unit
#' \eqn{i} in period \eqn{t} is the sign pattern of its running variables in
#' the OTHER periods, \eqn{\mathbf{V}_{i,-t} = (1\{R_{i,s} \ge c\})_{s \ne t}}.
#'
#' For each period \eqn{t} and each type value \eqn{v}, estimate the
#' local-linear RD jump of the type indicator \eqn{1\{\mathbf{V}_{i,-t} = v\}}
#' on the running variable \eqn{R_{i,t}}.  The jump
#' \eqn{\hat\pi_{t,(+)}(v) - \hat\pi_{t,(-)}(v)} is the output of [rd_period()]
#' with a binary outcome.  The joint Wald statistic across all (period, type)
#' pairs drops one reference type per period, because within each period the
#' type indicators sum to 1 (the full block is singular); df = number of kept
#' contrasts.  The covariance is built from the per-unit influence vectors
#' returned by [rd_period()], using the same within-period and cross-period
#' id-matching as the main estimator, scheme-aware.
#'
#' @param data a long data frame, one row per unit-period. A unit's type in
#'   period \eqn{t} is read from its running variable in the other period(s);
#'   units unobserved there are dropped from period \eqn{t}, so the panel need
#'   not be balanced.
#' @param x column name (string) for the running variable.
#' @param time column name (string) for the period.
#' @param id column name (string) for the unit identifier.
#' @param estimand `"att"` (default) or `"atu"`. Label only: the test is
#'   identical under either estimand, because the continuous-type-distribution
#'   assumption is symmetric in the two sides of the cutoff.
#' @param c cutoff (default 0).
#' @param h bandwidth.  If `NULL`, the bandwidth is determined by `bwselect`;
#'   an explicit numeric value overrides `bwselect` and is used directly.
#' @param bwselect bandwidth selection rule when `h = NULL`: `"cct"` (default)
#'   computes a per-cell CCT MSE-optimal bandwidth via [rd_bw_cct()] for each
#'   (period, type) RD; `"rot"` uses the `0.5 * IQR(x)` rule of thumb applied
#'   to the full sample (the previous default behaviour).  Ignored when `h` is
#'   supplied explicitly.
#' @param kernel kernel for the local-linear RD: `"triangular"` (default),
#'   `"epanechnikov"`, or `"uniform"`.
#' @param scheme covariance scheme for the joint Wald: `"auto"` detects from
#'   the data (same logic as [rddid()]), or one of `"cs"`, `"pc"`, `"pv"`.
#' @param bc use robust bias-corrected jumps and variances in the LL-Wald
#'   (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test
#'   with the bias-corrected [rddid()] estimator; `FALSE` uses the conventional
#'   local-linear jumps and variances.
#' @param ... currently unused.
#'
#' @return An object of class `"rd_typecont"`, a named list with:
#'   \item{ll_wald}{list with `stat` (chi-square), `df`, `p`.}
#'   \item{per_period}{named list (by period) of the per-period components the
#'     joint test aggregates; each entry has `ll_wald` (that period's own
#'     LL-Wald `stat`/`df`/`p`, restricted to its kept contrasts).}
#'   \item{meta}{list with `periods`, `type_values`, `h` (NA when
#'     `bwselect = "cct"`), `bwselect`, `scheme`, `bc`, `estimand`.}
#' @export
rd_typecont <- function(data, x, time, id,
                        estimand = c("att", "atu"),
                        c = 0,
                        h = NULL,
                        bwselect = c("cct", "rot"),
                        kernel = "triangular",
                        scheme = c("auto", "cs", "pc", "pv"),
                        bc = TRUE,
                        ...) {
  scheme   <- match.arg(scheme)
  bwselect <- match.arg(bwselect)
  estimand <- match.arg(estimand)

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
  # When bwselect = "rot" and h = NULL, h is set to the IQR-based pilot
  # 0.5*IQR(x) (current behaviour).
  # When bwselect = "cct" and h = NULL, h stays NULL; per-cell CCT is computed
  # inside the LL-Wald loop below.
  if (is.null(h) && bwselect == "rot") {
    all_x <- data[[x]]
    h     <- 0.5 * stats::IQR(all_x)
    if (h <= 0) h <- stats::sd(all_x)
  }

  # ----- build types --------------------------------------------------------
  # Per-period frames (id, R, type) from the shared canonical builder in
  # R/test_helpers.R.  `type` is the "+"/"-" sign-pattern string of the other
  # periods; units at the cutoff are treated as above it.
  pt <- .build_types(data, x, time, id, c = c)$period_types

  # All type values that appear anywhere across all periods
  all_type_values <- sort(unique(unlist(lapply(pt, `[[`, "type"))), method = "radix")  # locale-independent
  n_types <- length(all_type_values)

  # ----- detect scheme -----------------------------------------------------
  # Classify from every unit with a defined type in each period via the shared
  # primitive (side = 1{R >= c}, treated at the cutoff). The scheme is a
  # property of the design (repeated ids, side switching), not of a window;
  # restricting to a rule-of-thumb window under-detected switching whenever the
  # per-cell CCT bandwidths reached beyond it.
  if (scheme == "auto") {
    long <- do.call(rbind, lapply(plab, function(k) {
      df_k  <- pt[[k]]
      data.frame(period = k,
                 id     = df_k$id,
                 side   = as.integer(df_k$R >= c))
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
      # h is non-NULL (explicit or pre-set from 0.5*IQR for bwselect = "rot").
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
  # equivalent full-rank test at a common bandwidth (the dropped jump is minus
  # the sum of the rest; with per-cell CCT bandwidths the jumps do not sum
  # exactly to zero, so it is then a different, still valid, test). Types are
  # in radix order, so the dropped type is the all-below pattern and the kept
  # contrast at P = 2 is the paper's jump for the indicator 1{V_is = 1}.
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

  # Per-period components (the building blocks behind the joint test): each
  # period's own LL-Wald (chi^2 with its kept contrasts).
  per_period <- stats::setNames(lapply(seq_along(plab), function(ki) {
    list(ll_wald = per_period_wald[[ki]])
  }), plab)

  # ----- assemble output ---------------------------------------------------
  structure(
    list(
      ll_wald    = ll_result,
      per_period = per_period,
      meta = list(
        periods      = plab,
        type_values  = all_type_values,
        h            = if (!is.null(h)) h else NA_real_,
        bwselect     = bwselect,
        scheme       = use_scheme,
        bc           = bc,
        estimand     = estimand
      )
    ),
    class = "rd_typecont"
  )
}


#' @export
print.rd_typecont <- function(x, ...) {
  cat("Type-continuity test\n")
  h_str <- if (is.na(x$meta$h)) paste0("per-cell ", toupper(x$meta$bwselect)) else sprintf("%.4g", x$meta$h)
  cat(sprintf("  Periods: %s   Types: %s   h=%s   bwselect=%s   scheme=%s\n\n",
              paste(x$meta$periods, collapse = ", "),
              paste(x$meta$type_values, collapse = ", "),
              h_str, x$meta$bwselect, x$meta$scheme))
  est <- if (is.null(x$meta$estimand)) "att" else x$meta$estimand
  if (est == "atu")
    cat("  estimand: atu (test is unchanged; see ?rd_typecont)\n")

  cat("LL-Wald:\n")
  cat(sprintf("    chi2(%.0f) = %.4f   p = %.4f\n",
              x$ll_wald$df, x$ll_wald$stat, x$ll_wald$p))
  invisible(x)
}
