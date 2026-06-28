# Test of Assumption A9: type-homogeneous confounding.
# Manuscript ref: Leventer and Nevo, "Correcting Invalid RD Designs",
# paragraph "Type-homogeneous confounding (Assumption A9)".
#
# The test is run in COMPARISON PERIODS ONLY.  In comparison period t0 the
# outcome RD jump equals the pure confounding:
#   D_{t0}(v_{-t0}) = alpha_{t0,0}(v_{-t0}).
# We estimate that jump separately per type v_{-t0} using rd_period(), then
# form a Wald test that the jumps are equal across types.
#
# IMPORTANT: this test is NEITHER necessary NOR sufficient for homogeneity at
# the RD period; see documentation for rd_homog().

# ---------------------------------------------------------------------------
# This test has no local helpers: types come from the shared .build_types()
# (R/test_helpers.R) and cross-period covariances from .cov_scheme()/.cross_cov()
# (R/aggregate.R).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Test of type-homogeneous confounding (Assumption A9)
#'
#' Tests whether the outcome RD discontinuity is constant across types in
#' **comparison periods** (periods where the treatment of interest is absent).
#' In a comparison period \eqn{t_0} the discontinuity equals the pure
#' confounding: \eqn{D_{t_0}(\mathbf{v}_{-t_0}) = \alpha_{t_0,0}(\mathbf{v}_{-t_0})}.
#' The function estimates the jump by type (sign pattern of the OTHER periods'
#' running variables) via [rd_period()], then forms a joint Wald test that the
#' jumps are equal across types and comparison periods.
#'
#' ## Scope and interpretation
#'
#' **This test is run in comparison periods only.**  The type-homogeneous
#' confounding assumption (A9 in Leventer and Nevo) is needed at the *RD
#' period* \eqn{t_{\mathrm{RD}}}, but there the jump also contains the ATT and
#' \eqn{\alpha_{t_{\mathrm{RD}},0}(\mathbf{v}_{-t_{\mathrm{RD}}})} is not
#' separately observable.  Testing homogeneity in comparison periods is
#' therefore only **suggestive** evidence: it is **neither necessary nor
#' sufficient** for the assumption to hold at \eqn{t_{\mathrm{RD}}}.
#'
#' ## Null hypothesis
#'
#' \eqn{H_0 :} the outcome RD jump is equal across all type cells
#' \eqn{\mathbf{v}_{-t_0}}, jointly across all comparison periods:
#' \eqn{D_{t_0}(\mathbf{v}) = D_{t_0}(\mathbf{v}') \; \forall \, \mathbf{v}, \mathbf{v}', t_0}.
#'
#' The Wald statistic is formed from the per-type jump contrasts
#' \eqn{\Delta = D_{t_0,v} - D_{t_0,\mathrm{ref}}}, stacked across comparison
#' periods, with covariance estimated from the [rd_period()] influence vectors.
#' Within a period, cross-type covariance is zero (types partition the sample).
#' Across periods, id-matched covariance is used under the detected/requested
#' sampling scheme (see `scheme`).
#'
#' @param data A long data frame, one row per unit-period.
#' @param y,x,time Column names (character strings) for the outcome, running
#'   variable, and period indicator.
#' @param id Column name for the unit identifier.  Required (the test uses the
#'   other periods' running variables to assign types, which needs a panel).
#' @param comparisons Values of `time` to use as comparison periods.  Defaults
#'   to all periods except `t_rd` (if `t_rd` is supplied) or all periods (if
#'   `t_rd = NULL`).
#' @param t_rd Value of `time` for the RD period.  Used only to exclude it from
#'   comparison periods when `comparisons = NULL`; the RD period is **not**
#'   used in the test itself.
#' @param c Cutoff value for the running variable (default 0).
#' @param h Main bandwidth.  If `NULL` (default), a simple rule-of-thumb
#'   \eqn{h = 0.2 \times \mathrm{range}(x)} is applied within each
#'   comparison period.  For reproducible results, supply `h` explicitly.
#' @param kernel Kernel for the local-linear RD: `"triangular"` (default),
#'   `"epanechnikov"`, or `"uniform"`.
#' @param scheme Sampling scheme for the cross-period covariance:
#'   `"auto"` (detect from the id/side structure, default), `"cs"` (repeated
#'   cross-section, no cross-period terms), `"pc"` (panel, time-constant
#'   running variable), or `"pv"` (panel, time-varying running variable).
#' @param min_n Minimum number of observations per type-side before that type
#'   is included.  Types with fewer than `min_n` obs on either side in a given
#'   period are silently dropped (default 10).
#' @param bc Use robust bias-corrected per-type jumps and variances
#'   (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test
#'   with the bias-corrected [rddid()] estimator; `FALSE` uses the conventional
#'   local-linear jumps and variances.
#' @param type_by How a unit's type is defined. `"pattern"` (default) uses the
#'   full multi-period sign pattern of the other periods (the general case).
#'   `"rd_side"` uses only the unit's side of the cutoff in the RD period `t_rd`,
#'   a binary partition; this is the relevant partition for a single joint test
#'   across comparison periods that share a running variable, where it yields one
#'   contrast per period (a \eqn{\chi^2(P)} test for `P` comparison periods).
#'   Requires `t_rd`.
#' @param ... Further arguments passed to [rd_period()] (e.g. `p`, `q`, `b`).
#'
#' @return An object of class `"rd_homog"`, a list with:
#'   \describe{
#'     \item{`statistic`}{Joint Wald chi-squared statistic.}
#'     \item{`df`}{Degrees of freedom (total contrasts used).}
#'     \item{`p_value`}{p-value from the chi-squared distribution.}
#'     \item{`period_type_jumps`}{Data frame with one row per (comparison
#'       period, type) cell: period, type, jump estimate, SE, number of
#'       observations, and whether it was the reference type.}
#'     \item{`contrasts`}{Named numeric vector of jump contrasts
#'       (non-reference minus reference), stacked across periods.}
#'     \item{`cov_matrix`}{Estimated covariance matrix of the contrasts.}
#'     \item{`scheme`}{Sampling scheme used.}
#'     \item{`comparisons`}{Comparison periods actually used.}
#'     \item{`call`}{The matched call.}
#'   }
#'
#' @note
#' **Necessary and sufficient status:** This test is *neither necessary nor
#' sufficient* for Assumption A9 (type-homogeneous confounding) to hold at the
#' RD period.  Homogeneity in comparison periods is only suggestive because the
#' assumption is needed at \eqn{t_{\mathrm{RD}}}, where the confounding jump is
#' not separately identified.
#'
#' @seealso [rd_period()], [rddid()]
#'
#' @references
#' Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
#' Designs." Working paper.
#'
#' Imbens, G. W. and Lemieux, T. (2008). Regression discontinuity designs:
#' A guide to practice. *Journal of Econometrics*, 142(2), 615-635.
#'
#' Calonico, S., Cattaneo, M. D., and Titiunik, R. (2014). Robust
#' nonparametric confidence intervals for regression-discontinuity designs.
#' *Econometrica*, 82(6), 2295-2326.
#'
#' @examples
#' \dontrun{
#' # Two-period panel: period 1 is RD, period 2 is comparison.
#' # Type in period 2 = sign of R_{i,1} (the running variable of the other period).
#' set.seed(1)
#' n <- 500
#' r1 <- runif(n, -1, 1)
#' r2 <- runif(n, -1, 1)
#' # homogeneous confounding in comparison period (period 2)
#' y2 <- 0.3 * r2 + 0.4 * (r2 >= 0) + rnorm(n, 0, 0.3)
#' dat <- data.frame(
#'   id   = rep(seq_len(n), 2),
#'   time = rep(1:2, each = n),
#'   x    = c(r1, r2),
#'   y    = c(rnorm(n), y2)
#' )
#' rd_homog(dat, y = "y", x = "x", time = "time", id = "id",
#'          comparisons = 2, t_rd = 1, h = 0.3)
#' }
#' @export
rd_homog <- function(data, y, x, time, id,
                     comparisons = NULL, t_rd = NULL,
                     c = 0, h = NULL,
                     kernel = "triangular",
                     scheme = c("auto", "cs", "pc", "pv"),
                     min_n = 10L,
                     bc = TRUE,
                     type_by = c("pattern", "rd_side"),
                     ...) {
  cl      <- match.call()
  scheme  <- match.arg(scheme)
  type_by <- match.arg(type_by)

  # ---- validate columns ----
  for (nm in base::c(y, x, time, id))
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")

  times_all <- sort(unique(data[[time]]))

  # ---- comparison periods ----
  if (is.null(comparisons)) {
    comparisons <- if (!is.null(t_rd)) setdiff(times_all, t_rd) else times_all
  }
  if (length(comparisons) < 1L)
    stop("need at least one comparison period.")

  # ---- build types ----
  # shared canonical builder; types are "+"/"-" sign-pattern strings, units at
  # the cutoff treated as above it, partially-observed units dropped per period.
  bt        <- .build_types(data, x = x, time = time, id = id, c = c)
  type_list <- bt$period_types

  # type_by = "rd_side": override each comparison period's type with the unit's
  # side of the cutoff in the RD period t_rd (a binary partition), instead of the
  # full multi-period sign pattern. This is the relevant partition for a joint
  # homogeneity test across comparison periods when the comparison periods share
  # a running variable (so the other comparison period's side is collinear with
  # the running variable and carries no extra type information). Reduces each
  # period to one contrast, so P comparison periods give a chi^2(P) joint test.
  if (type_by == "rd_side") {
    if (is.null(t_rd))
      stop("type_by = \"rd_side\" requires `t_rd` (the RD period whose side defines the type).")
    sidecol <- paste0("side_", t_rd)
    if (!sidecol %in% names(bt$wide))
      stop("RD period '", t_rd, "' has no side column; cannot define rd_side types.")
    side_map <- stats::setNames(bt$wide[[sidecol]], as.character(bt$wide$id))
    for (tp in names(type_list)) {
      df0      <- type_list[[tp]]
      df0$type <- unname(side_map[as.character(df0$id)])
      type_list[[tp]] <- df0[!is.na(df0$type), , drop = FALSE]
    }
  }

  # ---- per-comparison-period, per-type fits ----
  # For each comparison period t0, for each type v:
  #   subset the data to period t0 AND units whose type label is v,
  #   run rd_period() on that subset.

  dots <- list(...)

  # Determine the sampling scheme for the cross-period covariance from the
  # id/side structure ACROSS the comparison periods, using the package-canonical
  # detector (shared with rddid() and rd_typecont()): units recurring across
  # periods with a constant side -> "pc", switching side -> "pv", none -> "cs".
  # (The earlier within-period heuristic could never observe cross-period
  # recurrence and silently collapsed to "cs", zeroing the joint test's
  # cross-period covariance.)
  comp_plist <- stats::setNames(lapply(as.character(comparisons), function(tp) {
    d_cp <- data[data[[time]] == tp, , drop = FALSE]
    list(id = d_cp[[id]], x = d_cp[[x]])
  }), as.character(comparisons))
  detected_scheme <- .detect_scheme(comp_plist, c = c)

  all_fits    <- list()   # key = "t0::type", value = rd_period object
  all_meta    <- list()   # key = "t0::type", value = (period, type, D, V_D, n)
  contrast_keys <- list() # per-period list of non-reference type keys

  for (tp in as.character(comparisons)) {
    # subset to this period
    rows  <- which(data[[time]] == tp)
    d_tp  <- data[rows, , drop = FALSE]

    # merge in type labels
    tdf   <- type_list[[tp]]   # data frame: id, type
    id_tp <- d_tp[[id]]
    m_idx <- match(id_tp, tdf$id)
    tvec  <- tdf$type[m_idx]   # type label per row

    valid_types <- sort(unique(tvec[!is.na(tvec)]))
    if (length(valid_types) < 2L) {
      # only one type (or no types): skip this period
      message("rd_homog: period ", tp, " has fewer than 2 types; skipping.")
      next
    }

    # reference type = first alphabetically
    ref_type <- valid_types[1L]
    for (vt in valid_types) {
      keep   <- !is.na(tvec) & tvec == vt
      y_vt   <- d_tp[[y]][keep]
      x_vt   <- d_tp[[x]][keep]
      id_vt  <- d_tp[[id]][keep]

      # default bandwidth
      bw_h   <- if (!is.null(h)) h else {
        rng <- diff(range(x_vt[is.finite(x_vt)], na.rm = TRUE))
        0.2 * rng
      }

      # skip if too few obs on either side
      n_pos <- sum(x_vt >= c, na.rm = TRUE)
      n_neg <- sum(x_vt <  c, na.rm = TRUE)
      if (n_pos < min_n || n_neg < min_n) next

      fit <- tryCatch(
        do.call(rd_period, base::c(list(y = y_vt, x = x_vt, h = bw_h,
                                        id = id_vt, c = c, kernel = kernel), dots)),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      key <- paste0(tp, "::", vt)
      all_fits[[key]] <- fit
      all_meta[[key]] <- list(period = tp, type = vt, ref = (vt == ref_type),
                               D = if (bc) fit$D_bc else fit$D,
                               V_D = if (bc) fit$V_D_bc else fit$V_D, n = fit$n)
    }

    # record contrast keys for this period (non-ref types that actually fitted)
    fitted_types_tp <- vapply(valid_types, function(vt) {
      paste0(tp, "::", vt) %in% names(all_fits)
    }, logical(1))
    fitted_types <- valid_types[fitted_types_tp]
    if (length(fitted_types) < 2L) next   # need at least ref + 1
    ref_fitted <- paste0(tp, "::", fitted_types[1L])
    if (!ref_fitted %in% names(all_fits)) next
    non_ref <- fitted_types[-1L]
    contrast_keys[[tp]] <- list(
      ref = ref_fitted,
      non_ref = paste0(tp, "::", non_ref)
    )
  }

  use_scheme <- if (scheme == "auto") detected_scheme else scheme

  # ---- build contrast vector and covariance matrix ----
  # stack non-ref types across periods:
  #   Delta[k] = D[non_ref[k]] - D[ref[k]]

  all_contrast_entries <- do.call(base::c, lapply(names(contrast_keys), function(tp) {
    contrast_keys[[tp]]$non_ref
  }))

  if (length(all_contrast_entries) == 0L)
    stop("rd_homog: no usable type contrasts found; check data, bandwidth, or min_n.")

  # map: for each contrast entry key, which ref key does it subtract?
  ref_map <- do.call(base::c, lapply(names(contrast_keys), function(tp) {
    stats::setNames(rep(contrast_keys[[tp]]$ref, length(contrast_keys[[tp]]$non_ref)),
                    contrast_keys[[tp]]$non_ref)
  }))

  K <- length(all_contrast_entries)

  # contrast vector Delta
  # all_meta$D already holds the bc-selected jump (D_bc when bc = TRUE).
  Delta <- vapply(all_contrast_entries, function(k) {
    all_meta[[k]]$D - all_meta[[ref_map[k]]]$D
  }, numeric(1))

  # covariance of Delta:
  # Var(D_A - D_refA - (D_B - D_refB))
  # = Var(D_A) + Var(D_refA) + Var(D_B) + Var(D_refB)
  #   - 2 Cov(D_A, D_refA) - 2 Cov(D_B, D_refB)
  #   + 2 Cov(D_A, D_B) - 2 Cov(D_A, D_refB) - 2 Cov(D_refA, D_B) + 2 Cov(D_refA, D_refB)
  # but within a period cross-type covariance = 0, so:
  #   If A and refA are in the same period: Cov(D_A, D_refA) = 0
  #   If A and B are in different periods: Cov(D_A, D_B) = .cov_scheme(...)
  #
  # Implementation: Sigma[i,j] = Cov(Delta_i, Delta_j)
  # Delta_i = D_{A_i} - D_{ref_i}
  # Cov(Delta_i, Delta_j) = Cov(D_{A_i}, D_{A_j}) - Cov(D_{A_i}, D_{ref_j})
  #                        - Cov(D_{ref_i}, D_{A_j}) + Cov(D_{ref_i}, D_{ref_j})

  # helper: Cov(D from fitA, D from fitB)
  cov_dd <- function(keyA, keyB) {
    fitA <- all_fits[[keyA]]; fitB <- all_fits[[keyB]]
    same_period <- (all_meta[[keyA]]$period == all_meta[[keyB]]$period)
    # same key = variance
    if (keyA == keyB) return(if (bc) fitA$V_D_bc else fitA$V_D)
    # same period, different type = 0 (disjoint samples)
    if (same_period) return(0)
    .cov_scheme(fitA, fitB, use_scheme, bc = bc)
  }

  Sigma <- matrix(NA_real_, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      Ai  <- all_contrast_entries[i]; rAi <- ref_map[Ai]
      Aj  <- all_contrast_entries[j]; rAj <- ref_map[Aj]
      Sigma[i, j] <- cov_dd(Ai, Aj) - cov_dd(Ai, rAj) -
                     cov_dd(rAi, Aj) + cov_dd(rAi, rAj)
    }
  }
  rownames(Sigma) <- colnames(Sigma) <- all_contrast_entries

  # ---- Wald statistic ----
  # Deliberately NOT .joint_wald(): this contrast covariance is a difference of
  # estimated covariances and can come back numerically indefinite.  We use an
  # eigen pseudo-inverse with a looser eps^0.5 tolerance that DROPS non-positive
  # eigen-directions (rather than inverting their signed magnitude as the SVD in
  # .joint_wald would), which is the more conservative choice here.
  ev   <- eigen(Sigma, symmetric = TRUE)
  tol  <- max(abs(ev$values)) * K * .Machine$double.eps^0.5
  pos  <- ev$values > tol
  df   <- sum(pos)
  if (df == 0L) stop("rd_homog: estimated covariance matrix is numerically zero.")

  Sigma_inv <- ev$vectors[, pos, drop = FALSE] %*%
               diag(1 / ev$values[pos], nrow = sum(pos)) %*%
               t(ev$vectors[, pos, drop = FALSE])

  W  <- as.numeric(t(Delta) %*% Sigma_inv %*% Delta)
  pv <- stats::pchisq(W, df = df, lower.tail = FALSE)

  # ---- period-type jump table ----
  jump_df <- do.call(rbind, lapply(names(all_meta), function(k) {
    m <- all_meta[[k]]
    data.frame(period    = m$period,
               type      = m$type,
               jump      = m$D,
               se        = sqrt(m$V_D),
               n         = m$n,
               reference = m$ref,
               stringsAsFactors = FALSE)
  }))
  rownames(jump_df) <- NULL

  # ---- output ----
  structure(
    list(
      statistic         = W,
      df                = df,
      p_value           = pv,
      period_type_jumps = jump_df,
      contrasts         = Delta,
      cov_matrix        = Sigma,
      scheme            = use_scheme,
      bc                = bc,
      comparisons       = names(contrast_keys),
      call              = cl
    ),
    class = "rd_homog"
  )
}

#' @export
print.rd_homog <- function(x, ...) {
  cat("Type-homogeneous confounding test (Assumption A9)\n")
  cat(sprintf("  Comparison periods: %s\n",
              paste(x$comparisons, collapse = ", ")))
  cat(sprintf("  Sampling scheme: %s\n", toupper(x$scheme)))
  cat(sprintf("  Wald statistic: %.4f   df: %d   p-value: %.4f\n",
              x$statistic, x$df, x$p_value))
  cat("\n  NOTE: This test is NEITHER necessary NOR sufficient for Assumption A9\n")
  cat("  at the RD period. It is run in comparison periods only and is only\n")
  cat("  suggestive of homogeneity where the assumption is needed (t_RD).\n")
  if (nrow(x$period_type_jumps) > 0L) {
    cat("\n  Per-type jumps (comparison periods):\n")
    df <- x$period_type_jumps
    df$ref_marker <- ifelse(df$reference, "[ref]", "")
    cat(sprintf("  %-10s  %-12s  %8s  %8s  %6s  %s\n",
                "Period", "Type", "Jump", "SE", "n", ""))
    for (k in seq_len(nrow(df)))
      cat(sprintf("  %-10s  %-12s  %8.4f  %8.4f  %6d  %s\n",
                  df$period[k], df$type[k], df$jump[k], df$se[k],
                  df$n[k], df$ref_marker[k]))
  }
  invisible(x)
}
