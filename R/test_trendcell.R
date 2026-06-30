# Test of Assumption A10: within-type confounding trend (ass:trend-cell).
# Manuscript ref: Leventer and Nevo, "Correcting Invalid RD Designs",
# paragraph "Within-cell confounding trend (Assumption A10)".
#
# The test is run in COMPARISON PERIODS ONLY.  In comparison period t0 the
# outcome RD jump equals the pure confounding:
#   D_{t0}(k) = alpha_{t0,0}(k),
# where k denotes the unit's cell (fixed across comparison periods).
# We estimate that jump separately per cell using rd_period(), then form a
# Wald test that the per-cell jumps are flat (constant) or linear in t0
# across comparison periods.
#
# IMPORTANT: this test is NEITHER necessary NOR sufficient for the trend-cell
# assumption to hold at the RD period; see documentation for rd_trendcell().

# ---------------------------------------------------------------------------
# Cell assignment for A10 must be FIXED across comparison periods.
# This differs from A9 (rd_homog) where the type is period-specific.
# We build a single cell_map once:
#   type_by = "rd_side"  → sign of the unit's running variable in t_rd.
#   type_by = "pattern"  → type from the t_rd perspective (sign pattern of all
#                          comparison periods); if t_rd is NULL, uses the first
#                          comparison period's perspective instead.
# Because the cell is fixed, the covariance structure of the stacked contrasts
# is block-diagonal by cell (cross-cell covariance is zero even across periods).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Test of within-cell confounding trend (Assumption A10)
#'
#' Tests whether the per-cell outcome RD discontinuity is constant (or linear)
#' **across comparison periods** within each cell.
#' In a comparison period \eqn{t_0} the discontinuity equals the pure
#' confounding: \eqn{D_{t_0}(k) = \alpha_{t_0,0}(k)}, where \eqn{k} denotes
#' the unit's cell (side of the cutoff in \eqn{t_{\mathrm{RD}}} under the
#' default \code{type_by = "rd_side"}).  The function estimates the jump per
#' cell via [rd_period()], then forms a joint Wald test that the within-cell
#' jumps conform to the hypothesised trend \eqn{g_0} across comparison periods.
#'
#' ## Scope and interpretation
#'
#' **This test is run in comparison periods only.**  Assumption A10
#' (within-cell confounding trend, \eqn{\text{ass:trend-cell}} in Leventer
#' and Nevo) is needed at the *RD period* \eqn{t_{\mathrm{RD}}}, but there
#' the jump also contains the ATT and the confounding is not separately
#' observable.  Testing the trend assumption in comparison periods is therefore
#' only **suggestive** evidence, analogous to a pre-trends check in
#' difference-in-differences: it is **neither necessary nor sufficient** for the
#' assumption to hold at \eqn{t_{\mathrm{RD}}}.
#'
#' ## Null hypothesis
#'
#' \deqn{H_0 : D_{t_0}(k) \text{ is } g_0 \text{ in } t_0,
#'   \text{ jointly for all cells } k,}
#' where \eqn{g_0} is `trend = "constant"` (equal jumps across all comparison
#' periods) or `trend = "linear"` (jumps lie on a line in \eqn{t_0}).
#'
#' Under `trend = "constant"`, the contrasts for cell \eqn{k} with
#' \eqn{|T_0|} comparison periods are
#' \eqn{D_{t_0}(k) - D_{t_1}(k)} for \eqn{t_0 \neq t_1} (reference =
#' first comparison period), yielding \eqn{|T_0| - 1} contrasts per cell.
#'
#' Under `trend = "linear"`, the testable contrasts are the **second differences**
#' of the time-ordered per-cell jumps:
#' \eqn{D_{t_{j+1}}(k) - 2 D_{t_j}(k) + D_{t_{j-1}}(k)},
#' yielding \eqn{|T_0| - 2} contrasts per cell.  **This requires at least 3
#' comparison periods per cell**; if no cell reaches this threshold the
#' function returns an object with `df = 0`, `statistic = NA`, and a message.
#'
#' @param data A long data frame, one row per unit-period.
#' @param y,x,time Column names (character strings) for the outcome, running
#'   variable, and period indicator.
#' @param id Column name for the unit identifier.  Required (the test needs a
#'   panel to define a fixed cell per unit).
#' @param comparisons Values of `time` to use as comparison periods.  Defaults
#'   to all periods except `t_rd` (if `t_rd` is supplied) or all periods (if
#'   `t_rd = NULL`).
#' @param t_rd Value of `time` for the RD period.  Required for
#'   `type_by = "rd_side"`.  Under `type_by = "pattern"`, if supplied, the
#'   cell is the sign pattern of all comparison periods from \eqn{t_{\mathrm{RD}}}'s
#'   perspective; if `NULL`, the pattern is taken from the first comparison
#'   period's perspective.
#' @param c Cutoff value for the running variable (default 0).
#' @param h Main bandwidth.  If `NULL` (default), bandwidth is chosen
#'   according to `bwselect`.  An explicit numeric value overrides `bwselect`
#'   and is used directly for every (cell, period) combination.
#' @param bwselect Bandwidth selection when `h = NULL`:
#'   `"cct"` (default) computes a per-(cell, period) CCT MSE-optimal bandwidth
#'   via [rd_bw_cct()] on that cell's outcome and running variable;
#'   `"rot"` uses the 0.2 × range rule of thumb (the original behavior).
#'   Ignored when `h` is supplied.
#' @param kernel Kernel for the local-linear RD: `"triangular"` (default),
#'   `"epanechnikov"`, or `"uniform"`.
#' @param scheme Sampling scheme for the cross-period covariance:
#'   `"auto"` (detect from the id/side structure, default), `"cs"` (repeated
#'   cross-section, no cross-period terms), `"pc"` (panel, time-constant
#'   running variable), or `"pv"` (panel, time-varying running variable).
#' @param min_n Minimum number of observations per cell-side before a
#'   (cell, period) pair is included.  Pairs with fewer than `min_n` obs on
#'   either side are silently dropped (default 10).
#' @param bc Use robust bias-corrected per-cell jumps and variances
#'   (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test
#'   with the bias-corrected [rddid()] estimator; `FALSE` uses the conventional
#'   local-linear jumps and variances.
#' @param type_by How a unit's cell is defined.  `"rd_side"` (default) uses
#'   the unit's side of the cutoff in the RD period `t_rd`, a binary partition;
#'   this is the canonical partition for a joint cross-period trend test.
#'   `"pattern"` uses the sign pattern of the other periods' running variables
#'   (fixed via the `t_rd` or first-comparison-period perspective).  The cell
#'   is **fixed** across comparison periods for all choices.
#' @param trend Trend form for the null hypothesis.  `"constant"` (default)
#'   tests equal per-cell jumps across all comparison periods; `"linear"` tests
#'   that the per-cell jumps lie on a line in time, using second-difference
#'   contrasts (requires \eqn{\geq 3} comparison periods per cell).
#' @param ... Further arguments passed to [rd_period()] (e.g. `p`, `q`, `b`).
#'
#' @return An object of class `"rd_trendcell"`, a list with:
#'   \describe{
#'     \item{`statistic`}{Joint Wald chi-squared statistic (`NA` when `df = 0`).}
#'     \item{`df`}{Degrees of freedom (number of positive eigenvalues used).
#'       `0` when `trend = "linear"` and no cell has \eqn{\geq 3} comparison
#'       periods.}
#'     \item{`p_value`}{p-value from the chi-squared distribution (`NA` when
#'       `df = 0`).}
#'     \item{`cell_period_jumps`}{Data frame with one row per (cell, comparison
#'       period): cell, period, jump estimate, SE, number of observations, and
#'       a flag indicating whether this is the reference period for that cell
#'       (under `trend = "constant"` only).}
#'     \item{`contrasts`}{Named numeric vector of jump contrasts stacked across
#'       cells.}
#'     \item{`cov_matrix`}{Estimated covariance matrix of the contrasts;
#'       block-diagonal by cell.}
#'     \item{`scheme`}{Sampling scheme used.}
#'     \item{`bc`}{Whether bias-corrected jumps were used.}
#'     \item{`trend`}{Trend form used (`"constant"` or `"linear"`).}
#'     \item{`comparisons`}{Comparison periods actually used (character).}
#'     \item{`call`}{The matched call.}
#'   }
#'
#' @note
#' **Necessary and sufficient status:** This test is *neither necessary nor
#' sufficient* for Assumption A10 (within-cell confounding trend) to hold at
#' the RD period.  Conformity to the trend in comparison periods is only
#' suggestive because the assumption is needed at \eqn{t_{\mathrm{RD}}},
#' where the confounding jump is not separately identified.  It is a pre-trend
#' check in the spirit of difference-in-differences and should be interpreted
#' as such.
#'
#' **`trend = "linear"`** requires at least 3 comparison periods per cell to
#' be informative.  With only 2 comparison periods the within-cell linear trend
#' is just-identified (any two points define a line), so no second-difference
#' contrast exists and the test returns `df = 0`.
#'
#' @seealso [rd_homog()], [rd_period()], [rddid()]
#'
#' @references
#' Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
#' Designs." Working paper.
#'
#' Calonico, S., Cattaneo, M. D., and Titiunik, R. (2014). Robust
#' nonparametric confidence intervals for regression-discontinuity designs.
#' *Econometrica*, 82(6), 2295-2326.
#'
#' @examples
#' \dontrun{
#' # Three-period panel: period 3 is RD, periods 1 and 2 are comparison.
#' # Constant within-cell confounding across periods 1 and 2 (H0 holds).
#' set.seed(1)
#' n <- 600
#' r_rd   <- runif(n, -1, 1)
#' r_comp <- runif(n, -1, 1)   # same running variable both comparison periods
#' cell   <- ifelse(r_rd >= 0, "+", "-")
#' # constant confounding: same jump in both comparison periods within each cell
#' alpha  <- ifelse(cell == "+", 0.5, 0.3)
#' y1 <- 0.3 * r_comp + alpha * (r_comp >= 0) + rnorm(n, 0, 0.3)
#' y2 <- 0.3 * r_comp + alpha * (r_comp >= 0) + rnorm(n, 0, 0.3)
#' y3 <- 0.3 * r_rd + 1.0 * (r_rd >= 0) + rnorm(n, 0, 0.3)  # RD period
#' dat <- data.frame(
#'   id   = rep(seq_len(n), 3),
#'   time = rep(1:3, each = n),
#'   x    = c(r_comp, r_comp, r_rd),
#'   y    = c(y1, y2, y3)
#' )
#' rd_trendcell(dat, y = "y", x = "x", time = "time", id = "id",
#'              comparisons = 1:2, t_rd = 3, h = 0.3)
#' }
#' @export
rd_trendcell <- function(data, y, x, time, id,
                         comparisons = NULL, t_rd = NULL,
                         c = 0, h = NULL,
                         bwselect = c("cct", "rot"),
                         kernel = "triangular",
                         scheme = c("auto", "cs", "pc", "pv"),
                         min_n = 10L,
                         bc = TRUE,
                         type_by = c("rd_side", "pattern"),
                         trend   = c("constant", "linear"),
                         ...) {
  cl       <- match.call()
  scheme   <- match.arg(scheme)
  type_by  <- match.arg(type_by)
  trend    <- match.arg(trend)
  bwselect <- match.arg(bwselect)

  # ---- validate columns ----
  for (nm in base::c(y, x, time, id))
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")

  times_all <- sort(unique(data[[time]]))

  # ---- comparison periods (sorted for consistent time ordering) ----
  if (is.null(comparisons)) {
    comparisons <- if (!is.null(t_rd)) setdiff(times_all, t_rd) else times_all
  }
  comparisons <- sort(comparisons)
  if (length(comparisons) < 1L)
    stop("need at least one comparison period.")

  # ---- build fixed cell assignment ----
  # A10 requires the cell to be FIXED across all comparison periods so that
  # cross-cell covariance is zero regardless of which periods we compare.
  bt <- .build_types(data, x = x, time = time, id = id, c = c)

  if (type_by == "rd_side") {
    if (is.null(t_rd))
      stop("type_by = \"rd_side\" requires `t_rd` (the RD period whose side defines the cell).")
    sidecol <- paste0("side_", t_rd)
    if (!sidecol %in% names(bt$wide))
      stop("RD period '", t_rd, "' has no side column in the data; cannot define rd_side cells.")
    cell_map <- stats::setNames(bt$wide[[sidecol]], as.character(bt$wide$id))
  } else {
    # pattern: use type from a fixed reference period so the cell is constant
    # across comparison periods.  Prefer t_rd (sign pattern of all comparison
    # periods); fall back to the first comparison period if t_rd is absent.
    ref_tp <- if (!is.null(t_rd) && as.character(t_rd) %in% names(bt$period_types))
      as.character(t_rd)
    else
      as.character(comparisons[1L])
    tdf_ref  <- bt$period_types[[ref_tp]]
    cell_map <- stats::setNames(tdf_ref$type, as.character(tdf_ref$id))
  }

  # ---- detect sampling scheme across comparison periods ----
  comp_plist <- stats::setNames(lapply(as.character(comparisons), function(tp) {
    d_cp <- data[data[[time]] == tp, , drop = FALSE]
    list(id = d_cp[[id]], x = d_cp[[x]])
  }), as.character(comparisons))
  detected_scheme <- .detect_scheme(comp_plist, c = c)
  use_scheme <- if (scheme == "auto") detected_scheme else scheme

  # ---- per-(cell, comparison-period) fits ----
  # key format: "cell::period"
  dots     <- list(...)
  all_fits <- list()
  all_meta <- list()

  for (tp in as.character(comparisons)) {
    rows  <- which(data[[time]] == tp)
    d_tp  <- data[rows, , drop = FALSE]
    id_tp <- d_tp[[id]]
    cvec  <- unname(cell_map[as.character(id_tp)])

    for (ck in sort(unique(cvec[!is.na(cvec)]))) {
      keep  <- !is.na(cvec) & cvec == ck
      y_ck  <- d_tp[[y]][keep]
      x_ck  <- d_tp[[x]][keep]
      id_ck <- d_tp[[id]][keep]

      # per-cell bandwidth (0.2*range per cell for rot; CCT per cell for cct)
      bw   <- .cell_bandwidth(y_ck, x_ck, c, kernel, h, bwselect)
      bw_h <- bw[["h"]]
      bw_b <- bw[["b"]]

      n_pos <- sum(x_ck >= c, na.rm = TRUE)
      n_neg <- sum(x_ck <  c, na.rm = TRUE)
      if (n_pos < min_n || n_neg < min_n) next

      fit <- tryCatch({
        call_args <- base::c(
          list(y = y_ck, x = x_ck, h = bw_h, b = bw_b, id = id_ck, c = c,
               kernel = kernel),
          dots)
        do.call(rd_period, call_args)
      }, error = function(e) NULL)
      if (is.null(fit)) next

      key            <- paste0(ck, "::", tp)
      all_fits[[key]] <- fit
      all_meta[[key]] <- list(cell   = ck,
                               period = tp,
                               D      = if (bc) fit$D_bc else fit$D,
                               V_D    = if (bc) fit$V_D_bc else fit$V_D,
                               n      = fit$n)
    }
  }

  # ---- per-cell: build contrast vector and covariance block ----
  # For each cell k, collect its valid (period, fit) pairs in time order.
  # Build:
  #   D_k       — m_k-vector of per-period jumps
  #   Sigma_k   — m_k x m_k covariance matrix of the D_k entries
  #   C_k       — n_contrasts_k x m_k contrast matrix
  #   Delta_k   = C_k %*% D_k
  #   Omega_k   = C_k %*% Sigma_k %*% t(C_k)   (covariance of Delta_k)
  #
  # The full contrast vector and covariance matrix are obtained by stacking
  # Delta_k and block-diagonally assembling Omega_k across cells.
  #
  # Within a cell, cov_dd returns the cross-period covariance from .cov_scheme()
  # (using the panel g-vectors from rd_period()).  Across cells, all covariances
  # are zero because cells partition the sample and a unit's cell is fixed across
  # periods.

  cov_dd <- function(keyA, keyB) {
    fitA <- all_fits[[keyA]]
    fitB <- all_fits[[keyB]]
    if (keyA == keyB) return(if (bc) fitA$V_D_bc else fitA$V_D)
    .cov_scheme(fitA, fitB, use_scheme, bc = bc)
  }

  valid_cells <- sort(unique(vapply(names(all_meta), function(k)
    all_meta[[k]]$cell, character(1))))

  Delta_all    <- numeric(0)
  Sigma_all    <- matrix(numeric(0), nrow = 0, ncol = 0)
  delta_labels <- character(0)
  jump_rows    <- list()   # for cell_period_jumps output table

  for (ck in valid_cells) {
    # ordered keys for this cell (time order = comparisons order, already sorted)
    cell_keys <- paste0(ck, "::", as.character(comparisons))
    cell_keys <- cell_keys[cell_keys %in% names(all_fits)]
    m_k <- length(cell_keys)

    # record per-row metadata (reference flag set below)
    for (k in cell_keys) {
      m <- all_meta[[k]]
      jump_rows[[k]] <- data.frame(
        cell      = m$cell,
        period    = m$period,
        jump      = m$D,
        se        = sqrt(m$V_D),
        n         = m$n,
        reference = FALSE,
        stringsAsFactors = FALSE
      )
    }

    n_contrasts_k <- switch(trend,
      constant = max(0L, m_k - 1L),
      linear   = max(0L, m_k - 2L)
    )
    if (n_contrasts_k == 0L) next

    # mark reference in the output table
    if (trend == "constant")
      jump_rows[[cell_keys[1L]]]$reference <- TRUE

    # D vector and jump covariance matrix for this cell
    D_k     <- vapply(cell_keys, function(k) all_meta[[k]]$D, numeric(1))
    Sigma_k <- matrix(NA_real_, nrow = m_k, ncol = m_k)
    for (i in seq_len(m_k))
      for (j in seq_len(m_k))
        Sigma_k[i, j] <- cov_dd(cell_keys[i], cell_keys[j])

    # contrast matrix C_k
    if (trend == "constant") {
      # rows: e_{i+1} - e_1  (difference from first-period reference)
      C_k <- matrix(0, nrow = n_contrasts_k, ncol = m_k)
      for (i in seq_len(n_contrasts_k)) {
        C_k[i, 1L]     <- -1
        C_k[i, i + 1L] <-  1
      }
      ref_period <- all_meta[[cell_keys[1L]]]$period
      labs_k <- paste0(ck, "::",
                       vapply(cell_keys[-1L], function(k) all_meta[[k]]$period,
                              character(1)),
                       "-", ref_period)
    } else {
      # second differences: rows e_{i} - 2*e_{i+1} + e_{i+2}
      C_k <- matrix(0, nrow = n_contrasts_k, ncol = m_k)
      for (i in seq_len(n_contrasts_k)) {
        C_k[i, i]       <-  1
        C_k[i, i + 1L]  <- -2
        C_k[i, i + 2L]  <-  1
      }
      labs_k <- paste0(ck, "::2nd_diff_", seq_len(n_contrasts_k))
    }

    Delta_k  <- as.numeric(C_k %*% D_k)
    Omega_k  <- C_k %*% Sigma_k %*% t(C_k)

    # append to overall contrast vector and block-expand covariance matrix
    n_prev   <- length(Delta_all)
    n_new    <- n_contrasts_k
    new_size <- n_prev + n_new

    Sigma_big <- matrix(0, nrow = new_size, ncol = new_size)
    if (n_prev > 0L)
      Sigma_big[seq_len(n_prev), seq_len(n_prev)] <- Sigma_all
    Sigma_big[(n_prev + 1L):new_size, (n_prev + 1L):new_size] <- Omega_k

    Delta_all    <- base::c(Delta_all, Delta_k)
    Sigma_all    <- Sigma_big
    delta_labels <- base::c(delta_labels, labs_k)
  }

  names(Delta_all) <- delta_labels
  if (length(delta_labels) > 0L)
    rownames(Sigma_all) <- colnames(Sigma_all) <- delta_labels

  # ---- build output jump table ----
  jump_df <- if (length(jump_rows) > 0L) {
    df0 <- do.call(rbind, jump_rows)
    rownames(df0) <- NULL
    df0
  } else {
    data.frame(cell = character(0), period = character(0),
               jump = numeric(0), se = numeric(0),
               n = integer(0), reference = logical(0),
               stringsAsFactors = FALSE)
  }

  # ---- handle no-contrast case ----
  if (length(Delta_all) == 0L) {
    if (trend == "linear") {
      message("rd_trendcell: linear trend is not testable -- no cell has ",
              "3 or more comparison periods (degrees of freedom = 0). ",
              "Returning an object with df = 0, statistic = NA, p_value = NA.")
      return(structure(
        list(statistic        = NA_real_,
             df               = 0L,
             p_value          = NA_real_,
             cell_period_jumps = jump_df,
             contrasts        = numeric(0),
             cov_matrix       = matrix(numeric(0), 0L, 0L),
             scheme           = use_scheme,
             bc               = bc,
             trend            = trend,
             comparisons      = as.character(comparisons),
             call             = cl),
        class = "rd_trendcell"
      ))
    }
    stop("rd_trendcell: no usable within-cell cross-period contrasts found; ",
         "check data, bandwidth, min_n, or number of comparison periods.")
  }

  # ---- Wald statistic via eigen pseudo-inverse ----
  # The contrast covariance is a linear combination of estimated covariances
  # (C_k %*% Sigma_k %*% C_k') and can be numerically indefinite.  Use the
  # conservative eigen pseudo-inverse that drops non-positive directions
  # (see .wald_eigen() in R/test_helpers.R for the implementation).
  ew <- .wald_eigen(Delta_all, Sigma_all)
  if (ew$df == 0L)
    stop("rd_trendcell: estimated covariance matrix is numerically zero.")
  W  <- ew$stat
  df <- ew$df
  pv <- ew$p

  # ---- output ----
  structure(
    list(statistic        = W,
         df               = df,
         p_value          = pv,
         cell_period_jumps = jump_df,
         contrasts        = Delta_all,
         cov_matrix       = Sigma_all,
         scheme           = use_scheme,
         bc               = bc,
         trend            = trend,
         comparisons      = as.character(comparisons),
         call             = cl),
    class = "rd_trendcell"
  )
}

#' @export
print.rd_trendcell <- function(x, ...) {
  cat("Within-cell confounding trend test (Assumption A10)\n")
  cat(sprintf("  Trend form: %s\n", x$trend))
  cat(sprintf("  Comparison periods: %s\n",
              paste(x$comparisons, collapse = ", ")))
  cat(sprintf("  Sampling scheme: %s\n", toupper(x$scheme)))
  if (is.na(x$statistic)) {
    cat(sprintf("  Wald statistic: NA   df: %d   p-value: NA\n", x$df))
    cat("\n  NOTE: df = 0; the linear within-cell trend is just-identified\n")
    cat("  with fewer than 3 comparison periods (no second-difference contrasts).\n")
  } else {
    cat(sprintf("  Wald statistic: %.4f   df: %d   p-value: %.4f\n",
                x$statistic, x$df, x$p_value))
  }
  cat("\n  NOTE: This test is NEITHER necessary NOR sufficient for Assumption A10\n")
  cat("  at the RD period. It is a pre-trend check run in comparison periods\n")
  cat("  only and is only suggestive of trend conformity where the assumption\n")
  cat("  is needed (t_RD).\n")
  if (!is.null(x$cell_period_jumps) && nrow(x$cell_period_jumps) > 0L) {
    cat("\n  Per-cell jumps (comparison periods):\n")
    df <- x$cell_period_jumps
    df$ref_marker <- ifelse(df$reference, "[ref]", "")
    cat(sprintf("  %-10s  %-12s  %8s  %8s  %6s  %s\n",
                "Cell", "Period", "Jump", "SE", "n", ""))
    for (k in seq_len(nrow(df)))
      cat(sprintf("  %-10s  %-12s  %8.4f  %8.4f  %6d  %s\n",
                  df$cell[k], df$period[k], df$jump[k], df$se[k],
                  df$n[k], df$ref_marker[k]))
  }
  invisible(x)
}
