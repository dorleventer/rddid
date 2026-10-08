# rd_trendcell.R -- rd_trendcell(): test of a constant (or linear) within-type confounding jump
# across the comparison periods, with its print method. Shared pieces: assumption_tests_helpers.R,
# sampling_scheme.R, cross_period_covariance.R.

# A unit's type is fixed across the comparison periods (in rd_homog it is period-specific), so
# different types are different units: the stacked contrasts have a block-diagonal covariance.

#' Each unit's type, named by unit id, fixed across the comparison periods.
#' @noRd
.trendcell_cell_map <- function(types, type_by, t_rd, comparisons) {
  if (type_by == "rd_side") {
    if (is.null(t_rd))
      stop("type_by = \"rd_side\" requires `t_rd` (the RD period whose side defines the cell).")
    side_col <- paste0("side_", t_rd)
    if (!side_col %in% names(types$wide))
      stop("RD period '", t_rd, "' has no side column in the data; cannot define rd_side cells.")
    cell_map <- stats::setNames(types$wide[[side_col]], as.character(types$wide$id))
  } else {
    # the types as seen from ONE period, so that they do not change across the comparison
    # periods: from t_rd (the sides in all other periods), else from the first comparison period
    ref_tp <- if (!is.null(t_rd) && as.character(t_rd) %in% names(types$period_types))
      as.character(t_rd)
    else
      as.character(comparisons[1L])
    ref_types <- types$period_types[[ref_tp]]
    cell_map  <- stats::setNames(ref_types$type, as.character(ref_types$id))
  }
  cell_map
}


#' The fit of every (type, comparison period) cell: list(fits, meta, skipped), keyed "type::period".
#' @noRd
.trendcell_fit_cells <- function(data, y, x, time, id, comparisons, cell_map, cutoff, kernel, h,
                                 bwselect, min_n, bc, p, q) {
  orders   <- list(p = p, q = q)   # polynomial orders of the per-cell fits
  skipped  <- character()
  all_fits <- list()
  all_meta <- list()

  for (tp in as.character(comparisons)) {
    rows    <- which(data[[time]] == tp)
    d_tp    <- data[rows, , drop = FALSE]
    id_tp   <- d_tp[[id]]
    cell_tp <- unname(cell_map[as.character(id_tp)])

    # radix sorts in byte order whatever the locale, so "+" comes before "-" on every machine
    for (cell_k in sort(unique(cell_tp[!is.na(cell_tp)]), method = "radix")) {
      keep <- !is.na(cell_tp) & cell_tp == cell_k
      y_k  <- d_tp[[y]][keep]
      x_k  <- d_tp[[x]][keep]
      id_k <- d_tp[[id]][keep]

      # the bandwidth comes before the min_n check, so CCT also runs on cells skipped below
      bw   <- .cell_bandwidth(y_k, x_k, cutoff, kernel, h, bwselect, p = p)
      bw_h <- bw[["h"]]
      bw_b <- bw[["b"]]

      n_pos <- sum(x_k >= cutoff, na.rm = TRUE)
      n_neg <- sum(x_k <  cutoff, na.rm = TRUE)
      if (n_pos < min_n || n_neg < min_n) {
        skipped <- c(skipped, sprintf("cell %s, period %s (n = %d below, %d above)",
                                      cell_k, tp, n_neg, n_pos))
        next
      }

      # rd_period() stops when a side has too few observations in the pilot window or a
      # singular design; such a cell is skipped and listed, not fatal to the test
      fit <- tryCatch({
        call_args <- base::c(
          list(y = y_k, x = x_k, h = bw_h, b = bw_b, id = id_k, c = cutoff,
               kernel = kernel),
          orders)
        do.call(rd_period, call_args)
      }, error = function(e) NULL)
      if (is.null(fit)) {
        skipped <- c(skipped, sprintf("cell %s, period %s (local-linear fit failed)", cell_k, tp))
        next
      }

      key             <- paste0(cell_k, "::", tp)
      all_fits[[key]] <- fit
      all_meta[[key]] <- list(cell   = cell_k,
                              period = tp,
                              D      = if (bc) fit$D_bc else fit$D,
                              V_D    = if (bc) fit$V_D_bc else fit$V_D,
                              n      = fit$n)
    }
  }
  list(fits = all_fits, meta = all_meta, skipped = skipped)
}

#' Covariance matrix of one type's jumps, rows and columns in `cell_keys` (time) order.
#' @noRd
.trendcell_jump_cov <- function(all_fits, cell_keys, use_scheme, bc) {
  # same period: the jump's variance. Two periods: .cov_scheme() -- 0 under "cs"; under "pc" the
  # covariance from units on the same side in both periods; under "pv" that, minus the
  # covariance from units that change side.
  cov_dd <- function(keyA, keyB) {
    fitA <- all_fits[[keyA]]
    fitB <- all_fits[[keyB]]
    if (keyA == keyB) return(if (bc) fitA$V_D_bc else fitA$V_D)
    .cov_scheme(fitA, fitB, use_scheme, bc = bc)
  }

  m_k     <- length(cell_keys)
  Sigma_k <- matrix(NA_real_, nrow = m_k, ncol = m_k)
  for (i in seq_len(m_k))
    for (j in seq_len(m_k))
      Sigma_k[i, j] <- cov_dd(cell_keys[i], cell_keys[j])
  Sigma_k
}

#' One type's contrast matrix (rows = contrasts, columns = its periods) and labels: list(C, labels).
#' @noRd
.trendcell_contrast_matrix <- function(trend, cell_k, cell_keys, all_meta, n_contrasts_k, m_k) {
  if (trend == "constant") {
    # rows e_{i+1} - e_1: m_k jumps give m_k - 1 free differences, all taken from the first period
    C_k <- matrix(0, nrow = n_contrasts_k, ncol = m_k)
    for (i in seq_len(n_contrasts_k)) {
      C_k[i, 1L]     <- -1
      C_k[i, i + 1L] <-  1
    }
    ref_period <- all_meta[[cell_keys[1L]]]$period
    labels_k <- paste0(cell_k, "::",
                       vapply(cell_keys[-1L], function(key) all_meta[[key]]$period,
                              character(1)),
                       "-", ref_period)
  } else {
    # Second differences IN TIME of the period-ordered jumps: for consecutive periods
    # (t0, t1, t2) the weights (t2 - t1, -(t2 - t0), t1 - t0) annihilate any jump that is
    # linear in t, whatever the spacing; scaled by 2 / (t2 - t0) they are exactly (1, -2, 1)
    # for equally spaced periods (2026-10-08 fix: the rows used to be index-based, which rejects
    # a true linear trend when the comparison periods are unequally spaced).
    tv <- suppressWarnings(as.numeric(vapply(cell_keys, function(key) all_meta[[key]]$period,
                                             character(1))))
    if (anyNA(tv))
      stop("rd_trendcell: trend = \"linear\" needs numeric period values.")
    C_k <- matrix(0, nrow = n_contrasts_k, ncol = m_k)
    for (i in seq_len(n_contrasts_k)) {
      t0 <- tv[i]; t1 <- tv[i + 1L]; t2 <- tv[i + 2L]
      C_k[i, i]       <-  2 * (t2 - t1) / (t2 - t0)
      C_k[i, i + 1L]  <- -2 * (t2 - t0) / (t2 - t0)
      C_k[i, i + 2L]  <-  2 * (t1 - t0) / (t2 - t0)
    }
    labels_k <- paste0(cell_k, "::2nd_diff_", seq_len(n_contrasts_k))
  }
  list(C = C_k, labels = labels_k)
}

#' `Sigma_all` grown by the diagonal block `Omega_k`, with zeros between the blocks.
#' @noRd
.trendcell_block_append <- function(Sigma_all, Omega_k, n_prev, n_new) {
  new_size <- n_prev + n_new

  # different types are different units, so contrasts of different types do not covary
  Sigma_big <- matrix(0, nrow = new_size, ncol = new_size)
  if (n_prev > 0L)
    Sigma_big[seq_len(n_prev), seq_len(n_prev)] <- Sigma_all
  Sigma_big[(n_prev + 1L):new_size, (n_prev + 1L):new_size] <- Omega_k
  Sigma_big
}

#' Stacked contrasts, their covariance and the jump-table rows: list(Delta, Sigma, jump_rows).
#' @noRd
.trendcell_stack_contrasts <- function(all_fits, all_meta, comparisons, trend, use_scheme, bc) {
  # radix: the same locale-free type order as in .trendcell_fit_cells()
  valid_cells <- sort(unique(vapply(names(all_meta), function(key)
    all_meta[[key]]$cell, character(1))), method = "radix")

  Delta_all    <- numeric(0)
  Sigma_all    <- matrix(numeric(0), nrow = 0, ncol = 0)
  delta_labels <- character(0)
  jump_rows    <- list()

  for (cell_k in valid_cells) {
    # this type's fitted periods, in time order (comparisons is sorted)
    cell_keys <- paste0(cell_k, "::", as.character(comparisons))
    cell_keys <- cell_keys[cell_keys %in% names(all_fits)]
    m_k <- length(cell_keys)

    for (key in cell_keys) {
      meta <- all_meta[[key]]
      jump_rows[[key]] <- data.frame(
        cell      = meta$cell,
        period    = meta$period,
        jump      = meta$D,
        se        = sqrt(meta$V_D),
        n         = meta$n,
        reference = FALSE,
        stringsAsFactors = FALSE
      )
    }

    n_contrasts_k <- switch(trend,
      constant = max(0L, m_k - 1L),
      linear   = max(0L, m_k - 2L)
    )
    if (n_contrasts_k == 0L) next

    # under "constant" the first period is the one every other period is compared with
    if (trend == "constant")
      jump_rows[[cell_keys[1L]]]$reference <- TRUE

    D_k     <- vapply(cell_keys, function(key) all_meta[[key]]$D, numeric(1))
    Sigma_k <- .trendcell_jump_cov(all_fits, cell_keys, use_scheme, bc)

    contrast <- .trendcell_contrast_matrix(trend, cell_k, cell_keys, all_meta, n_contrasts_k, m_k)
    C_k      <- contrast$C
    labels_k <- contrast$labels

    Delta_k <- as.numeric(C_k %*% D_k)
    Omega_k <- C_k %*% Sigma_k %*% t(C_k)

    Sigma_big    <- .trendcell_block_append(Sigma_all, Omega_k, n_prev = length(Delta_all),
                                            n_new = n_contrasts_k)
    Delta_all    <- base::c(Delta_all, Delta_k)
    Sigma_all    <- Sigma_big
    delta_labels <- base::c(delta_labels, labels_k)
  }

  names(Delta_all) <- delta_labels
  if (length(delta_labels) > 0L)
    rownames(Sigma_all) <- colnames(Sigma_all) <- delta_labels
  list(Delta = Delta_all, Sigma = Sigma_all, jump_rows = jump_rows)
}

#' The `cell_period_jumps` data frame; with no fitted cell, an empty one with the same columns.
#' @noRd
.trendcell_jump_table <- function(jump_rows) {
  if (length(jump_rows) > 0L) {
    jumps <- do.call(rbind, jump_rows)
    rownames(jumps) <- NULL
    jumps
  } else {
    data.frame(cell = character(0), period = character(0),
               jump = numeric(0), se = numeric(0),
               n = integer(0), reference = logical(0),
               stringsAsFactors = FALSE)
  }
}

#' The result when no type has three comparison periods under `trend = "linear"`: NA, df 0, NA.
#' @noRd
.trendcell_untestable <- function(jump_df, use_scheme, bc, estimand, trend, comparisons, cl) {
  message("rd_trendcell: linear trend is not testable -- no cell has ",
          "3 or more comparison periods (degrees of freedom = 0). ",
          "Returning an object with df = 0, statistic = NA, p_value = NA.")
  structure(
    list(statistic         = NA_real_,
         df                = 0L,
         p_value           = NA_real_,
         cell_period_jumps = jump_df,
         contrasts         = numeric(0),
         cov_matrix        = matrix(numeric(0), 0L, 0L),
         scheme            = use_scheme,
         bc                = bc,
         estimand          = estimand,
         trend             = trend,
         comparisons       = as.character(comparisons),
         call              = cl),
    class = "rd_trendcell"
  )
}

#' Test of constant within-type confounding
#'
#' When the running variable moves over time, some units are above the cutoff
#' in one period and below it in another. A unit's **type** is the side of the
#' cutoff it is on in the other period(s); by default here, its side in the RD
#' period, which stays fixed across the comparison periods. The
#' confounding-trend assumption of [rddid()] must then hold within each type.
#' It concerns the RD period, where the confounding jump is not observed
#' separately, so, like a pre-trends check in difference-in-differences,
#' `rd_trendcell()` tests it across the comparison periods: the null is that
#' **within each type, the confounding jump is the same in every comparison
#' period** (with `trend = "linear"`: moves linearly in time). A rejection
#' means the comparison periods do not support the trend assumption, and the
#' estimate of [rddid()] under that assumption can be biased; if the jumps
#' move linearly, consider `rddid(trend = "linear")`. With a running variable
#' fixed over time (as in [rddid_sim]) the types are degenerate and the test is
#' not informative.
#'
#' @details
#' ## What is estimated
#'
#' Each unit gets one type, fixed across the comparison periods. In each
#' comparison period and for each type, a local-linear RD of the outcome on the
#' running variable gives that type's confounding jump. Under
#' `trend = "constant"` each type's jump in every comparison period is compared
#' with its jump in the first one (one difference fewer than the number of
#' comparison periods, per type). Under `trend = "linear"` the test uses the
#' second differences of the time-ordered jumps (two fewer per type), so it
#' needs at least three comparison periods; if no type has jumps in three, the
#' function returns `statistic = NA` and `df = 0`, with a message. All
#' differences are tested jointly by a Wald statistic. Their covariance is
#' estimated and can be numerically indefinite, so the statistic uses only its
#' positive directions, and `df` counts them.
#'
#' ## Shared units and the sampling scheme
#'
#' A unit's type is fixed, so different types are different units and their
#' jumps are independent. Within a type the same units can appear in several
#' comparison periods, and the covariance follows `scheme`, matching units on
#' `id`: none under `"cs"`; from the units on the same side of the cutoff in
#' both periods under `"pc"`; under `"pv"` also from the units that change
#' side, with the opposite sign. `"auto"` reads the scheme off the comparison
#' periods, by the rule of [rddid()].
#'
#' ## Options
#'
#' `type_by = "rd_side"` (default) types each unit by its side of the cutoff in
#' the RD period, so `t_rd` is required; units not observed in the RD period
#' are left out. `type_by = "pattern"` types units by their sides in all
#' periods in `data` other than the RD period (other than the first comparison
#' period if `t_rd` is `NULL`). A (type, period) cell with fewer than `min_n`
#' observations on either side of the cutoff is dropped, and a message lists
#' the dropped cells. `bc = TRUE` (default) tests the bias-corrected jumps with
#' their robust variance, as in the `Robust` row of [rddid()]; `bc = FALSE`
#' uses the conventional jumps and variances. With `bwselect = "cct"` (default)
#' each cell gets its own CCT bandwidths from [rd_bw_cct()], computed on that
#' cell's outcome and running variable. With `bwselect = "rot"` the rule of
#' thumb is `h = b = 0.2` times the range of the running variable within each
#' cell. A numeric `h` is used as both bandwidths in every cell.
#'
#' ## ATU designs
#'
#' When the comparison periods are uniformly treated, the comparison-period
#' jumps are the confounding jumps among treated units. The computation is the
#' same, so `estimand` only labels the output.
#'
#' @param data a data frame in long format, one row per unit and period, from a
#'   panel (it need not be balanced). A unit's type is read from the periods in
#'   which it is observed; a unit missing from a period that its type needs is
#'   left out of the cells that use that type.
#' @param id name of the unit-identifier column (a string). Required: types are
#'   read across periods.
#' @param t_rd the RD period. Required with `type_by = "rd_side"` (the
#'   default), where a unit's side of the cutoff in the RD period is its type.
#'   The test does not use the RD period's outcomes.
#' @param comparisons the comparison periods in which the test runs, taken in
#'   time order. `NULL` (default) uses every period other than `t_rd` (every
#'   period if `t_rd` is `NULL`).
#' @param estimand `"att"` (default) or `"atu"`, as in [rddid()]. Label only:
#'   the test is the same either way.
#' @param trend the trend assumption tested within each type: `"constant"`
#'   (default; the confounding jump is the same in every comparison period) or
#'   `"linear"` (it moves linearly in time; needs at least three comparison
#'   periods). Use the `trend` of the [rddid()] call being checked.
#'   With `"linear"` the second differences are taken in time (the period values), so
#'   unequally spaced comparison periods are handled.
#' @param h a bandwidth to use, as both main and pilot bandwidth, in every cell.
#'   If given, `bwselect` is ignored.
#' @param bwselect the bandwidth rule when `h` is not given: `"cct"` (default;
#'   each cell's own CCT bandwidths from [rd_bw_cct()]) or `"rot"` (the rule of
#'   thumb `0.2` times the range of the running variable within the cell).
#' @param min_n the minimum number of observations on each side of the cutoff
#'   for a (type, period) cell to enter the test (default 10). Smaller cells
#'   are dropped, with a message listing them.
#' @param scheme the sampling scheme, which sets the covariance across
#'   comparison periods in the test: `"cs"`, `"pc"` or `"pv"` (as in
#'   [rddid()]), or `"auto"` (default), which reads it off the comparison
#'   periods by the rule of [rddid()]. See Details.
#' @param bc logical. `TRUE` (default): test the bias-corrected jumps with their
#'   robust variance, as in the `Robust` row of [rddid()]; `FALSE`: the
#'   conventional jumps and variances.
#' @param type_by how types are defined: `"rd_side"` (default; the unit's side
#'   of the cutoff in the RD period, which needs `t_rd`) or `"pattern"` (its
#'   sides in all periods other than the RD period). Either way the type is
#'   fixed across the comparison periods.
#' @param p,q orders of the local polynomials in every cell, for the point
#'   estimate and the bias correction (defaults 1 and 2; `q` must exceed `p`).
#'   The per-cell CCT bandwidths are chosen for the order-`p` fit.
#' @inheritParams rddid
#'
#' @return An object of class `"rd_trendcell"`, a list with:
#'   \describe{
#'     \item{`statistic`, `df`, `p_value`}{the Wald statistic, its degrees of
#'       freedom (the number of positive directions of the covariance used) and
#'       its chi-squared p-value; `NA`, `0`, `NA` when `trend = "linear"` and
#'       no type has jumps in three comparison periods.}
#'     \item{`scheme`}{the sampling scheme used.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`call`}{the matched call.}
#'     \item{`cell_period_jumps`}{a data frame with one row per (type, comparison
#'       period) cell that was fitted: `cell` (the type: the unit's side(s),
#'       `"+"` above and `"-"` below the cutoff), `period`, `jump` (the cell's
#'       confounding jump,
#'       bias-corrected when `bc = TRUE`), `se`, `n` (observations in the cell)
#'       and `reference` (`TRUE` for each type's first comparison period under
#'       `trend = "constant"`).}
#'     \item{`contrasts`}{named numeric vector of the tested differences,
#'       stacked across types.}
#'     \item{`cov_matrix`}{the estimated covariance matrix of `contrasts`
#'       (block-diagonal by type).}
#'     \item{`bc`, `trend`}{as passed.}
#'     \item{`comparisons`}{the comparison periods, in time order.}
#'   }
#'
#' @references
#' Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression Discontinuity
#' Designs Using Multiple Time-Period Data.* arXiv:2408.05847.
#' \url{https://arxiv.org/abs/2408.05847}
#'
#' Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust nonparametric
#' confidence intervals for regression-discontinuity designs. *Econometrica*
#' 82(6), 2295-2326.
#'
#' @seealso [rddid()] for the estimate; [rddid_sim_pv] for example data;
#'   `tidy()` in [rddid-tidiers] for a one-row summary.
#' @family tests of the assumptions
#'
#' @examples
#' # rddid_sim_pv: the running variable moves, so some units change side
#' tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' tr          # the Wald test, then each type's jump in each comparison period
#' # trend = "linear" needs three comparison periods; with two it is not testable
#' rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
#'              trend = "linear")
#' @export
rd_trendcell <- function(data, y, x, time, id,
                         t_rd = NULL, comparisons = NULL,
                         estimand = c("att", "atu"),
                         trend   = c("constant", "linear"),
                         c = 0, h = NULL,
                         bwselect = c("cct", "rot"),
                         kernel = "triangular",
                         scheme = c("auto", "cs", "pc", "pv"),
                         min_n = 10L,
                         bc = TRUE,
                         type_by = c("rd_side", "pattern"),
                         p = 1L, q = 2L) {
  cl       <- match.call()
  scheme   <- match.arg(scheme)
  type_by  <- match.arg(type_by)
  trend    <- match.arg(trend)
  bwselect <- match.arg(bwselect)
  estimand <- match.arg(estimand)
  kernel   <- match.arg(kernel, c("triangular", "epanechnikov", "uniform"))
  cutoff   <- c   # the cutoff; `c` stays the argument name for rdrobust users

  .check_columns(data, c(y, x, time, id))

  times_all <- sort(unique(data[[time]]))
  if (is.null(comparisons))
    comparisons <- if (!is.null(t_rd)) setdiff(times_all, t_rd) else times_all
  comparisons <- sort(comparisons)   # the contrasts below are in time order
  if (length(comparisons) < 1L)
    stop("need at least one comparison period.")

  types    <- .build_types(data, x = x, time = time, id = id, c = cutoff)
  cell_map <- .trendcell_cell_map(types, type_by, t_rd, comparisons)

  detected_scheme <- .detect_scheme_comparisons(data, x, time, id, comparisons, cutoff)
  use_scheme      <- if (scheme == "auto") detected_scheme else scheme

  cells   <- .trendcell_fit_cells(data, y, x, time, id, comparisons, cell_map, cutoff, kernel, h,
                                  bwselect, min_n, bc, p, q)
  stacked <- .trendcell_stack_contrasts(cells$fits, cells$meta, comparisons, trend, use_scheme, bc)
  Delta_all <- stacked$Delta
  Sigma_all <- stacked$Sigma
  jump_df   <- .trendcell_jump_table(stacked$jump_rows)

  if (length(cells$skipped))
    message("rd_trendcell: skipped ", length(cells$skipped), " cell(s) with fewer than min_n = ",
            min_n, " observations on a side or a failed fit: ",
            paste(cells$skipped, collapse = "; "))

  if (length(Delta_all) == 0L) {
    if (trend == "linear")
      return(.trendcell_untestable(jump_df, use_scheme, bc, estimand, trend, comparisons, cl))
    stop("rd_trendcell: no usable within-cell cross-period contrasts found; ",
         "check data, bandwidth, min_n, or number of comparison periods.")
  }

  # Sigma_all is built from estimated covariances and can be numerically indefinite, so the Wald
  # statistic uses only its positive eigen-directions (.wald_eigen(), assumption_tests_helpers.R)
  wald <- .wald_eigen(Delta_all, Sigma_all, scale = max(abs(data[[y]]), na.rm = TRUE))
  if (isTRUE(wald$degenerate))
    stop("rd_trendcell: the covariance of the contrasts is zero to working precision (the ",
         "outcome has no residual variation within the cells); the test is undefined.")
  if (wald$df == 0L)
    stop("rd_trendcell: estimated covariance matrix is numerically zero.")

  structure(
    list(statistic         = wald$stat,
         df                = wald$df,
         p_value           = wald$p,
         cell_period_jumps = jump_df,
         contrasts         = Delta_all,
         cov_matrix        = Sigma_all,
         scheme            = use_scheme,
         bc                = bc,
         estimand          = estimand,
         trend             = trend,
         comparisons       = as.character(comparisons),
         call              = cl),
    class = "rd_trendcell"
  )
}

#' @export
print.rd_trendcell <- function(x, ...) {
  null_text <- if (x$trend == "linear")
    "within each type, the confounding jump moves linearly across the comparison periods"
  else
    "within each type, the confounding jump is the same in every comparison period"
  .print_test_header("a constant within-type confounding discontinuity", "rd_trendcell",
                     null_text, x$scheme, TRUE, x$estimand)
  cat(sprintf("  Comparison periods: %s   Trend: %s\n\n", paste(x$comparisons, collapse = ", "),
              x$trend))
  .print_wald(x$statistic, x$df, x$p_value, label = "Wald")
  if (is.na(x$statistic))
    cat("    (a linear within-type trend needs at least 3 comparison periods to be testable)\n")
  .print_jump_table(x$cell_period_jumps, c("cell", "period"), c("Type", "Period"))
  invisible(x)
}
