# rd_homog.R -- rd_homog(): test of homogeneous confounding (equal confounding jump across types
# within each comparison period), with its print method. Shared pieces: assumption_tests_helpers.R,
# sampling_scheme.R, cross_period_covariance.R.

# ---------------------------------------------------------------------------
# Steps of rd_homog(), in the order it calls them
# ---------------------------------------------------------------------------

#' Per-period type tables in which a unit's type is its side of the cutoff in the RD period.
#' @noRd
.homog_rd_side_types <- function(type_list, wide, side_col) {
  # One binary partition (above / below in the RD period) instead of the sign pattern of all
  # other periods: when the comparison periods share a running variable, a unit's side in
  # another comparison period is collinear with the running variable and carries no extra type
  # information. Each comparison period then gives one contrast.
  side_map <- stats::setNames(wide[[side_col]], as.character(wide$id))
  for (tp in names(type_list)) {
    types_tp      <- type_list[[tp]]
    types_tp$type <- unname(side_map[as.character(types_tp$id)])
    # units not observed in the RD period have no type
    type_list[[tp]] <- types_tp[!is.na(types_tp$type), , drop = FALSE]
  }
  type_list
}

#' Sampling scheme ("cs", "pc" or "pv") read off the comparison periods by .detect_scheme().
#' @noRd
.homog_detect_scheme <- function(data, x, time, id, comparisons, cutoff) {
  # Only the comparison periods enter: the scheme here sets the covariance between
  # comparison-period jumps. rddid() reads it off every period, so the two can differ.
  comp_plist <- stats::setNames(lapply(as.character(comparisons), function(tp) {
    d_cp <- data[data[[time]] == tp, , drop = FALSE]
    list(id = d_cp[[id]], x = d_cp[[x]])
  }), as.character(comparisons))
  .detect_scheme(comp_plist, c = cutoff)
}

#' Local-linear fit of every (comparison period, type) cell: list(fits, meta, contrast_keys,
#' skipped), with `fits` and `meta` keyed "period::type".
#' @noRd
.homog_fit_cells <- function(data, y, x, time, id, comparisons, type_list, cutoff, kernel,
                             h, bwselect, min_n, bc, p, q) {
  poly_orders   <- list(p = p, q = q)
  skipped       <- character()
  fits          <- list()   # "period::type" -> rd_period fit
  meta          <- list()   # "period::type" -> (period, type, ref, D, V_D, n)
  contrast_keys <- list()   # period -> (reference key, non-reference keys)

  for (tp in as.character(comparisons)) {
    rows     <- which(data[[time]] == tp)
    d_tp     <- data[rows, , drop = FALSE]
    types_tp <- type_list[[tp]]
    id_tp    <- d_tp[[id]]
    m_idx    <- match(id_tp, types_tp$id)
    row_type <- types_tp$type[m_idx]

    # method = "radix" sorts by byte value on every machine ("+" is 43, "-" is 45), so the
    # decreasing order puts the all-below pattern ("-", "--", ...) first whatever the locale
    valid_types <- sort(unique(row_type[!is.na(row_type)]), method = "radix",
                        decreasing = TRUE)
    if (length(valid_types) < 2L) {
      message("rd_homog: period ", tp, " has fewer than 2 types; skipping.")
      next
    }

    # Reference type = the all-below pattern (below the cutoff in every other period, or in
    # t_rd under type_by = "rd_side"). Equal jumps across J types are J - 1 restrictions, so a
    # period contributes (type - reference) for its other types only.
    ref_type <- valid_types[1L]
    for (type_v in valid_types) {
      in_cell <- !is.na(row_type) & row_type == type_v
      y_cell  <- d_tp[[y]][in_cell]
      x_cell  <- d_tp[[x]][in_cell]
      id_cell <- d_tp[[id]][in_cell]

      # the bandwidth comes before the min_n check, so a CCT fallback message can come from a
      # cell that is then skipped
      bw   <- .cell_bandwidth(y_cell, x_cell, cutoff, kernel, h, bwselect)
      bw_h <- bw[["h"]]
      bw_b <- bw[["b"]]

      n_above <- sum(x_cell >= cutoff, na.rm = TRUE)
      n_below <- sum(x_cell <  cutoff, na.rm = TRUE)
      if (n_above < min_n || n_below < min_n) {
        skipped <- c(skipped, sprintf("period %s, type %s (n = %d below, %d above)",
                                      tp, type_v, n_below, n_above))
        next
      }

      # Any error in the fit (too few points in a bandwidth window, a singular local design, an
      # unusable bandwidth) skips the cell instead of stopping the test; the cell is then listed
      # in the skipped-cells message.
      fit <- tryCatch({
        call_args <- base::c(list(y = y_cell, x = x_cell, h = bw_h, b = bw_b,
                                  id = id_cell, c = cutoff, kernel = kernel), poly_orders)
        do.call(rd_period, call_args)
      }, error = function(e) NULL)
      if (is.null(fit)) {
        skipped <- c(skipped, sprintf("period %s, type %s (local-linear fit failed)",
                                      tp, type_v))
        next
      }

      key <- paste0(tp, "::", type_v)
      fits[[key]] <- fit
      meta[[key]] <- list(period = tp, type = type_v, ref = (type_v == ref_type),
                          D = if (bc) fit$D_bc else fit$D,
                          V_D = if (bc) fit$V_D_bc else fit$V_D, n = fit$n)
    }

    keys_tp <- .homog_contrast_keys(tp, valid_types, fits)
    if (!is.null(keys_tp)) contrast_keys[[tp]] <- keys_tp
  }

  list(fits = fits, meta = meta, contrast_keys = contrast_keys, skipped = skipped)
}

#' Keys of one period's contrasts, list(ref, non_ref), or NULL when fewer than two of its types
#' were fitted.
#' @noRd
.homog_contrast_keys <- function(tp, valid_types, fits) {
  is_fitted <- vapply(valid_types, function(type_v) {
    paste0(tp, "::", type_v) %in% names(fits)
  }, logical(1))
  fitted_types <- valid_types[is_fitted]
  if (length(fitted_types) < 2L) return(NULL)
  # The reference is the first FITTED type: if the all-below cell was skipped, the next type
  # stands in, while the jump table's `reference` column marks only the all-below type.
  ref_fitted <- paste0(tp, "::", fitted_types[1L])
  if (!ref_fitted %in% names(fits)) return(NULL)
  non_ref <- fitted_types[-1L]
  list(
    ref = ref_fitted,
    non_ref = paste0(tp, "::", non_ref)
  )
}

#' Keys of the non-reference cells of every period, stacked: the order of the contrast vector.
#' @noRd
.homog_contrast_entries <- function(contrast_keys) {
  do.call(base::c, lapply(names(contrast_keys), function(tp) {
    contrast_keys[[tp]]$non_ref
  }))
}

#' Map from each contrast key to the reference key it is differenced against.
#' @noRd
.homog_ref_map <- function(contrast_keys) {
  do.call(base::c, lapply(names(contrast_keys), function(tp) {
    stats::setNames(rep(contrast_keys[[tp]]$ref, length(contrast_keys[[tp]]$non_ref)),
                    contrast_keys[[tp]]$non_ref)
  }))
}

#' Covariance matrix of the stacked contrasts, rows and columns named by contrast key.
#' @noRd
.homog_sigma <- function(contrast_entries, ref_map, fits, meta, use_scheme, bc) {
  K <- length(contrast_entries)

  # covariance of the jumps of two cells
  cov_dd <- function(keyA, keyB) {
    fitA <- fits[[keyA]]
    fitB <- fits[[keyB]]
    same_period <- (meta[[keyA]]$period == meta[[keyB]]$period)
    if (keyA == keyB) return(if (bc) fitA$V_D_bc else fitA$V_D)
    # two types of one period are disjoint sets of units
    if (same_period) return(0)
    # across periods the scheme decides: 0 under "cs"; under "pc" the same-side covariance of
    # the units in both periods; under "pv" that minus the opposite-side covariance (switchers)
    .cov_scheme(fitA, fitB, use_scheme, bc = bc)
  }

  # contrast i is D[key_i] - D[ref_i], so
  #   Sigma[i, j] = Cov(D[key_i], D[key_j]) - Cov(D[key_i], D[ref_j])
  #               - Cov(D[ref_i], D[key_j]) + Cov(D[ref_i], D[ref_j])
  Sigma <- matrix(NA_real_, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      key_i <- contrast_entries[i]
      ref_i <- ref_map[key_i]
      key_j <- contrast_entries[j]
      ref_j <- ref_map[key_j]
      Sigma[i, j] <- cov_dd(key_i, key_j) - cov_dd(key_i, ref_j) -
                     cov_dd(ref_i, key_j) + cov_dd(ref_i, ref_j)
    }
  }
  rownames(Sigma) <- colnames(Sigma) <- contrast_entries
  Sigma
}

#' Jump table: one row per fitted (period, type) cell with period, type, jump, se, n, reference.
#' @noRd
.homog_jump_table <- function(meta) {
  jump_df <- do.call(rbind, lapply(names(meta), function(key) {
    cell <- meta[[key]]
    data.frame(period    = cell$period,
               type      = cell$type,
               jump      = cell$D,
               se        = sqrt(cell$V_D),
               n         = cell$n,
               reference = cell$ref,
               stringsAsFactors = FALSE)
  }))
  rownames(jump_df) <- NULL
  jump_df
}

# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Test of homogeneous confounding
#'
#' When the running variable moves over time, some units are above the cutoff
#' in one period and below it in another. A unit's **type** is the side of the
#' cutoff it is on in the other period(s); by default here, its side in the RD
#' period. `rd_homog()` tests the null that **in each comparison period the
#' confounding jump is the same for every type**, using the comparison periods
#' only, where the jump in the outcome is the confounding jump. Homogeneous
#' confounding and composition stability ([rd_compstable()]) are alternatives:
#' the estimate of [rddid()] needs one of the two (together with type
#' continuity, [rd_typecont()]). A rejection here alone therefore does not
#' invalidate the estimate; if composition stability is rejected as well, the
#' comparison periods mix the types differently from the RD period and the
#' estimate can be biased. With a running variable fixed over time (as in
#' [rddid_sim]) the types are degenerate and the test is not informative.
#'
#' @details
#' ## What is estimated
#'
#' In each comparison period the units are split by type, and a local-linear
#' RD of the outcome on the running variable within each type gives that
#' type's confounding jump. Within each period every type's jump is compared
#' with the jump of a reference type (the type below the cutoff in the other
#' period(s)), and these differences are tested jointly across the comparison
#' periods by a Wald statistic. The null is equality across types within each
#' period, not equality across periods (that is [rd_trendcell()]). The
#' covariance of the differences is estimated and can be numerically
#' indefinite, so the Wald statistic uses only its positive directions, and
#' `df` counts them.
#'
#' ## Shared units and the sampling scheme
#'
#' Within a period the types are different units, so their jumps are
#' independent. Across comparison periods the same units can appear in both,
#' and the covariance follows `scheme`, matching units on `id`: none under
#' `"cs"`; from the units on the same side of the cutoff in both periods under
#' `"pc"`; under `"pv"` also from the units that change side, with the opposite
#' sign. `"auto"` reads the scheme off the comparison periods, by the rule of
#' [rddid()].
#'
#' ## Options
#'
#' `type_by = "rd_side"` (default) types each unit by its side of the cutoff in
#' the RD period, so `t_rd` is required and each comparison period gives one
#' difference; units not observed in the RD period are left out.
#' `type_by = "pattern"` types units by their sides in all other periods in
#' `data`. A (period, type) cell with fewer than `min_n` observations on either
#' side of the cutoff is dropped, and a message lists the dropped cells; a
#' period with fewer than two types is skipped. `bc = TRUE` (default) tests the
#' bias-corrected jumps with their robust variance, as in the `Robust` row of
#' [rddid()]; `bc = FALSE` uses the conventional jumps and variances. With
#' `bwselect = "cct"` (default) each cell gets its own CCT bandwidths from
#' [rd_bw_cct()], computed on that cell's outcome and running variable. With
#' `bwselect = "rot"` the rule of thumb is `h = b = 0.2` times the range of the
#' running variable within each cell. A numeric `h` is used as both bandwidths
#' in every cell.
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
#' @param comparisons the comparison periods in which the test runs. `NULL`
#'   (default) uses every period other than `t_rd` (every period if `t_rd` is
#'   `NULL`).
#' @param estimand `"att"` (default) or `"atu"`, as in [rddid()]. Label only:
#'   the test is the same either way.
#' @param h a bandwidth to use, as both main and pilot bandwidth, in every cell.
#'   If given, `bwselect` is ignored.
#' @param bwselect the bandwidth rule when `h` is not given: `"cct"` (default;
#'   each cell's own CCT bandwidths from [rd_bw_cct()]) or `"rot"` (the rule of
#'   thumb `0.2` times the range of the running variable within the cell).
#' @param min_n the minimum number of observations on each side of the cutoff
#'   for a (period, type) cell to enter the test (default 10). Smaller cells
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
#'   sides in all other periods in `data`).
#' @param p,q orders of the local polynomials in every cell, for the point
#'   estimate and the bias correction (defaults 1 and 2; `q` must exceed `p`).
#'   The CCT bandwidths are always chosen for a local-linear fit.
#' @inheritParams rddid
#'
#' @return An object of class `"rd_homog"`, a list with:
#'   \describe{
#'     \item{`statistic`, `df`, `p_value`}{the Wald statistic, its degrees of
#'       freedom (the number of positive directions of the covariance used) and
#'       its chi-squared p-value.}
#'     \item{`scheme`}{the sampling scheme used.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`call`}{the matched call.}
#'     \item{`period_type_jumps`}{a data frame with one row per (comparison
#'       period, type) cell that was fitted: `period`, `type` (the unit's
#'       side(s), `"+"` above and `"-"` below the cutoff), `jump` (the cell's
#'       confounding jump,
#'       bias-corrected when `bc = TRUE`), `se`, `n` (observations in the cell)
#'       and `reference` (`TRUE` for the reference type).}
#'     \item{`contrasts`}{named numeric vector of the tested differences (type
#'       minus reference), stacked across periods.}
#'     \item{`cov_matrix`}{the estimated covariance matrix of `contrasts`.}
#'     \item{`bc`}{as passed.}
#'     \item{`comparisons`}{the comparison periods that contributed a
#'       difference.}
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
#' hc <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' hc          # the Wald test, then each type's jump in each comparison period
#' hc$period_type_jumps
#' @export
rd_homog <- function(data, y, x, time, id,
                     t_rd = NULL, comparisons = NULL,
                     estimand = c("att", "atu"),
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
  bwselect <- match.arg(bwselect)
  estimand <- match.arg(estimand)
  kernel   <- match.arg(kernel, c("triangular", "epanechnikov", "uniform"))
  cutoff   <- c   # the cutoff; `c` stays the argument name for rdrobust users

  for (nm in base::c(y, x, time, id))
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")

  times_all <- sort(unique(data[[time]]))
  if (is.null(comparisons)) {
    comparisons <- if (!is.null(t_rd)) setdiff(times_all, t_rd) else times_all
  }
  if (length(comparisons) < 1L)
    stop("need at least one comparison period.")

  # types are "+"/"-" sign-pattern strings (a unit at the cutoff counts as above); a unit
  # unobserved in a period its type needs is dropped from that period
  types     <- .build_types(data, x = x, time = time, id = id, c = cutoff)
  type_list <- types$period_types
  if (type_by == "rd_side") {
    if (is.null(t_rd))
      stop("type_by = \"rd_side\" requires `t_rd` (the RD period whose side defines the type).")
    side_col <- paste0("side_", t_rd)
    if (!side_col %in% names(types$wide))
      stop("RD period '", t_rd, "' has no side column; cannot define rd_side types.")
    type_list <- .homog_rd_side_types(type_list, types$wide, side_col)
  }

  detected_scheme <- .homog_detect_scheme(data, x, time, id, comparisons, cutoff)
  cells <- .homog_fit_cells(data, y, x, time, id, comparisons, type_list, cutoff, kernel,
                            h, bwselect, min_n, bc, p, q)
  use_scheme <- if (scheme == "auto") detected_scheme else scheme

  contrast_entries <- .homog_contrast_entries(cells$contrast_keys)
  if (length(contrast_entries) == 0L)
    stop("rd_homog: no usable type contrasts found; check data, bandwidth, or min_n.")
  ref_map <- .homog_ref_map(cells$contrast_keys)
  # meta$D is already the bias-corrected jump when bc = TRUE
  Delta <- vapply(contrast_entries, function(key) {
    cells$meta[[key]]$D - cells$meta[[ref_map[key]]]$D
  }, numeric(1))
  Sigma <- .homog_sigma(contrast_entries, ref_map, cells$fits, cells$meta, use_scheme, bc)

  # Deliberately not .joint_wald(): Sigma is a difference of estimated covariances and can
  # come back numerically indefinite, so .wald_eigen() keeps only its positive directions
  # (eigenvalues above its relative tolerance); df = 0 means none is left.
  wald <- .wald_eigen(Delta, Sigma)
  if (wald$df == 0L) stop("rd_homog: estimated covariance matrix is numerically zero.")
  wald_stat <- wald$stat
  wald_df   <- wald$df
  p_value   <- wald$p

  jump_df <- .homog_jump_table(cells$meta)

  if (length(cells$skipped))
    message("rd_homog: skipped ", length(cells$skipped), " cell(s) with fewer than min_n = ",
            min_n, " observations on a side or a failed fit: ",
            paste(cells$skipped, collapse = "; "))

  structure(
    list(
      statistic         = wald_stat,
      df                = wald_df,
      p_value           = p_value,
      period_type_jumps = jump_df,
      contrasts         = Delta,
      cov_matrix        = Sigma,
      scheme            = use_scheme,
      bc                = bc,
      estimand          = estimand,
      comparisons       = names(cells$contrast_keys),
      call              = cl
    ),
    class = "rd_homog"
  )
}

#' @export
print.rd_homog <- function(x, ...) {
  .print_test_header("homogeneous confounding", "rd_homog",
                     "in each comparison period the confounding jump is the same for every type",
                     x$scheme, TRUE, x$estimand)
  cat(sprintf("  Comparison periods: %s\n\n", paste(x$comparisons, collapse = ", ")))
  .print_wald(x$statistic, x$df, x$p_value, label = "Wald")
  .print_jump_table(x$period_type_jumps, c("period", "type"), c("Period", "Type"))
  invisible(x)
}
