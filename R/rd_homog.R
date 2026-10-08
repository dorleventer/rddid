# rd_homog.R -- rd_homog(): test of homogeneous confounding (equal confounding jump across types
# within each comparison period), with its print method. Shared pieces: assumption_tests_helpers.R,
# sampling_scheme.R, cross_period_covariance.R.

# ---------------------------------------------------------------------------
# This test has no local helpers: types come from the shared .build_types()
# (R/test_helpers.R) and cross-period covariances from .cov_scheme()/.cross_cov()
# (R/aggregate.R).
# ---------------------------------------------------------------------------

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

  dots     <- list(p = p, q = q)   # polynomial orders of the per-cell fits
  skipped  <- character()

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

    # locale-independent order, all-below pattern ("-", "--", ...) first
    valid_types <- sort(unique(tvec[!is.na(tvec)]), method = "radix", decreasing = TRUE)
    if (length(valid_types) < 2L) {
      # only one type (or no types): skip this period
      message("rd_homog: period ", tp, " has fewer than 2 types; skipping.")
      next
    }

    # reference type = the all-below pattern (V = 0 in every other period, or
    # V_{t_rd} = 0 under type_by = "rd_side"); contrasts are (type) - (reference),
    # so with binary types contrasts = D(1) - D(0)
    ref_type <- valid_types[1L]
    for (vt in valid_types) {
      keep   <- !is.na(tvec) & tvec == vt
      y_vt   <- d_tp[[y]][keep]
      x_vt   <- d_tp[[x]][keep]
      id_vt  <- d_tp[[id]][keep]

      # per-cell bandwidth (0.2*range per cell for rot; CCT per cell for cct)
      bw   <- .cell_bandwidth(y_vt, x_vt, c, kernel, h, bwselect)
      bw_h <- bw[["h"]]
      bw_b <- bw[["b"]]

      # skip if too few obs on either side
      n_pos <- sum(x_vt >= c, na.rm = TRUE)
      n_neg <- sum(x_vt <  c, na.rm = TRUE)
      if (n_pos < min_n || n_neg < min_n) {
        skipped <- c(skipped, sprintf("period %s, type %s (n = %d below, %d above)", tp, vt, n_neg, n_pos))
        next
      }

      fit <- tryCatch({
        call_args <- base::c(list(y = y_vt, x = x_vt, h = bw_h, b = bw_b,
                                   id = id_vt, c = c, kernel = kernel), dots)
        do.call(rd_period, call_args)
      }, error = function(e) NULL)
      if (is.null(fit)) {
        skipped <- c(skipped, sprintf("period %s, type %s (local-linear fit failed)", tp, vt))
        next
      }

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
  # estimated covariances and can come back numerically indefinite.  Use the
  # conservative eigen pseudo-inverse that drops non-positive directions
  # (see .wald_eigen() in R/test_helpers.R for the implementation).
  ew <- .wald_eigen(Delta, Sigma)
  if (ew$df == 0L) stop("rd_homog: estimated covariance matrix is numerically zero.")
  W  <- ew$stat
  df <- ew$df
  pv <- ew$p

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

  if (length(skipped))
    message("rd_homog: skipped ", length(skipped), " cell(s) with fewer than min_n = ", min_n,
            " observations on a side or a failed fit: ", paste(skipped, collapse = "; "))

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
      estimand          = estimand,
      comparisons       = names(contrast_keys),
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
