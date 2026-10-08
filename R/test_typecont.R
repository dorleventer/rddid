#' Test of type continuity
#'
#' When the running variable moves over time, some units are above the cutoff
#' in one period and below it in another. A unit's **type** is the side of the
#' cutoff it is on in the other period(s); with two periods, "above in the
#' other period" or "below in the other period". `rd_typecont()` tests the
#' null that **the share of each type jumps by zero at the cutoff, in every
#' period**. A rejection means that units sort across the cutoff by type, so
#' the jump in the outcome at the cutoff can reflect who is on each side as
#' well as the policies, and the estimate of [rddid()] can be biased. With a
#' running variable fixed over time (as in [rddid_sim]) every unit's type is
#' its own side, the shares jump from 0 to 1 by construction, and the test is
#' not informative.
#'
#' @details
#' ## What is estimated
#'
#' In each period and for each type, a local-linear RD of the indicator "the
#' unit is of this type" on the running variable estimates the jump in that
#' type's share at the cutoff. The shares sum to one within a period, so one
#' reference type per period is dropped (the type below the cutoff in every
#' other period), and the remaining jumps are tested jointly by a Wald
#' statistic, chi-squared with as many degrees of freedom as independent jumps
#' tested. Each period's own Wald test is reported as well. A unit's type in a
#' period needs its side in every other period, so a unit missing from some
#' period is left out of the regressions of the other periods. The test treats
#' all periods alike: `t_rd` and `comparisons` only select which periods enter.
#'
#' ## Shared units and the sampling scheme
#'
#' Within a period the type-indicator regressions use the same units, so their
#' covariance always enters. Across periods it follows `scheme`, matching units
#' on `id`: none under `"cs"`; from the units on the same side of the cutoff in
#' both periods under `"pc"`; under `"pv"` also from the units that change
#' side, with the opposite sign. `"auto"` reads the scheme off the data as
#' [rddid()] does.
#'
#' ## Options
#'
#' `bc = TRUE` (default) tests the bias-corrected jumps with their robust
#' variance, as in the `Robust` row of [rddid()]; `bc = FALSE` uses the
#' conventional jumps and variances. With `bwselect = "cct"` (default) each
#' (period, type) regression gets its own CCT bandwidths from [rd_bw_cct()].
#' With `bwselect = "rot"` the rule of thumb is `h = b = 0.5 * IQR(x)`, the
#' interquartile range of the running variable over all periods used (`sd(x)`
#' if that is zero), the same in every regression. A numeric `h` is used as
#' both bandwidths in every regression.
#'
#' ## ATU designs
#'
#' The null treats the two sides of the cutoff alike, so the test is the same
#' under `estimand = "att"` and `"atu"`; `estimand` only labels the output.
#'
#' @param data a data frame in long format, one row per unit and period, from a
#'   panel (it need not be balanced). A unit's type is read from the periods in
#'   which it is observed; a unit missing from a period that its type needs is
#'   left out of the cells that use that type.
#' @param id name of the unit-identifier column (a string). Required: types are
#'   read across periods.
#' @param t_rd,comparisons optional: the RD period and the comparison periods,
#'   to restrict the test to these periods. With `comparisons = NULL` (default)
#'   every period in `data` enters, whatever `t_rd`; since the test treats all
#'   periods alike, `t_rd` alone changes nothing and only lets you write the
#'   same call as for [rddid()].
#' @param estimand `"att"` (default) or `"atu"`, as in [rddid()]. Label only:
#'   the test is the same either way.
#' @param h a bandwidth to use, as both main and pilot bandwidth, in every
#'   regression of the test. If given, `bwselect` is ignored.
#' @param bwselect the bandwidth rule when `h` is not given: `"cct"` (default;
#'   each regression's own CCT bandwidths from [rd_bw_cct()]) or `"rot"` (the
#'   rule of thumb `0.5 * IQR(x)`, the same in every regression).
#' @param scheme the sampling scheme, which sets the covariance across periods
#'   in the test: `"cs"`, `"pc"` or `"pv"` (as in [rddid()]), or `"auto"`
#'   (default), which reads it off the data as [rddid()] does. See Details.
#' @param bc logical. `TRUE` (default): test the bias-corrected jumps with their
#'   robust variance, as in the `Robust` row of [rddid()]; `FALSE`: the
#'   conventional jumps and variances.
#' @inheritParams rddid
#'
#' @return An object of class `"rd_typecont"`, a list with:
#'   \describe{
#'     \item{`statistic`, `df`, `p_value`}{the joint Wald statistic over all
#'       periods, its degrees of freedom and its chi-squared p-value.}
#'     \item{`scheme`}{the sampling scheme used; `scheme_requested` is the
#'       argument as passed.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`call`}{the matched call.}
#'     \item{`per_period`}{a list by period; each element holds `ll_wald`, that
#'       period's own Wald test (`stat`, `df`, `p`).}
#'     \item{`ll_wald`}{the joint test again, as a list (`stat`, `df`, `p`).}
#'     \item{`meta`}{a list with `periods`, `type_values` (the types, written
#'       as the sides in the other periods in time order, e.g. `"+-"`), `h`
#'       (the common bandwidth, `NA` with `bwselect = "cct"`), `bwselect`,
#'       `scheme`, `bc` and `estimand`.}
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
#' tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id")
#' tc          # the null, the joint Wald test, then each period's own test
#' tc$p_value
#' @export
rd_typecont <- function(data, x, time, id,
                        t_rd = NULL, comparisons = NULL,
                        estimand = c("att", "atu"),
                        c = 0,
                        h = NULL,
                        bwselect = c("cct", "rot"),
                        kernel = "triangular",
                        scheme = c("auto", "cs", "pc", "pv"),
                        bc = TRUE) {
  cl       <- match.call()
  scheme   <- match.arg(scheme)
  bwselect <- match.arg(bwselect)
  estimand <- match.arg(estimand)
  kernel   <- match.arg(kernel, c("triangular", "epanechnikov", "uniform"))

  # ----- input checks -------------------------------------------------------
  for (nm in c(x, time, id)) {
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")
  }
  # `t_rd`/`comparisons` only select which periods enter (the test treats every period
  # alike); with both NULL every period in `data` is used
  if (!is.null(comparisons)) {
    use_periods <- c(t_rd, comparisons)
    if (!all(use_periods %in% data[[time]]))
      stop("periods not in `data`: ", paste(setdiff(use_periods, data[[time]]), collapse = ", "))
    data <- data[data[[time]] %in% use_periods, , drop = FALSE]
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
      statistic  = ll_result$stat,
      df         = ll_result$df,
      p_value    = ll_result$p,
      scheme     = use_scheme,
      scheme_requested = scheme,
      estimand   = estimand,
      ll_wald    = ll_result,
      per_period = per_period,
      call       = cl,
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
  .print_test_header("a continuous type distribution", "rd_typecont",
                     "the share of each type jumps by zero at the cutoff, in every period",
                     x$scheme, identical(x$scheme_requested, "auto"), x$estimand)
  cat(sprintf("  Periods: %s   Types: %s   Bandwidth: %s\n\n",
              paste(x$meta$periods, collapse = ", "), paste(x$meta$type_values, collapse = ", "),
              .bw_label_test(x$meta$h, x$meta$bwselect)))
  .print_wald(x$statistic, x$df, x$p_value, label = "Joint Wald")
  for (k in names(x$per_period)) {
    pp <- x$per_period[[k]]$ll_wald
    if (!is.null(pp)) .print_wald(pp$stat, pp$df, pp$p, label = sprintf("Period %s:", k), indent = "    ")
  }
  invisible(x)
}
