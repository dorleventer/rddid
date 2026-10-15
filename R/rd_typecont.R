# rd_typecont.R -- rd_typecont(): test of a continuous type distribution at the cutoff, with its
# print method and its internal steps .typecont_*(). Shared pieces: assumption_tests_helpers.R
# (types, Wald, per-cell bandwidth), sampling_scheme.R, cross_period_covariance.R.
#
# Layout of the stacked jumps: cell (type idx_v, period idx_t) is entry
# (idx_t - 1) * n_types + idx_v of `theta` and row/column of `Sigma` (the type index runs fastest).


#' The sampling scheme used: `scheme` itself, or under "auto" the one ("cs", "pc" or "pv") read
#' off the units that have a type in each period.
#' @noRd
.typecont_scheme <- function(scheme, period_types, period_labels, cutoff) {
  if (scheme != "auto") return(scheme)
  # The scheme is a property of the design (repeated ids, side switching), so it is read from
  # every unit with a type, not from a bandwidth window: a window misses the switchers outside
  # it, and per-cell CCT bandwidths can reach beyond any common window. Only units with a type
  # enter (observed in every period), unlike rddid(), which uses every row. A unit on the
  # cutoff counts as above it.
  long <- do.call(rbind, lapply(period_labels, function(k) {
    units_k <- period_types[[k]]
    data.frame(period = k,
               id     = units_k$id,
               side   = as.integer(units_k$R >= cutoff))
  }))
  .scheme_from_long(long)
}

#' Fits of every (period, type) cell: list(fits = type x period list-matrix of rd_period fits,
#' NULL where the fit failed; theta = the stacked jumps, NA where it failed).
#' @noRd
.typecont_fit_cells <- function(period_types, period_labels, type_values, n_periods, n_types,
                                cutoff, kernel, h, bwselect, bc) {
  fits <- vector("list", n_periods * n_types)
  dim(fits) <- c(n_types, n_periods)
  dimnames(fits) <- list(type   = as.character(type_values),
                         period = period_labels)

  theta <- numeric(n_periods * n_types)
  idx   <- 0L

  for (idx_t in seq_along(period_labels)) {
    units_t <- period_types[[period_labels[idx_t]]]
    for (idx_v in seq_along(type_values)) {
      idx    <- idx + 1L
      type_v <- type_values[idx_v]
      y_v    <- as.numeric(units_t$type == type_v)   # the RD outcome: "the unit is of type v"
      # h is NULL only under bwselect = "cct", where each cell gets its own CCT bandwidths;
      # otherwise it is the user's h or the rule of thumb, the same in every cell
      bw     <- .cell_bandwidth(y_v, units_t$R, cutoff, kernel, h, bwselect)
      h_cell <- bw[["h"]]
      b_cell <- bw[["b"]]
      # A cell whose fit fails is left out of the tests (NULL fit, NA jump). The usual cause is
      # too few units on one side within the bandwidths (a type that is rare near the cutoff);
      # any other rd_period error in the cell is hidden the same way.
      fit <- tryCatch(
        rd_period(y = y_v, x = units_t$R, h = h_cell, b = b_cell, id = units_t$id,
                  c = cutoff, p = 1L, q = 2L, kernel = kernel),
        error = function(e) NULL
      )
      fits[idx_v, idx_t] <- list(fit)   # `[<-` with list() can store a NULL fit; `[[<-` cannot
      theta[idx] <- if (is.null(fit)) NA_real_ else if (bc) fit$D_bc else fit$D
    }
  }
  list(fits = fits, theta = theta)
}

#' Covariance matrix of the stacked jumps theta, (n_periods * n_types) square, with zero rows
#' and columns for the cells whose fit failed.
#' @noRd
.typecont_sigma <- function(fits, period_labels, type_values, n_periods, n_types, scheme, bc) {
  n_cells <- n_periods * n_types
  Sigma   <- matrix(0, n_cells, n_cells)

  for (idx_t1 in seq_along(period_labels)) {
    for (idx_v1 in seq_along(type_values)) {
      row_idx <- (idx_t1 - 1L) * n_types + idx_v1
      fit_row <- fits[[idx_v1, idx_t1]]
      if (is.null(fit_row)) next

      for (idx_t2 in seq_along(period_labels)) {
        for (idx_v2 in seq_along(type_values)) {
          col_idx <- (idx_t2 - 1L) * n_types + idx_v2
          fit_col <- fits[[idx_v2, idx_t2]]
          if (is.null(fit_col)) next

          if (idx_t1 == idx_t2) {
            # Same period: both regressions use the same units, each on one side of the
            # cutoff, so only the same-side term enters, whatever the scheme.
            Sigma[row_idx, col_idx] <- .cross_cov(fit_row, fit_col, bc = bc)$pc
          } else {
            # Across periods the scheme decides: "cs" none; "pc" the units on the same side in
            # both periods; "pv" also the units that change side, with the opposite sign.
            Sigma[row_idx, col_idx] <- .cov_scheme(fit_row, fit_col, scheme, bc = bc)
          }
        }
      }
    }
  }
  Sigma
}

#' Logical vector over the entries of theta: TRUE for the jumps the tests use (in each period,
#' every present type but the last, the reference).
#' @noRd
.typecont_keep_rows <- function(theta, period_labels, n_types) {
  # The type indicators sum to one within a period, so their jumps sum to zero and Sigma is
  # singular by construction. A pseudo-inverse of it is fragile (a structural-zero singular
  # value survives the cut on some LAPACK builds and its 1/sv inflates the statistic), so one
  # reference type per period is dropped instead and the other k - 1 present types are kept.
  # At a common bandwidth this is the same test, the dropped jump being minus the sum of the
  # rest; with per-cell CCT bandwidths the jumps do not sum exactly to zero, so it is a
  # different, still valid, test. Types are in radix order, so the dropped (last present) type
  # is the all-below pattern when present, and with two periods the kept contrast is the
  # paper's jump in the indicator 1{V_is = 1}.
  keep <- logical(length(theta))
  for (idx_t in seq_along(period_labels)) {
    rows_t  <- (idx_t - 1L) * n_types + seq_len(n_types)
    present <- rows_t[!is.na(theta[rows_t])]
    if (length(present) >= 2L) keep[present[-length(present)]] <- TRUE
  }
  keep
}

#' Wald test of the jumps in `rows` of theta: list(stat, df, p); stat 0, df 0, p 1 if none.
#' @noRd
.typecont_wald <- function(theta, Sigma, rows) {
  if (length(rows) == 0L) {
    return(list(stat = NA_real_, df = 0L, p = NA_real_))   # nothing testable
  }
  .joint_wald(theta[rows], Sigma[rows, rows, drop = FALSE])
}

#' Each period's own tests on its kept jumps: a list by period of list(ll_wald = the Wald test,
#' test = the test under `joint`).
#' @noRd
.typecont_per_period <- function(theta, Sigma, keep, period_labels, n_types, joint) {
  # With two types this is the single type-share jump of the period, chi^2(1) (sup-t and sWald
  # are then the same two-sided normal test). The joint test is not the sum of these: it also
  # carries the covariance across periods.
  stats::setNames(lapply(seq_along(period_labels), function(idx_t) {
    rows_t <- (idx_t - 1L) * n_types + seq_len(n_types)
    keep_t <- rows_t[keep[rows_t]]
    tests  <- .joint_tests(theta[keep_t], Sigma[keep_t, keep_t, drop = FALSE])
    list(ll_wald = .typecont_wald(theta, Sigma, keep_t), test = tests[[joint]])
  }), period_labels)
}

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
#' other period), and the remaining jumps are tested jointly (see "The joint
#' statistic"). Each period's own test is reported as well. A unit's type in a
#' period needs its side in every other period, so a unit missing from some
#' period is left out of the regressions of the other periods. The test treats
#' all periods alike: `t_rd` and `comparisons` only select which periods enter.
#'
#' ## The joint statistic
#'
#' `joint = "supt"` (default) standardizes each jump by its standard error and
#' rejects when the largest |z| exceeds the sup-t critical value: the quantile
#' of max |Z_k| for Z normal with the estimated correlation of the jumps. The
#' printout names the (period, type) cell with the largest |z|. `"swald"` sums
#' the squared z's and compares the sum with its distribution under the same
#' normal (the sWald test of Fusejima, Ishihara and Sawada, 2025). `"wald"` is
#' the chi-squared Wald statistic, which inverts the estimated covariance. Both
#' sup-t and sWald use the covariance only through its correlations, so they
#' avoid the error of inverting it, which grows fastest with the number of
#' jumps. They still divide each jump by its estimated standard error, so a
#' variance that is too small, or noisy, inflates them too, but less: at
#' n = 1,000 the Wald rejects a true null 7-9% of the time at 9 jumps (three
#' periods) and 15% at 28 (four periods), sup-t 6-7% and 4-5%, in the package
#' authors' simulations. The critical values of sup-t and sWald are simulated with 100,000
#' normal draws under a fixed seed, so a result is reproducible, and the
#' caller's random-number state is left as it was. With one jump all three are
#' the same two-sided normal test. All three are always returned in
#' `joint_tests`.
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
#' @param joint the joint statistic: `"supt"` (default; the largest |z| against
#'   a sup-t critical value), `"swald"` (the sum of squared z's) or `"wald"` (the
#'   chi-squared Wald statistic). See "The joint statistic" in Details.
#' @inheritParams rddid
#'
#' @return An object of class `"rd_typecont"`, a list with:
#'   \describe{
#'     \item{`statistic`, `df`, `p_value`}{the joint test over all periods under
#'       `joint`: the statistic (the largest |z| for `"supt"`, the sum of squared
#'       z's for `"swald"`, the chi-squared statistic for `"wald"`), the number of
#'       jumps tested (`"supt"`, `"swald"`) or the degrees of freedom (`"wald"`),
#'       and the p-value.}
#'     \item{`joint_method`}{the joint statistic used; `joint_tests` has all three
#'       (a data frame with `method`, `statistic`, `df`, `p_value`).}
#'     \item{`z`}{the z-statistics of the jumps tested, named `"<period>:<type>"`.}
#'     \item{`scheme`}{the sampling scheme used; `scheme_requested` is the
#'       argument as passed.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`fits`}{the per-cell [rd_period()] fits, a list-matrix indexed by type and
#'       period (`NULL` where a cell could not be fitted).}
#'     \item{`data`}{the typed data by period: for each period a data frame with `id`, `R`
#'       (the running variable) and `type`.}
#'     \item{`sides`}{one row per unit with its running variable (`R_<period>`) and side of
#'       the cutoff (`side_<period>`, `"+"`/`"-"`) in every period; [plot.rd_typecont()]
#'       reads it.}
#'     \item{`call`}{the matched call.}
#'     \item{`per_period`}{a list by period; each element holds `test`, that
#'       period's own test under `joint` (`stat`, `df`, `p`), and `ll_wald`, its
#'       Wald test.}
#'     \item{`ll_wald`}{the joint Wald test, as a list (`stat`, `df`, `p`). A type
#'       share that is constant near the cutoff on both sides with no jump carries no
#'       information and is left out (`dropped` counts them); one that jumps
#'       deterministically makes the test undefined (`NA`).}
#'     \item{`meta`}{a list with `periods`, `type_values` (the types, written
#'       as the sides in the other periods in time order, e.g. `"+-"`), `h`
#'       (the common bandwidth, `NA` with `bwselect = "cct"`), `bwselect`,
#'       `scheme`, `bc`, `estimand` and `joint`.}
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
#' tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
#' tc          # the null, the joint sup-t test, then each period's own test
#' tc$p_value
#' tc$joint_tests   # sup-t, sWald and the chi-squared Wald side by side
#' @export
rd_typecont <- function(data, x, time, id,
                        t_rd = NULL, comparisons = NULL,
                        estimand = c("att", "atu"),
                        c = 0,
                        h = NULL,
                        bwselect = c("cct", "rot"),
                        kernel = "triangular",
                        scheme = c("auto", "cs", "pc", "pv"),
                        bc = TRUE,
                        joint = c("supt", "swald", "wald")) {
  cl       <- match.call()
  scheme   <- match.arg(scheme)
  joint    <- match.arg(joint)
  bwselect <- match.arg(bwselect)
  estimand <- match.arg(estimand)
  kernel   <- match.arg(kernel, c("triangular", "epanechnikov", "uniform"))
  cutoff   <- c   # the cutoff; `c` stays the argument name for rdrobust users

  # ----- inputs and periods (the stops stay here so that errors name rd_typecont()) -----
  .check_columns(data, c(x, time, id))
  # `t_rd`/`comparisons` only select which periods enter: the test treats every period alike
  if (!is.null(t_rd) && !t_rd %in% data[[time]])
    stop("`t_rd` (", t_rd, ") is not a period in `data`.")
  if (!is.null(comparisons)) {
    use_periods <- c(t_rd, comparisons)
    if (!all(use_periods %in% data[[time]]))
      stop("periods not in `data`: ", paste(setdiff(use_periods, data[[time]]), collapse = ", "))
    data <- data[data[[time]] %in% use_periods, , drop = FALSE]
  }
  data          <- data[stats::complete.cases(data[, c(x, time, id)]), , drop = FALSE]
  periods       <- sort(unique(data[[time]]))
  period_labels <- as.character(periods)
  n_periods     <- length(periods)
  if (n_periods < 2L) stop("need at least 2 periods to define a type.")

  # bwselect = "cct" leaves h NULL: each cell then gets its own CCT bandwidths
  if (is.null(h) && bwselect == "rot") h <- .rot_bandwidth_iqr(data[[x]])

  # ----- types, scheme, fits, covariance -----
  types        <- .build_types(data, x, time, id, c = cutoff)
  .stop_if_no_switchers(types$wide, period_labels, "rd_typecont")
  period_types <- types$period_types
  # radix: a locale-independent order, so the reference type dropped is the same on every machine
  type_values <- sort(unique(unlist(lapply(period_types, `[[`, "type"))), method = "radix")
  n_types     <- length(type_values)
  use_scheme  <- .typecont_scheme(scheme, period_types, period_labels, cutoff)

  cells <- .typecont_fit_cells(period_types, period_labels, type_values, n_periods, n_types,
                               cutoff, kernel, h, bwselect, bc)
  theta <- cells$theta
  Sigma <- .typecont_sigma(cells$fits, period_labels, type_values, n_periods, n_types,
                           use_scheme, bc)

  # ----- tests: joint over all periods, then each period's own -----
  keep        <- .typecont_keep_rows(theta, period_labels, n_types)
  rows        <- which(keep)
  wald_joint  <- .typecont_wald(theta, Sigma, rows)
  joint_tests <- .joint_tests(theta[rows], Sigma[rows, rows, drop = FALSE])
  chosen      <- joint_tests[[joint]]
  per_period  <- .typecont_per_period(theta, Sigma, keep, period_labels, n_types, joint)
  # z of each tested jump, named period:type (the type index runs fastest in theta)
  cell_names  <- paste0(period_labels[(rows - 1L) %/% n_types + 1L], ":",
                        type_values[(rows - 1L) %% n_types + 1L])
  z           <- stats::setNames(theta[rows] / sqrt(pmax(diag(Sigma)[rows], 0)), cell_names)

  structure(
    list(
      statistic  = chosen$stat,
      df         = chosen$df,
      p_value    = chosen$p,
      joint_method = joint,
      joint_tests = .joint_tests_df(joint_tests),
      z          = z,
      scheme     = use_scheme,
      scheme_requested = scheme,
      estimand   = estimand,
      ll_wald    = wald_joint,
      per_period = per_period,
      fits       = cells$fits,
      data       = period_types,
      sides      = types$wide,
      call       = cl,
      meta = list(
        periods      = period_labels,
        type_values  = type_values,
        c            = cutoff,
        t_rd         = t_rd,
        kernel       = kernel,
        h            = if (!is.null(h)) h else NA_real_,
        bwselect     = bwselect,
        scheme       = use_scheme,
        bc           = bc,
        estimand     = estimand,
        joint        = joint
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
              paste(x$meta$periods, collapse = ", "),
              paste(x$meta$type_values, collapse = ", "),
              .bw_label_test(x$meta$h, x$meta$bwselect)))
  method <- x$joint_method %||% "wald"
  largest <- if (method == "supt" && length(x$z) && any(is.finite(x$z)))
    names(x$z)[which.max(abs(x$z))] else NULL
  .print_joint(list(stat = x$statistic, df = x$df, p = x$p_value), method, "Joint",
               cell_label = if (!is.null(largest)) sprintf("period:type %s", largest))
  for (k in names(x$per_period)) {
    period_test <- x$per_period[[k]]$test %||% x$per_period[[k]]$ll_wald
    if (!is.null(period_test)) {
      .print_joint(period_test, if (is.null(x$per_period[[k]]$test)) "wald" else method,
                   sprintf("Period %s:", k), indent = "    ")
    }
  }
  invisible(x)
}
