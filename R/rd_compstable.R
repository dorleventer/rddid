# rd_compstable.R -- rd_compstable(): test of composition stability between the RD period and
# each comparison period (reflected-sample construction), with its print method. Shared pieces:
# assumption_tests_helpers.R, cross_period_covariance.R.


#' Wide table of the panel, one row per unit (column `id`): `R_<period>`, the running variable
#' (mirrored under "atu"), and `side_<period>`, 1 if the unit is above the original cutoff in that
#' period and 0 if below (NA if unobserved); returns the data frame.
#' @noRd
.compstable_wide <- function(data, x, time, id, all_periods, period_labels, cutoff, estimand) {
  wide <- data.frame(id = unique(data[[id]]), stringsAsFactors = FALSE)
  for (lab in period_labels) {
    rows_t <- data[data[[time]] == all_periods[match(lab, period_labels)], , drop = FALSE]
    pos    <- match(wide$id, rows_t[[id]])
    wide[[paste0("R_", lab)]]    <- rows_t[[x]][pos]
    wide[[paste0("side_", lab)]] <- as.integer(rows_t[[x]][pos] >= cutoff)
    # Under "atu" the sample is selected on the mirrored x (below the original cutoff), but the
    # TYPE keeps the original orientation, "above the cutoff in the other period": with x
    # mirrored and no ties at the cutoff, original above == mirrored x < 0. Same Wald either way
    # (pi(0) = 1 - pi(1)); this convention makes the reported jump comparable across estimands.
    if (estimand == "atu")
      wide[[paste0("side_", lab)]] <- 1L - wide[[paste0("side_", lab)]]
  }
  wide
}

#' Reflected sample of one pair. The units above the cutoff in the RD period enter at
#' x' = x - cutoff >= 0, those above it in the comparison period t0 at x' = -(x - cutoff) <= 0,
#' each with its type; units without a type (unobserved in a period the type needs) are dropped.
#' Returns list(x_trd, id_trd, type_trd, x_t0, id_t0, type_t0, n_trd, n_t0, n_both).
#' @noRd
.compstable_reflect <- function(wide, t_rd, t0, period_labels, cutoff) {
  t_rd_str <- as.character(t_rd)
  t0_str   <- as.character(t0)
  # A unit's type is its sides in the other periods (period order), then its side in the
  # partner period of the pair (t0 for an RD-period unit, t_rd for a comparison-period unit).
  other_periods <- setdiff(period_labels, base::c(t_rd_str, t0_str))

  above_trd   <- wide[!is.na(wide[[paste0("R_", t_rd_str)]]) &
                        wide[[paste0("R_", t_rd_str)]] >= cutoff, , drop = FALSE]
  partner_trd <- above_trd[[paste0("side_", t0_str)]]

  above_t0   <- wide[!is.na(wide[[paste0("R_", t0_str)]]) &
                       wide[[paste0("R_", t0_str)]] >= cutoff, , drop = FALSE]
  partner_t0 <- above_t0[[paste0("side_", t_rd_str)]]

  if (length(other_periods) > 0L) {
    # a unit unobserved in an other period has no type (as in .build_types)
    other_trd <- apply(above_trd[, paste0("side_", other_periods), drop = FALSE], 1,
                       function(r) if (anyNA(r)) NA_character_ else paste(r, collapse = ""))
    other_t0  <- apply(above_t0[,  paste0("side_", other_periods), drop = FALSE], 1,
                       function(r) if (anyNA(r)) NA_character_ else paste(r, collapse = ""))
  } else {
    # two periods: the type is the partner side alone
    other_trd <- rep("", nrow(above_trd))
    other_t0  <- rep("", nrow(above_t0))
  }

  partner_trd[is.na(partner_trd)] <- NA_integer_
  partner_t0[is.na(partner_t0)]   <- NA_integer_

  type_trd <- ifelse(is.na(partner_trd) | is.na(other_trd), NA_character_,
                     paste0(other_trd, as.character(partner_trd)))
  type_t0  <- ifelse(is.na(partner_t0) | is.na(other_t0),  NA_character_,
                     paste0(other_t0,  as.character(partner_t0)))

  xref_trd <- above_trd[[paste0("R_", t_rd_str)]] - cutoff
  xref_t0  <-  -(above_t0[[paste0("R_", t0_str)]] - cutoff)
  # A t0-above unit exactly at the cutoff reflects to 0 and would be assigned to the right
  # (t_rd) group by rd_period's `x >= c` split; keep it on the reflected (left) side with a
  # negative value that carries full kernel weight.
  xref_t0[xref_t0 == 0] <- -.Machine$double.xmin

  id_trd <- above_trd$id
  id_t0  <- above_t0$id

  keep_trd <- !is.na(type_trd)
  keep_t0  <- !is.na(type_t0)

  xref_trd <- xref_trd[keep_trd]
  id_trd   <- id_trd[keep_trd]
  type_trd <- type_trd[keep_trd]
  xref_t0  <- xref_t0[keep_t0]
  id_t0    <- id_t0[keep_t0]
  type_t0  <- type_t0[keep_t0]

  n_trd  <- length(id_trd)
  n_t0   <- length(id_t0)
  n_both <- length(intersect(id_trd, id_t0))

  list(x_trd = xref_trd, id_trd = id_trd, type_trd = type_trd,
       x_t0 = xref_t0, id_t0 = id_t0, type_t0 = type_t0,
       n_trd = n_trd, n_t0 = n_t0, n_both = n_both)
}

#' Local-linear RD of each type indicator on the reflected sample of one pair, at the artificial
#' cutoff 0; returns list(theta = the jumps (bias-corrected if `bc`, NA for a failed fit),
#' fits = the rd_period fits).
#' @noRd
.compstable_fit_types <- function(refl, type_values, n_types, kernel, h, bwselect, bc) {
  # rd_period splits the stacked sample at 0: the RD-period units form the "+" side, the
  # comparison-period units the "-" side
  x_all    <- base::c(refl$x_trd, refl$x_t0)
  id_all   <- base::c(refl$id_trd, refl$id_t0)
  type_all <- base::c(refl$type_trd, refl$type_t0)

  theta <- numeric(n_types)
  fits  <- vector("list", n_types)
  names(fits) <- type_values

  for (idx_v in seq_along(type_values)) {
    type_v <- type_values[idx_v]
    y_v    <- as.numeric(type_all == type_v)
    # h = NULL (bwselect = "cct"): this regression's own CCT bandwidths; otherwise h is used as
    # both bandwidths
    bw     <- .cell_bandwidth(y_v, x_all, 0, kernel, h, bwselect)
    h_cell <- bw[["h"]]
    b_cell <- bw[["b"]]
    # Any rd_period error (typically: at most q + 1 = 3 units in a side's window) skips the
    # type: its jump is NA and it leaves the Wald test. Errors of the CCT bandwidth above are
    # not caught.
    fit <- tryCatch(
      rd_period(y = y_v, x = x_all, h = h_cell, b = b_cell, id = id_all,
                c = 0, p = 1L, q = 2L, kernel = kernel),
      error = function(e) NULL
    )
    # `fits[i] <- list(NULL)` keeps a NULL in place; `fits[[i]] <- NULL` would delete the
    # element and misalign the types (2026-10-08 fix: the last type's failure then crashed)
    fits[idx_v]   <- list(fit)
    theta[idx_v]  <- if (is.null(fit)) NA_real_ else if (bc) fit$D_bc else fit$D
  }
  list(theta = theta, fits = fits)
}

#' Covariance matrix of the type jumps of one pair from the rd_period influence vectors g;
#' returns Sigma (row and column k: type k; zero rows and columns for failed fits).
#' @noRd
.compstable_sigma <- function(fits, n_types, use_scheme, bc) {
  # Every type indicator is fitted on the same stacked sample, so the same-side terms enter
  # every entry. Under "pv" a unit above the cutoff in both periods sits on both sides of the
  # artificial cutoff (matched on id), and the opposite-side terms are subtracted as well: on
  # the diagonal twice the cross-side sum, off it the `pv` term of .cross_cov(). "cs" and "pc"
  # treat the two sides as independent. With binary types only one jump is tested (see
  # .compstable_kept_types()), so the off-diagonal matters with three or more periods only.
  Sigma <- matrix(0, n_types, n_types)

  for (type_idx in seq_len(n_types)) {
    fit1 <- fits[[type_idx]]
    if (is.null(fit1)) next
    for (type_idx2 in seq_len(n_types)) {
      fit2 <- fits[[type_idx2]]
      if (is.null(fit2)) next

      if (type_idx == type_idx2) {
        g_plus   <- if (bc) fit1$sides$`+`$g_bc else fit1$sides$`+`$g
        g_minus  <- if (bc) fit1$sides$`-`$g_bc else fit1$sides$`-`$g
        var_jump <- if (bc) fit1$V_D_bc else fit1$V_D
        if (use_scheme == "pv") {
          cross_side <- .match_sum(fit1$sides$`+`$id, g_plus,
                                   fit1$sides$`-`$id, g_minus)
          Sigma[type_idx, type_idx] <- var_jump - 2 * cross_side
        } else {
          Sigma[type_idx, type_idx] <- var_jump
        }
        next
      }

      xcov <- .cross_cov(fit1, fit2, bc = bc)
      Sigma[type_idx, type_idx2] <- xcov$pc - (if (use_scheme == "pv") xcov$pv else 0)
    }
  }
  Sigma
}

#' Indices of the type jumps that enter the Wald test of one pair; returns an integer vector.
#' @noRd
.compstable_kept_types <- function(theta, n_types) {
  # The type indicators sum to 1 on each side of the artificial cutoff, so when every type is
  # fitted at a common bandwidth the n_types jumps sum to zero and their covariance is
  # rank-deficient by construction. Drop one reference type (the first in radix order: partner
  # side 0 / the all-below pattern) and test the rest: with binary types the kept jump is the
  # paper's pi_{tRD,(+)}(1) - pi_{t0,(+)}(1), chi-square with 1 df. With per-type CCT
  # bandwidths the jumps do not sum exactly to zero, so this is then a different (still valid)
  # full-rank test rather than an equivalent one. If a fit failed there is no exact redundancy
  # among the survivors: keep them all.
  present <- which(!is.na(theta))
  if (length(present) == n_types && n_types >= 2L) present[-1L] else present
}

#' Test of one pair on its reflected sample (`refl`, from .compstable_reflect()): one jump per
#' type, their covariance, the reference-type drop and the Wald test; returns the pair's element
#' of `pairs` (ll_wald, jumps, jump_se, type_values, scheme, n_trd, n_t0, n_both).
#' @noRd
.compstable_pair_test <- function(refl, scheme, kernel, h, bwselect, bc) {
  # radix: a locale-independent order, so the reference type dropped is the same on every machine
  type_values <- sort(unique(base::c(refl$type_trd, refl$type_t0)), method = "radix")
  n_types     <- length(type_values)
  # "auto": "pv" when some unit is above the cutoff in both periods (it then sits on both sides
  # of the artificial cutoff), "cs" otherwise
  use_scheme <- if (scheme != "auto") scheme else {
    if (refl$n_both > 0L) "pv" else "cs"
  }

  fit   <- .compstable_fit_types(refl, type_values, n_types, kernel, h, bwselect, bc)
  theta <- fit$theta
  Sigma <- .compstable_sigma(fit$fits, n_types, use_scheme, bc)

  ok_idx <- .compstable_kept_types(theta, n_types)
  wald   <- if (length(ok_idx) == 0L) list(stat = 0, df = 0L, p = 1) else
    .joint_wald(theta[ok_idx], Sigma[ok_idx, ok_idx, drop = FALSE])

  list(
    ll_wald     = wald,
    # binary types: the single share jump pi_{tRD,(+)}(1) - pi_{t0,(+)}(1); the standard errors
    # are the square roots of the diagonal of Sigma (dependence-adjusted under "pv")
    jumps       = if (length(ok_idx)) {
      stats::setNames(theta[ok_idx], type_values[ok_idx])
    } else {
      numeric(0)
    },
    jump_se     = if (length(ok_idx)) {
      stats::setNames(sqrt(diag(Sigma)[ok_idx]), type_values[ok_idx])
    } else {
      numeric(0)
    },
    type_values = type_values,
    scheme      = use_scheme,
    n_trd       = refl$n_trd,
    n_t0        = refl$n_t0,
    n_both      = refl$n_both,
    # for plot.rd_compstable(): the reflected sample and each type's fit
    fits        = fit$fits,
    sample      = refl[c("x_trd", "type_trd", "x_t0", "type_t0")]
  )
}

#' Joint test over pairs: the sums of the pair Wald statistics and of their degrees of freedom,
#' with the chi-squared p-value; returns list(stat, df, p).
#' @noRd
.compstable_joint <- function(pairs_out) {
  # Summing assumes independent pairs, which is only approximate: every pair's "+" group is the
  # same RD-period set (and a unit above the cutoff in two comparison periods enters two "-"
  # groups). The paper's test is per pair; the joint test is a convenience summary.
  stat    <- 0
  wald_df <- 0L
  for (key in names(pairs_out)) {
    pair    <- pairs_out[[key]]
    stat    <- stat + pair$ll_wald$stat
    wald_df <- wald_df + pair$ll_wald$df
  }
  p_value <- if (wald_df == 0L) 1 else
    stats::pchisq(stat, df = wald_df, lower.tail = FALSE)
  list(stat = stat, df = wald_df, p = p_value)
}

#' Test of composition stability
#'
#' When the running variable moves over time, some units are above the cutoff
#' in one period and below it in another. A unit's **type** is the side of the
#' cutoff it is on in the other period(s); with two periods, "above in the
#' other period" or "below in the other period". `rd_compstable()` tests the
#' null that **the share of each type among the units just above the cutoff is the
#' same in the RD period and in each comparison period**. Composition
#' stability and homogeneous confounding ([rd_homog()]) are alternatives: the
#' estimate of [rddid()] needs one of the two (together with type continuity,
#' [rd_typecont()]). A rejection here alone therefore does not invalidate the
#' estimate; if homogeneous confounding is rejected as well, the comparison
#' periods mix the types differently from the RD period and the estimate can
#' be biased. With a running variable fixed over time (as in [rddid_sim]) the
#' types are degenerate and the test is not informative.
#'
#' @details
#' ## What is estimated
#'
#' The test runs on each pair of the RD period and one comparison period. Take
#' the units above the cutoff in each of the two periods and stack them into
#' one artificial sample: the RD-period units at their distance above the
#' cutoff, the comparison-period units reflected to the same distance below an
#' artificial cutoff at zero. In a local-linear RD of a type indicator on this
#' reflected running variable, the jump at zero is the type's share among the
#' RD-period units just above the cutoff minus its share among the
#' comparison-period units just above the cutoff. With two periods the type is
#' the unit's side in the other period of the pair, and there is one jump per
#' pair; with more periods the type also records the unit's sides in the
#' remaining periods. The shares sum to one, so one reference type is dropped
#' and the remaining jumps are tested jointly by a Wald statistic, one test per
#' pair. The joint test over pairs adds up the pair statistics and degrees of
#' freedom, which treats the pairs as independent although they share the
#' RD-period units; read it as approximate (the paper's test is per pair). A
#' pair with fewer than three units in either group is skipped with a warning.
#'
#' ## Shared units and the sampling scheme
#'
#' A unit above the cutoff in both periods of a pair appears on both sides of
#' the artificial cutoff. Under `"pv"` the covariance between its two
#' appearances, matched on `id`, is subtracted from the variance of the jump;
#' `"cs"` and `"pc"` treat the two groups as independent and give the same
#' test. `"auto"` (default) uses `"pv"` for a pair in which some unit is above
#' the cutoff in both periods and `"cs"` otherwise; the scheme is reported per
#' pair.
#'
#' ## Options
#'
#' `bc = TRUE` (default) tests the bias-corrected jumps with their robust
#' variance, as in the `Robust` row of [rddid()]; `bc = FALSE` uses the
#' conventional jumps and variances. With `bwselect = "cct"` (default) each
#' type-indicator regression in the reflected sample gets its own CCT
#' bandwidths from [rd_bw_cct()]. With `bwselect = "rot"` the rule of thumb is
#' `h = b = 0.5 * IQR(x)`, the interquartile range of the running variable in
#' `data` (`sd(x)` if that is zero), the same in every regression. A numeric
#' `h` is used as both bandwidths in every regression.
#'
#' ## ATU designs
#'
#' With `estimand = "atu"` (comparison periods uniformly treated) the running
#' variable is mirrored around the cutoff before the construction above, so
#' the test is on the units *below* the cutoff: the null becomes that the share
#' of each type among the units below the cutoff is the same in the RD period
#' and in each comparison period. The type keeps its meaning, so the reported
#' jump is the change in the share of below-cutoff units that are above the
#' cutoff in the other period. This is the only one of the four tests whose
#' computation changes with `estimand`. Units exactly at the cutoff count as
#' above it in the original design and cannot be placed in the mirrored one,
#' so `"atu"` stops with an error if any `x == c`; put the cutoff between
#' support points (e.g. `c = 4999.5` for integer populations).
#'
#' @param data a data frame in long format, one row per unit and period, from a
#'   panel (it need not be balanced). A unit's type is read from the periods in
#'   which it is observed; a unit missing from a period that its type needs is
#'   left out of the cells that use that type.
#' @param id name of the unit-identifier column (a string). Required: types are
#'   read across periods.
#' @param estimand `"att"` (default) or `"atu"`, as in [rddid()]. Under `"atu"`
#'   the test is on the shares among the units below the cutoff; see "ATU
#'   designs" in Details.
#' @param h a bandwidth to use, as both main and pilot bandwidth, in every
#'   regression of the test. If given, `bwselect` is ignored.
#' @param bwselect the bandwidth rule when `h` is not given: `"cct"` (default;
#'   each regression's own CCT bandwidths from [rd_bw_cct()]) or `"rot"` (the
#'   rule of thumb `0.5 * IQR(x)`, the same in every regression).
#' @param scheme the sampling scheme for the covariance between the two groups
#'   of a pair: `"auto"` (default; `"pv"` for a pair in which some unit is above
#'   the cutoff in both periods, `"cs"` otherwise), `"cs"`, `"pc"` or `"pv"`
#'   (as in [rddid()]). See Details.
#' @param bc logical. `TRUE` (default): test the bias-corrected jumps with their
#'   robust variance, as in the `Robust` row of [rddid()]; `FALSE`: the
#'   conventional jumps and variances.
#' @inheritParams rddid
#'
#' @return An object of class `"rd_compstable"`, a list with:
#'   \describe{
#'     \item{`statistic`, `df`, `p_value`}{the joint test over pairs: the sum
#'       of the pair Wald statistics, the sum of their degrees of freedom, and
#'       the chi-squared p-value (approximate; see Details).}
#'     \item{`scheme`}{the sampling scheme used: one value if every pair used
#'       the same, `"mixed"` otherwise; `scheme_requested` is the argument as
#'       passed.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`call`}{the matched call.}
#'     \item{`t_rd`, `comparisons`}{the RD period and the comparison periods
#'       used.}
#'     \item{`pairs`}{a list with one element per pair, named
#'       `"<RD period>::<comparison period>"`, each holding `ll_wald` (the
#'       pair's Wald test: `stat`, `df`, `p`); `jumps` and `jump_se` (the
#'       tested share jumps and their standard errors, named by type: the
#'       unit's sides, `1` above and `0` below the cutoff, with its side in the
#'       other period of the pair last);
#'       `type_values` (the types present); `scheme` (the scheme used for the
#'       pair); and `n_trd`, `n_t0`, `n_both` (the number of units above the
#'       cutoff in the RD period, in the comparison period, and in both; below
#'       the cutoff under `estimand = "atu"`); `fits` (each type's [rd_period()]
#'       fit on the reflected sample, `NULL` where it failed) and `sample` (the
#'       reflected sample: `x_trd`, `type_trd`, `x_t0`, `type_t0`), which feed
#'       [plot.rd_compstable()].}
#'     \item{`joint`}{the joint test again, as `ll_wald` (`stat`, `df`, `p`).}
#'     \item{`meta`}{a list with `t_rd`, `comparisons`, `h` (the common
#'       bandwidth, `NA` with `bwselect = "cct"`), `bwselect`, `c` (the cutoff
#'       as passed, before any mirroring), `bc` and `estimand`.}
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
#' cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
#' cs          # one Wald test per (RD period, comparison period) pair, then their sum
#' cs$pairs[["3::1"]]$jumps
#' # comparison periods uniformly treated: the shares below the cutoff
#' rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3,
#'               estimand = "atu")
#' @export
rd_compstable <- function(data, x, time, id, t_rd,
                          comparisons = NULL,
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
  cutoff   <- c   # the cutoff; `c` stays the argument name for rdrobust users

  # ----- inputs and periods (the stops stay here so that errors name rd_compstable()) -----
  .check_columns(data, c(x, time, id))
  c_orig <- cutoff
  if (estimand == "atu") {
    if (any(data[[x]] == cutoff, na.rm = TRUE))
      stop("estimand = \"atu\": ", sum(data[[x]] == cutoff, na.rm = TRUE),
           " observation(s) have x == c. Units at the cutoff are treated in the original ",
           "design but cannot be placed on the treated side of the mirrored design; ",
           "set the cutoff between support points (e.g. c = 4999.5 for integer populations) ",
           "so that no unit sits on it.")
    # mirror x around the cutoff: the units below it are now the ones above the cutoff 0
    data[[x]] <- cutoff - data[[x]]
    cutoff    <- 0
  }
  data <- data[stats::complete.cases(data[, base::c(x, time, id)]), , drop = FALSE]

  all_periods <- sort(unique(data[[time]]))
  if (!t_rd %in% all_periods) stop("`t_rd` (", t_rd, ") is not a period in `data`.")
  if (is.null(comparisons)) comparisons <- setdiff(all_periods, t_rd)
  if (length(comparisons) == 0L) stop("no comparison periods found.")
  missing_comp <- setdiff(comparisons, all_periods)
  if (length(missing_comp) > 0L)
    stop("comparison periods not in data: ", paste(missing_comp, collapse = ", "))

  # bwselect = "cct" leaves h NULL: each type regression then gets its own CCT bandwidths
  if (is.null(h) && bwselect == "rot") h <- .rot_bandwidth_iqr(data[[x]])

  # ----- one Wald test per (RD period, comparison period) pair -----
  period_labels <- as.character(all_periods)
  wide <- .compstable_wide(data, x, time, id, all_periods, period_labels, cutoff, estimand)
  pairs_out <- list()
  for (t0 in comparisons) {
    pair_key <- paste0(as.character(t_rd), "::", as.character(t0))
    refl     <- .compstable_reflect(wide, t_rd, t0, period_labels, cutoff)
    if (refl$n_trd < 3L || refl$n_t0 < 3L) {
      warning("rd_compstable: pair ", pair_key,
              " has too few above-cutoff observations; skipping.")
      next
    }
    pairs_out[[pair_key]] <- .compstable_pair_test(refl, scheme, kernel, h, bwselect, bc)
  }

  joint        <- .compstable_joint(pairs_out)
  pair_schemes <- vapply(pairs_out, function(pair) pair$scheme, character(1))
  structure(
    list(
      statistic        = joint$stat,
      df               = joint$df,
      p_value          = joint$p,
      scheme           = if (length(unique(pair_schemes)) == 1L) unique(pair_schemes) else "mixed",
      scheme_requested = scheme,
      estimand         = estimand,
      t_rd             = t_rd,
      comparisons      = comparisons,
      pairs            = pairs_out,
      joint            = list(ll_wald = joint),
      call             = cl,
      meta = list(
        t_rd        = t_rd,
        comparisons = comparisons,
        h           = if (!is.null(h)) h else NA_real_,
        bwselect    = bwselect,
        kernel      = kernel,
        c           = c_orig,
        bc          = bc,
        estimand    = estimand
      )
    ),
    class = "rd_compstable"
  )
}


#' @export
print.rd_compstable <- function(x, ...) {
  side <- if (identical(x$estimand, "atu")) "below" else "above"
  h0 <- sprintf(paste0("the share of each type among the units just %s the cutoff is the same in ",
                       "the RD period and in each comparison period"), side)
  atu_note <- paste0("the units below the cutoff are the ones untreated in the RD period, ",
                     "so the test is on their shares (mirrored design)")
  .print_test_header("composition stability", "rd_compstable", h0,
                     x$scheme, identical(x$scheme_requested, "auto"), x$estimand,
                     atu_note = atu_note)
  cat(sprintf("  RD period: %s   Comparison periods: %s   Bandwidth: %s\n\n",
              x$t_rd, paste(x$comparisons, collapse = ", "),
              .bw_label_test(x$meta$h, x$meta$bwselect)))
  for (key in names(x$pairs)) {
    pair <- x$pairs[[key]]
    .print_wald(pair$ll_wald$stat, pair$ll_wald$df, pair$ll_wald$p,
                label = sprintf("Pair %s:", key))
    cat(sprintf("    n %s the cutoff: %d (RD period), %d (comparison), %d in both\n",
                side, pair$n_trd, pair$n_t0, pair$n_both))
  }
  if (length(x$pairs) > 1L) {
    cat("\n")
    .print_wald(x$statistic, x$df, x$p_value,
                label = "Joint over pairs (sum of chi-squared):")
  }
  invisible(x)
}
