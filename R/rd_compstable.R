# rd_compstable.R -- rd_compstable(): test of composition stability between the RD period and
# each comparison period (reflected-sample construction), with its print method. Shared pieces:
# assumption_tests_helpers.R, cross_period_covariance.R.

#' Test of composition stability
#'
#' When the running variable moves over time, some units are above the cutoff
#' in one period and below it in another. A unit's **type** is the side of the
#' cutoff it is on in the other period(s); with two periods, "above in the
#' other period" or "below in the other period". `rd_compstable()` tests the
#' null that **the share of each type among the units above the cutoff is the
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
#'       the cutoff under `estimand = "atu"`).}
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

  # ---- input validation -------------------------------------------------------
  for (nm in base::c(x, time, id)) {
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")
  }

  c_orig <- c
  if (estimand == "atu") {
    if (any(data[[x]] == c, na.rm = TRUE))
      stop("estimand = \"atu\": ", sum(data[[x]] == c, na.rm = TRUE),
           " observation(s) have x == c. Units at the cutoff are treated in the original ",
           "design but cannot be placed on the treated side of the mirrored design; ",
           "set the cutoff between support points (e.g. c = 4999.5 for integer populations) ",
           "so that no unit sits on it.")
    data[[x]] <- c - data[[x]]
    c <- 0
  }

  data <- data[stats::complete.cases(data[, base::c(x, time, id)]), , drop = FALSE]

  all_periods <- sort(unique(data[[time]]))
  if (!t_rd %in% all_periods)
    stop("`t_rd` (", t_rd, ") is not a period in `data`.")

  if (is.null(comparisons)) {
    comparisons <- setdiff(all_periods, t_rd)
  }
  if (length(comparisons) == 0L)
    stop("no comparison periods found.")
  missing_comp <- setdiff(comparisons, all_periods)
  if (length(missing_comp) > 0L)
    stop("comparison periods not in data: ",
         paste(missing_comp, collapse = ", "))

  # ---- default bandwidth -------------------------------------------------------
  # When bwselect = "rot" and h = NULL, use 0.5*IQR (current behaviour).
  # When bwselect = "cct" and h = NULL, h stays NULL; per-cell CCT is computed
  # inside the per-pair type loop below in the reflected space (c = 0).
  if (is.null(h) && bwselect == "rot") {
    all_x <- data[[x]]
    h     <- 0.5 * stats::IQR(all_x)
    if (h <= 0) h <- stats::sd(all_x)
  }

  # ---- wide pivot (id x period running variable and side) ----------------------
  # We need, for each unit, its running variable in each period.
  plab_all <- as.character(all_periods)
  wide <- data.frame(id = unique(data[[id]]), stringsAsFactors = FALSE)
  for (k in plab_all) {
    sub <- data[data[[time]] == all_periods[match(k, plab_all)], , drop = FALSE]
    m   <- match(wide$id, sub[[id]])
    wide[[paste0("R_", k)]]    <- sub[[x]][m]
    wide[[paste0("side_", k)]] <- as.integer(sub[[x]][m] >= c)
    # Under "atu" the sample is selected on the mirrored x (below the original
    # cutoff), but the TYPE keeps the original orientation, "above the cutoff in
    # the other period": with x mirrored and no ties at the cutoff, original
    # above == mirrored x < 0. Same Wald either way (pi(0) = 1 - pi(1)); this
    # convention makes the reported jump comparable across estimands.
    if (estimand == "atu")
      wide[[paste0("side_", k)]] <- 1L - wide[[paste0("side_", k)]]
  }

  # ---- per-pair analysis -------------------------------------------------------
  pairs_out <- list()

  for (t0 in comparisons) {
    pair_key <- paste0(as.character(t_rd), "::", as.character(t0))

    t_rd_str <- as.character(t_rd)
    t0_str   <- as.character(t0)

    # Shared "other" periods u = periods excluding both t_rd and t0
    u_periods <- setdiff(plab_all, base::c(t_rd_str, t0_str))

    # ---------- build the reflected cross-section --------------------------------
    # Take above-cutoff units from t_rd:
    #   reflected x' = R_{i,t_rd} - c   (positive, above artificial 0)
    #   "type" (u,b): u = sides of u_periods, b = side of t0 (the partner)
    above_trd <- wide[!is.na(wide[[paste0("R_", t_rd_str)]]) &
                        wide[[paste0("R_", t_rd_str)]] >= c, , drop = FALSE]
    # partner side for t_rd rows = t0's side (period-t0 side)
    b_trd <- above_trd[[paste0("side_", t0_str)]]

    # Take above-cutoff units from t0:
    #   reflected x' = -(R_{i,t0} - c)  (negative, below artificial 0)
    #   "type" (u,b): u = sides of u_periods, b = side of t_rd (the partner)
    above_t0  <- wide[!is.na(wide[[paste0("R_", t0_str)]]) &
                        wide[[paste0("R_", t0_str)]] >= c, , drop = FALSE]
    b_t0 <- above_t0[[paste0("side_", t_rd_str)]]

    # Build u-string (sides of shared other periods)
    if (length(u_periods) > 0L) {
      # a unit unobserved in a shared period has no type (as in .build_types)
      u_trd <- apply(above_trd[, paste0("side_", u_periods), drop = FALSE], 1,
                     function(r) if (anyNA(r)) NA_character_ else paste(r, collapse = ""))
      u_t0  <- apply(above_t0[,  paste0("side_", u_periods), drop = FALSE], 1,
                     function(r) if (anyNA(r)) NA_character_ else paste(r, collapse = ""))
    } else {
      # P = 2: no shared other periods; u is empty
      u_trd <- rep("", nrow(above_trd))
      u_t0  <- rep("", nrow(above_t0))
    }

    # Handle NAs in partner sides (units not observed in a period)
    b_trd[is.na(b_trd)] <- NA_integer_
    b_t0[is.na(b_t0)]   <- NA_integer_

    # type string = paste(u, b)
    type_trd <- ifelse(is.na(b_trd) | is.na(u_trd), NA_character_,
                       paste0(u_trd, as.character(b_trd)))
    type_t0  <- ifelse(is.na(b_t0) | is.na(u_t0),  NA_character_,
                       paste0(u_t0,  as.character(b_t0)))

    # Reflected running variable
    xref_trd <- above_trd[[paste0("R_", t_rd_str)]] - c   # >= 0
    xref_t0  <-  -(above_t0[[paste0("R_", t0_str)]] - c)  # <= 0
    # A t0-above unit exactly at the cutoff reflects to 0 and would be assigned
    # to the right (t_rd) group by rd_period's `x >= c` split; keep it on the
    # reflected (left) side with a negative value that carries full kernel weight.
    xref_t0[xref_t0 == 0] <- -.Machine$double.xmin

    id_trd <- above_trd$id
    id_t0  <- above_t0$id

    # Drop rows with NA type (unit not observed in one of the periods)
    keep_trd <- !is.na(type_trd)
    keep_t0  <- !is.na(type_t0)

    xref_trd  <- xref_trd[keep_trd];  id_trd   <- id_trd[keep_trd]
    type_trd  <- type_trd[keep_trd]
    xref_t0   <- xref_t0[keep_t0];    id_t0    <- id_t0[keep_t0]
    type_t0   <- type_t0[keep_t0]

    n_trd_obs <- length(id_trd)
    n_t0_obs  <- length(id_t0)
    n_both    <- length(intersect(id_trd, id_t0))

    if (n_trd_obs < 3L || n_t0_obs < 3L) {
      warning("rd_compstable: pair ", pair_key,
              " has too few above-cutoff observations; skipping.")
      next
    }

    # All type values present in this pair
    all_type_vals <- sort(unique(base::c(type_trd, type_t0)), method = "radix")  # locale-independent
    n_types <- length(all_type_vals)

    # ---- auto-detect scheme ----------------------------------------------------
    use_scheme <- if (scheme != "auto") scheme else {
      if (n_both > 0L) "pv" else "cs"
    }

    # ---- LL-Wald -----------------------------------------------------------
    # For each type value v, run rd_period on the type indicator
    # y = 1{type == v}, x = xref, on the reflected data (trd above → "+", t0 above → "-")
    # The "+" side uses (xref_trd, id_trd, type_trd)
    # The "-" side uses (xref_t0,  id_t0,  type_t0)
    # We call rd_period on the combined reflected data; rd_period itself
    # splits by sign of x (>= c = 0 or < 0).

    x_all    <- base::c(xref_trd, xref_t0)
    id_all   <- base::c(id_trd,   id_t0)
    type_all <- base::c(type_trd, type_t0)

    theta <- numeric(n_types)
    fits  <- vector("list", n_types)
    names(fits) <- all_type_vals

    for (vi in seq_along(all_type_vals)) {
      v   <- all_type_vals[vi]
      y_v <- as.numeric(type_all == v)
      # Per-cell bandwidth in the reflected space (c = 0): CCT when h = NULL
      # and bwselect = "cct"; otherwise h is non-NULL (explicit or pre-set
      # from 0.5*IQR for bwselect = "rot").
      bw     <- .cell_bandwidth(y_v, x_all, 0, kernel, h, bwselect)
      h_cell <- bw[["h"]]
      b_cell <- bw[["b"]]
      fit <- tryCatch(
        rd_period(y = y_v, x = x_all, h = h_cell, b = b_cell, id = id_all,
                  c = 0, p = 1L, q = 2L, kernel = kernel),
        error = function(e) NULL
      )
      fits[[vi]] <- fit
      theta[vi]  <- if (is.null(fit)) NA_real_ else if (bc) fit$D_bc else fit$D
    }

    # Build covariance matrix
    # The "+" side of the artificial cutoff = t_rd-above units
    # The "-" side = t_0-above units
    # Units in both appear in fits[[v]]$sides$`+`$id AND fits[[v]]$sides$`-`$id
    # Diagonal = Var(D_v) = sum(g_+^2) + sum(g_-^2), minus 2 x the shared-unit
    #   cross-side term under "pv" (a unit above in both periods sits on both
    #   sides of the artificial cutoff).
    # Off-diagonal (types v != v'): both indicator fits run on the SAME reflected
    #   sample, so every unit enters both with different 0/1 outcomes and the
    #   same-side term sum(g_v g_v') is nonzero (as in rd_typecont's within-period
    #   block); under "pv" the opposite-side term is subtracted as well. With
    #   binary types only one jump is kept (below), so the off-diagonal matters
    #   for P >= 3 only.

    N     <- n_types
    Sigma <- matrix(0, N, N)

    for (a in seq_len(N)) {
      fit_a <- fits[[a]]
      if (is.null(fit_a)) next
      for (b_idx in seq_len(N)) {
        fit_b <- fits[[b_idx]]
        if (is.null(fit_b)) next

        # Within-type variance (a == b): standard HC1 formula.
        # Under the pv scheme the two artificial-cutoff sides share units (a
        # unit above in both periods contributes a g on BOTH sides), so the
        # diagonal must subtract the id-matched cross-side term — the same
        # "same-side minus opposite-side" logic used for the off-diagonal.
        if (a == b_idx) {
          gp <- if (bc) fit_a$sides$`+`$g_bc else fit_a$sides$`+`$g
          gm <- if (bc) fit_a$sides$`-`$g_bc else fit_a$sides$`-`$g
          V_aa <- if (bc) fit_a$V_D_bc else fit_a$V_D
          if (use_scheme == "pv") {
            cross <- .match_sum(fit_a$sides$`+`$id, gp,
                                fit_a$sides$`-`$id, gm)
            Sigma[a, a] <- V_aa - 2 * cross
          } else {
            Sigma[a, a] <- V_aa
          }
          next
        }

        # Cross-type off-diagonal: same-side term always (same sample, different
        # indicator outcomes); opposite-side term only when the two artificial
        # sides share units ("pv").
        cc <- .cross_cov(fit_a, fit_b, bc = bc)
        Sigma[a, b_idx] <- cc$pc - (if (use_scheme == "pv") cc$pv else 0)
      }
    }

    # The type indicators sum to 1 on each side of the artificial cutoff, so when
    # every type is fitted at a common bandwidth the n_types jumps sum to zero
    # and their covariance is rank-deficient by construction. Drop one reference
    # type (the first in radix order: partner side 0 / the all-below pattern) and
    # test the rest: with binary types the kept jump is the paper's
    # pi_{tRD,(+)}(1) - pi_{t0,(+)}(1), chi-square with 1 df. With per-type CCT
    # bandwidths the jumps do not sum exactly to zero, so this is then a
    # different (still valid) full-rank test rather than an equivalent one. If a
    # fit failed there is no exact redundancy among the survivors: keep them all.
    present   <- which(!is.na(theta))
    ok_idx    <- if (length(present) == n_types && n_types >= 2L) present[-1L] else present
    ll_result <- if (length(ok_idx) == 0L) list(stat = 0, df = 0L, p = 1) else
      .joint_wald(theta[ok_idx], Sigma[ok_idx, ok_idx, drop = FALSE])

    pairs_out[[pair_key]] <- list(
      ll_wald    = ll_result,
      # the tested jumps (binary types: the single share jump pi_{tRD,(+)}(1) - pi_{t0,(+)}(1))
      # and their dependence-adjusted standard errors (diagonal of Sigma)
      jumps      = if (length(ok_idx)) stats::setNames(theta[ok_idx], all_type_vals[ok_idx]) else numeric(0),
      jump_se    = if (length(ok_idx)) stats::setNames(sqrt(diag(Sigma)[ok_idx]), all_type_vals[ok_idx]) else numeric(0),
      type_values = all_type_vals,
      scheme     = use_scheme,
      n_trd      = n_trd_obs,
      n_t0       = n_t0_obs,
      n_both     = n_both
    )
  }

  # ---- joint result across all pairs ------------------------------------------
  # Sum the per-pair statistics and df. This assumes independent pairs, which
  # is only approximate: every pair's "+" group is the same t_rd-above set (and
  # a unit above the cutoff in two comparison periods enters two "-" groups).
  # The paper's test is per pair; the joint is a convenience summary.
  # For the Wald: sum chi-sq statistics with summed df.
  joint_ll_stat <- 0
  joint_ll_df   <- 0L

  for (pk in names(pairs_out)) {
    pr <- pairs_out[[pk]]
    joint_ll_stat <- joint_ll_stat + pr$ll_wald$stat
    joint_ll_df   <- joint_ll_df   + pr$ll_wald$df
  }
  joint_ll_p <- if (joint_ll_df == 0L) 1 else
    stats::pchisq(joint_ll_stat, df = joint_ll_df, lower.tail = FALSE)

  pair_schemes <- vapply(pairs_out, function(pr) pr$scheme, character(1))
  structure(
    list(
      statistic   = joint_ll_stat,
      df          = joint_ll_df,
      p_value     = joint_ll_p,
      scheme      = if (length(unique(pair_schemes)) == 1L) unique(pair_schemes) else "mixed",
      scheme_requested = scheme,
      estimand    = estimand,
      t_rd        = t_rd,
      comparisons = comparisons,
      pairs = pairs_out,
      joint = list(
        ll_wald = list(stat = joint_ll_stat, df = joint_ll_df, p = joint_ll_p)
      ),
      call = cl,
      meta = list(
        t_rd        = t_rd,
        comparisons = comparisons,
        h           = if (!is.null(h)) h else NA_real_,
        bwselect    = bwselect,
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
  .print_test_header("composition stability", "rd_compstable",
                     sprintf("the share of each type among the units %s the cutoff is the same in the RD period and in each comparison period", side),
                     x$scheme, identical(x$scheme_requested, "auto"), x$estimand,
                     atu_note = "the units below the cutoff are the ones untreated in the RD period, so the test is on their shares (mirrored design)")
  cat(sprintf("  RD period: %s   Comparison periods: %s   Bandwidth: %s\n\n",
              x$t_rd, paste(x$comparisons, collapse = ", "),
              .bw_label_test(x$meta$h, x$meta$bwselect)))
  for (pk in names(x$pairs)) {
    pr <- x$pairs[[pk]]
    .print_wald(pr$ll_wald$stat, pr$ll_wald$df, pr$ll_wald$p,
                label = sprintf("Pair %s:", pk))
    cat(sprintf("    n %s the cutoff: %d (RD period), %d (comparison), %d in both\n",
                side, pr$n_trd, pr$n_t0, pr$n_both))
  }
  if (length(x$pairs) > 1L) {
    cat("\n")
    .print_wald(x$statistic, x$df, x$p_value, label = "Joint over pairs (sum of chi-squared):")
  }
  invisible(x)
}
