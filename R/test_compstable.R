#' Test composition stability across periods
#'
#' Wald test of the composition-stability assumption (Section 4.4 of Leventer
#' and Nevo): for each RD-period / comparison-period pair
#' \eqn{(t_{\mathrm{RD}}, t_0)}, the share of each type among the units just
#' above the cutoff is the same in the two periods,
#' \eqn{\pi_{t_{\mathrm{RD}},(+)}(v) = \pi_{t_0,(+)}(v)}. With two periods
#' the type is binary (the unit's side in the other period) and this is the
#' single share-jump test of the paper's Section 4.4; with more periods the
#' type \eqn{(\mathbf{u}, b)} collects the unit's sides in the other
#' comparison periods, \eqn{\mathbf{u}}, and its side \eqn{b} in the partner
#' period of the pair.
#'
#' ## Reflection construction
#'
#' For each RD-period / comparison-period pair \eqn{(t_{\mathrm{RD}}, t_0)},
#' take the above-cutoff units of each period.  For \eqn{t_0}-above units flip
#' the centred running variable: \eqn{x' = -(R_{i,t_0} - c)}, placing them
#' just BELOW an artificial cutoff at 0.  For \eqn{t_{\mathrm{RD}}}-above
#' units set \eqn{x' = R_{i,t_{\mathrm{RD}}} - c}, keeping them just ABOVE 0.
#' Stack the two groups into one artificial cross-section.  At the artificial
#' cutoff the left/right limits of the \eqn{(\mathbf{u},b)} type share are
#' then \eqn{\pi_{t_0,(+)}(\mathbf{u},b)} and
#' \eqn{\pi_{t_{\mathrm{RD}},(+)}(\mathbf{u},b)}, so a jump at 0 equals the
#' composition difference.  A local-linear RD of each \eqn{(\mathbf{u},b)}
#' indicator on the reflected running variable, with a joint Wald that all
#' jumps are zero (dropping one reference type, since the type shares sum to
#' 1 and the full set of jumps is rank-deficient), is then this test.
#'
#' ## Unit-level wrinkle
#'
#' A unit that is above the cutoff in BOTH periods \eqn{t_{\mathrm{RD}}} and
#' \eqn{t_0} appears on BOTH sides of the artificial cutoff (as a
#' \eqn{t_{\mathrm{RD}}}-above observation above the artificial 0 and a
#' \eqn{t_0}-above observation below it).  The covariance matrix between the
#' left-side and right-side intercept estimates must therefore include the
#' id-matched cross-side covariance term (the two g-vectors can share unit
#' ids).  The function computes \eqn{(\text{cov}_{++} + \text{cov}_{--} -
#' \text{cov}_{+-} - \text{cov}_{-+})} — the same formula as the PV scheme in
#' the main estimator — rather than assuming the two sides are independent.
#'
#' ## ATU designs
#'
#' With `estimand = "atu"` the running variable is mirrored,
#' \eqn{x \to c - x} (and the cutoff reset to 0), before the construction
#' above runs. This takes the units BELOW the original cutoff in each period.
#' The type indicator keeps its original orientation, "above the cutoff in the
#' other period", so the jump estimates
#' \eqn{\pi_{t_{\mathrm{RD}},(-)}(1) - \pi_{t_0,(-)}(1)}, the change across
#' periods in the share of below-cutoff units that are above the cutoff in
#' the other period -- the composition-stability condition the ATU requires
#' (Leventer and Nevo, Section 6). (Stating the jump for the complementary
#' type, "below in the other period", would flip its sign and leave the Wald
#' test unchanged.) Units with `x == c` are treated in the
#' original design but cannot be placed on the treated side of the mirrored
#' design, so `estimand = "atu"` errors if any are present; place the cutoff
#' between support points (e.g. `c = 4999.5` for integer populations) so that
#' no unit sits on it.
#'
#' @param data a long data frame, one row per unit-period. A unit's type in
#'   period \eqn{t} is read from its running variable in the other period(s);
#'   units unobserved there are dropped from period \eqn{t}, so the panel need
#'   not be balanced.
#' @param x Column name (string) for the running variable.
#' @param time Column name (string) for the period indicator.
#' @param id Column name (string) for the unit identifier.
#' @param t_rd Value of `time` identifying the RD period.
#' @param comparisons Values of `time` to use as comparison periods.  If
#'   `NULL` (default), all periods except `t_rd` are used.
#' @param estimand `"att"` (default) or `"atu"`. Under `"atu"` the running
#'   variable is mirrored before the test runs, so the test is on the
#'   below-cutoff shares instead of the above-cutoff shares; this is the ONE
#'   function among the five with `estimand` where the computation actually
#'   differs. See "ATU designs" above.
#' @param c Cutoff for the running variable (default 0).
#' @param h Bandwidth.  If `NULL` (default), the bandwidth is determined by
#'   `bwselect`; an explicit numeric value overrides `bwselect` and is used
#'   directly.
#' @param bwselect Bandwidth selection rule when `h = NULL`: `"cct"` (default)
#'   computes a per-cell CCT MSE-optimal bandwidth via [rd_bw_cct()] for each
#'   type indicator RD in the reflected space; `"rot"` uses \eqn{0.5 \times
#'   \mathrm{IQR}(x)} as a rule-of-thumb applied to the full sample.  Ignored
#'   when `h` is supplied explicitly.
#' @param kernel Kernel for the local-linear RD: `"triangular"` (default),
#'   `"epanechnikov"`, or `"uniform"`.
#' @param scheme Covariance scheme for the Wald test:
#'   \describe{
#'     \item{`"auto"`}{Detects whether any unit appears on both sides of the
#'       artificial cutoff (i.e., above the true cutoff in both periods).
#'       If yes, uses `"pv"` (time-varying panel); otherwise `"cs"`.}
#'     \item{`"cs"`}{Treats the two sides as independent.}
#'     \item{`"pc"`}{Includes same-side cross-period covariance only.}
#'     \item{`"pv"`}{Full panel with time-varying running variable: includes
#'       same-side minus opposite-side cross-period covariance.}
#'   }
#' @param bc Use robust bias-corrected jumps and variances in the LL-Wald
#'   (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the test
#'   with the bias-corrected [rddid()] estimator; `FALSE` uses the conventional
#'   local-linear jumps and variances.
#' @param ... Currently unused.
#'
#' @return An object of class `"rd_compstable"`, a named list with:
#'   \describe{
#'     \item{`pairs`}{A list, one element per \eqn{(t_{\mathrm{RD}}, t_0)}
#'       pair (named `"trd::t0"`), each containing:
#'       \describe{
#'         \item{`ll_wald`}{list with `stat`, `df`, `p`.}
#'         \item{`jumps`, `jump_se`}{The tested type-share jump(s) and their
#'           standard errors (below-cutoff units when `estimand = "atu"`).}
#'         \item{`type_values`}{Character vector of \eqn{(\mathbf{u},b)} type
#'           labels present in this pair.}
#'         \item{`scheme`}{Scheme actually used.}
#'         \item{`n_trd`}{Number of above-cutoff units from \eqn{t_RD}
#'           (below-cutoff units when `estimand = "atu"`).}
#'         \item{`n_t0`}{Number of above-cutoff units from \eqn{t_0}
#'           (below-cutoff units when `estimand = "atu"`).}
#'         \item{`n_both`}{Number of units above the cutoff in both periods
#'           (below-cutoff units when `estimand = "atu"`).}
#'       }
#'     }
#'     \item{`joint`}{Joint result over all pairs (stacked Wald):
#'       \describe{
#'         \item{`ll_wald`}{list with `stat`, `df`, `p` (sum of the per-pair
#'           statistics and df, which assumes independent pairs; the pairs
#'           share the RD-period above-cutoff group, so treat the joint as
#'           approximate — the paper's test is per pair).}
#'       }
#'     }
#'     \item{`meta`}{list with `t_rd`, `comparisons`, `h` (NA when
#'       `bwselect = "cct"`), `bwselect`, `c` (the original, unmirrored
#'       cutoff, as passed), `bc`, `estimand`.}
#'   }
#'
#' @references
#' Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
#' Designs." Working paper.
#'
#' @seealso [rd_typecont()], [rd_period()], [rddid()]
#'
#' @examples
#' \dontrun{
#' # Two-period panel with no composition shift (null DGP).
#' set.seed(1)
#' n <- 500
#' eta <- rnorm(n)
#' dat <- data.frame(
#'   id   = rep(seq_len(n), 2),
#'   time = rep(1:2, each = n),
#'   R    = c(eta + rnorm(n), eta + rnorm(n))
#' )
#' rd_compstable(dat, x = "R", time = "time", id = "id", t_rd = 2,
#'               comparisons = 1, h = 0.5)
#' }
#' @export
rd_compstable <- function(data, x, time, id, t_rd,
                          comparisons = NULL,
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

  structure(
    list(
      pairs = pairs_out,
      joint = list(
        ll_wald = list(stat = joint_ll_stat, df = joint_ll_df, p = joint_ll_p)
      ),
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
  est <- if (is.null(x$meta$estimand)) "att" else x$meta$estimand
  cat("Composition-stability test\n")
  if (est == "atu")
    cat("  estimand: atu -- test is on the below-cutoff shares (mirrored design)\n")
  cat(sprintf("  RD period: %s   Comparison periods: %s\n",
              x$meta$t_rd,
              paste(x$meta$comparisons, collapse = ", ")))
  h_str <- if (is.na(x$meta$h)) paste0("per-cell ", toupper(x$meta$bwselect)) else sprintf("%.4g", x$meta$h)
  cat(sprintf("  h=%s   bwselect=%s\n\n", h_str, x$meta$bwselect))

  for (pk in names(x$pairs)) {
    pr <- x$pairs[[pk]]
    cat(sprintf("Pair %s  [scheme=%s  n_trd=%d  n_t0=%d  n_both=%d]\n",
                pk, pr$scheme, pr$n_trd, pr$n_t0, pr$n_both))
    cat(sprintf("  LL-Wald:  chi2(%.0f) = %.4f   p = %.4f\n\n",
                pr$ll_wald$df, pr$ll_wald$stat, pr$ll_wald$p))
  }

  if (length(x$pairs) > 1L) {
    cat("Joint (across all pairs):\n")
    cat(sprintf("  LL-Wald:  chi2(%.0f) = %.4f   p = %.4f\n",
                x$joint$ll_wald$df, x$joint$ll_wald$stat, x$joint$ll_wald$p))
  }
  invisible(x)
}
