#' Test composition stability across periods (Assumption A8)
#'
#' Tests Assumption A8 ("composition stability") from Leventer and Nevo:
#' \deqn{\pi_{t_{\mathrm{RD}},(+)}(\mathbf{u}, b) =
#'        \pi_{t_0,(+)}(\mathbf{u}, b)
#'        \quad \forall\,(\mathbf{u}, b),}
#' where \eqn{\mathbf{u} = \mathbf{v}_{-\{t_0, t_{\mathrm{RD}}\}}} are the
#' sides of the OTHER periods (shared between the two confounding objects) and
#' \eqn{b \in \{0,1\}} is the "partner" side.  This is a cross-period
#' covariate-continuity statement: the above-cutoff \eqn{(\mathbf{u},b)} type
#' mix must be the same at the RD period and at each comparison period.
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
#' composition difference.  The same two tests from the type-continuity
#' assessment (A7) then apply:
#'
#' \enumerate{
#'   \item **LL-Wald** (necessary AND sufficient): local-linear RD of each
#'     \eqn{(\mathbf{u},b)} indicator on the reflected running variable; joint
#'     Wald that all jumps are zero, with a Moore-Penrose pseudo-inverse
#'     handling the singularity (type shares sum to 1).
#'   \item **Canay-Kamat permutation** (necessary AND sufficient): approximate
#'     sign randomisation test comparing the two sides of the artificial
#'     cutoff, permuted at the **unit** level (see "Unit-level wrinkle" below).
#' }
#'
#' Both tests are **necessary AND sufficient** for Assumption A8.
#'
#' ## Unit-level wrinkle
#'
#' A unit that is above the cutoff in BOTH periods \eqn{t_{\mathrm{RD}}} and
#' \eqn{t_0} appears on BOTH sides of the artificial cutoff (as a
#' \eqn{t_{\mathrm{RD}}}-above observation above the artificial 0 and a
#' \eqn{t_0}-above observation below it).  This has two consequences handled
#' by the function:
#'
#' (i) **Wald test**: the covariance matrix between the left-side and
#'   right-side intercept estimates must include the id-matched cross-side
#'   covariance term (the two g-vectors can share unit ids).  The function
#'   computes \eqn{(\text{cov}_{++} + \text{cov}_{--} - \text{cov}_{+-} -
#'   \text{cov}_{-+})} — the same formula as the PV scheme in the main
#'   estimator — rather than assuming the two sides are independent.
#'
#' (ii) **Permutation test**: permutes at the unit level.  Each unit
#'   contributes its observations (possibly >1) as a block; the side label is
#'   permuted across units, not across individual rows.  This follows
#'   Amro and Pauly (2017) and Derrick et al. (2022) (see References).
#'
#' @param data A long data frame, one row per unit-period (balanced or
#'   unbalanced panel).
#' @param x Column name (string) for the running variable.
#' @param time Column name (string) for the period indicator.
#' @param id Column name (string) for the unit identifier.
#' @param t_rd Value of `time` identifying the RD period.
#' @param comparisons Values of `time` to use as comparison periods.  If
#'   `NULL` (default), all periods except `t_rd` are used.
#' @param c Cutoff for the running variable (default 0).
#' @param h Bandwidth.  If `NULL` (default), uses \eqn{0.5 \times
#'   \mathrm{IQR}(x)} across the full sample as a simple default; supply an
#'   explicit bandwidth for reproducible results.
#' @param q Number of observations nearest the artificial cutoff on each side
#'   for the Canay-Kamat permutation test. `NULL` (default) selects `q` per
#'   \eqn{(t_{\mathrm{RD}}, t_0)} pair by the Canay & Kamat (2018) rule of thumb
#'   (see [rd_typecont()]); this is the recommended choice, since a fixed `q`
#'   over-rejects in finite samples when the type distribution varies steeply in
#'   the running variable at the cutoff. Supply an integer to force a fixed `q`
#'   on every pair. The per-pair `q` actually used is returned in `meta$q_used`.
#' @param S Number of permutation replications (default 499).
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
#' @param ... Currently unused.
#'
#' @return An object of class `"rd_compstable"`, a named list with:
#'   \describe{
#'     \item{`pairs`}{A list, one element per \eqn{(t_{\mathrm{RD}}, t_0)}
#'       pair (named `"trd::t0"`), each containing:
#'       \describe{
#'         \item{`ll_wald`}{list with `stat`, `df`, `p`.}
#'         \item{`ck_perm`}{list with `stat` (observed sum of |mean diffs|)
#'           and `p`.}
#'         \item{`type_values`}{Character vector of \eqn{(\mathbf{u},b)} type
#'           labels present in this pair.}
#'         \item{`scheme`}{Scheme actually used.}
#'         \item{`q`}{Number of nearest observations per side actually used.}
#'         \item{`n_trd`}{Number of above-cutoff units from \eqn{t_RD}.}
#'         \item{`n_t0`}{Number of above-cutoff units from \eqn{t_0}.}
#'         \item{`n_both`}{Number of units above the cutoff in both periods.}
#'       }
#'     }
#'     \item{`joint`}{Joint result over all pairs (stacked Wald + minimum-p
#'       permutation envelope):
#'       \describe{
#'         \item{`ll_wald`}{list with `stat`, `df`, `p` (stacked across all
#'           pairs, assuming independence across pairs).}
#'         \item{`ck_perm`}{list with `stat` (sum of per-pair stats) and `p`.}
#'       }
#'     }
#'     \item{`meta`}{list with `t_rd`, `comparisons`, `h`, `q` (`"rot"` when the
#'       rule of thumb is used), `q_used` (per-pair `q`), `S`, `c`.}
#'   }
#'
#' @note
#' **Necessary and sufficient status:** Both the LL-Wald and the Canay-Kamat
#' permutation test are necessary AND sufficient for Assumption A8 (composition
#' stability).  See Leventer and Nevo for the proof.
#'
#' @references
#' Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
#' Designs." Working paper.
#'
#' Amro, L. and Pauly, M. (2017). Permuting longitudinal data in spite of the
#' dependencies. *Journal of Statistical Computation and Simulation*, 87(15),
#' 3033-3044.
#'
#' Derrick, B., Broad, A., Ruck, A., and White, P. (2022). The impact of
#' repeated measures on the permutation test. *Journal of Applied Quantitative
#' Methods*, 17(1).
#'
#' Canay, I. A. and Kamat, V. (2018). Approximate permutation tests and
#' induced order statistics in the regression discontinuity design. *Review of
#' Economic Studies*, 85(3), 1577-1608.
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
#'               comparisons = 1, h = 0.5, S = 99)
#' }
#' @export
rd_compstable <- function(data, x, time, id, t_rd,
                          comparisons = NULL,
                          c = 0,
                          h = NULL,
                          q = NULL,
                          S  = 499L,
                          kernel = "triangular",
                          scheme = c("auto", "cs", "pc", "pv"),
                          ...) {
  scheme <- match.arg(scheme)

  # ---- input validation -------------------------------------------------------
  for (nm in base::c(x, time, id)) {
    if (!nm %in% names(data))
      stop("column '", nm, "' not found in `data`.")
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
  if (is.null(h)) {
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
  }

  # ---- per-pair analysis -------------------------------------------------------
  pairs_out <- list()
  q_used    <- list()   # per-pair q actually used (rule of thumb or fixed)

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
      u_trd <- apply(above_trd[, paste0("side_", u_periods), drop = FALSE],
                     1, paste, collapse = "")
      u_t0  <- apply(above_t0[,  paste0("side_", u_periods), drop = FALSE],
                     1, paste, collapse = "")
    } else {
      # P = 2: no shared other periods; u is empty
      u_trd <- rep("", nrow(above_trd))
      u_t0  <- rep("", nrow(above_t0))
    }

    # Handle NAs in partner sides (units not observed in a period)
    b_trd[is.na(b_trd)] <- NA_integer_
    b_t0[is.na(b_t0)]   <- NA_integer_

    # type string = paste(u, b)
    type_trd <- ifelse(is.na(b_trd), NA_character_,
                       paste0(u_trd, as.character(b_trd)))
    type_t0  <- ifelse(is.na(b_t0),  NA_character_,
                       paste0(u_t0,  as.character(b_t0)))

    # Reflected running variable
    xref_trd <- above_trd[[paste0("R_", t_rd_str)]] - c   # > 0
    xref_t0  <-  -(above_t0[[paste0("R_", t0_str)]] - c)  # < 0

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
    all_type_vals <- sort(unique(base::c(type_trd, type_t0)))
    n_types <- length(all_type_vals)

    # ---- auto-detect scheme ----------------------------------------------------
    use_scheme <- if (scheme != "auto") scheme else {
      if (n_both > 0L) "pv" else "cs"
    }

    # ---- (1) LL-Wald -----------------------------------------------------------
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
      fit <- tryCatch(
        rd_period(y = y_v, x = x_all, h = h, b = h, id = id_all,
                  c = 0, p = 1L, q = 2L, kernel = kernel),
        error = function(e) NULL
      )
      fits[[vi]] <- fit
      theta[vi]  <- if (is.null(fit)) NA_real_ else fit$D
    }

    # Build covariance matrix
    # The "+" side of the artificial cutoff = t_rd-above units
    # The "-" side = t_0-above units
    # Units in both appear in fits[[v]]$sides$`+`$id AND fits[[v]]$sides$`-`$id
    # Cross-type covariance on the same side = 0 (types partition units on each side)
    # Diagonal = within-type Var(D_v) = sum(g_+^2) + sum(g_-^2)
    # Off-diagonal (different types, same pair): 0 within side; cross-side =
    #   .match_sum on id (for pv scheme: cov_{+,+} + cov_{-,-} - cov_{+,-} - cov_{-,+})
    # But since types partition units: cross-type, same-side cov = 0.
    # So Sigma[vi, vj] for vi != vj is purely from cross-side terms when
    # a unit switches type between periods (only if type definitions differ).
    # In practice this CAN be non-zero because the partner side b differs:
    # a unit above in both periods has b_trd = its t0-side and b_t0 = its t_rd-side.
    # These are generally equal (both encode the same pair) so the types MATCH,
    # meaning a unit lands in the SAME type on both sides.
    # We still code the general formula.

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
          if (use_scheme == "pv") {
            cross <- .match_sum(fit_a$sides$`+`$id, fit_a$sides$`+`$g,
                                fit_a$sides$`-`$id, fit_a$sides$`-`$g)
            Sigma[a, a] <- fit_a$V_D - 2 * cross
          } else {
            Sigma[a, a] <- fit_a$V_D
          }
          next
        }

        # Cross-type off-diagonal: on the same side types partition units, so a
        # nonzero entry needs a unit landing in type a on one side and type b on
        # the other.  This is the scheme-combined cross covariance of the two
        # reflected fits (cs → 0, pc → same-side, pv → same-side − opposite).
        Sigma[a, b_idx] <- .cov_scheme(fit_a, fit_b, use_scheme)
      }
    }

    # Drop NA rows/cols
    ok_idx    <- which(!is.na(theta))
    theta_ok  <- theta[ok_idx]
    Sigma_ok  <- Sigma[ok_idx, ok_idx, drop = FALSE]
    ll_result <- .joint_wald(theta_ok, Sigma_ok)

    # ---- (2) Canay-Kamat permutation at the unit level -------------------------
    # Construct the 2q nearest observations at the unit level.
    # "+" side: t_rd-above units, sorted by x' = R_{i,t_rd} - c (ascending)
    # "-" side: t_0-above units, sorted by |x'| = R_{i,t_0} - c (ascending)

    # Sort each side by distance from artificial cutoff
    ord_trd <- order(xref_trd)
    ord_t0  <- order(-xref_t0)     # most negative (closest to 0) first

    # Number of nearest observations per side: the Canay-Kamat rule of thumb
    # on the pooled reflected sample (x_all is the signed artificial-cutoff
    # distance, type_all the partner (u,b) type) when q is NULL, otherwise the
    # user-supplied fixed q. Same fixed-q over-rejection exposure as rd_typecont.
    q_pair <- if (is.null(q)) .q_rot(x_all, type_all, 0) else as.integer(q)
    q_used[[pair_key]] <- q_pair

    q_trd <- min(q_pair, n_trd_obs)
    q_t0  <- min(q_pair, n_t0_obs)

    sel_trd  <- ord_trd[seq_len(q_trd)]
    sel_t0   <- ord_t0[seq_len(q_t0)]

    id_near_trd  <- id_trd[sel_trd]
    id_near_t0   <- id_t0[sel_t0]
    type_near_trd <- type_trd[sel_trd]
    type_near_t0  <- type_t0[sel_t0]

    # Build the pooled observation table for the permutation:
    # Each row is a unit-observation: (id, side [1=trd, 0=t0], type)
    perm_df <- data.frame(
      id   = base::c(id_near_trd,  id_near_t0),
      side = base::c(rep(1L, q_trd), rep(0L, q_t0)),
      type = base::c(type_near_trd, type_near_t0),
      stringsAsFactors = FALSE
    )

    # Unique units present in the permutation set
    unique_ids <- unique(perm_df$id)
    n_units    <- length(unique_ids)
    # The "+" side has q_trd observations; the "-" side has q_t0 observations
    # For unit-level permutation: permute which units go to which side.
    # Each unit is assigned entirely to one side OR keeps its actual side
    # (units appearing twice get both rows permuted together).
    # We implement: build a unit → side mapping; permute that mapping;
    # re-assign each row's "side" from its unit's permuted side.
    # For units that appear on both sides, they get ONE randomly chosen side.

    # Observed statistic: sum over types of |mean(type==v | side==1) - mean(type==v | side==0)|
    ck_obs <- 0
    for (v in all_type_vals) {
      g_trd_near <- as.numeric(type_near_trd == v)
      g_t0_near  <- as.numeric(type_near_t0  == v)
      if (length(g_trd_near) == 0L || length(g_t0_near) == 0L) next
      ck_obs <- ck_obs + abs(mean(g_trd_near) - mean(g_t0_near))
    }

    # Unit-level (partially-overlapping samples) permutation null distribution.
    #
    # The observed statistic compares the type distribution of the t_rd-near
    # group against the t0-near group. Under H0 (composition stable) the two
    # one-sided type distributions are equal, so the group labels are
    # exchangeable. A unit above the cutoff in BOTH periods contributes one
    # observation to each group (its t_rd-side type and its t0-side type), so the
    # two samples partially overlap. We use the partially-overlapping-samples
    # permutation \citep{amro2017permuting,derrick2022review}: units present in
    # both groups are PAIRED and their two labels swapped with probability 1/2;
    # units in only one group are freely relabelled while preserving group sizes.
    # With distinct units per period (cs scheme) there is no overlap and this
    # reduces to the standard two-sample label permutation.
    #
    # Each period contributes one row per unit, so a unit appears at most once
    # per side in the near set; map id -> its single type on each side.
    trd_type_of <- stats::setNames(type_near_trd, as.character(id_near_trd))
    t0_type_of  <- stats::setNames(type_near_t0,  as.character(id_near_t0))

    both_ids     <- intersect(id_near_trd, id_near_t0)
    only_trd_ids <- setdiff(id_near_trd, both_ids)
    only_t0_ids  <- setdiff(id_near_t0,  both_ids)

    both_trd_types <- unname(trd_type_of[as.character(both_ids)])   # paired, t_rd side
    both_t0_types  <- unname(t0_type_of[as.character(both_ids)])    # paired, t0 side
    free_trd_types <- unname(trd_type_of[as.character(only_trd_ids)])
    free_t0_types  <- unname(t0_type_of[as.character(only_t0_ids)])
    free_pool      <- base::c(free_trd_types, free_t0_types)
    n_free         <- length(free_pool)
    n_free_trd     <- length(free_trd_types)
    n_both_near    <- length(both_ids)   # paired units within the near window

    # statistic for a (trd-types, t0-types) split: sum_v |p_trd(v) - p_t0(v)|
    .ck_stat <- function(tt, t0) {
      if (length(tt) == 0L || length(t0) == 0L) return(NA_real_)
      s <- 0
      for (v in all_type_vals) s <- s + abs(mean(tt == v) - mean(t0 == v))
      s
    }

    ck_perm_dist <- replicate(S, {
      if (n_both_near > 0L) {
        swap <- stats::runif(n_both_near) < 0.5
        p_trd_both <- ifelse(swap, both_t0_types, both_trd_types)
        p_t0_both  <- ifelse(swap, both_trd_types, both_t0_types)
      } else {
        p_trd_both <- character(0); p_t0_both <- character(0)
      }
      if (n_free > 0L) {
        idx <- sample.int(n_free, n_free_trd)
        p_trd_free <- free_pool[idx]
        p_t0_free  <- free_pool[-idx]
      } else {
        p_trd_free <- character(0); p_t0_free <- character(0)
      }
      .ck_stat(base::c(p_trd_both, p_trd_free), base::c(p_t0_both, p_t0_free))
    })

    ck_perm_dist_clean <- ck_perm_dist[!is.na(ck_perm_dist)]
    ck_p <- if (length(ck_perm_dist_clean) == 0L) NA_real_ else {
      (1 + sum(ck_perm_dist_clean >= ck_obs)) / (length(ck_perm_dist_clean) + 1)
    }

    pairs_out[[pair_key]] <- list(
      ll_wald    = ll_result,
      ck_perm    = list(stat = ck_obs, p = ck_p),
      type_values = all_type_vals,
      scheme     = use_scheme,
      q          = q_pair,
      n_trd      = n_trd_obs,
      n_t0       = n_t0_obs,
      n_both     = n_both
    )
  }

  # ---- joint result across all pairs ------------------------------------------
  # Stack the theta/Sigma across pairs (assuming independence across pairs,
  # which holds exactly when each pair uses a distinct t0).
  # For the Wald: sum chi-sq statistics with summed df.
  joint_ll_stat <- 0
  joint_ll_df   <- 0L
  joint_ck_stat <- 0
  joint_ck_perms <- NULL

  for (pk in names(pairs_out)) {
    pr <- pairs_out[[pk]]
    joint_ll_stat <- joint_ll_stat + pr$ll_wald$stat
    joint_ll_df   <- joint_ll_df   + pr$ll_wald$df
    joint_ck_stat <- joint_ck_stat + pr$ck_perm$stat
  }
  joint_ll_p <- if (joint_ll_df == 0L) 1 else
    stats::pchisq(joint_ll_stat, df = joint_ll_df, lower.tail = FALSE)

  # Joint permutation p: re-run permutations jointly (independence across pairs
  # means we can add the stats from each pair's permutation distribution, but
  # that requires synchronised replication indices). Simpler: report the
  # Fisher combined p-value if there are multiple pairs; otherwise use the
  # single pair's p-value.
  all_ck_ps <- vapply(pairs_out, function(pr) pr$ck_perm$p, numeric(1))
  joint_ck_p <- if (length(all_ck_ps) == 0L) NA_real_ else if (length(all_ck_ps) == 1L) {
    all_ck_ps[[1L]]
  } else {
    # Fisher's combined test: -2 sum log(p), chi-sq with 2*K df
    valid_ps <- all_ck_ps[!is.na(all_ck_ps) & all_ck_ps > 0]
    if (length(valid_ps) == 0L) NA_real_ else {
      fisher_stat <- -2 * sum(log(valid_ps))
      stats::pchisq(fisher_stat, df = 2L * length(valid_ps), lower.tail = FALSE)
    }
  }

  structure(
    list(
      pairs = pairs_out,
      joint = list(
        ll_wald = list(stat = joint_ll_stat, df = joint_ll_df, p = joint_ll_p),
        ck_perm = list(stat = joint_ck_stat, p = joint_ck_p)
      ),
      meta = list(
        t_rd        = t_rd,
        comparisons = comparisons,
        h           = h,
        q           = if (is.null(q)) "rot" else as.integer(q),
        q_used      = q_used,
        S           = S,
        c           = c
      )
    ),
    class = "rd_compstable"
  )
}


#' @export
print.rd_compstable <- function(x, ...) {
  cat("Composition-stability test (Assumption A8)\n")
  cat(sprintf("  RD period: %s   Comparison periods: %s\n",
              x$meta$t_rd,
              paste(x$meta$comparisons, collapse = ", ")))
  q_lbl <- if (identical(x$meta$q, "rot")) "rule-of-thumb" else x$meta$q
  cat(sprintf("  h=%.4g   q=%s   S=%d\n\n",
              x$meta$h, q_lbl, x$meta$S))

  for (pk in names(x$pairs)) {
    pr <- x$pairs[[pk]]
    cat(sprintf("Pair %s  [scheme=%s  q=%d  n_trd=%d  n_t0=%d  n_both=%d]\n",
                pk, toupper(pr$scheme), pr$q, pr$n_trd, pr$n_t0, pr$n_both))
    cat(sprintf("  (1) LL-Wald [nec & suff]:  chi2(%.0f) = %.4f   p = %.4f\n",
                pr$ll_wald$df, pr$ll_wald$stat, pr$ll_wald$p))
    cat(sprintf("  (2) CK permutation [nec & suff]:  stat = %.4f   p = %.4f\n\n",
                pr$ck_perm$stat, pr$ck_perm$p))
  }

  if (length(x$pairs) > 1L) {
    cat("Joint (across all pairs):\n")
    cat(sprintf("  LL-Wald:  chi2(%.0f) = %.4f   p = %.4f\n",
                x$joint$ll_wald$df, x$joint$ll_wald$stat, x$joint$ll_wald$p))
    cat(sprintf("  CK perm (Fisher):  p = %.4f\n", x$joint$ck_perm$p))
  }

  cat("\nNOTE: Both tests are necessary AND sufficient for Assumption A8.\n")
  invisible(x)
}
