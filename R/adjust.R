# Composition-adjusted RD-DID estimator (Section 4.4 / Theorem 3, "thm:adjust").
#
# When composition stability (A8) is rejected but the within-type confounding
# trend (A10, ass:trend-cell) is credible, the composition effect is estimated
# and removed instead of assumed away. We use the WITHIN-TYPE (reweighting) form,
# which is the more interpretable read and is equivalent to the correction form
# under A7:
#
#   ATT_adj(t_rd) = sum_a pi_{t_rd,(+)}(v_comp = a) * ATT(t_rd | v_comp = a),
#   ATT(t_rd | v_comp = a) = D_{t_rd}(v_comp = a) - g0({ D_{t0}(v_rd = a) : t0 }).
#
# THE CROSSING (verified against code/simulations/s35_comp_adjust_sim.R):
#   * D_{t_rd}(v_comp = a) : RD-period jump among units grouped by their
#     COMPARISON-regime side (v_comp).
#   * D_{t0}(v_rd = a)     : comparison-period jump among units grouped by their
#     RD-period side (v_rd).  (Same label a, different grouping variable.)
#   * pi_{t_rd,(+)}(v_comp = a) : RD-period above-cutoff share of units with
#     v_comp = a, i.e. the above-cutoff intercept of a local-linear RD of the
#     indicator 1{v_comp = a}.
#
# The reweighting aggregate equals D_{t_rd} - sum_a pi(a) * g0({D_{t0}(v_rd=a)})
# up to the A7 adding-up gap, which we report.

#' Composition-adjusted RD-DID estimator (Theorem 3)
#'
#' Estimates the composition-adjusted ATT for a time-varying running variable,
#' the route of Section 4.4 / Theorem `thm:adjust` to take when composition
#' stability (A8, [rd_compstable()]) is rejected but the within-type confounding
#' trend (A10, [rd_trendcell()]) is credible. Reports the per-type RD-DID effects
#' `ATT(t_rd | v_comp = a)`, the RD-period composition shares, and their
#' share-weighted aggregate, alongside the unadjusted estimator.
#'
#' Standard errors are not computed here; a unit-level cluster bootstrap wrapper
#' is added separately (`se = "bootstrap"`).
#'
#' @param data long data frame, one row per unit-period.
#' @param y,x,time,id column names (strings).
#' @param t_rd RD-period value of `time`.
#' @param comparisons comparison-period values of `time`; `NULL` uses all others.
#' @param c cutoff (default 0).
#' @param h common bandwidth applied to every block (jumps and shares); `b` is
#'   the pilot for bias correction (defaults to `h`).
#' @param b pilot bandwidth (defaults to `h`).
#' @param weights `"constant"`, `"linear"`, or a numeric vector over
#'   `comparisons` (the trend `g0`, as in [rddid()]).
#' @param kernel `"triangular"` (default), `"epanechnikov"`, `"uniform"`.
#' @param min_n minimum observations per side of any block (default 10).
#' @param p,q point / bias-correction polynomial orders (default 1, 2).
#'
#' @return object of class `"rd_adjust"`.
#' @export
rd_adjust <- function(data, y, x, time, id, t_rd, comparisons = NULL,
                      c = 0, h, b = h, weights = "constant",
                      kernel = "triangular", min_n = 10L, p = 1L, q = 2L) {
  cl <- match.call()
  for (nm in c(y, x, time, id))
    if (!nm %in% names(data)) stop("column '", nm, "' not found in `data`.")
  if (missing(h) || is.null(h)) stop("supply a bandwidth `h`.")

  times <- data[[time]]
  if (is.null(comparisons)) comparisons <- sort(setdiff(unique(times), t_rd))
  comparisons <- sort(comparisons)
  w <- .rddid_weights(weights, comparisons, t_rd)

  # ---- per-unit sides: v_rd (side at t_rd) and v_comp (comparison-regime side) ----
  side_at <- function(period) {
    d <- data[times == period, , drop = FALSE]
    stats::setNames(as.integer(d[[x]] >= c), as.character(d[[id]]))
  }
  v_rd_map <- side_at(t_rd)
  # comparison side: must be constant across the supplied comparison periods
  # (the canonical panel shares one running variable across them). Take the
  # earliest comparison period and check the others agree where observed.
  v_comp_map <- side_at(comparisons[1L])
  if (length(comparisons) > 1L) {
    for (t0 in comparisons[-1L]) {
      m <- side_at(t0)
      shared <- intersect(names(m), names(v_comp_map))
      if (any(m[shared] != v_comp_map[shared], na.rm = TRUE))
        warning("comparison periods disagree on a unit's side; using the earliest ",
                "comparison period for v_comp. The within-type form assumes the ",
                "comparison regime shares one running variable.")
    }
  }

  ids_chr <- as.character(data[[id]])
  v_rd   <- v_rd_map[ids_chr]
  v_comp <- v_comp_map[ids_chr]

  # ---- one local-linear RD jump on a (period, group) subsample -> D, D_bc ----
  jump <- function(period, group_side, a) {
    keep <- times == period & !is.na(group_side) & group_side == a
    yy <- data[[y]][keep]; xx <- data[[x]][keep]; ii <- data[[id]][keep]
    if (sum(xx >= c, na.rm = TRUE) < min_n || sum(xx < c, na.rm = TRUE) < min_n)
      return(c(D = NA_real_, D_bc = NA_real_))
    f <- tryCatch(rd_period(yy, xx, h = h, b = b, id = ii, c = c, p = p, q = q,
                            kernel = kernel),
                  error = function(e) NULL)
    if (is.null(f)) return(c(D = NA_real_, D_bc = NA_real_))
    c(D = f$D, D_bc = f$D_bc)
  }

  # ---- RD-period above-cutoff share of v_comp = a (indicator-RD intercept) ----
  share_plus <- function(a) {
    keep <- times == t_rd & !is.na(v_comp)
    z  <- as.numeric(v_comp[keep] == a)
    xx <- data[[x]][keep]; ii <- data[[id]][keep]
    f <- tryCatch(rd_period(z, xx, h = h, b = b, id = ii, c = c, p = p, q = q,
                            kernel = kernel),
                  error = function(e) NULL)
    if (is.null(f)) return(c(pi = NA_real_, pi_bc = NA_real_))
    c(pi = f$sides[["+"]]$beta0, pi_bc = f$sides[["+"]]$beta0_bc)
  }

  labels <- c(0L, 1L)

  # ---- blocks ----
  # comparison jumps grouped by v_rd: D_{t0}(v_rd = a)
  Dt0 <- lapply(as.character(comparisons), function(t0)
    vapply(labels, function(a) jump(as.numeric(t0), v_rd, a), numeric(2)))
  names(Dt0) <- as.character(comparisons)            # each: 2 x 2 (rows D/D_bc, cols a)
  # RD-period jumps grouped by v_comp: D_{t_rd}(v_comp = a)
  Dtrd_a <- vapply(labels, function(a) jump(t_rd, v_comp, a), numeric(2))
  # shares pi_{t_rd,(+)}(v_comp = a)
  Pi <- vapply(labels, share_plus, numeric(2))       # rows pi/pi_bc, cols a
  # aggregate RD-period and comparison jumps (for unadjusted + cross-check)
  agg_jump <- function(period) {
    keep <- times == period
    yy <- data[[y]][keep]; xx <- data[[x]][keep]; ii <- data[[id]][keep]
    f <- rd_period(yy, xx, h = h, b = b, id = ii, c = c, p = p, q = q, kernel = kernel)
    c(D = f$D, D_bc = f$D_bc)
  }
  Dtrd <- agg_jump(t_rd)
  Dt0_agg <- vapply(as.character(comparisons), function(t0) agg_jump(as.numeric(t0)),
                    numeric(2))                       # rows D/D_bc, cols t0

  # ---- assemble, for each method m in {D (conv), D_bc} ----
  assemble <- function(row) {  # row = 1 (conv) or 2 (bc)
    pi_a   <- Pi[row, ]
    pi_a   <- pmin(pmax(pi_a, 0), 1)                  # clip shares to [0,1]
    dtrd_a <- Dtrd_a[row, ]
    # g0 of the comparison cell jumps, per label a: sum_t0 w_t0 D_{t0}(v_rd=a)
    g0_a <- vapply(labels, function(a) {
      ja <- vapply(as.character(comparisons), function(t0) Dt0[[t0]][row, which(labels == a)],
                   numeric(1))
      sum(w * ja)
    }, numeric(1))
    att_by_type <- dtrd_a - g0_a                      # ATT(t_rd | v_comp = a)
    att_adj     <- sum(pi_a * att_by_type)            # within-type aggregate
    d_trd       <- as.numeric(Dtrd[row])              # strip the "D"/"D_bc" name
    # reweighting cross-check: D_trd - sum_a pi(a) g0({D_t0(v_rd=a)})
    att_reweight <- d_trd - sum(pi_a * g0_a)
    # unadjusted: D_trd - sum_t0 w_t0 D_t0 (aggregate jumps)
    att_unadj   <- d_trd - sum(w * Dt0_agg[row, ])
    list(att_adj = att_adj, att_reweight = att_reweight, att_unadj = att_unadj,
         addingup_gap = att_adj - att_reweight,
         pi = pi_a, att_by_type = att_by_type, dtrd_a = dtrd_a, g0_a = g0_a)
  }
  conv <- assemble(1L); bc <- assemble(2L)

  within_type <- data.frame(
    label        = labels,
    pi_plus      = conv$pi,          pi_plus_bc   = bc$pi,
    D_trd_a      = conv$dtrd_a,      D_trd_a_bc   = bc$dtrd_a,
    g0_comp_a    = conv$g0_a,        g0_comp_a_bc = bc$g0_a,
    att_by_type  = conv$att_by_type, att_by_type_bc = bc$att_by_type
  )
  comp_jumps <- do.call(rbind, lapply(as.character(comparisons), function(t0)
    data.frame(period = as.numeric(t0), label = labels,
               D = Dt0[[t0]][1, ], D_bc = Dt0[[t0]][2, ])))

  structure(list(
    att_adj      = c(conv = conv$att_adj,      bc = bc$att_adj),
    att_unadj    = c(conv = conv$att_unadj,    bc = bc$att_unadj),
    att_reweight = c(conv = conv$att_reweight, bc = bc$att_reweight),
    addingup_gap = c(conv = conv$addingup_gap, bc = bc$addingup_gap),
    within_type  = within_type,
    comp_jumps   = comp_jumps,
    weights = w, weights_type = if (is.numeric(weights)) "custom" else weights,
    t_rd = t_rd, comparisons = comparisons, h = h, b = b, c = c, call = cl
  ), class = "rd_adjust")
}

#' @export
print.rd_adjust <- function(x, ...) {
  cat("Composition-adjusted RD-DID estimate (Theorem 3)\n")
  cat(sprintf("  RD period: %s   comparison periods: %s   trend: %s\n",
              x$t_rd, paste(x$comparisons, collapse = ", "), x$weights_type))
  cat(sprintf("  ATT_adj (conv) = %+.4g   ATT_adj (bc) = %+.4g\n",
              x$att_adj["conv"], x$att_adj["bc"]))
  cat(sprintf("  ATT_unadj (conv) = %+.4g   ATT_unadj (bc) = %+.4g\n",
              x$att_unadj["conv"], x$att_unadj["bc"]))
  cat(sprintf("  adding-up gap (within-type vs reweighting): conv %.2g, bc %.2g\n",
              x$addingup_gap["conv"], x$addingup_gap["bc"]))
  cat("\n  Per-type effects ATT(t_rd | v_comp = a):\n")
  wt <- x$within_type
  for (i in seq_len(nrow(wt)))
    cat(sprintf("    a=%d:  pi(+)=%.3f   D_trd(a)=%+.3g   g0_comp(a)=%+.3g   ATT(a)=%+.3g\n",
                wt$label[i], wt$pi_plus[i], wt$D_trd_a[i], wt$g0_comp_a[i], wt$att_by_type[i]))
  invisible(x)
}
