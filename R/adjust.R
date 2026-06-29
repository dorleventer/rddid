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
#     v_comp = a (the above-cutoff intercept of a local-linear RD of the
#     indicator 1{v_comp = a}).
#
# Standard errors are a unit-level (cluster) nonparametric bootstrap: resample
# units with replacement, rebuild the panel, recompute the whole plug-in. This is
# the paper's sanctioned equivalent to the closed-form a'Sigma a variance.

# ---- internal: one composition-adjusted fit on a given long data frame --------
# Returns the point quantities plus a flat `boot_vec` of the scalars we bootstrap.
.rd_adjust_fit <- function(data, y, x, time, id, t_rd, comparisons, w,
                           c, h, b, kernel, min_n, p, q) {
  times <- data[[time]]
  labels <- c(0L, 1L)

  side_at <- function(period) {
    d <- data[times == period, , drop = FALSE]
    stats::setNames(as.integer(d[[x]] >= c), as.character(d[[id]]))
  }
  v_rd_map   <- side_at(t_rd)
  v_comp_map <- side_at(comparisons[1L])
  ids_chr <- as.character(data[[id]])
  v_rd   <- v_rd_map[ids_chr]
  v_comp <- v_comp_map[ids_chr]

  jump <- function(period, group_side, a) {
    keep <- times == period & !is.na(group_side) & group_side == a
    yy <- data[[y]][keep]; xx <- data[[x]][keep]; ii <- data[[id]][keep]
    if (sum(xx >= c, na.rm = TRUE) < min_n || sum(xx < c, na.rm = TRUE) < min_n)
      return(c(D = NA_real_, D_bc = NA_real_))
    f <- tryCatch(rd_period(yy, xx, h = h, b = b, id = ii, c = c, p = p, q = q,
                            kernel = kernel), error = function(e) NULL)
    if (is.null(f)) return(c(D = NA_real_, D_bc = NA_real_))
    c(D = f$D, D_bc = f$D_bc)
  }
  share_plus <- function(a) {
    keep <- times == t_rd & !is.na(v_comp)
    z  <- as.numeric(v_comp[keep] == a)
    xx <- data[[x]][keep]; ii <- data[[id]][keep]
    f <- tryCatch(rd_period(z, xx, h = h, b = b, id = ii, c = c, p = p, q = q,
                            kernel = kernel), error = function(e) NULL)
    if (is.null(f)) return(c(pi = NA_real_, pi_bc = NA_real_))
    c(pi = f$sides[["+"]]$beta0, pi_bc = f$sides[["+"]]$beta0_bc)
  }
  agg_jump <- function(period) {
    keep <- times == period
    f <- tryCatch(rd_period(data[[y]][keep], data[[x]][keep], h = h, b = b,
                            id = data[[id]][keep], c = c, p = p, q = q, kernel = kernel),
                  error = function(e) NULL)
    if (is.null(f)) return(c(D = NA_real_, D_bc = NA_real_))
    c(D = f$D, D_bc = f$D_bc)
  }

  Dt0    <- lapply(as.character(comparisons), function(t0)
              vapply(labels, function(a) jump(as.numeric(t0), v_rd, a), numeric(2)))
  names(Dt0) <- as.character(comparisons)
  Dtrd_a <- vapply(labels, function(a) jump(t_rd, v_comp, a), numeric(2))
  Pi     <- vapply(labels, share_plus, numeric(2))
  Dtrd   <- agg_jump(t_rd)
  Dt0_agg <- vapply(as.character(comparisons), function(t0) agg_jump(as.numeric(t0)),
                    numeric(2))

  assemble <- function(row) {
    pi_a   <- pmin(pmax(Pi[row, ], 0), 1)
    dtrd_a <- Dtrd_a[row, ]
    g0_a   <- vapply(labels, function(a)
                sum(w * vapply(as.character(comparisons),
                               function(t0) Dt0[[t0]][row, which(labels == a)], numeric(1))),
                numeric(1))
    att_by_type <- dtrd_a - g0_a
    d_trd  <- as.numeric(Dtrd[row])
    list(att_adj = sum(pi_a * att_by_type),
         att_reweight = d_trd - sum(pi_a * g0_a),
         att_unadj = d_trd - sum(w * Dt0_agg[row, ]),
         pi = pi_a, att_by_type = att_by_type, dtrd_a = dtrd_a, g0_a = g0_a)
  }
  conv <- assemble(1L); bc <- assemble(2L)

  boot_vec <- c(
    att_adj.conv = conv$att_adj,   att_adj.bc = bc$att_adj,
    att_unadj.conv = conv$att_unadj, att_unadj.bc = bc$att_unadj,
    att_reweight.conv = conv$att_reweight,
    stats::setNames(conv$att_by_type, paste0("att.a", labels, ".conv")),
    stats::setNames(bc$att_by_type,   paste0("att.a", labels, ".bc")),
    stats::setNames(conv$pi,          paste0("pi.a", labels))
  )

  list(conv = conv, bc = bc, labels = labels, Dt0 = Dt0, boot_vec = boot_vec)
}

#' Composition-adjusted RD-DID estimator (Theorem 3)
#'
#' Estimates the composition-adjusted ATT for a time-varying running variable,
#' the route of Section 4.4 / Theorem `thm:adjust` to take when composition
#' stability (A8, [rd_compstable()]) is rejected but the within-type confounding
#' trend (A10, [rd_trendcell()]) is credible. Reports the per-type RD-DID effects
#' `ATT(t_rd | v_comp = a)`, the RD-period composition shares, and their
#' share-weighted aggregate, alongside the unadjusted estimator. Standard errors
#' are a unit-level cluster bootstrap.
#'
#' @param data long data frame, one row per unit-period.
#' @param y,x,time,id column names (strings).
#' @param t_rd RD-period value of `time`.
#' @param comparisons comparison-period values of `time`; `NULL` uses all others.
#' @param c cutoff (default 0).
#' @param h common bandwidth for every block (jumps and shares); cells are thin,
#'   so a bandwidth wider than the aggregate is appropriate.
#' @param b pilot bandwidth for bias correction (defaults to `h`).
#' @param weights `"constant"`, `"linear"`, or a numeric vector over
#'   `comparisons` (the trend `g0`, as in [rddid()]).
#' @param kernel `"triangular"` (default), `"epanechnikov"`, `"uniform"`.
#' @param min_n minimum observations per side of any block (default 10).
#' @param p,q point / bias-correction polynomial orders (default 1, 2).
#' @param se `"bootstrap"` (default) or `"none"`.
#' @param B bootstrap replications (default 500).
#'
#' @return object of class `"rd_adjust"`.
#' @export
rd_adjust <- function(data, y, x, time, id, t_rd, comparisons = NULL,
                      c = 0, h, b = h, weights = "constant",
                      kernel = "triangular", min_n = 10L, p = 1L, q = 2L,
                      se = c("bootstrap", "none"), B = 500L) {
  cl <- match.call()
  se <- match.arg(se)
  for (nm in c(y, x, time, id))
    if (!nm %in% names(data)) stop("column '", nm, "' not found in `data`.")
  if (missing(h) || is.null(h)) stop("supply a bandwidth `h`.")

  times <- data[[time]]
  if (is.null(comparisons)) comparisons <- sort(setdiff(unique(times), t_rd))
  comparisons <- sort(comparisons)
  w <- .rddid_weights(weights, comparisons, t_rd)

  # comparison-side consistency check (the within-type form assumes the
  # comparison regime shares one running variable across its periods)
  if (length(comparisons) > 1L) {
    s1 <- stats::setNames(as.integer(data[[x]][times == comparisons[1L]] >= c),
                          as.character(data[[id]][times == comparisons[1L]]))
    for (t0 in comparisons[-1L]) {
      m <- stats::setNames(as.integer(data[[x]][times == t0] >= c),
                           as.character(data[[id]][times == t0]))
      sh <- intersect(names(m), names(s1))
      if (length(sh) && any(m[sh] != s1[sh], na.rm = TRUE)) {
        warning("comparison periods disagree on a unit's side; using the earliest ",
                "comparison period for v_comp.")
        break
      }
    }
  }

  fit <- .rd_adjust_fit(data, y, x, time, id, t_rd, comparisons, w,
                        c, h, b, kernel, min_n, p, q)
  conv <- fit$conv; bc <- fit$bc; labels <- fit$labels

  # ---- bootstrap SE (resample units, rebuild panel, recompute) ----
  se_vec <- NULL
  if (se == "bootstrap") {
    idx_by_id <- split(seq_len(nrow(data)), as.character(data[[id]]))
    U <- names(idx_by_id)
    bmat <- matrix(NA_real_, nrow = B, ncol = length(fit$boot_vec),
                   dimnames = list(NULL, names(fit$boot_vec)))
    for (rb in seq_len(B)) {
      drawn <- sample(U, length(U), replace = TRUE)
      rows_list <- idx_by_id[drawn]
      bd <- data[unlist(rows_list, use.names = FALSE), , drop = FALSE]
      bd[[id]] <- rep(seq_along(drawn), lengths(rows_list))   # fresh cluster ids
      bf <- tryCatch(.rd_adjust_fit(bd, y, x, time, id, t_rd, comparisons, w,
                                    c, h, b, kernel, min_n, p, q),
                     error = function(e) NULL)
      if (!is.null(bf)) bmat[rb, ] <- bf$boot_vec
    }
    se_vec <- apply(bmat, 2L, stats::sd, na.rm = TRUE)
    attr(se_vec, "n_ok") <- sum(stats::complete.cases(bmat[, "att_adj.conv", drop = FALSE]))
  }
  getse <- function(nm) if (is.null(se_vec)) NA_real_ else unname(se_vec[nm])

  within_type <- data.frame(
    label          = labels,
    pi_plus        = conv$pi,          pi_plus_bc     = bc$pi,
    D_trd_a        = conv$dtrd_a,      D_trd_a_bc     = bc$dtrd_a,
    g0_comp_a      = conv$g0_a,        g0_comp_a_bc   = bc$g0_a,
    att_by_type    = conv$att_by_type, att_by_type_bc = bc$att_by_type,
    att_by_type_se = getse(paste0("att.a", labels, ".conv")),
    att_by_type_bc_se = getse(paste0("att.a", labels, ".bc")),
    pi_plus_se     = getse(paste0("pi.a", labels))
  )
  comp_jumps <- do.call(rbind, lapply(as.character(comparisons), function(t0)
    data.frame(period = as.numeric(t0), label = labels,
               D = fit$Dt0[[t0]][1, ], D_bc = fit$Dt0[[t0]][2, ])))

  structure(list(
    att_adj      = c(conv = conv$att_adj,      bc = bc$att_adj),
    att_adj_se   = c(conv = getse("att_adj.conv"),   bc = getse("att_adj.bc")),
    att_unadj    = c(conv = conv$att_unadj,    bc = bc$att_unadj),
    att_unadj_se = c(conv = getse("att_unadj.conv"), bc = getse("att_unadj.bc")),
    att_reweight = c(conv = conv$att_reweight, bc = bc$att_reweight),
    addingup_gap = c(conv = conv$att_adj - conv$att_reweight,
                     bc   = bc$att_adj - bc$att_reweight),
    within_type  = within_type, comp_jumps = comp_jumps,
    weights = w, weights_type = if (is.numeric(weights)) "custom" else weights,
    t_rd = t_rd, comparisons = comparisons, h = h, b = b, c = c,
    se = se, B = if (se == "bootstrap") B else 0L,
    n_boot_ok = if (is.null(se_vec)) NA_integer_ else attr(se_vec, "n_ok"),
    call = cl
  ), class = "rd_adjust")
}

#' @export
print.rd_adjust <- function(x, ...) {
  cat("Composition-adjusted RD-DID estimate (Theorem 3)\n")
  cat(sprintf("  RD period: %s   comparison periods: %s   trend: %s\n",
              x$t_rd, paste(x$comparisons, collapse = ", "), x$weights_type))
  sefmt <- function(v, s) if (is.na(s)) sprintf("%+.4g", v) else sprintf("%+.4g (%.3g)", v, s)
  cat(sprintf("  ATT_adj   (conv) = %s   (bc) = %s\n",
              sefmt(x$att_adj["conv"], x$att_adj_se["conv"]),
              sefmt(x$att_adj["bc"],   x$att_adj_se["bc"])))
  cat(sprintf("  ATT_unadj (conv) = %s   (bc) = %s\n",
              sefmt(x$att_unadj["conv"], x$att_unadj_se["conv"]),
              sefmt(x$att_unadj["bc"],   x$att_unadj_se["bc"])))
  if (x$se == "bootstrap")
    cat(sprintf("  SE: unit cluster bootstrap, B=%d (%d ok)\n", x$B, x$n_boot_ok))
  cat(sprintf("  adding-up gap (within-type vs reweighting): conv %.2g, bc %.2g\n",
              x$addingup_gap["conv"], x$addingup_gap["bc"]))
  cat("\n  Per-type effects ATT(t_rd | v_comp = a):\n")
  wt <- x$within_type
  for (i in seq_len(nrow(wt)))
    cat(sprintf("    a=%d:  pi(+)=%.3f   D_trd(a)=%+.3g   g0_comp(a)=%+.3g   ATT(a)=%s\n",
                wt$label[i], wt$pi_plus[i], wt$D_trd_a[i], wt$g0_comp_a[i],
                sefmt(wt$att_by_type[i], wt$att_by_type_se[i])))
  invisible(x)
}
