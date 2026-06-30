# Correction-form family of composition-adjusted RD-DID ATT estimators.
#
# This is the "combine layer": it assembles the family of adjusted ATT
# estimators from the three per-block primitives the package already provides --
#   D_t : the period discontinuity            (rd_period, R/rd_period.R)
#   S_t : within-period composition term      (Prop sdisc, R/sadjust.R)
#   C   : cross-period composition term        (Prop cdisc, R/cterm.R)
# -- and supplies ALL standard errors from ONE shared unit-level cluster
# bootstrap. Selecting the bandwidths once and recomputing every primitive on the
# SAME resampled draws carries the cross-primitive covariances automatically, so
# the standard errors of the assembled combinations (which mix D, S and C) are
# valid. Bootstrapping the pieces in separate loops would lose those covariances.
#
# Headline form is the CORRECTION form of Theorem `thm:adjust`:
#   ATT_unadj = D_{t_rd} - g0({D_{t0}})
#   ATT_s     = (D_{t_rd} - S_{t_rd}) - g0({D_{t0} - S_{t0}})
#   ATT_c     = D_{t_rd} - g0({D_{t0}}) - Cagg
#   ATT_sc    = (D_{t_rd} - S_{t_rd}) - g0({D_{t0} - S_{t0}}) - Cagg
# with g0({x_{t0}}) = sum_{t0} w_{t0} x_{t0} the trend extrapolation (the package
# trend weights of .rddid_weights) and Cagg = sum_{t0} w_{t0} C_{t0,t_rd}.
#
# Pairwise (P = 2) binary-type design, mirroring R/sadjust.R and R/cterm.R: each
# comparison period t0 forms the pair {t0, t_rd}. S_{t0} uses partner = t_rd;
# S_{t_rd} is the trend-weighted average over partners, S_{t_rd} = sum_{t0} w_{t0}
# S_{t_rd}(partner = t0), so the single-comparison case (w = 1) reduces to exactly
# the two-period hand formulas
#   adj-c  = D_{t_rd} - D_{t0} - C_{t0,t_rd},
#   adj-sc = (D_{t_rd} - S_{t_rd}) - (D_{t0} - S_{t0}) - C_{t0,t_rd}.
#
# This is the correction form and is independent of rd_adjust's within-type
# att_adj (the reweighting form), which is left as-is (legacy).

# ---- internal: period-t discontinuity D_t (conv + bc) at a FIXED bandwidth ----
.rd_att_Dt <- function(data, y, x, time, period, bw, c, kernel) {
  keep <- data[[time]] == period
  f <- tryCatch(
    rd_period(y = data[[y]][keep], x = data[[x]][keep],
              h = unname(bw["h"]), b = unname(bw["b"]), c = c, kernel = kernel),
    error = function(e) NULL)
  if (is.null(f)) c(NA_real_, NA_real_) else c(f$D, f$D_bc)
}

# ---- internal: one full assembly on a long data frame at FIXED bandwidths ------
# Reusable point helper: rebuilds every primitive from `data`, so the shared
# bootstrap calls it once per replication on the resampled panel. Returns a flat
# named numeric vector of every reported quantity (point estimate). The bootstrap
# SE of each quantity is the SD of its column over the replications, which is what
# makes the combination SEs carry the right covariances.
.rd_att_point <- function(data, y, x, time, id, comparisons, t_rd, periods,
                          w, c, bws, kernel, labels) {
  cp  <- as.character(periods)
  ct0 <- as.character(comparisons)
  rdc <- as.character(t_rd)

  # period discontinuities D_t (conv + bc) for every period
  Dmat <- vapply(periods, function(p)
    .rd_att_Dt(data, y, x, time, p, bws$D[[as.character(p)]], c, kernel),
    numeric(2))
  dimnames(Dmat) <- list(c("conv", "bc"), cp)

  # S_{t0} (focal = t0, partner = t_rd) for each comparison
  Scomp <- vapply(comparisons, function(t0) {
    fr <- .rd_sadjust_frame(data, y, x, time, id, t0, t_rd, c)
    o  <- tryCatch(.rd_sadjust_one(fr, bws$S_comp[[as.character(t0)]], c, kernel),
                   error = function(e) rep(NA_real_, 4L))
    c(o[["S"]], o[["S_bc"]])
  }, numeric(2))
  dimnames(Scomp) <- list(c("conv", "bc"), ct0)

  # S_{t_rd}(partner = t0) for each comparison, then trend-weighted aggregate
  Srd <- vapply(comparisons, function(t0) {
    fr <- .rd_sadjust_frame(data, y, x, time, id, t_rd, t0, c)
    o  <- tryCatch(.rd_sadjust_one(fr, bws$S_rd[[as.character(t0)]], c, kernel),
                   error = function(e) rep(NA_real_, 4L))
    c(o[["S"]], o[["S_bc"]])
  }, numeric(2))
  dimnames(Srd) <- list(c("conv", "bc"), ct0)
  S_rd <- c(conv = sum(w * Srd["conv", ]), bc = sum(w * Srd["bc", ]))

  # C_{t0,t_rd} (t = t0 comparison role, s = t_rd RD role) for each comparison
  Cmat <- vapply(comparisons, function(t0) {
    o <- tryCatch(
      .rd_c_point(data, y, x, time, id, t0, t_rd, c,
                  bws$C[[as.character(t0)]], kernel, labels),
      error = function(e) c(NA_real_, NA_real_))
    c(o[["C"]], o[["C_bc"]])
  }, numeric(2))
  dimnames(Cmat) <- list(c("conv", "bc"), ct0)
  Cagg <- c(conv = sum(w * Cmat["conv", ]), bc = sum(w * Cmat["bc", ]))

  # per-period S aligned over all periods (comparisons carry S_{t0}; t_rd the agg)
  Sconv <- stats::setNames(numeric(length(periods)), cp)
  Sbc   <- stats::setNames(numeric(length(periods)), cp)
  Sconv[ct0] <- Scomp["conv", ]; Sbc[ct0] <- Scomp["bc", ]
  Sconv[rdc] <- S_rd[["conv"]];  Sbc[rdc] <- S_rd[["bc"]]

  adjsD_conv <- Dmat["conv", ] - Sconv
  adjsD_bc   <- Dmat["bc", ]   - Sbc

  # ATT assembly (correction form)
  Drd_conv <- Dmat["conv", rdc]; Drd_bc <- Dmat["bc", rdc]
  g0D_conv <- sum(w * Dmat["conv", ct0]); g0D_bc <- sum(w * Dmat["bc", ct0])
  g0adjs_conv <- sum(w * (Dmat["conv", ct0] - Scomp["conv", ]))
  g0adjs_bc   <- sum(w * (Dmat["bc", ct0]   - Scomp["bc", ]))
  adjsrd_conv <- Drd_conv - S_rd[["conv"]]; adjsrd_bc <- Drd_bc - S_rd[["bc"]]

  att_unadj_conv <- Drd_conv - g0D_conv;     att_unadj_bc <- Drd_bc - g0D_bc
  att_s_conv     <- adjsrd_conv - g0adjs_conv; att_s_bc   <- adjsrd_bc - g0adjs_bc
  att_c_conv     <- att_unadj_conv - Cagg[["conv"]]; att_c_bc <- att_unadj_bc - Cagg[["bc"]]
  att_sc_conv    <- att_s_conv - Cagg[["conv"]];     att_sc_bc <- att_s_bc - Cagg[["bc"]]

  c(
    stats::setNames(Dmat["conv", ], paste0("D.", cp, ".conv")),
    stats::setNames(Dmat["bc", ],   paste0("D.", cp, ".bc")),
    stats::setNames(Sconv,          paste0("S.", cp, ".conv")),
    stats::setNames(Sbc,            paste0("S.", cp, ".bc")),
    stats::setNames(adjsD_conv,     paste0("adj_s_D.", cp, ".conv")),
    stats::setNames(adjsD_bc,       paste0("adj_s_D.", cp, ".bc")),
    stats::setNames(Cmat["conv", ], paste0("C.", ct0, ".conv")),
    stats::setNames(Cmat["bc", ],   paste0("C.", ct0, ".bc")),
    Cagg.conv = unname(Cagg[["conv"]]), Cagg.bc = unname(Cagg[["bc"]]),
    ATT.unadj.conv = att_unadj_conv, ATT.unadj.bc = att_unadj_bc,
    ATT.s.conv     = att_s_conv,     ATT.s.bc     = att_s_bc,
    ATT.c.conv     = att_c_conv,     ATT.c.bc     = att_c_bc,
    ATT.sc.conv    = att_sc_conv,    ATT.sc.bc    = att_sc_bc
  )
}

#' Correction-form composition-adjusted RD-DID ATT family (Theorem 3)
#'
#' Assembles the family of composition-adjusted average-treatment-on-the-treated
#' (ATT) estimators in their CORRECTION form from the three per-block primitives
#' of the package -- the period discontinuity \eqn{D_t} ([rd_period()]), the
#' within-period composition term \eqn{S_t} (Proposition `sdisc`, [rd_sadjust()]),
#' and the cross-period composition term \eqn{C} (Proposition `cdisc`,
#' [rd_c()]) -- and exposes the full \eqn{D}/\eqn{S}/\eqn{C} decomposition. The
#' headline quantities are
#' \deqn{\mathrm{ATT}_{\mathrm{unadj}} = D_{t_{RD}} - g_0(\{D_{t_0}\}),}
#' \deqn{\mathrm{ATT}_{s} = (D_{t_{RD}}-S_{t_{RD}}) - g_0(\{D_{t_0}-S_{t_0}\}),}
#' \deqn{\mathrm{ATT}_{c} = D_{t_{RD}} - g_0(\{D_{t_0}\}) - C_{\mathrm{agg}},}
#' \deqn{\mathrm{ATT}_{sc} = (D_{t_{RD}}-S_{t_{RD}}) - g_0(\{D_{t_0}-S_{t_0}\}) -
#'   C_{\mathrm{agg}},}
#' where \eqn{g_0(\{x_{t_0}\}) = \sum_{t_0} w_{t_0} x_{t_0}} is the trend
#' extrapolation with the package trend weights and
#' \eqn{C_{\mathrm{agg}} = \sum_{t_0} w_{t_0} C_{t_0,t_{RD}}}.
#'
#' The estimator uses the pairwise (\eqn{P = 2}) binary-type design: each
#' comparison period `t0` forms the pair `{t0, t_rd}`. \eqn{S_{t_0}} uses
#' partner `t_rd`; \eqn{S_{t_{RD}}} is the trend-weighted average over partners,
#' so the single-comparison case (\eqn{w = 1}) reduces exactly to
#' \eqn{\mathrm{adj}\text{-}c = D_{t_{RD}} - D_{t_0} - C_{t_0,t_{RD}}} and
#' \eqn{\mathrm{adj}\text{-}sc = (D_{t_{RD}}-S_{t_{RD}}) - (D_{t_0}-S_{t_0}) -
#' C_{t_0,t_{RD}}}.
#'
#' All standard errors come from ONE shared unit-level cluster bootstrap: the
#' per-block CCT bandwidths are chosen once on the original sample and held fixed;
#' each replication resamples units with replacement, rebuilds the panel with
#' fresh sequential cluster ids, and recomputes every primitive at the fixed
#' bandwidths before re-forming every ATT. Because all quantities are recomputed
#' on the same draws, the SD-over-replications standard errors carry the
#' cross-primitive covariances, so the combination SEs are valid.
#'
#' This is the correction form and does not reimplement or depend on the
#' within-type (reweighting) `att_adj` of [rd_adjust()], which is retained as-is.
#'
#' @param data long data frame, one row per unit-period.
#' @param y,x,time,id column names (strings) for outcome, running variable,
#'   period, and unit id.
#' @param comparisons comparison-period values of `time`; the pair `{t0, t_rd}`
#'   is formed for each.
#' @param t_rd RD-period value of `time`.
#' @param c cutoff (default 0).
#' @param h optional fixed numeric bandwidth used for every block; if `NULL`
#'   (default) a per-block CCT bandwidth is selected once on the original sample.
#' @param bwselect bandwidth selector label when `h` is `NULL` (default
#'   `"cct"`); retained for forward compatibility.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param trend trend extrapolation weights: `"constant"` (default), `"linear"`,
#'   or a numeric vector with one entry per comparison period (as in
#'   [rd_adjust()]).
#' @param se `"bootstrap"` (default) or `"none"`.
#' @param B bootstrap replications (default 499).
#'
#' @return An object of class `"rd_att"`: a list with
#'   \itemize{
#'     \item `att` -- the ATT table (one row per `estimator` in
#'       `{unadj, s, c, sc}`, columns `conv`, `bc`, `se`, `bc_se`);
#'     \item `per_period` -- the per-period decomposition table (`period`,
#'       `role`, `D`, `D_bc`, `D_se`, `D_bc_se`, `S`, `S_bc`, `S_se`, `S_bc_se`,
#'       `adj_s_D`, `adj_s_D_bc`, `adj_s_D_se`, `adj_s_D_bc_se`);
#'     \item `cterm` -- the cross-period composition table (`comparison`, `C`,
#'       `C_bc`, `C_se`, `C_bc_se`);
#'     \item `cagg` -- the aggregate `Cagg` (conv + bc + SE);
#'     \item meta fields (`t_rd`, `comparisons`, `weights`, `trend`, `c`,
#'       `kernel`, `bwselect`, `se`, `B`, `n_boot_ok`, `call`).
#'   }
#' @export
rd_att <- function(data, y, x, time, id, comparisons, t_rd, c = 0,
                   h = NULL, bwselect = "cct", kernel = "triangular",
                   trend = "constant", se = c("bootstrap", "none"), B = 499L) {
  cl <- match.call()
  se <- match.arg(se)
  for (nm in c(y, x, time, id))
    if (!nm %in% names(data)) stop("column '", nm, "' not found in `data`.")
  comparisons <- sort(unique(comparisons))
  if (length(comparisons) < 1L) stop("supply at least one comparison period.")
  if (t_rd %in% comparisons) stop("`t_rd` must not be in `comparisons`.")
  labels  <- c(0L, 1L)
  periods <- sort(unique(c(comparisons, t_rd)))
  w <- .rddid_weights(trend, comparisons, t_rd)

  # ---- bandwidths fixed once on the original sample --------------------------
  fixed_pair <- if (!is.null(h)) c(h = h, b = h) else NULL
  fixed_set  <- function(keys)
    stats::setNames(rep(list(fixed_pair), length(keys)), keys)

  bw_D <- stats::setNames(lapply(periods, function(p) {
    if (!is.null(h)) return(fixed_pair)
    keep <- data[[time]] == p
    rd_bw_cct(as.numeric(data[[y]][keep]), as.numeric(data[[x]][keep]),
              c = c, kernel = kernel)
  }), as.character(periods))

  bw_S_comp <- stats::setNames(lapply(comparisons, function(t0) {
    fr <- .rd_sadjust_frame(data, y, x, time, id, t0, t_rd, c)
    if (!is.null(h)) list(pi = fixed_pair, ab = fixed_pair, be = fixed_pair,
                          D = fixed_pair)
    else .rd_sadjust_pointbw(fr, c, kernel)
  }), as.character(comparisons))

  bw_S_rd <- stats::setNames(lapply(comparisons, function(t0) {
    fr <- .rd_sadjust_frame(data, y, x, time, id, t_rd, t0, c)
    if (!is.null(h)) list(pi = fixed_pair, ab = fixed_pair, be = fixed_pair,
                          D = fixed_pair)
    else .rd_sadjust_pointbw(fr, c, kernel)
  }), as.character(comparisons))

  bw_C <- stats::setNames(lapply(comparisons, function(t0) {
    frames0 <- .rd_c_frames(data, y, x, time, id, t0, t_rd, c)
    if (!is.null(h)) {
      per <- fixed_set(as.character(labels))
      list(D = per, pis = per, pit = per)
    } else .rd_c_pointbw(frames0, labels, c, kernel)
  }), as.character(comparisons))

  bws <- list(D = bw_D, S_comp = bw_S_comp, S_rd = bw_S_rd, C = bw_C)

  # ---- point estimate --------------------------------------------------------
  pt <- .rd_att_point(data, y, x, time, id, comparisons, t_rd, periods,
                      w, c, bws, kernel, labels)

  # ---- shared unit-level cluster bootstrap -----------------------------------
  se_vec <- stats::setNames(rep(NA_real_, length(pt)), names(pt))
  n_boot_ok <- NA_integer_
  if (se == "bootstrap") {
    idx_by_id <- split(seq_len(nrow(data)), as.character(data[[id]]))
    U <- names(idx_by_id)
    bmat <- matrix(NA_real_, nrow = B, ncol = length(pt),
                   dimnames = list(NULL, names(pt)))
    for (rb in seq_len(B)) {
      drawn <- sample(U, length(U), replace = TRUE)
      rows_list <- idx_by_id[drawn]
      bd <- data[unlist(rows_list, use.names = FALSE), , drop = FALSE]
      bd[[id]] <- rep(seq_along(drawn), lengths(rows_list))   # fresh cluster ids
      bmat[rb, ] <- tryCatch(
        .rd_att_point(bd, y, x, time, id, comparisons, t_rd, periods,
                      w, c, bws, kernel, labels),
        error = function(e) rep(NA_real_, length(pt)))
    }
    se_vec <- apply(bmat, 2L, stats::sd, na.rm = TRUE)
    n_boot_ok <- sum(stats::complete.cases(bmat[, "ATT.unadj.conv", drop = FALSE]))
  }

  g <- function(nm) unname(pt[nm])
  gs <- function(nm) unname(se_vec[nm])
  cp  <- as.character(periods)
  ct0 <- as.character(comparisons)
  role <- ifelse(periods == t_rd, "rd", "comparison")

  per_period <- data.frame(
    period        = periods,
    role          = role,
    D             = g(paste0("D.", cp, ".conv")),
    D_bc          = g(paste0("D.", cp, ".bc")),
    D_se          = gs(paste0("D.", cp, ".conv")),
    D_bc_se       = gs(paste0("D.", cp, ".bc")),
    S             = g(paste0("S.", cp, ".conv")),
    S_bc          = g(paste0("S.", cp, ".bc")),
    S_se          = gs(paste0("S.", cp, ".conv")),
    S_bc_se       = gs(paste0("S.", cp, ".bc")),
    adj_s_D       = g(paste0("adj_s_D.", cp, ".conv")),
    adj_s_D_bc    = g(paste0("adj_s_D.", cp, ".bc")),
    adj_s_D_se    = gs(paste0("adj_s_D.", cp, ".conv")),
    adj_s_D_bc_se = gs(paste0("adj_s_D.", cp, ".bc")),
    stringsAsFactors = FALSE, row.names = NULL
  )

  cterm <- data.frame(
    comparison = comparisons,
    C          = g(paste0("C.", ct0, ".conv")),
    C_bc       = g(paste0("C.", ct0, ".bc")),
    C_se       = gs(paste0("C.", ct0, ".conv")),
    C_bc_se    = gs(paste0("C.", ct0, ".bc")),
    stringsAsFactors = FALSE, row.names = NULL
  )

  ests <- c("unadj", "s", "c", "sc")
  att <- data.frame(
    estimator = ests,
    conv      = g(paste0("ATT.", ests, ".conv")),
    bc        = g(paste0("ATT.", ests, ".bc")),
    se        = gs(paste0("ATT.", ests, ".conv")),
    bc_se     = gs(paste0("ATT.", ests, ".bc")),
    stringsAsFactors = FALSE, row.names = NULL
  )

  cagg <- c(conv = g("Cagg.conv"), bc = g("Cagg.bc"),
            se = gs("Cagg.conv"), bc_se = gs("Cagg.bc"))

  structure(list(
    att          = att,
    per_period   = per_period,
    cterm        = cterm,
    cagg         = cagg,
    t_rd         = t_rd,
    comparisons  = comparisons,
    weights      = w,
    trend        = if (is.numeric(trend)) "custom" else trend,
    c            = c,
    h            = h,
    bwselect     = if (is.null(h)) bwselect else "fixed",
    kernel       = kernel,
    se           = se,
    B            = if (se == "bootstrap") B else 0L,
    n_boot_ok    = n_boot_ok,
    call         = cl
  ), class = "rd_att")
}

#' @export
print.rd_att <- function(x, ...) {
  cat("Composition-adjusted RD-DID ATT family (correction form, Theorem 3)\n")
  cat(sprintf("  RD period: %s   comparison periods: %s   trend: %s   bw: %s\n",
              x$t_rd, paste(x$comparisons, collapse = ", "), x$trend, x$bwselect))
  sefmt <- function(v, s) if (is.na(s)) sprintf("%+.4g", v) else
    sprintf("%+.4g (%.3g)", v, s)
  lab <- c(unadj = "ATT_unadj", s = "ATT_s    ", c = "ATT_c    ", sc = "ATT_sc   ")
  cat("\n  ATT estimators (conv | bias-corrected):\n")
  for (i in seq_len(nrow(x$att))) {
    e <- x$att$estimator[i]
    cat(sprintf("    %s  conv = %s   bc = %s\n", lab[[e]],
                sefmt(x$att$conv[i], x$att$se[i]),
                sefmt(x$att$bc[i],   x$att$bc_se[i])))
  }
  cat("\n  Per-period decomposition (D, S, adj_s_D = D - S; conv):\n")
  pp <- x$per_period
  for (i in seq_len(nrow(pp)))
    cat(sprintf("    t=%g [%-10s]: D=%s  S=%s  D-S=%+.4g\n",
                pp$period[i], pp$role[i], sefmt(pp$D[i], pp$D_se[i]),
                sefmt(pp$S[i], pp$S_se[i]), pp$adj_s_D[i]))
  cat("\n  Cross-period composition C_{t0,t_rd} (conv):\n")
  ct <- x$cterm
  for (i in seq_len(nrow(ct)))
    cat(sprintf("    t0=%g:  C=%s\n", ct$comparison[i],
                sefmt(ct$C[i], ct$C_se[i])))
  cat(sprintf("    Cagg = %s\n", sefmt(x$cagg[["conv"]], x$cagg[["se"]])))
  if (x$se == "bootstrap")
    cat(sprintf("\n  SE: shared unit cluster bootstrap, B=%d (%d ok)\n",
                x$B, x$n_boot_ok))
  invisible(x)
}
