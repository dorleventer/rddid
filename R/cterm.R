# Cross-period composition term C of Proposition `cdisc` (Section 4.4 / 4.5).
#
# For an ORDERED pair of periods (t, s) -- t in the comparison role, s in the
# RD-period role -- the cross-period composition term is the Kitagawa
# (between-group) part of the change in the confounding discontinuity:
#
#   C_{t,s} = sum_a ( pi_{s,(+)}(v_t = a) - pi_{t,(+)}(v_s = a) ) * D_t(v_s = a),
#
# the linear-g0 estimator of eq. (Chat) in Appendix `app:est-adjust`, built from
# three single-period RD blocks (each a local-linear RD of eq. wls_rd):
#
#   * D_t(v_s = a)        : period-t jump among units grouped by their side a in
#                           the PARTNER period s (the within-cell comparison jump,
#                           alpha_{t,0}(a) = the observed period-t discontinuity).
#   * pi_{s,(+)}(v_t = a) : period-s above-cutoff share of units with side a in
#                           period t (the above-cutoff intercept of a local-linear
#                           RD of the indicator 1{v_t = a}).
#   * pi_{t,(+)}(v_s = a) : period-t above-cutoff share of units with side a in
#                           period s (the above-cutoff intercept of a local-linear
#                           RD of the indicator 1{v_s = a}).
#
# Their difference Delta_pi(a) = pi_{s,(+)}(v_t=a) - pi_{t,(+)}(v_s=a) is the
# cell-share change; weighting the comparison jump D_t(v_s=a) by it and summing
# over a gives C. In the composition-adjusted ATT of Theorem `thm:adjust`, C is
# used as C_{t0,t_rd} (t = t0 the comparison period, s = t_rd the RD period): the
# correction-form ATT is D_{t_rd} - g0({D_{t0}}) - sum_{t0} w_{t0} C_{t0,t_rd}.
#
# C is the PURE cross-period composition term: it is NOT the within-period term S
# of Proposition `sdisc` (which weights by the control-side outcome limit, not the
# RD-period confounding). The two are combined only at the ATT level.
#
# C_{t,s} is NOT symmetric in (t, s): the comparison jump D is taken at the FIRST
# period t, while the share weight pi_{s,(+)} is taken at the SECOND period s.
#
# Bandwidths default to per-block CCT, selected once on the original sample and
# held FIXED across bootstrap replications (recomputing CCT on the thin cells is
# unstable). Standard errors are a unit-level cluster bootstrap, the sanctioned
# equivalent of the closed-form a'Sigma a variance (App `app:est-adjust`).

# ---- internal: each unit's side (1 above, 0 below) of the cutoff in `period` ---
.rd_c_sidemap <- function(data, x, time, id, period, c) {
  keep <- data[[time]] == period
  stats::setNames(as.integer(data[[x]][keep] >= c),
                  as.character(data[[id]][keep]))
}

# ---- internal: the two working frames for an ordered pair (t, s) ---------------
# frame_t: period-t rows with the partner side `part` = each unit's side in s.
# frame_s: period-s rows with the partner side `part` = each unit's side in t.
# Units absent from the partner period get part = NA and drop from cell fits.
.rd_c_frames <- function(data, y, x, time, id, t, s, c) {
  side_s <- .rd_c_sidemap(data, x, time, id, s, c)
  side_t <- .rd_c_sidemap(data, x, time, id, t, c)
  kt <- data[[time]] == t
  ks <- data[[time]] == s
  idt <- as.character(data[[id]][kt])
  ids <- as.character(data[[id]][ks])
  list(
    frame_t = data.frame(
      y    = as.numeric(data[[y]][kt]),
      x    = as.numeric(data[[x]][kt]),
      part = unname(side_s[idt]),
      stringsAsFactors = FALSE),
    frame_s = data.frame(
      x    = as.numeric(data[[x]][ks]),
      part = unname(side_t[ids]),
      stringsAsFactors = FALSE)
  )
}

# ---- internal: point-estimate CCT bandwidths for the C blocks ------------------
# Returns list(D = list(per a), pis = list(per a), pit = list(per a)), each a
# named c(h=, b=). labels are the binary sides c(0L, 1L).
.rd_c_pointbw <- function(frames, labels, c, kernel) {
  ft <- frames$frame_t[stats::complete.cases(frames$frame_t$y, frames$frame_t$x,
                                             frames$frame_t$part), , drop = FALSE]
  fs <- frames$frame_s[stats::complete.cases(frames$frame_s$x,
                                             frames$frame_s$part), , drop = FALSE]
  cct <- function(yy, xx) rd_bw_cct(yy, xx, c = c, kernel = kernel)
  per_a <- function(fun) stats::setNames(lapply(labels, fun), as.character(labels))
  list(
    D   = per_a(function(a) { keep <- ft$part == a; cct(ft$y[keep], ft$x[keep]) }),
    pis = per_a(function(a) cct(as.numeric(fs$part == a), fs$x)),
    pit = per_a(function(a) cct(as.numeric(ft$part == a), ft$x))
  )
}

# ---- internal: one C fit (conv + bc) at FIXED bandwidths -----------------------
# Reusable point helper: rebuilds frames from `data`, so a shared bootstrap (or a
# combined D/S/C bootstrap) can call it with a resampled data frame. Returns the
# scalar pair c(C, C_bc) plus a per-cell detail matrix as attribute "cells".
.rd_c_point <- function(data, y, x, time, id, t, s, c, bws, kernel, labels) {
  frames <- .rd_c_frames(data, y, x, time, id, t, s, c)
  ft <- frames$frame_t[stats::complete.cases(frames$frame_t$y, frames$frame_t$x,
                                             frames$frame_t$part), , drop = FALSE]
  fs <- frames$frame_s[stats::complete.cases(frames$frame_s$x,
                                             frames$frame_s$part), , drop = FALSE]

  rp <- function(yy, xx, bw) tryCatch(
    rd_period(y = yy, x = xx, h = unname(bw["h"]), b = unname(bw["b"]),
              c = c, kernel = kernel),
    error = function(e) NULL)

  cells <- vapply(labels, function(a) {
    ai <- as.character(a)
    keep <- ft$part == a
    fD   <- rp(ft$y[keep], ft$x[keep], bws$D[[ai]])
    fpis <- rp(as.numeric(fs$part == a), fs$x, bws$pis[[ai]])
    fpit <- rp(as.numeric(ft$part == a), ft$x, bws$pit[[ai]])
    if (is.null(fD) || is.null(fpis) || is.null(fpit))
      return(rep(NA_real_, 8L))
    D       <- fD$D;            D_bc       <- fD$D_bc
    pis     <- fpis$sides[["+"]]$beta0;  pis_bc <- fpis$sides[["+"]]$beta0_bc
    pit     <- fpit$sides[["+"]]$beta0;  pit_bc <- fpit$sides[["+"]]$beta0_bc
    dpi     <- pis    - pit
    dpi_bc  <- pis_bc - pit_bc
    c(D = D, D_bc = D_bc, dpi = dpi, dpi_bc = dpi_bc,
      pis = pis, pit = pit, contrib = dpi * D, contrib_bc = dpi_bc * D_bc)
  }, numeric(8L))
  rownames(cells) <- c("D", "D_bc", "dpi", "dpi_bc", "pis", "pit",
                       "contrib", "contrib_bc")
  colnames(cells) <- as.character(labels)

  out <- c(C = sum(cells["contrib", ]), C_bc = sum(cells["contrib_bc", ]))
  attr(out, "cells") <- cells
  out
}

#' Cross-period composition term C (Proposition `cdisc`)
#'
#' Estimates the cross-period composition term \eqn{C_{t,s}} between two time
#' periods as a first-class primitive: the Kitagawa (between-group) part of the
#' change in the confounding discontinuity from period `t` to period `s`. With
#' binary sides \eqn{a\in\{0,1\}} (the unit's side of the cutoff) it is
#' \deqn{C_{t,s}=\sum_a\bigl(\pi_{s,(+)}(v_t{=}a)-\pi_{t,(+)}(v_s{=}a)\bigr)\,
#'   D_t(v_s{=}a),}
#' the linear-\eqn{g_0} estimator of Appendix `app:est-adjust`, where
#' \eqn{D_t(v_s{=}a)} is the period-`t` discontinuity among units with side `a`
#' in the partner period `s`, and \eqn{\pi_{\cdot,(+)}} are above-cutoff cell
#' shares (above-cutoff intercepts of a local-linear RD of the side indicator).
#'
#' In the composition-adjusted ATT (Theorem `thm:adjust`) the ordered pair is
#' `(t, s) = (t0, t_rd)`: `t` is the comparison period, `s` the RD period, and
#' the correction-form ATT subtracts \eqn{C_{t_0,t_{RD}}} from the unadjusted
#' aggregate. \eqn{C_{t,s}} is the pure CROSS-period term and is distinct from
#' the within-period term \eqn{S_t} of [rd_sadjust()]; the two are combined only
#' at the ATT level. \eqn{C_{t,s}} is not symmetric in `(t, s)`: the comparison
#' jump \eqn{D} is taken at `t`, the share weight \eqn{\pi_{s,(+)}} at `s`.
#'
#' Bandwidths default to a per-block CCT selection, made once on the original
#' sample and held fixed across the unit-level cluster bootstrap that supplies
#' the standard error (recomputing CCT on the thin cells is unstable).
#'
#' @param data long data frame, one row per unit-period.
#' @param y,x,time,id column names (strings) for outcome, running variable,
#'   period, and unit id.
#' @param t first (comparison-role) period value of `time`.
#' @param s second (RD-role) period value of `time`.
#' @param c cutoff (default 0).
#' @param h optional fixed numeric bandwidth used for every block; if `NULL`
#'   (default) a per-block CCT bandwidth is selected.
#' @param bwselect bandwidth selector label when `h` is `NULL` (default
#'   `"cct"`); retained for forward compatibility.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param se `"bootstrap"` (default) or `"none"`.
#' @param B bootstrap replications (default 499).
#'
#' @return An object of class `"rd_c"`: a list with the conventional and
#'   bias-corrected composition term (`C`, `C_bc`), their cluster-bootstrap
#'   standard errors (`C_se`, `C_bc_se`), the ordered periods (`t`, `s`), a
#'   per-cell `cells` data frame (`label`, `dpi`, `dpi_bc`, `D`, `D_bc`,
#'   `pi_s_plus`, `pi_t_plus`, `contrib`, `contrib_bc`), and meta fields.
#' @export
rd_c <- function(data, y, x, time, id, t, s, c = 0,
                 h = NULL, bwselect = "cct", kernel = "triangular",
                 se = c("bootstrap", "none"), B = 499L) {
  cl <- match.call()
  se <- match.arg(se)
  for (nm in c(y, x, time, id))
    if (!nm %in% names(data)) stop("column '", nm, "' not found in `data`.")
  if (t == s) stop("`t` and `s` must be different periods.")
  labels <- c(0L, 1L)

  # bandwidths fixed at the point estimate (per-block CCT or a single fixed `h`)
  frames0 <- .rd_c_frames(data, y, x, time, id, t, s, c)
  bws <- if (!is.null(h)) {
    hp <- c(h = h, b = h)
    per <- stats::setNames(rep(list(hp), length(labels)), as.character(labels))
    list(D = per, pis = per, pit = per)
  } else {
    .rd_c_pointbw(frames0, labels, c, kernel)
  }

  pt    <- .rd_c_point(data, y, x, time, id, t, s, c, bws, kernel, labels)
  cells <- attr(pt, "cells")

  # ---- bootstrap SE (resample units, rebuild frames, recompute at fixed bw) ----
  C_se <- NA_real_; C_bc_se <- NA_real_; n_boot_ok <- NA_integer_
  if (se == "bootstrap") {
    idx_by_id <- split(seq_len(nrow(data)), as.character(data[[id]]))
    U <- names(idx_by_id)
    bmat <- matrix(NA_real_, nrow = B, ncol = 2L,
                   dimnames = list(NULL, c("C", "C_bc")))
    for (rb in seq_len(B)) {
      drawn <- sample(U, length(U), replace = TRUE)
      rows_list <- idx_by_id[drawn]
      bd <- data[unlist(rows_list, use.names = FALSE), , drop = FALSE]
      bd[[id]] <- rep(seq_along(drawn), lengths(rows_list))   # fresh cluster ids
      bmat[rb, ] <- tryCatch(
        .rd_c_point(bd, y, x, time, id, t, s, c, bws, kernel, labels)[c("C", "C_bc")],
        error = function(e) c(NA_real_, NA_real_))
    }
    C_se      <- stats::sd(bmat[, "C"],    na.rm = TRUE)
    C_bc_se   <- stats::sd(bmat[, "C_bc"], na.rm = TRUE)
    n_boot_ok <- sum(stats::complete.cases(bmat))
  }

  cells_df <- data.frame(
    label      = labels,
    dpi        = cells["dpi", ],
    dpi_bc     = cells["dpi_bc", ],
    D          = cells["D", ],
    D_bc       = cells["D_bc", ],
    pi_s_plus  = cells["pis", ],
    pi_t_plus  = cells["pit", ],
    contrib    = cells["contrib", ],
    contrib_bc = cells["contrib_bc", ],
    stringsAsFactors = FALSE,
    row.names  = NULL
  )

  structure(list(
    C         = unname(pt["C"]),
    C_bc      = unname(pt["C_bc"]),
    C_se      = C_se,
    C_bc_se   = C_bc_se,
    t         = t,
    s         = s,
    cells     = cells_df,
    c         = c,
    h         = h,
    bwselect  = if (is.null(h)) bwselect else "fixed",
    kernel    = kernel,
    se        = se,
    B         = if (se == "bootstrap") B else 0L,
    n_boot_ok = n_boot_ok,
    call      = cl
  ), class = "rd_c")
}

#' @export
print.rd_c <- function(x, ...) {
  cat("Cross-period composition term C (Proposition cdisc)\n")
  cat(sprintf("  periods: t (comparison role) = %s   s (RD role) = %s   bw: %s\n",
              x$t, x$s, x$bwselect))
  sefmt <- function(v, sd) if (is.na(sd)) sprintf("%+.4g", v) else
    sprintf("%+.4g (%.3g)", v, sd)
  cat(sprintf("  C (conv) = %s\n", sefmt(x$C, x$C_se)))
  cat(sprintf("  C (bc)   = %s\n", sefmt(x$C_bc, x$C_bc_se)))
  if (x$se == "bootstrap")
    cat(sprintf("  SE: unit cluster bootstrap, B=%d (%d ok)\n", x$B, x$n_boot_ok))
  cat("\n  Per-cell contributions  Delta_pi(a) * D_t(v_s = a):\n")
  cc <- x$cells
  for (i in seq_len(nrow(cc)))
    cat(sprintf("    a=%d:  Delta_pi=%+.3f   D_t(a)=%+.4g   contrib=%+.4g\n",
                cc$label[i], cc$dpi[i], cc$D[i], cc$contrib[i]))
  invisible(x)
}
