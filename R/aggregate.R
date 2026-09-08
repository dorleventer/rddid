# Cross-period covariance and aggregation of per-period RD fits.
# All variances flow from the per-unit influence-times-residual vectors `g`
# returned by rd_period(): the discontinuity variance is sum(g^2) per side, and
# a cross-period covariance is the sum of g_t * g_s over units shared between
# the two periods' windows. Same-side matches build C^same, opposite-side matches
# build C^opp (eq:cross-decomp), and the sampling scheme selects which enter the
# aggregate variance (eq:var-cs, eq:var-pc, eq:var-pv). Map: dev/appB_map.md §2.3.

#' Sum of g_t * g_s over units present in both windows, matched on id
#' @keywords internal
#' @noRd
.match_sum <- function(idA, gA, idB, gB) {
  m <- match(idA, idB)
  ok <- !is.na(m)
  if (!any(ok)) return(0)
  sum(gA[ok] * gB[m[ok]])
}

#' Same-side (PC) and opposite-side (PV) cross-period covariance of two D-hats
#'
#' @param ft,fs `rd_period` fits for two distinct periods.
#' @param bc use bias-corrected `g` vectors.
#' @return list with `pc` = Cov(b+_t,b+_s)+Cov(b-_t,b-_s) (same-side, the
#'   paper's C^same), `pv` = Cov(b+_t,b-_s)+Cov(b-_t,b+_s) (opposite-side,
#'   C^opp), and the two same-side parts `pc_p` = Cov(b+_t,b+_s) and
#'   `pc_m` = Cov(b-_t,b-_s), which the period-specific bandwidth rule scales
#'   separately by side (lem:cov-pc).
#' @keywords internal
#' @noRd
.cross_cov <- function(ft, fs, bc = FALSE) {
  gt_p <- if (bc) ft$sides$`+`$g_bc else ft$sides$`+`$g
  gt_m <- if (bc) ft$sides$`-`$g_bc else ft$sides$`-`$g
  gs_p <- if (bc) fs$sides$`+`$g_bc else fs$sides$`+`$g
  gs_m <- if (bc) fs$sides$`-`$g_bc else fs$sides$`-`$g
  it_p <- ft$sides$`+`$id; it_m <- ft$sides$`-`$id
  is_p <- fs$sides$`+`$id; is_m <- fs$sides$`-`$id
  pc_p <- .match_sum(it_p, gt_p, is_p, gs_p)      # Cov(b+_t, b+_s)
  pc_m <- .match_sum(it_m, gt_m, is_m, gs_m)      # Cov(b-_t, b-_s)
  list(
    pc = pc_p + pc_m,                                                            # C^same
    pv = .match_sum(it_p, gt_p, is_m, gs_m) + .match_sum(it_m, gt_m, is_p, gs_p), # C^opp
    pc_p = pc_p, pc_m = pc_m
  )
}

#' Scheme-combined cross-period covariance Cov(D_t, D_s) of two D-hats
#'
#' Collapses `.cross_cov()` to a single number under a sampling scheme:
#' `cs` → 0 (independent), `pc` → same-side only, `pv` → same-side minus
#' opposite-side (`C^same - C^opp`, eq:cross-decomp). Used wherever a pair of
#' `rd_period` fits in different periods needs its scalar covariance.
#'
#' @param ft,fs `rd_period` fits for two distinct periods.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`.
#' @param bc use bias-corrected `g` vectors.
#' @keywords internal
#' @noRd
.cov_scheme <- function(ft, fs, scheme, bc = FALSE) {
  if (scheme == "cs") return(0)
  cc <- .cross_cov(ft, fs, bc = bc)
  if (scheme == "pc") cc$pc else cc$pc - cc$pv
}

#' Aggregate per-period fits into the RD-DID estimator and its variances
#'
#' Forms ATT-hat = sum_tau coef_tau * D-hat_tau and the variance under the three
#' sampling schemes (eq:var-cs, eq:var-pc, eq:var-pv). With
#' Cov(D_t,D_s) = C^same - C^opp (eq:cross-decomp), the general quadratic form
#' specialises to: CS drops all cross terms, PC keeps same-side (C^same) only, PV
#' keeps both (C^same - C^opp).
#'
#' @param fits named list of `rd_period` objects (names = period labels).
#' @param coef named numeric vector of period coefficients (same names as fits).
#' @param bc use bias-corrected estimates / variances.
#' @return named numeric vector: `est`, `V_cs`, `V_pc`, `V_pv`.
#' @keywords internal
#' @noRd
.aggregate_fits <- function(fits, coef, bc = FALSE) {
  keys <- names(coef)
  Dget <- function(k) if (bc) fits[[k]]$D_bc else fits[[k]]$D
  Vget <- function(k) if (bc) fits[[k]]$V_D_bc else fits[[k]]$V_D

  est   <- sum(coef * vapply(keys, Dget, numeric(1)))
  Vdiag <- sum(coef^2 * vapply(keys, Vget, numeric(1)))

  add_pc <- 0
  add_pv <- 0
  if (length(keys) >= 2L) {
    for (i in seq_along(keys)) {
      for (j in seq_along(keys)) {
        if (j <= i) next
        cc <- .cross_cov(fits[[keys[i]]], fits[[keys[j]]], bc = bc)
        wij <- 2 * coef[[i]] * coef[[j]]
        add_pc <- add_pc + wij * cc$pc
        add_pv <- add_pv + wij * (cc$pc - cc$pv)
      }
    }
  }
  c(est = est, V_cs = Vdiag, V_pc = Vdiag + add_pc, V_pv = Vdiag + add_pv)
}
