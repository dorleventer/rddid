# bandwidth_joint.R -- the bandwidth rules that target the diff-in-disc estimate itself:
# .bw_constants() estimates the curvature at each period's own CCT pilot, and the variance
# constants and the bias design moments at one rule-of-thumb bandwidth shared by every period
# (.bw_rot());
# .bw_joint() minimises the aggregate AMSE over one common h (bwselect = "joint", the default).
# Labels such as eq:common_h_opt, lem:agg-var or "B.4 P3" in the internal notes below point into the
# paper's Appendix B; the code-to-paper map is dev/appB_map.md (developer material, not shipped).

# Original file header:
# Bandwidth selection for the aggregate diff-in-disc estimator.
#
# Two modes, matching Section 5.4 (sec:est-bw) and Appendix B.4 (app:est-bw) of
# Leventer and Nevo; the code <-> equation map is dev/appB_map.md (section 2.5).
#   * "cct"   - per-period MSE-optimal bandwidth (Imbens-Kalyanaraman / CCT) applied to
#               each D-hat_t on its own (B.4 P1). Delegated to rdrobust::rdbwselect.
#               Leaves the cross-period bias cancellation unexploited.
#   * "joint" - one common h* minimising the asymptotic MSE of the AGGREGATE estimator
#               (B.4 P3, eq:common_h_opt):
#                 h* = ( ((p+1)!)^2/(2(p+1)) * V^S / B^2 )^{1/(2p+3)} n^{-1/(2p+3)},
#               with B = sum_tau w~_tau B_tau (signed: comparison biases can offset the
#               RD-period bias) and V^S the scheme-specific aggregate variance constant
#               (lem:agg-var at a common h).
# Both use the conventional variance (CCT's mserd convention); the regularization
# of B.4 ("Regularization") is applied in "joint". (The period-specific "iter" rule was
# removed after tag v0.4.0.9000-iter; it was dropped from the paper, T1633.)


#' Rule-of-thumb preliminary bandwidth shared by every period (App. B.4, "Bandwidth estimation")
#'
#' rdrobust's rule of thumb for the preliminary bandwidth of `rdbwselect()`,
#' `C_K * min(sd, IQR / 1.349) * N^{-1/5}` (`C_K` = 2.576 triangular, 1.843 uniform, 2.34
#' Epanechnikov), computed once from the running variables of all periods pooled, with `N` the
#' average number of rows per period. Every period's variance constant is estimated at it, so
#' the cross-period covariances are estimated at a bandwidth ratio of one.
#' @keywords internal
#' @noRd
.bw_rot <- function(plist, kernel = "triangular") {
  x <- unlist(lapply(plist, function(d) d$x), use.names = FALSE)
  x <- x[is.finite(x)]
  s <- min(stats::sd(x), stats::IQR(x) / 1.349)
  if (!is.finite(s) || s <= 0) s <- stats::sd(x)
  ck <- switch(kernel, triangular = 2.576, uniform = 1.843, epanechnikov = 2.34)
  nbar <- mean(vapply(plist, nrow, numeric(1)))
  ck * s * nbar^(-1 / 5)
}

#' Per-period asymptotic constants for the aggregate AMSE (App. B.4)
#'
#' Used by the joint selector; the paper's "Bandwidth estimation" paragraph of App. B.4.
#' Two sets of fits, as in CCT's own selector:
#'   * curvature: every period at its OWN pilot pair (default: its CCT pair from
#'     [rd_bw_cct()]); gives the order-q coefficient on `(x - c)^{p+1}` at `b_t` on each side
#'     (`curv`) and the pilot `b/h` ratios;
#'   * design side: every period at ONE rule-of-thumb bandwidth `h_rot` ([.bw_rot()]); gives
#'     the bias design moments `bmom` (rdrobust's `BConst`), so that
#'     `B-hat_t = (p+1)! [bmom_+(h_rot) curv_+(b_t) - bmom_-(h_rot) curv_-(b_t)]`, the
#'     covariance matrix of the `B-hat_t` for the regularization, with the cross-period
#'     covariances the scheme allows (none under `"cs"`, same side under `"pc"`, same minus
#'     opposite side under `"pv"`, as for the jumps), and
#'     the variance constant `V-hat_t = n_t h_rot V-hat(D_t(h_rot))` (`v_const`) and, under
#'     `"pc"`, the same-side cross-period term `h_rot Cov-hat(D_t(h_rot), D_s(h_rot))`, estimated
#'     directly at a bandwidth ratio of one (the ratio the common bandwidth needs). Under
#'     `"cs"` the covariance is zero and under `"pv"` it is O(1/n) = o(1/(nh)) (`lem:cov-pv`),
#'     so the limit has no cross term and none is used.
#' Nothing here depends on which period is the RD period, so the selector built on it is
#' invariant to that label; this matters when the estimator aggregates several RD periods.
#'
#' @param plist named per-period data list; `coef` named period coefficients.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`.
#' @param pilot_bws optional named list of per-period `c(h, b)` for the bias side; defaults
#'   to CCT.
#' @param h_rot optional bandwidth for the variance side; defaults to [.bw_rot()].
#' @return list with `keys`, `coef` (in `keys` order), `bias_const`, `var_const`, `n_obs`,
#'   `var_bias` (K x K covariance of the `B-hat_t`), `var_cross` (K x K same-side cross term,
#'   zero off PC), `h_rot`, `h_pilot`, `ratio`, `pilot_fits`, `var_fits`, `pilot_bws`.
#' @keywords internal
#' @noRd
.bw_constants <- function(plist, coef, scheme = "cs", pilot_bws = NULL, h_rot = NULL,
                          c = 0, p = 1L, q = 2L, kernel = "triangular") {
  cutoff <- c
  keys <- names(coef)
  K    <- length(keys)
  fp1  <- factorial(p + 1)
  if (is.null(pilot_bws))
    pilot_bws <- .bw_cct(plist, c = cutoff, p = p, q = q, kernel = kernel)
  missing_k <- setdiff(keys, names(pilot_bws))
  if (length(missing_k) > 0L)
    .abort("`pilot_bws` is missing entries for period(s): ",
           paste(missing_k, collapse = ", "), ".")
  if (is.null(h_rot)) h_rot <- .bw_rot(plist[keys], kernel = kernel)

  fit_at <- function(k, h, b)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = h, b = b, id = plist[[k]]$id,
              c = cutoff, p = p, q = q, kernel = kernel)

  # The estimated bias of D_t on a side is h^{p+1} * bmom * curv: a design moment at h and the
  # curvature (the order-q coefficient on (x - c)^{p+1}) from the pilot fit at b. As in
  # rdrobust's selector, the curvature needs its own pilot bandwidth, but the design moment,
  # like the variance, uses only the running variable, so both are taken at the shared h_rot.
  # pilot fits: every period at its own CCT pair (the curvature, and the b/h ratios)
  pilot_fits <- stats::setNames(lapply(keys, function(k)
    fit_at(k, pilot_bws[[k]][["h"]], pilot_bws[[k]][["b"]])), keys)
  # variance-side fits: every period at the common rule-of-thumb bandwidth (the variance
  # constants and the design moments)
  var_fits <- stats::setNames(lapply(keys, function(k) fit_at(k, h_rot, h_rot)), keys)
  period_coef <- unname(vapply(keys, function(k) coef[[k]], numeric(1)))
  ratio       <- vapply(keys, function(k)
    unname(pilot_bws[[k]][["b"]] / pilot_bws[[k]][["h"]]), numeric(1))
  h_pilot     <- vapply(keys, function(k) unname(pilot_bws[[k]][["h"]]), numeric(1))

  # B-hat_t = (p+1)! [bmom_+(h_rot) curv_+(b_t) - bmom_-(h_rot) curv_-(b_t)]. psi is the
  # per-unit influence of one side's term, (p+1)! bmom(h_rot) g_curv(b_t), unsigned: the
  # minus sign of the side below is applied through "same - opp" in the covariances.
  psi <- function(k, side)
    fp1 * var_fits[[k]]$sides[[side]]$bmom * pilot_fits[[k]]$sides[[side]]$g_curv
  bias_const <- vapply(keys, function(k)
    fp1 * (var_fits[[k]]$sides[["+"]]$bmom * pilot_fits[[k]]$sides[["+"]]$curv -
           var_fits[[k]]$sides[["-"]]$bmom * pilot_fits[[k]]$sides[["-"]]$curv), numeric(1))

  # Cov(B-hat_t, B-hat_s) from the psi vectors, matched by id across periods like the jumps'
  # g: this is what the regularization term adds to the squared bias of the aggregate
  ids <- function(k, side) pilot_fits[[k]]$sides[[side]]$id
  ms  <- function(a, sa, b, sb) .match_sum(ids(a, sa), psi(a, sa), ids(b, sb), psi(b, sb))
  var_bias <- matrix(0, K, K, dimnames = list(keys, keys))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    a <- keys[i]; b <- keys[j]
    var_bias[i, j] <- if (i == j) {
      sum(psi(a, "+")^2) + sum(psi(a, "-")^2)
    } else if (scheme == "cs") {
      0
    } else {
      same <- ms(a, "+", b, "+") + ms(a, "-", b, "-")
      opp  <- if (scheme == "pv") ms(a, "+", b, "-") + ms(a, "-", b, "+") else 0
      same - opp
    }
  }

  var_const <- vapply(keys, function(k) var_fits[[k]]$v_const, numeric(1))
  n_obs     <- vapply(keys, function(k) var_fits[[k]]$n,       numeric(1))
  var_cross <- matrix(0, K, K, dimnames = list(keys, keys))
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K)) for (j in seq_len(K)) if (i != j)
      var_cross[i, j] <- h_rot * .cross_cov(var_fits[[keys[i]]], var_fits[[keys[j]]], bc = FALSE)$pc
    if (all(var_cross == 0))
      .warn("scheme = \"pc\" but no unit is active in two periods' windows: ",
            "the cross-period term is identically zero (same as scheme = \"cs\").")
  }
  list(keys = keys, coef = period_coef, bias_const = bias_const, var_const = var_const,
       n_obs = n_obs, var_bias = var_bias, var_cross = var_cross, h_rot = h_rot,
       h_pilot = h_pilot, ratio = ratio, pilot_fits = pilot_fits, var_fits = var_fits,
       pilot_bws = pilot_bws)
}

#' Joint AMSE-optimal common bandwidth for the aggregate estimator
#'
#' Feasible plug-in for eq:common_h_opt (App. B.4, "Bandwidth estimation"). With the
#' constants of [.bw_constants()] (curvature at each period's own CCT pilot; variance and bias
#' design moments at the shared rule-of-thumb `h_rot`), the aggregate objective at a common `h` is
#'   AMSE^S(h) = (h^{p+1}/(p+1)!)^2 (B^2 + reg) + Veff / h,
#'   B    = sum_tau w~_tau B-hat_tau                       (signed; can cancel),
#'   Veff = sum_tau w~_tau^2 V-hat_tau / n_tau
#'          + 1{S = PC} sum_{t != s} w~_t w~_s h_rot Cov-hat(D_t(h_rot), D_s(h_rot)),
#'   reg  = lambda Var-hat(sum_tau w~_tau B-hat_tau)   (with the cross-period covariances),
#' (`Veff` is `V^S(t_RD) / n` of `lem:agg-var` at a common h, so the `n^{-1/(2p+3)}`
#' factor is absorbed), whose minimizer is
#'   h* = ( ((p+1)!)^2/(2(p+1)) * Veff / (B^2 + reg) )^{1/(2p+3)}.
#' The pilot
#' bandwidth keeps each period's CCT ratio, `b_t = h* b_t^CCT / h_t^CCT` (B.4 after
#' eq:common_h_opt; consistent with Assumption R(f)).
#'
#'
#' @param plist named per-period data list; `coef` named period coefficients
#'   (RD period = +1, comparisons = -w_t); `t_rd` the RD-period label (kept for
#'   call compatibility; the selector does not use it).
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`; sets which variance scales the
#'   bandwidth. Under `"pc"` the same-side cross-period covariances enter (they are
#'   O(1/(nh)), `lem:cov-pc`); under `"pv"` they are O(1/n) = o(1/(nh)) (`lem:cov-pv`)
#'   and drop, so the covariance-free CS form is used (`lem:agg-var`).
#' @param pilot_bws optional named list of per-period `c(h, b)` pilots; defaults to
#'   each period's CCT/IK pair.
#' @param regularize add the App. B.4 regularization term to the squared bias so a
#'   small, noisy estimated bias constant cannot inflate `h*` (mirrors the
#'   `regularize = TRUE` default of [rdrobust::rdbwselect]). Default `TRUE`.
#' @param reg_const multiple of the bias-estimate variance used as the
#'   regularization term (the paper's `lambda`; default 3, the CCT convention).
#' @param constants optional precomputed [.bw_constants()] output (reuses the
#'   constants without refitting).
#' @return list with the common `h`, the per-period pilot bandwidths `b` (named by
#'   period), the bias constant `B`, the variance scale `Veff`, the regularization
#'   term `reg`, the rule-of-thumb `h_rot`, the `pilot_bws` used, and the `constants`.
#' @keywords internal
#' @noRd
.bw_joint <- function(plist, coef, t_rd = NULL, scheme = "cs", pilot_bws = NULL,
                      c = 0, p = 1L, q = 2L, kernel = "triangular",
                      regularize = TRUE, reg_const = 3, constants = NULL) {
  cutoff <- c
  if (is.null(constants))
    constants <- .bw_constants(plist, coef, scheme = scheme, pilot_bws = pilot_bws,
                               c = cutoff, p = p, q = q, kernel = kernel)
  const <- constants
  fp1   <- factorial(p + 1)
  K     <- length(const$keys)

  # aggregate bias constant: sum_t coef_t B-hat_t
  B_agg <- sum(const$coef * const$bias_const)
  # variance scale: own-period terms sum_t coef_t^2 V-hat_t / n_t, plus, under "pc", the
  # same-side cross-period terms (estimated at h_rot in every period, so at rho = 1)
  cf   <- const$coef
  Veff <- sum(cf^2 * const$var_const / const$n_obs) + drop(crossprod(cf, const$var_cross %*% cf))
  # regularization: reg_const x Var(sum_t coef_t B-hat_t), with the cross-period covariances of
  # the B-hat_t. At a common h it carries the same h^{2(p+1)} factor as B_agg^2, so it simply
  # joins B_agg^2 in the denominator of h*.
  reg   <- if (regularize) reg_const * drop(crossprod(cf, const$var_bias %*% cf)) else 0
  denom <- B_agg^2 + reg

  if (!is.finite(denom) || denom <= 0)
    .abort("aggregate bias constant B is ~0 and regularization is off: the ",
           "AMSE-optimal bandwidth diverges. Use bwselect = \"cct\" or pass h.")
  if (!is.finite(Veff) || Veff <= 0)
    .abort("aggregate variance scale is not positive; check the sampling scheme and ",
           "the per-period fits.")
  # the closed-form minimiser of the AMSE in h: for p = 1 the leading constant is 1 and the
  # exponent 1/5
  h_star <- (fp1^2 / (2 * (p + 1)) * Veff / denom)^(1 / (2 * p + 3))
  radius <- max(vapply(plist, function(d) max(abs(d$x - cutoff), na.rm = TRUE), numeric(1)))
  if (is.finite(radius) && h_star > radius)
    .warn("joint AMSE-optimal bandwidth h* = ", signif(h_star, 4),
          " exceeds the running-variable radius ", signif(radius, 4),
          " (the period biases nearly cancel); consider bwselect = \"cct\" or a fixed h.")
  b <- stats::setNames(h_star * const$ratio, const$keys)
  list(h = h_star, b = b, B = B_agg, Veff = Veff, reg = reg, h_rot = const$h_rot,
       pilot_bws = const$pilot_bws, constants = const)
}
