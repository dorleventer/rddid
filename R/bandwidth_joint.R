# bandwidth_joint.R -- the bandwidth rules that target the RD-DID estimate itself:
# .bw_constants() estimates each period's bias and variance constants at its own CCT pilot;
# .bw_joint() minimises the aggregate AMSE over one common h (bwselect = "joint", the default).
# The kernel constants for the PC cross-period term come from kernel_constants.R.
# Labels such as eq:common_h_opt, lem:agg-var or "B.4 P3" in the internal notes below point into the
# paper's Appendix B; the code-to-paper map is dev/appB_map.md (developer material, not shipped).

# Original file header:
# Bandwidth selection for the aggregate RD-DID estimator.
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


#' Per-period asymptotic constants for the aggregate AMSE (App. B.4)
#'
#' Used by the joint selector. Fits every period at its OWN pilot pair
#' (default: its CCT/IK pair from [rd_bw_cct()]) and returns the plug-ins of
#' eq:amse-att: the curvature constant `B-hat_t` (`rd_period`'s `b_const`), the
#' variance constant `V-hat_t` (`v_const`), the period sample size `n_t`, the
#' variance of `B-hat_t` for the regularization term, the pilot `h`'s and `b/h`
#' ratios, and, under `"pc"`, the per-side h-free scales `kappa` of the same-side
#' cross-period covariance (`lem:cov-pc`). Nothing here depends on which period is
#' the RD period, so the selector built on it is invariant to that label; this
#' matters when the estimator aggregates several RD periods.
#'
#' PC cross-period term (`lem:cov-pc`). Per side, the same-side covariance of two
#' intercepts satisfies `n h_s Cov -> sigma_{t,s,(side)} c_side(p, h_t/h_s) / f(c)`,
#' with `c_side(p, rho)` the kernel constant of App. B.3 (the paper's
#' omega_(side),p(rho); `.kc_c` here). The pilot fits give the plug-in
#' `P_side = Cov-hat(h0_t, h0_s)`, from which the h-free scale
#'   kappa_side = P_side * h0_s / c_side(p, h0_t/h0_s)   (~ sigma / (f(c) n))
#' is read off; at a common bandwidth h the term is `kappa_side * c_side(p, 1) / h`.
#' Under `"cs"` the covariance is zero and under `"pv"` it is O(1/n) = o(1/(nh))
#' (`lem:cov-pv`), so the term is omitted for both.
#'
#' @param plist named per-period data list; `coef` named period coefficients.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`.
#' @param pilot_bws optional named list of per-period `c(h, b)`; defaults to CCT.
#' @return list with `keys`, `coef` (in `keys` order), `bias_const`, `var_const`, `n_obs`,
#'   `var_bias_const`, `h_pilot`, `ratio`, `kappa_plus`, `kappa_minus`, `pilot_fits`,
#'   `pilot_bws`.
#' @keywords internal
#' @noRd
.bw_constants <- function(plist, coef, scheme = "cs", pilot_bws = NULL,
                          c = 0, p = 1L, q = 2L, kernel = "triangular") {
  cutoff <- c
  keys <- names(coef)
  fp1  <- factorial(p + 1)
  if (is.null(pilot_bws))
    pilot_bws <- .bw_cct(plist, c = cutoff, p = p, q = q, kernel = kernel)
  missing_k <- setdiff(keys, names(pilot_bws))
  if (length(missing_k) > 0L)
    .abort("`pilot_bws` is missing entries for period(s): ",
           paste(missing_k, collapse = ", "), ".")

  # one pilot fit per period at its own CCT pair; rd_period() returns the plug-in constants
  # B-hat_t (bias, via the conventional / bias-corrected gap) and V-hat_t (variance)
  pilot_fits <- lapply(keys, function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = pilot_bws[[k]][["h"]],
              b = pilot_bws[[k]][["b"]], id = plist[[k]]$id,
              c = cutoff, p = p, q = q, kernel = kernel))
  names(pilot_fits) <- keys
  period_coef <- unname(vapply(keys, function(k) coef[[k]], numeric(1)))
  bias_const  <- vapply(keys, function(k) pilot_fits[[k]]$b_const, numeric(1))
  var_const   <- vapply(keys, function(k) pilot_fits[[k]]$v_const, numeric(1))
  n_obs       <- vapply(keys, function(k) pilot_fits[[k]]$n,       numeric(1))
  ratio       <- vapply(keys, function(k)
    unname(pilot_bws[[k]][["b"]] / pilot_bws[[k]][["h"]]), numeric(1))
  h_pilot     <- vapply(keys, function(k) unname(pilot_bws[[k]][["h"]]), numeric(1))

  # Var(B-hat_t) = ((p+1)! / h0^{p+1})^2 Var(D - D_bc), from the g_diff influence vectors;
  # this is what the regularization term adds to the squared bias
  var_bias_const <- vapply(keys, function(k) {
    s <- pilot_fits[[k]]$sides
    (fp1 / h_pilot[[k]]^(p + 1))^2 * (sum(s[["+"]]$g_diff^2) + sum(s[["-"]]$g_diff^2))
  }, numeric(1))

  # Same-side cross-period covariance under scheme "pc", precomputed per side so the objective
  # can rescale it to any pair of bandwidths: kappa_side[i, j] = pilot covariance x h0_j /
  # c_side(p, h0_i / h0_j), where c_side is the kernel constant of kernel_constants.R.
  K <- length(keys)
  kappa_plus  <- matrix(0, K, K)
  kappa_minus <- matrix(0, K, K)
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        cc   <- .cross_cov(pilot_fits[[keys[i]]], pilot_fits[[keys[j]]], bc = FALSE)
        rho0 <- h_pilot[i] / h_pilot[j]
        kappa_plus[i, j]  <- cc$pc_p * h_pilot[j] / .kc_c(p, "+", rho0, kernel)
        kappa_minus[i, j] <- cc$pc_m * h_pilot[j] / .kc_c(p, "-", rho0, kernel)
      }
    }
    if (all(kappa_plus == 0) && all(kappa_minus == 0))
      .warn("scheme = \"pc\" but no unit is active in two periods' windows: ",
            "the cross-period term is identically zero (same as scheme = \"cs\").")
  }
  list(keys = keys, coef = period_coef, bias_const = bias_const, var_const = var_const,
       n_obs = n_obs, var_bias_const = var_bias_const, h_pilot = h_pilot, ratio = ratio,
       kappa_plus = kappa_plus, kappa_minus = kappa_minus,
       pilot_fits = pilot_fits, pilot_bws = pilot_bws)
}

#' Joint AMSE-optimal common bandwidth for the aggregate estimator
#'
#' Feasible plug-in for eq:common_h_opt (App. B.4 P3). With the per-period constants
#' of [.bw_constants()] (each period at its own CCT pilot), the aggregate objective
#' at a common `h` is
#'   AMSE^S(h) = (h^{p+1}/(p+1)!)^2 (B^2 + reg) + Veff / h,
#'   B    = sum_tau w~_tau B-hat_tau                       (signed; can cancel),
#'   Veff = sum_tau w~_tau^2 V-hat_tau / n_tau
#'          + 1{S = PC} 2 sum_{t<s} w~_t w~_s [kappa_+ c_+(p,1) + kappa_- c_-(p,1)],
#'   reg  = lambda sum_tau w~_tau^2 Var(B-hat_tau),
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
#'   term `reg`, the `pilot_bws` used, and the `constants`.
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
  # variance scale: own-period terms sum_t coef_t^2 V-hat_t / n_t ...
  Veff <- sum(const$coef^2 * const$var_const / const$n_obs)
  # ... plus, under "pc", the same-side cross-period terms at a common bandwidth (rho = 1)
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        Veff <- Veff + 2 * const$coef[i] * const$coef[j] *
          (const$kappa_plus[i, j]  * .kc_c(p, "+", 1, kernel) +
           const$kappa_minus[i, j] * .kc_c(p, "-", 1, kernel))
      }
    }
  }
  # regularization: reg_const x sum_t coef_t^2 Var(B-hat_t). At a common h it carries the same
  # h^{2(p+1)} factor as B_agg^2, so it simply joins B_agg^2 in the denominator of h*.
  reg   <- if (regularize) reg_const * sum(const$coef^2 * const$var_bias_const) else 0
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
  list(h = h_star, b = b, B = B_agg, Veff = Veff, reg = reg,
       pilot_bws = const$pilot_bws, constants = const)
}
