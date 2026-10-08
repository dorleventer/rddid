# bandwidth_joint.R -- the bandwidth rules that target the RD-DID estimate itself:
# .bw_constants() estimates each period's bias and variance constants at its own CCT pilot;
# .bw_joint() minimises the aggregate AMSE over one common h (bwselect = "joint", the default);
# .bw_joint_iter() minimises it over one h per period by coordinate descent (bwselect = "iter").
# The kernel constants for the PC cross-period term come from kernel_constants.R.
# Labels such as eq:common_h_opt, lem:agg-var or "B.4 P3" in the internal notes below point into the
# paper's Appendix B; the code-to-paper map is dev/appB_map.md (developer material, not shipped).

# Original file header:
# Bandwidth selection for the aggregate RD-DID estimator.
#
# Three modes, matching Section 5.3 (sec:est-bw) and Appendix B.4 (app:est-bw) of
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
#   * "iter"  - period-specific bandwidths by coordinate descent on the aggregate AMSE
#               (B.4 P4, eq:update, Algorithm alg:coorddesc).
# All three use the conventional variance (CCT's mserd convention); the regularization
# of B.4 ("Regularization") is applied in "joint" and "iter".


#' Per-period asymptotic constants for the aggregate AMSE (App. B.4)
#'
#' Shared by the two joint selectors. Fits every period at its OWN pilot pair
#' (default: its CCT/IK pair from [rd_bw_cct()]) and returns the plug-ins of
#' eq:amse-att: the curvature constant `B-hat_t` (`rd_period`'s `b_const`), the
#' variance constant `V-hat_t` (`v_const`), the period sample size `n_t`, the
#' variance of `B-hat_t` for the regularization term, the pilot `h`'s and `b/h`
#' ratios, and, under `"pc"`, the per-side h-free scales `kappa` of the same-side
#' cross-period covariance (`lem:cov-pc`; see the `.bw_joint_iter` header for the
#' construction). Nothing here depends on which period is the RD period, so the
#' selectors built on it are invariant to that label; this matters when the
#' estimator aggregates several RD periods.
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
    stop("`pilot_bws` is missing entries for period(s): ",
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
      warning("scheme = \"pc\" but no unit is active in two periods' windows: ",
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
#' This is exactly the scalar restriction of the `.bw_joint_iter` objective, so the
#' two selectors optimize one function (pinned in test-appB-bandwidth.R). The pilot
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
#' @param constants optional precomputed [.bw_constants()] output (used by
#'   `.bw_joint_iter` to seed without refitting).
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
    stop("aggregate bias constant B is ~0 and regularization is off: the ",
         "AMSE-optimal bandwidth diverges. Use bwselect = \"cct\" or pass h.")
  if (!is.finite(Veff) || Veff <= 0)
    stop("aggregate variance scale is not positive; check the sampling scheme and ",
         "the per-period fits.")
  # the closed-form minimiser of the AMSE in h: for p = 1 the leading constant is 1 and the
  # exponent 1/5
  h_star <- (fp1^2 / (2 * (p + 1)) * Veff / denom)^(1 / (2 * p + 3))
  radius <- max(vapply(plist, function(d) max(abs(d$x - cutoff), na.rm = TRUE), numeric(1)))
  if (is.finite(radius) && h_star > radius)
    warning("joint AMSE-optimal bandwidth h* = ", signif(h_star, 4),
            " exceeds the running-variable radius ", signif(radius, 4),
            " (the period biases nearly cancel); consider bwselect = \"cct\" or a fixed h.")
  b <- stats::setNames(h_star * const$ratio, const$keys)
  list(h = h_star, b = b, B = B_agg, Veff = Veff, reg = reg,
       pilot_bws = const$pilot_bws, constants = const)
}

#' Period-specific joint AMSE bandwidths by coordinate descent
#'
#' Minimises the aggregate AMSE of ATT-hat (eq:amse-att) over a SEPARATE bandwidth per
#' period, cycling one period at a time (eq:update, Algorithm alg:coorddesc): holding
#' the others fixed, each update is a scalar minimisation of
#'   AMSE^S(h_tau | rest) = ( sum_t w~_t h_t^{p+1} B_t / (p+1)! )^2
#'                        + sum_t w~_t^2 V_t / (n_t h_t)
#'                        + 1{S = PC} sum_t sum_{s != t} w~_t w~_s C_{t,s}(h_t/h_s) / (n h_s)
#'                        + regularization,
#' over `[lo, hmax]`. Each update weakly lowers the objective, so the returned
#' bandwidths weakly improve on their starting point (B.4 P4). The starting point is
#' controlled by `start` (default: the common joint-optimal h* of [.bw_joint()],
#' computed from the same constants).
#'
#' PC cross-period term (`lem:cov-pc`, B.4 P4). Per side, the same-side covariance of
#' two intercepts satisfies `n h_s Cov -> sigma_{t,s,(side)} c_side(p, h_t/h_s) / f(c)`,
#' with `c_side(p, rho)` the kernel constant of App. B.3 (the paper's
#' omega_(side),p(rho); `.kc_c` here). The pilot fits give
#' the plug-in `P_side = Cov-hat(h0_t, h0_s)`, from which the h-free scale
#'   kappa_side = P_side * h0_s / c_side(p, h0_t/h0_s)   (~ sigma / (f(c) n))
#' is read off; at candidate bandwidths the term is `kappa_side * c_side(p, h_t/h_s) / h_s`.
#' This is symmetric in (t, s) because `c(1/rho) = rho c(rho)`. Under `"cs"` the
#' covariance is zero and under `"pv"` it is O(1/n) = o(1/(nh)) (`lem:cov-pv`), so the
#' term is omitted for both.
#'
#' @param plist named per-period data; `coef` named period coefficients;
#'   `t_rd` RD-period label (kept for call compatibility; not used by the selector).
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`; determines whether the
#'   cross-period covariance term enters the AMSE (only for `"pc"`).
#' @param pilot_bws optional named list of per-period `c(h, b)` used to
#'   estimate the per-period constants and the regularization penalty; defaults
#'   to CCT. The role of `pilot_bws` is to supply the asymptotic constants
#'   (curvature, variance, regularization, PC covariance) - it does NOT control
#'   the descent starting point; that is controlled by `start`.
#' @param start seed for the coordinate descent. Three modes:
#'   \describe{
#'     \item{`"hstar"` (default)}{Start all periods at the common joint-optimal
#'       h* (the same scalar h* for every period).}
#'     \item{`"cct"`}{Start each period at its own CCT/IK pilot h (i.e.
#'       `pilot_bws[[k]]["h"]` for period `k`).}
#'     \item{numeric / named list}{Supply a manual per-period seed. Must cover
#'       every period in `names(coef)`. A numeric vector is matched positionally
#'       to `names(coef)`; a named vector/list is matched by name.}
#'   }
#' @param hmax cap on each bandwidth (defaults to the running-variable radius).
#' @param maxit,tol coordinate-descent controls.
#' @return list with per-period `bws` (named list of `c(h, b)`), the constants
#'   `b_const`/`v_const`, the iteration count `niter`, the objective value at the
#'   returned bandwidths `objective`, and the objective function `amse_fun` (a
#'   function of the vector of bandwidths in `names(coef)` order) for diagnostics
#'   and tests.
#' @keywords internal
#' @noRd
.bw_joint_iter <- function(plist, coef, t_rd = NULL, scheme = "cs", pilot_bws = NULL,
                           start = "hstar",
                           c = 0, p = 1L, q = 2L, kernel = "triangular",
                           regularize = TRUE, reg_const = 3,
                           hmax = NULL, maxit = 50L, tol = 1e-4) {
  cutoff <- c
  const <- .bw_constants(plist, coef, scheme = scheme, pilot_bws = pilot_bws,
                         c = cutoff, p = p, q = q, kernel = kernel)
  keys <- const$keys
  if (is.null(hmax))
    hmax <- max(vapply(plist, function(d) max(abs(d$x - cutoff), na.rm = TRUE), numeric(1)))
  amse <- .bw_iter_objective(const, scheme, p, kernel, regularize, reg_const)

  # search box for every period: from 3% of the running-variable radius to the radius (a
  # package choice; the paper's update step is unconstrained, see the boundary warning below)
  lo <- 0.03 * hmax
  h0 <- .bw_iter_start(start, const, plist, coef, t_rd, scheme, cutoff, p, q, kernel,
                       regularize, reg_const)
  h  <- pmin(pmax(h0, lo), hmax)

  # coordinate descent: sweep the periods, minimising the objective in each bandwidth with the
  # others fixed; stop when no bandwidth moved by more than tol x radius, or after maxit sweeps
  it <- 0L
  repeat {
    it <- it + 1L
    h_old <- h
    for (j in seq_along(keys)) {
      obj <- function(hj) { hv <- h; hv[j] <- hj; amse(hv) }
      h[j] <- stats::optimize(obj, interval = c(lo, hmax))$minimum
    }
    if (max(abs(h - h_old)) < tol * hmax || it >= maxit) break
  }
  # a boundary solution means the objective has no interior minimizer (e.g. the period biases
  # nearly cancel with regularize = FALSE): say so rather than return the box edge silently
  at_edge <- h <= lo * (1 + 1e-8) | h >= hmax * (1 - 1e-8)
  if (any(at_edge))
    warning("period-specific bandwidth for period(s) ",
            paste(keys[at_edge], collapse = ", "),
            " hit the search boundary [", signif(lo, 3), ", ", signif(hmax, 3),
            "]; the AMSE has no interior minimizer there. Consider regularize = TRUE, ",
            "bwselect = \"cct\", or a fixed h.")
  bws <- stats::setNames(lapply(seq_along(keys), function(j)
    c(h = unname(h[j]), b = unname(h[j] * const$ratio[j]))), keys)
  list(bws = bws, b_const = const$bias_const, v_const = const$var_const, niter = it,
       objective = amse(h), amse_fun = amse, pilot_bws = const$pilot_bws)
}

#' The aggregate AMSE as a function of one bandwidth per period
#'
#' Returns a closure `amse(hv)` over the pilot constants: squared leading bias of the aggregate,
#' the regularization term, the own-period variances, and (under "pc") the same-side
#' cross-period covariances rescaled to the bandwidth pair (hv[i], hv[j]).
#' @keywords internal
#' @noRd
.bw_iter_objective <- function(const, scheme, p, kernel, regularize, reg_const) {
  fp1  <- factorial(p + 1)
  cf   <- const$coef
  bias <- const$bias_const
  vari <- const$var_const
  nobs <- const$n_obs
  var_bias    <- const$var_bias_const
  kappa_plus  <- const$kappa_plus
  kappa_minus <- const$kappa_minus
  K <- length(const$keys)
  function(hv) {
    # leading bias of the aggregate: sum_t coef_t h_t^{p+1} B_t / (p+1)!
    Bbar <- sum(cf * hv^(p + 1) * bias) / fp1
    # regularization: reg_const x sum_t coef_t^2 (h_t^{p+1}/(p+1)!)^2 Var(B-hat_t)
    pen <- if (regularize) reg_const * sum(cf^2 * (hv^(p + 1) / fp1)^2 * var_bias) else 0
    # own-period variances: sum_t coef_t^2 V_t / (n_t h_t)
    var_term <- sum(cf^2 * vari / (nobs * hv))
    # same-side cross-period term under "pc" (each unordered pair once, factor 2; the kernel
    # constant c_side(rho) rescales the pilot covariance to the pair (h_i, h_j))
    cov_term <- 0
    if (scheme == "pc" && K >= 2L) {
      for (i in seq_len(K - 1L)) {
        for (j in (i + 1L):K) {
          rho <- hv[i] / hv[j]
          cov_term <- cov_term + 2 * cf[i] * cf[j] *
            (kappa_plus[i, j]  * .kc_c(p, "+", rho, kernel) +
             kappa_minus[i, j] * .kc_c(p, "-", rho, kernel)) / hv[j]
        }
      }
    }
    unname(Bbar^2 + pen + var_term + cov_term)
  }
}

#' Starting bandwidths of the coordinate descent, from `start`
#'
#' `"hstar"`: the common AMSE-optimal h in every period (from the same constants, no refit);
#' `"cct"`: each period's own CCT pilot h; or a numeric vector / named list, one per period.
#' @keywords internal
#' @noRd
.bw_iter_start <- function(start, const, plist, coef, t_rd, scheme, cutoff, p, q, kernel,
                           regularize, reg_const) {
  keys <- const$keys
  if (is.character(start) && length(start) == 1L) {
    start <- match.arg(start, c("hstar", "cct"))
    if (start == "hstar") {
      jb <- .bw_joint(plist, coef, t_rd, scheme = scheme, c = cutoff, p = p, q = q,
                      kernel = kernel, regularize = regularize,
                      reg_const = reg_const, constants = const)
      return(rep(jb$h, length(keys)))
    }
    return(const$h_pilot)
  }
  if (is.list(start)) start <- unlist(start)
  if (!is.numeric(start))
    stop("`start` must be \"hstar\", \"cct\", or a numeric vector/list of per-period bandwidths.")
  if (!is.null(names(start))) {
    missing_k <- setdiff(keys, names(start))
    if (length(missing_k) > 0L)
      stop("`start` is missing entries for period(s): ",
           paste(missing_k, collapse = ", "), ".")
    h0 <- unname(start[keys])
  } else {
    if (length(start) != length(keys))
      stop("`start` has length ", length(start), " but there are ", length(keys),
           " periods; supply a named vector/list or one value per period in order.")
    h0 <- unname(start)
  }
  if (any(!is.finite(h0) | h0 <= 0))
    stop("all values in `start` must be finite and positive.")
  h0
}

