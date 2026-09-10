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

#' Per-period CCT/IK MSE-optimal bandwidths
#'
#' @param plist named list of per-period data frames with columns `y`, `x`, `id`.
#' @keywords internal
#' @noRd
.bw_cct <- function(plist, c = 0, p = 1L, kernel = "triangular") {
  # rd_bw_cct() wraps rdrobust::rdbwselect(bwselect = "mserd") with the
  # finiteness / failure guards and the documented 0.5*IQR fallback.
  lapply(plist, function(d) rd_bw_cct(d$y, d$x, c = c, p = p, kernel = kernel))
}

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
#' estimator aggregates several RD periods (drafts/agg_att_note.tex, 2026-09-10).
#'
#' @param plist named per-period data list; `coef` named period coefficients.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`.
#' @param pilot_bws optional named list of per-period `c(h, b)`; defaults to CCT.
#' @return list with `keys`, `cf` (coefficients in `keys` order), `bt`, `vt`, `nt`,
#'   `var_b`, `h0v`, `ratio`, `kap_p`, `kap_m`, `fitp` (the pilot fits), `pilot_bws`.
#' @keywords internal
#' @noRd
.bw_constants <- function(plist, coef, scheme = "cs", pilot_bws = NULL,
                          c = 0, p = 1L, q = 2L, kernel = "triangular") {
  keys <- names(coef)
  fp1 <- factorial(p + 1)
  if (is.null(pilot_bws))
    pilot_bws <- .bw_cct(plist, c = c, p = p, kernel = kernel)
  missing_k <- setdiff(keys, names(pilot_bws))
  if (length(missing_k) > 0L)
    stop("`pilot_bws` is missing entries for period(s): ",
         paste(missing_k, collapse = ", "), ".")

  # per-period asymptotic constants from a pilot fit: B-hat_t (curvature, via the
  # conventional/BC gap) and V-hat_t (variance constant); rd_period returns both
  # (eq:per-period-orders).
  fitp <- lapply(keys, function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = pilot_bws[[k]][["h"]],
              b = pilot_bws[[k]][["b"]], id = plist[[k]]$id,
              c = c, p = p, q = q, kernel = kernel))
  names(fitp) <- keys
  cf <- unname(vapply(keys, function(k) coef[[k]], numeric(1)))
  bt <- vapply(keys, function(k) fitp[[k]]$b_const, numeric(1))
  vt <- vapply(keys, function(k) fitp[[k]]$v_const, numeric(1))
  nt <- vapply(keys, function(k) fitp[[k]]$n,       numeric(1))
  ratio <- vapply(keys, function(k)
    unname(pilot_bws[[k]][["b"]] / pilot_bws[[k]][["h"]]), numeric(1))
  h0v <- vapply(keys, function(k) unname(pilot_bws[[k]][["h"]]), numeric(1))

  # per-period variance of the bias-constant estimate, Var(B-hat_t) =
  # ((p+1)!/h0^{p+1})^2 Var(D - D_bc), for the B.4 regularization term.
  var_b <- vapply(keys, function(k) {
    s <- fitp[[k]]$sides
    (fp1 / h0v[[k]]^(p + 1))^2 * (sum(s[["+"]]$g_diff^2) + sum(s[["-"]]$g_diff^2))
  }, numeric(1))

  # PC covariance precomputation (see the .bw_joint_iter header): per side,
  # kappa_side[i, j] = pilot same-side covariance * h0_j / c_side(p, h0_i / h0_j).
  # Only under scheme "pc".
  K <- length(keys)
  kap_p <- matrix(0, K, K)
  kap_m <- matrix(0, K, K)
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        cc   <- .cross_cov(fitp[[keys[i]]], fitp[[keys[j]]], bc = FALSE)
        rho0 <- h0v[i] / h0v[j]
        kap_p[i, j] <- cc$pc_p * h0v[j] / .kc_c(p, "+", rho0, kernel)
        kap_m[i, j] <- cc$pc_m * h0v[j] / .kc_c(p, "-", rho0, kernel)
      }
    }
    if (all(kap_p == 0) && all(kap_m == 0))
      warning("scheme = \"pc\" but no unit is active in two periods' windows: ",
              "the cross-period term is identically zero (same as scheme = \"cs\").")
  }
  list(keys = keys, cf = cf, bt = bt, vt = vt, nt = nt, var_b = var_b,
       h0v = h0v, ratio = ratio, kap_p = kap_p, kap_m = kap_m,
       fitp = fitp, pilot_bws = pilot_bws)
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
#' Before 2026-09-10 every period was fitted at the RD period's CCT pilot; the
#' selector then depended on which period carried the `t_rd` label, which is
#' arbitrary for an aggregate over several RD periods. The label no longer enters.
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
  if (is.null(constants))
    constants <- .bw_constants(plist, coef, scheme = scheme, pilot_bws = pilot_bws,
                               c = c, p = p, q = q, kernel = kernel)
  cs <- constants
  fp1 <- factorial(p + 1)
  K <- length(cs$keys)

  # estimated aggregate bias constant B(t_RD) = sum_tau coef_tau B-hat_tau (B.4 P3)
  B <- sum(cs$cf * cs$bt)
  # variance scale: own-period terms sum_tau coef_tau^2 V-hat_tau / n_tau (lem:agg-var)
  Veff <- sum(cs$cf^2 * cs$vt / cs$nt)
  # PC same-side cross-period term at a common h (rho = 1), lem:agg-var / B.4 P3
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        Veff <- Veff + 2 * cs$cf[i] * cs$cf[j] *
          (cs$kap_p[i, j] * .kc_c(p, "+", 1, kernel) +
           cs$kap_m[i, j] * .kc_c(p, "-", 1, kernel))
      }
    }
  }
  # Regularization (B.4): lambda * sum_tau coef_tau^2 (h^{p+1}/(p+1)!)^2 Var(B-hat_tau).
  # At a common h this carries the same h^{2(p+1)} factor as B^2, so it enters the
  # denominator of h*.
  reg <- if (regularize) reg_const * sum(cs$cf^2 * cs$var_b) else 0
  denom <- B^2 + reg

  if (!is.finite(denom) || denom <= 0)
    stop("aggregate bias constant B is ~0 and regularization is off: the ",
         "AMSE-optimal bandwidth diverges. Use bwselect = \"cct\" or pass h.")
  if (!is.finite(Veff) || Veff <= 0)
    stop("aggregate variance scale is not positive; check the sampling scheme and ",
         "the per-period fits.")
  # eq:common_h_opt; for p = 1 the leading constant is ((2)!)^2/(2*2) = 1 and the
  # exponent 1/5.
  h_star <- (fp1^2 / (2 * (p + 1)) * Veff / denom)^(1 / (2 * p + 3))
  radius <- max(vapply(plist, function(d) max(abs(d$x - c), na.rm = TRUE), numeric(1)))
  if (is.finite(radius) && h_star > radius)
    warning("joint AMSE-optimal bandwidth h* = ", signif(h_star, 4),
            " exceeds the running-variable radius ", signif(radius, 4),
            " (the period biases nearly cancel); consider bwselect = \"cct\" or a fixed h.")
  b <- stats::setNames(h_star * cs$ratio, cs$keys)
  list(h = h_star, b = b, B = B, Veff = Veff, reg = reg,
       pilot_bws = cs$pilot_bws, constants = cs)
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
  cs <- .bw_constants(plist, coef, scheme = scheme, pilot_bws = pilot_bws,
                      c = c, p = p, q = q, kernel = kernel)
  keys <- cs$keys
  fp1 <- factorial(p + 1)
  cf <- cs$cf; bt <- cs$bt; vt <- cs$vt; nt <- cs$nt
  ratio <- cs$ratio; h0v <- cs$h0v; var_b <- cs$var_b
  kap_p <- cs$kap_p; kap_m <- cs$kap_m
  K <- length(keys)
  if (is.null(hmax))
    hmax <- max(vapply(plist, function(d) max(abs(d$x - c), na.rm = TRUE), numeric(1)))

  amse <- function(hv) {
    # leading bias of the aggregate, B.4 P2: sum_t w~_t h_t^{p+1} B_t / (p+1)!
    Bbar <- sum(cf * hv^(p + 1) * bt) / fp1
    # regularization: lambda * sum_t w~_t^2 (h_t^{p+1}/(p+1)!)^2 Var(B-hat_t)
    pen <- if (regularize) reg_const * sum(cf^2 * (hv^(p + 1) / fp1)^2 * var_b) else 0
    # own-period variances, lem:agg-var: sum_t w~_t^2 V_t / (n_t h_t)
    var_term <- sum(cf^2 * vt / (nt * hv))
    # PC same-side cross-period term, lem:agg-var / B.4 P4 (ordered double sum =
    # 2 x the unordered sum, by c(1/rho) = rho c(rho))
    cov_term <- 0
    if (scheme == "pc" && K >= 2L) {
      for (i in seq_len(K - 1L)) {
        for (j in (i + 1L):K) {
          rho <- hv[i] / hv[j]
          cov_term <- cov_term + 2 * cf[i] * cf[j] *
            (kap_p[i, j] * .kc_c(p, "+", rho, kernel) +
             kap_m[i, j] * .kc_c(p, "-", rho, kernel)) / hv[j]
        }
      }
    }
    unname(Bbar^2 + pen + var_term + cov_term)
  }
  # Seed coordinate descent -- mode controlled by `start`.
  lo <- 0.03 * hmax
  if (is.character(start) && length(start) == 1L) {
    start <- match.arg(start, c("hstar", "cct"))
    if (start == "hstar") {
      # Default: start all periods at the common joint-optimal h*, from the same
      # constants (no refit).
      jb <- .bw_joint(plist, coef, t_rd, scheme = scheme, c = c, p = p, q = q,
                      kernel = kernel, regularize = regularize,
                      reg_const = reg_const, constants = cs)
      h0 <- rep(jb$h, length(keys))
    } else {
      # "cct": start each period at its own CCT pilot h.
      h0 <- h0v
    }
  } else {
    # Manual: numeric vector or named list supplied by the user.
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
  }
  h <- pmin(pmax(h0, lo), hmax)
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
  # eq:update is an unconstrained argmin; the search box [lo, hmax] is the package's.
  # A boundary solution means the objective has no interior minimizer (e.g. the
  # period biases nearly cancel with regularize = FALSE): say so rather than
  # return the box edge silently.
  at_edge <- h <= lo * (1 + 1e-8) | h >= hmax * (1 - 1e-8)
  if (any(at_edge))
    warning("period-specific bandwidth for period(s) ",
            paste(keys[at_edge], collapse = ", "),
            " hit the search boundary [", signif(lo, 3), ", ", signif(hmax, 3),
            "]; the AMSE has no interior minimizer there. Consider regularize = TRUE, ",
            "bwselect = \"cct\", or a fixed h.")
  bws <- stats::setNames(lapply(seq_along(keys), function(j)
    c(h = unname(h[j]), b = unname(h[j] * ratio[j]))), keys)
  list(bws = bws, b_const = bt, v_const = vt, niter = it,
       objective = amse(h), amse_fun = amse, pilot_bws = cs$pilot_bws)
}

#' CCT (MSE-optimal) bandwidth for a single local-linear RD
#'
#' Returns the Calonico–Cattaneo–Titiunik MSE-optimal bandwidths `h` and `b`
#' for a single local-linear RD using \pkg{rdrobust}. Falls back gracefully when
#' \pkg{rdrobust} is unavailable, the call fails, or the returned bandwidth is
#' non-positive/non-finite.
#'
#' @param y Outcome vector.
#' @param x Running variable vector.
#' @param c Cutoff (default 0).
#' @param p Polynomial order (default 1L, local linear).
#' @param kernel Kernel type: `"triangular"` (default), `"epanechnikov"`, or
#'   `"uniform"`.
#'
#' @return A named numeric vector `c(h = ..., b = ...)` with the main and pilot
#'   bandwidths.  When the CCT computation is unavailable, both equal
#'   `0.5 * IQR(x)` (or `sd(x)` if IQR is zero), and a message is emitted
#'   naming the reason.
#' @export
rd_bw_cct <- function(y, x, c = 0, p = 1L, kernel = "triangular") {
  fallback <- function(reason) {
    h0 <- 0.5 * stats::IQR(x)
    if (!is.finite(h0) || h0 <= 0) h0 <- stats::sd(x)
    message("rd_bw_cct: CCT unavailable (", reason,
            "); falling back to 0.5*IQR h=", round(h0, 4))
    c(h = h0, b = h0)
  }
  if (!requireNamespace("rdrobust", quietly = TRUE))
    return(fallback("rdrobust not installed"))
  bw <- tryCatch(
    rdrobust::rdbwselect(y = y, x = x, c = c, p = p, kernel = kernel,
                         bwselect = "mserd"),
    error = function(e) e
  )
  if (inherits(bw, "error"))
    return(fallback(conditionMessage(bw)))
  h_val <- as.numeric(bw$bws[1, 1])
  b_val <- as.numeric(bw$bws[1, 3])
  if (!is.finite(h_val) || h_val <= 0)
    return(fallback(paste0("rdbwselect returned h=", h_val)))
  c(h = h_val, b = b_val)
}
