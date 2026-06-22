# Bandwidth selection for the aggregate RD-DID estimator.
#
# Two modes, matching Section 4.3 of Leventer and Nevo:
#   * "cct"   - within-period MSE-optimal bandwidth (Imbens-Kalyanaraman / CCT)
#               applied to each D-hat_t on its own. Never degenerates; leaves
#               the cross-period bias cancellation unexploited. Each period uses
#               its own (h_t, b_t).
#   * "joint" - one common h* minimising the asymptotic MSE of the AGGREGATE
#               estimator, h* = (V / B^2)^(1/5) n^(-1/5), with
#               B = b_{tRD} - sum_t w_t b_t (signed: comparison biases can offset
#               the RD-period bias) and V the aggregate variance constant.

#' Per-period CCT/IK MSE-optimal bandwidths
#'
#' @param plist named list of per-period data frames with columns `y`, `x`, `id`.
#' @keywords internal
#' @noRd
.bw_cct <- function(plist, c = 0, p = 1L, kernel = "triangular") {
  if (!requireNamespace("rdrobust", quietly = TRUE))
    stop("bwselect = \"cct\" needs the rdrobust package; install it or pass h.")
  lapply(plist, function(d) {
    bw <- rdrobust::rdbwselect(y = d$y, x = d$x, c = c, p = p,
                               kernel = kernel, bwselect = "mserd")
    # bws columns: h(left), h(right), b(left), b(right)
    c(h = as.numeric(bw$bws[1, 1]), b = as.numeric(bw$bws[1, 3]))
  })
}

#' Joint AMSE-optimal common bandwidth for the aggregate estimator
#'
#' Feasible plug-in: fit every period at a pilot bandwidth, read the aggregate
#' bias constant off the bias correction (ATT_conv - ATT_bc = (h0^2/2) B), and
#' the variance scale off h0 * V_agg(h0); then h* = (Veff / B^2)^(1/5).
#'
#' @param plist named per-period data list; `coef` named period coefficients
#'   (RD period = +1, comparisons = -w_t); `t_rd` the RD-period label.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`; sets which variance scales the
#'   bandwidth (under `"pv"` the cross-period covariances are o(1/(nh)) and drop,
#'   so the covariance-free CS form is used, per Section 4.3).
#' @param pilot optional `c(h, b)` pilot bandwidth; defaults to the RD period's
#'   CCT/IK bandwidth.
#' @param regularize add a regularization term to the denominator so a small,
#'   noisy estimated bias constant cannot inflate `h*` (mirrors the
#'   `regularize = TRUE` default of [rdrobust::rdbwselect]). Default `TRUE`.
#' @param reg_const multiple of the bias-estimate variance used as the
#'   regularization term (default 3, the CCT convention).
#' @return list with `h`, `b`, the bias constant `B`, the variance scale `Veff`,
#'   the regularization term `reg`, and the `pilot` used.
#' @keywords internal
#' @noRd
.bw_joint <- function(plist, coef, t_rd, scheme = "cs", pilot = NULL,
                      c = 0, p = 1L, q = 2L, kernel = "triangular",
                      regularize = TRUE, reg_const = 3) {
  if (is.null(pilot)) {
    if (!requireNamespace("rdrobust", quietly = TRUE))
      stop("joint bandwidth needs a pilot: install rdrobust or pass pilot = c(h, b).")
    bw <- rdrobust::rdbwselect(y = plist[[t_rd]]$y, x = plist[[t_rd]]$x, c = c,
                               p = p, kernel = kernel, bwselect = "mserd")
    pilot <- c(h = as.numeric(bw$bws[1, 1]), b = as.numeric(bw$bws[1, 3]))
  }
  h0 <- pilot[["h"]]; b0 <- pilot[["b"]]

  fits <- lapply(names(coef), function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = h0, b = b0, id = plist[[k]]$id,
              c = c, p = p, q = q, kernel = kernel))
  names(fits) <- names(coef)

  aggc <- .aggregate_fits(fits, coef, bc = FALSE)
  aggb <- .aggregate_fits(fits, coef, bc = TRUE)
  Vfld <- switch(scheme, cs = "V_cs", pc = "V_pc", pv = "V_cs")

  # estimated aggregate bias constant: ATT_conv - ATT_bc = (h0^2/2) B
  B <- 2 * (aggc[["est"]] - aggb[["est"]]) / h0^2
  Veff <- h0 * aggc[[Vfld]]

  # CCT-style regularization: R = reg_const * Var(B-hat), B-hat = 2*(conv-bc)/h0^2,
  # so Var(B-hat) = (2/h0^2)^2 Var(conv-bc). Var(conv-bc) is small (the two are
  # nearly collinear); estimate it from the per-unit influence difference g_diff.
  var_delta <- sum(vapply(names(coef), function(k) {
    s <- fits[[k]]$sides
    coef[[k]]^2 * (sum(s[["+"]]$g_diff^2) + sum(s[["-"]]$g_diff^2))
  }, numeric(1)))
  reg <- if (regularize) reg_const * (2 / h0^2)^2 * var_delta else 0
  denom <- B^2 + reg

  if (!is.finite(denom) || denom <= 0)
    stop("aggregate bias constant B is ~0 and regularization is off: the ",
         "AMSE-optimal bandwidth diverges. Use bwselect = \"cct\" or pass h.")
  h_star <- (Veff / denom)^(1 / 5)
  list(h = h_star, b = h_star * (b0 / h0), B = B, Veff = Veff, reg = reg,
       pilot = pilot)
}

#' Period-specific joint AMSE bandwidths by coordinate descent
#'
#' Minimises the aggregate AMSE of ATT-hat over a SEPARATE bandwidth per period,
#' rather than a single common h. The simultaneous optimum can corner (Appendix
#' C), so we cycle one period at a time: holding the others fixed, each
#' coordinate problem AMSE(h_tau | rest) = Bbar(h_tau)^2 + omega_tau v_tau /
#' (n h_tau) is coercive (-> Inf as h_tau -> 0 and as h_tau -> Inf, since
#' Bbar grows like h_tau^2), so it always has an interior minimum -- no
#' degeneracy. Initialised at per-period CCT.
#'
#' AMSE uses the (CS / diagonal) form Bbar = (1/2) sum_tau coef_tau h_tau^2 b_tau,
#' V = sum_tau coef_tau^2 v_tau / (n_tau h_tau); cross-period covariances are
#' o(1/(nh)) (time-varying R) or folded into the per-period v at the bandwidth
#' scale, matching the App C period-specific treatment.
#'
#' @param plist named per-period data; `coef` named period coefficients;
#'   `t_rd` RD-period label.
#' @param pilot_bws optional named list of per-period `c(h, b)` used both to
#'   initialise and to estimate the per-period constants; defaults to CCT.
#' @param hmax cap on each bandwidth (defaults to the running-variable radius).
#' @param maxit,tol coordinate-descent controls.
#' @return list with per-period `bws` (named list of `c(h, b)`), the constants
#'   `b_const`/`v_const`, and iteration count `niter`.
#' @keywords internal
#' @noRd
.bw_joint_iter <- function(plist, coef, t_rd, pilot_bws = NULL,
                           c = 0, p = 1L, q = 2L, kernel = "triangular",
                           regularize = TRUE, reg_const = 3,
                           hmax = NULL, maxit = 50L, tol = 1e-4) {
  keys <- names(coef)
  if (is.null(pilot_bws))
    pilot_bws <- .bw_cct(plist, c = c, p = p, kernel = kernel)

  # per-period asymptotic constants from a pilot fit (b_t = curvature gap,
  # v_t = variance constant); rd_period already returns these.
  fitp <- lapply(keys, function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = pilot_bws[[k]]["h"],
              b = pilot_bws[[k]]["b"], id = plist[[k]]$id,
              c = c, p = p, q = q, kernel = kernel))
  names(fitp) <- keys
  cf <- vapply(keys, function(k) coef[[k]],        numeric(1))
  bt <- vapply(keys, function(k) fitp[[k]]$b_const, numeric(1))
  vt <- vapply(keys, function(k) fitp[[k]]$v_const, numeric(1))
  nt <- vapply(keys, function(k) fitp[[k]]$n,       numeric(1))
  ratio <- vapply(keys, function(k)
    unname(pilot_bws[[k]]["b"] / pilot_bws[[k]]["h"]), numeric(1))
  if (is.null(hmax))
    hmax <- max(vapply(plist, function(d) max(abs(d$x - c)), numeric(1)))

  # per-period variance of the bias-constant estimate, Var(b_hat_t) =
  # (2/h0^2)^2 Var(conv-bc), for regularizing wide bandwidths on noisy curvature.
  var_b <- vapply(keys, function(k) {
    s <- fitp[[k]]$sides
    (2 / pilot_bws[[k]]["h"]^2)^2 * (sum(s[["+"]]$g_diff^2) + sum(s[["-"]]$g_diff^2))
  }, numeric(1))

  amse <- function(hv) {
    Bbar <- 0.5 * sum(cf * hv^2 * bt)
    # regularization: each period's bias-estimate variance penalises wide h_t,
    # so a poorly-estimated curvature cannot push its bandwidth to an extreme.
    pen <- if (regularize) reg_const * sum(cf^2 * (hv^4 / 4) * var_b) else 0
    Bbar^2 + pen + sum(cf^2 * vt / (nt * hv))
  }
  h <- vapply(keys, function(k) unname(pilot_bws[[k]]["h"]), numeric(1))
  lo <- 0.03 * hmax
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
  bws <- stats::setNames(lapply(seq_along(keys), function(j)
    c(h = unname(h[j]), b = unname(h[j] * ratio[j]))), keys)
  list(bws = bws, b_const = bt, v_const = vt, niter = it)
}
