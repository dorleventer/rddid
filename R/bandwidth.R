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
#' coordinate problem AMSE(h_tau | rest) is coercive (Appendix C, Lemma
#' lem:coercive), so it always has an interior minimum -- no degeneracy.
#' The starting point of the descent is controlled by `start` (default: the
#' common joint-optimal h*).
#'
#' Under `scheme = "pc"` the same-side cross-period covariance term is included
#' in the AMSE (it is O(1/(nh)), same order as the per-period variance).  Under
#' `"cs"` and `"pv"` it is omitted (CS covariance is zero; PV covariances are
#' o(1/(nh)) and dropped for bandwidth selection per Appendix C).
#'
#' The PC covariance term is
#'   `2 * sum_{tau<rho} cf_tau*cf_rho*P_tr*M0_tr / max(hv_tau,hv_rho)`
#' where `P_tr` is the pilot same-side cross-period covariance (from
#' `.cross_cov()$pc` at the per-period CCT pilot fits) and `M0_tr` = max of
#' the two pilot bandwidths.  This is the pilot-scaling form: since
#' `Cov = c/(n*max(h))` and `c/n = P*max(h0)`, n cancels.
#'
#' @param plist named per-period data; `coef` named period coefficients;
#'   `t_rd` RD-period label.
#' @param scheme one of `"cs"`, `"pc"`, `"pv"`; determines whether the
#'   cross-period covariance term enters the AMSE (only for `"pc"`).
#' @param pilot_bws optional named list of per-period `c(h, b)` used to
#'   estimate the per-period constants and the regularization penalty; defaults
#'   to CCT.  The role of `pilot_bws` is to supply the asymptotic constants
#'   (curvature, variance, regularization, PC covariance) — it does NOT control
#'   the descent starting point; that is controlled by `start`.
#' @param start seed for the coordinate descent.  Three modes:
#'   \describe{
#'     \item{`"hstar"` (default)}{Start all periods at the common joint-optimal
#'       h* (the same scalar h* for every period).  This is the original behaviour.}
#'     \item{`"cct"`}{Start each period at its own CCT/IK pilot h (i.e.
#'       `pilot_bws[[k]]["h"]` for period `k`).  Reproduces the pre-h* seed.}
#'     \item{numeric / named list}{Supply a manual per-period seed.  Must cover
#'       every period in `names(coef)`.  A numeric vector is matched positionally
#'       to `names(coef)`; a named vector/list is matched by name.}
#'   }
#' @param hmax cap on each bandwidth (defaults to the running-variable radius).
#' @param maxit,tol coordinate-descent controls.
#' @return list with per-period `bws` (named list of `c(h, b)`), the constants
#'   `b_const`/`v_const`, and iteration count `niter`.
#' @keywords internal
#' @noRd
.bw_joint_iter <- function(plist, coef, t_rd, scheme = "cs", pilot_bws = NULL,
                           start = "hstar",
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

  # PC covariance precomputation ------------------------------------------
  # P_mat[i,j] = pilot same-side cross-period covariance (from .cross_cov)
  # M0_mat[i,j] = max of the two pilot bandwidths
  # Both are only computed (and used) when scheme == "pc".
  K <- length(keys)
  P_mat  <- matrix(0, K, K)
  M0_mat <- matrix(0, K, K)
  h0v <- vapply(keys, function(k) unname(pilot_bws[[k]]["h"]), numeric(1))
  if (scheme == "pc" && K >= 2L) {
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        P_mat[i, j]  <- .cross_cov(fitp[[keys[i]]], fitp[[keys[j]]], bc = FALSE)$pc
        P_mat[j, i]  <- P_mat[i, j]
        M0_mat[i, j] <- max(h0v[i], h0v[j])
        M0_mat[j, i] <- M0_mat[i, j]
      }
    }
  }

  amse <- function(hv) {
    Bbar <- 0.5 * sum(cf * hv^2 * bt)
    # regularization: each period's bias-estimate variance penalises wide h_t,
    # so a poorly-estimated curvature cannot push its bandwidth to an extreme.
    pen <- if (regularize) reg_const * sum(cf^2 * (hv^4 / 4) * var_b) else 0
    var_term <- sum(cf^2 * vt / (nt * hv))
    # PC same-side cross-period covariance term (Appendix C, eq. amse-ps).
    # For CS and PV this is zero (zero covariance / lower order respectively).
    cov_term <- 0
    if (scheme == "pc" && K >= 2L) {
      for (i in seq_len(K - 1L)) {
        for (j in (i + 1L):K) {
          cov_term <- cov_term +
            2 * cf[i] * cf[j] * P_mat[i, j] * M0_mat[i, j] / max(hv[i], hv[j])
        }
      }
    }
    Bbar^2 + pen + var_term + cov_term
  }
  # Seed coordinate descent -- mode controlled by `start`.
  lo <- 0.03 * hmax
  if (is.character(start) && length(start) == 1L) {
    start <- match.arg(start, c("hstar", "cct"))
    if (start == "hstar") {
      # Default: start all periods at the common joint-optimal h*.
      jb <- .bw_joint(plist, coef, t_rd, scheme = scheme,
                      pilot = pilot_bws[[t_rd]], c = c, p = p, q = q,
                      kernel = kernel, regularize = regularize,
                      reg_const = reg_const)
      h0 <- rep(jb$h, length(keys))
    } else {
      # "cct": start each period at its own CCT pilot h.
      h0 <- vapply(keys, function(k) unname(pilot_bws[[k]]["h"]), numeric(1))
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
  bws <- stats::setNames(lapply(seq_along(keys), function(j)
    c(h = unname(h[j]), b = unname(h[j] * ratio[j]))), keys)
  list(bws = bws, b_const = bt, v_const = vt, niter = it)
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
