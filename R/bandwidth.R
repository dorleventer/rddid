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
