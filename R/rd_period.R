# rd_period.R -- the numerical core: one period's local-linear RD (conventional and
# bias-corrected jump, per-unit influence vectors g that the cross-period covariances use).
# Called by rddid() for every period, by the bandwidth rules (pilot fits) and by the four
# assumption tests (one fit per cell). Nothing here depends on other periods.

#' Local-linear RD in one period (building block)
#'
#' The building block that [rddid()] runs in every period: the jump in the
#' outcome at the cutoff in one period, estimated by local-linear RD, both
#' conventional and bias-corrected, with standard errors. Unlike [rddid()] it
#' takes vectors (`y`, `x`, `id`), not a data frame and column names, and the
#' bandwidths must be given, for example from [rd_bw_cct()]. At a given pair
#' (`h`, `b`) it reproduces the conventional and bias-corrected estimates and
#' the conventional and robust standard errors of [rdrobust::rdrobust()] with
#' `vce = "hc1"`.
#'
#' @details
#' On each side of the cutoff the function keeps, for every observation within
#' the main or pilot bandwidth, its influence on the intercept times its
#' residual (`g`). These vectors are what the rest of the package reuses: the
#' variance of the jump is the sum of the squared `g` on both sides, and in a
#' panel the covariance between two periods' jumps sums the products of `g`
#' over the units present in both periods, matched on `id`. Variances use the
#' HC1 convention of \pkg{rdrobust}: residuals are scaled by
#' \eqn{\sqrt{n_s / (n_s - k)}}{sqrt(n_s / (n_s - k))}, with \eqn{n_s} the
#' observations used on that side and \eqn{k} the number of fitted
#' coefficients; the bias-corrected variance uses the residuals of the order-`q`
#' pilot fit at `b`.
#'
#' @param y the outcome (a numeric vector).
#' @param x the running variable (a numeric vector, same length as `y`).
#' @param h the main bandwidth (point estimate).
#' @param b the pilot bandwidth (bias correction); defaults to `h`.
#' @param id optional unit identifiers (a vector, same length as `y`), needed
#'   only to combine this period with others in a panel. With `NULL` (default)
#'   the observations are numbered `1, 2, ...`, so they cannot be matched
#'   across periods.
#' @inheritParams rddid
#'
#' @return An object of class `"rd_period"`, a list with:
#'   \describe{
#'     \item{`D`, `V_D`}{the conventional jump and its variance.}
#'     \item{`D_bc`, `V_D_bc`}{the bias-corrected jump and its robust
#'       variance.}
#'     \item{`b_const`, `v_const`}{plug-in constants used by the bandwidth
#'       rules (an estimate of the bias constant, from the gap between the
#'       conventional and bias-corrected jumps, and `n * h * V_D`).}
#'     \item{`n`}{the number of observations supplied.}
#'     \item{`h`, `b`, `c`, `p`, `q`, `kernel`}{as passed.}
#'     \item{`sides`}{a list with one element per side of the cutoff, `"+"`
#'       (above) and `"-"` (below), each holding the `id` of the observations
#'       used, their conventional and bias-corrected `g` vectors (`g`, `g_bc`)
#'       and `g_diff` (for the variance of the estimated bias), the
#'       conventional intercept at the cutoff `beta0`, its bias-corrected
#'       version `beta0_bc`, and the conventional `slope`. The fitted line on a
#'       side is `beta0 + slope * (x - c)`.}
#'   }
#'
#' @references
#' Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust nonparametric
#' confidence intervals for regression-discontinuity designs. *Econometrica*
#' 82(6), 2295-2326.
#'
#' @family RD-DID estimation
#'
#' @examples
#' # the RD period (year 3) of rddid_sim, where the jump is 0.5 + 1 = 1.5
#' d3 <- rddid_sim[rddid_sim$year == 3, ]
#' bw <- rd_bw_cct(y = d3$Y, x = d3$R)
#' fit <- rd_period(y = d3$Y, x = d3$R, h = bw[["h"]], b = bw[["b"]], id = d3$id)
#' fit
#' c(jump = fit$D_bc, se = sqrt(fit$V_D_bc))
#' @export
rd_period <- function(y, x, h, b = h, id = NULL, c = 0, p = 1L, q = 2L,
                      kernel = "triangular") {
  # one main and one pilot bandwidth per call: a vector would be recycled across the
  # observations by the kernel weights below, silently giving each unit its own bandwidth
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0)
    stop("`h` must be a single positive number.", call. = FALSE)
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= 0)
    stop("`b` must be a single positive number (one pilot bandwidth for the period).",
         call. = FALSE)
  stopifnot(length(y) == length(x), q > p, p >= 1L)
  cutoff <- c                              # `c` stays the argument name, as in rdrobust
  y <- as.numeric(y)
  x <- as.numeric(x)
  # n is the number of rows passed in, BEFORE incomplete rows are dropped: it is the n_t that
  # the bandwidth constants below are scaled by (v_const = n h V_D). Changing it would move
  # every joint/iter bandwidth, so it is kept as it is.
  n <- length(y)
  if (is.null(id)) id <- seq_len(n)
  ok <- stats::complete.cases(y, x, id)
  y  <- y[ok]
  x  <- x[ok]
  id <- id[ok]

  # one weighted local-polynomial fit per side of the cutoff
  above <- x >= cutoff
  R_side <- .rd_side_fit(y[above],  x[above],  id[above],  cutoff, h, b, p, q, kernel)   # (+)
  L_side <- .rd_side_fit(y[!above], x[!above], id[!above], cutoff, h, b, p, q, kernel)   # (-)

  D    <- R_side$beta0    - L_side$beta0
  D_bc <- R_side$beta0_bc - L_side$beta0_bc
  V_D    <- sum(R_side$g^2)    + sum(L_side$g^2)
  V_D_bc <- sum(R_side$g_bc^2) + sum(L_side$g_bc^2)

  # plug-in constants for the bandwidth rules, from the asymptotic orders
  #   Var(D_t(h)) = V_t / (n_t h)            ->  V-hat_t = n h V_D
  #   Bias(D_t(h)) = h^{p+1}/(p+1)! * B_t    ->  B-hat_t = (p+1)! (D - D_bc) / h^{p+1},
  # since D - D_bc is the estimated bias of the conventional estimator at (h, b).
  v_const <- n * h * V_D
  b_const <- factorial(p + 1L) * (D - D_bc) / h^(p + 1L)

  structure(
    list(D = D, D_bc = D_bc, V_D = V_D, V_D_bc = V_D_bc,
         b_const = b_const, v_const = v_const,
         n = n, h = h, b = b, c = cutoff, p = p, q = q, kernel = kernel,
         sides = list(`+` = R_side, `-` = L_side)),
    class = "rd_period")
}

#' Local-polynomial fit on one side of the cutoff: conventional and bias-corrected
#' intercept, and the per-unit influence vectors behind every variance in the package
#'
#' Follows Calonico, Cattaneo and Titiunik (2014): the order-p fit at bandwidth h gives the
#' conventional intercept; an order-q fit at the pilot bandwidth b estimates the (p+1)-th
#' derivative, whose contribution is subtracted to give the bias-corrected intercept. Both
#' estimators are linear in y, so each has a per-unit weight (its "influence"), and
#' g = weight x residual is what the variances and cross-period covariances sum.
#'
#' @return list: `beta0`, `beta0_bc` (intercepts), `slope`, `id` (units in the active set),
#'   `g`, `g_bc` (influence x residual, conventional and bias-corrected), `g_diff` (influence of
#'   conventional minus bias-corrected, paired with the pilot residuals; feeds Var(B-hat) in the
#'   bandwidth regularization).
#' @keywords internal
#' @noRd
.rd_side_fit <- function(ys, xs, ids, cutoff, h, b, p, q, kernel) {
  u_h <- (xs - cutoff) / h
  u_b <- (xs - cutoff) / b
  w_h <- .rd_kweight(u_h, kernel)
  w_b <- .rd_kweight(u_b, kernel)
  # the active set is the union of the two windows; units outside both have zero weight in
  # every quantity below, and the pilot window (b >= h usually) is the larger one
  active <- (w_h > 0) | (w_b > 0)
  # the order-q fit needs more than q + 1 points, or its HC1 factor n_s / (n_s - (q + 1)) is
  # infinite
  if (sum(active) <= q + 1L || length(unique(xs[w_h > 0])) <= p + 1L ||
      length(unique(xs[w_b > 0])) <= q + 1L)
    stop("too few distinct running-variable values within the bandwidth on one side of the ",
         "cutoff (need more than p + 1 in the main window and more than q + 1 in the pilot ",
         "window); use a wider bandwidth or check the data near the cutoff.")
  xs  <- xs[active]
  ys  <- ys[active]
  ids <- ids[active]
  w_h <- w_h[active]
  w_b <- w_b[active]

  # design matrices: order q (pilot fit), with the order-p design as its first p + 1 columns
  Rq <- outer(xs - cutoff, 0:q, `^`)
  Rp <- Rq[, 1:(p + 1L), drop = FALSE]
  invG_p <- .qrXXinv(sqrt(w_h) * Rp)     # (X' A(h) X)^{-1}
  invG_q <- .qrXXinv(sqrt(w_b) * Rq)     # (X' A(b) X)^{-1}, order q

  # bias-correction weight matrix (CCT eq. for the bias-corrected estimator, without the 1/h
  # and 1/n scalings, which cancel): Q = X'A(h) - h^{p+1} * theta * e_{p+1}' Gq^{-1} X'A(b).
  # e_p1 picks the (p+1)-th power, which sits at position p + 2 of the length-(q+1) vector.
  # Keep the grouping of the matrix products exactly as written: a regrouping changes the last
  # bits of every bias-corrected number.
  e_p1 <- numeric(q + 1L)
  e_p1[p + 2L] <- 1
  theta <- crossprod(Rp * w_h, ((xs - cutoff) / h)^(p + 1L))  # X' A(h) S_{p+1}
  Aq_b  <- t(Rq * w_b)                                        # X' A(b)
  Qmat  <- t(Rp * w_h) - h^(p + 1L) * (theta %*% (t(e_p1) %*% invG_q %*% Aq_b))

  beta_p  <- invG_p %*% crossprod(Rp * w_h, ys)   # conventional coefficients
  beta_q  <- invG_q %*% crossprod(Rq * w_b, ys)   # order-q coefficients (pilot residuals)
  beta_bc <- invG_p %*% (Qmat %*% ys)             # bias-corrected coefficients

  # influence of each unit on the intercept (first row of the hat-type matrix), then
  # g = influence x residual, with rdrobust's HC1 factor sqrt(n_s / (n_s - k)) on the residuals
  # of the fit they come from: order p for the conventional, order q for the bias-corrected
  a_c   <- as.numeric(invG_p[1, ] %*% t(Rp * w_h))
  a_bc  <- as.numeric(invG_p[1, ] %*% Qmat)
  res_c <- sqrt(length(ys) / (length(ys) - (p + 1L))) * (ys - Rp %*% beta_p)
  res_b <- sqrt(length(ys) / (length(ys) - (q + 1L))) * (ys - Rq %*% beta_q)

  list(
    beta0    = beta_p[1L],
    beta0_bc = beta_bc[1L],
    slope    = as.numeric(beta_p[2L]),
    id       = ids,
    g        = a_c * as.numeric(res_c),
    g_bc     = a_bc * as.numeric(res_b),
    # influence on (conventional - bias-corrected): both are linear in the same y, so the
    # difference has per-unit weight (a_c - a_bc), paired with the pilot residuals like g_bc.
    # Used for Var(B-hat) in the bandwidth regularization.
    g_diff   = (a_c - a_bc) * as.numeric(res_b)
  )
}

#' @export
print.rd_period <- function(x, ...) {
  cat(sprintf("Single-period RD (p=%d, h=%.4g, b=%.4g, kernel=%s, n=%d)\n",
              x$p, x$h, x$b, x$kernel, x$n))
  cat(sprintf("  D (conventional)   = %+.5g  (se %.4g)\n", x$D, sqrt(x$V_D)))
  cat(sprintf("  D (bias-corrected) = %+.5g  (se %.4g)\n", x$D_bc, sqrt(x$V_D_bc)))
  invisible(x)
}

# ---- Cholesky inverse of a Gram matrix (used only by .rd_side_fit) --------------------------
# Columns are rescaled to unit norm before the Cholesky factorisation and the inverse is scaled
# back, which keeps the inverse accurate when the polynomial columns differ in magnitude.
#' Inverse of a weighted Gram matrix via Cholesky, given the square-root design
#'
#' The columns of `x` are powers of the centred running variable, so
#' `crossprod(x)` has a condition number of order 1e4–1e6 at orders 2–3 for
#' typical bandwidths. The inverse is computed after scaling each column of `x`
#' to unit norm and undoing the scaling afterwards (exact algebra, `G = D G* D`),
#' which keeps the result stable to ~1e-13 across BLAS implementations instead
#' of ~1e-10.
#' @keywords internal
#' @noRd
.qrXXinv <- function(x) {
  s  <- sqrt(colSums(x^2))
  Gi <- chol2inv(chol(crossprod(x / rep(s, each = nrow(x)))))
  Gi / outer(s, s)
}
