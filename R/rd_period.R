#' Single-period local-linear RD discontinuity
#'
#' Estimates the period-\eqn{t} outcome discontinuity
#' \eqn{D_t = \beta^{(0)}_{(+)} - \beta^{(0)}_{(-)}} by a standard local-linear
#' regression discontinuity, conventional and robust-bias-corrected
#' (Calonico, Cattaneo and Titiunik 2014). This is the per-period engine that
#' [rddid()] aggregates across periods; it implements the matrix form of
#' Appendix B (`app:est-bias`, `app:est-var`) of Leventer and Nevo. The
#' object-by-object map between this code and the paper's equations is
#' `dev/appB_map.md`.
#'
#' The function returns, for each side of the cutoff, the per-unit influence
#' weight on the intercept times its local-linear residual, `g`. These `g`
#' vectors are the single primitive the rest of the package reuses:
#' \itemize{
#'   \item the discontinuity variance is \eqn{V(\hat D_t) = \sum g_{(+)}^2 + \sum g_{(-)}^2};
#'   \item any cross-period covariance is a merge of two periods' `g` vectors on
#'     shared unit `id` (see the internal covariance helpers);
#'   \item the bandwidth constants are read off `V_D` and the conventional /
#'     bias-corrected gap.
#' }
#' Variances use the HC1 finite-sample convention of \pkg{rdrobust}
#' (`vce = "hc1"`): residuals are scaled by `sqrt(n_s / (n_s - k))` with `n_s`
#' the side's active sample and `k` the number of fitted coefficients; the BC
#' variance uses the residuals of the order-`q` pilot fit at `b`, as in
#' \pkg{rdrobust}. The paper's Appendix B states neither convention (see
#' `dev/appB_map.md`, section 0).
#'
#' @param y outcome vector.
#' @param x running variable.
#' @param h main (point-estimate) bandwidth.
#' @param b pilot (bias-correction) bandwidth; defaults to `h`.
#' @param id optional unit identifiers, needed only when this period will be
#'   combined with others under panel sampling. If `NULL`, sequential ids are
#'   assigned and no cross-period matching is possible.
#' @param c cutoff (default 0).
#' @param p point-estimate polynomial order (default 1, local linear).
#' @param q bias-correction polynomial order (default 2); must exceed `p`.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#'
#' @return An object of class `"rd_period"`: a list with the conventional and
#'   bias-corrected discontinuity (`D`, `D_bc`) and their variances (`V_D`,
#'   `V_D_bc`), the per-period plug-in bias and variance constants used for
#'   joint bandwidth selection (`b_const = (p+1)! (D - D_bc) / h^{p+1}` and
#'   `v_const = n h V_D`, the plug-ins for the constants `B_t` and `V_t` of
#'   `eq:per-period-orders` in Appendix B.3), the effective
#'   sample size `n`, the bandwidths, and a per-side list `sides` holding, for
#'   each side `"+"` and `"-"`, the active units' `id`, their conventional /
#'   bias-corrected `g` vectors, the conventional intercept at the cutoff
#'   (`beta0`), and the conventional local-linear `slope`.  The fitted line on
#'   a side is `beta0 + slope * (x - c)`.
#' @export
rd_period <- function(y, x, h, b = h, id = NULL, c = 0, p = 1L, q = 2L,
                      kernel = "triangular") {
  stopifnot(length(y) == length(x), h > 0, b > 0, q > p, p >= 1L)
  y <- as.numeric(y)
  x <- as.numeric(x)
  n <- length(y)
  if (is.null(id)) id <- seq_len(n)
  ok <- stats::complete.cases(y, x, id)
  y <- y[ok]; x <- x[ok]; id <- id[ok]

  side_fit <- function(keep) {
    xs <- x[keep]; ys <- y[keep]; ids <- id[keep]
    u_h <- (xs - c) / h
    u_b <- (xs - c) / b
    w_h <- .rd_kweight(u_h, kernel)        # main-bandwidth kernel weights
    w_b <- .rd_kweight(u_b, kernel)        # pilot-bandwidth kernel weights
    active <- (w_h > 0) | (w_b > 0)        # active set = main and pilot windows
    if (sum(active) <= q + 1L)
      stop("too few observations in the bias-correction window on one side; ",
           "widen the bandwidth.")
    xs <- xs[active]; ys <- ys[active]; ids <- ids[active]
    w_h <- w_h[active]; w_b <- w_b[active]

    # design matrices: order q (for bias), order p nested inside
    Rq <- outer(xs - c, 0:q, `^`)          # n_s x (q+1)
    Rp <- Rq[, 1:(p + 1L), drop = FALSE]

    invG_p <- .qrXXinv(sqrt(w_h) * Rp)     # (X' A(h) X)^{-1}
    invG_q <- .qrXXinv(sqrt(w_b) * Rq)     # (X' A(b) X)^{-1}, order q

    # bias-correction weight matrix of eq:bc_Q (without the 1/h, 1/n scalings,
    # which cancel): Q = X'A(h) - h^{p+1} * theta * e_{p+1,q}' Gq^{-1} X'A(b)
    e_p1 <- numeric(q + 1L); e_p1[p + 2L] <- 1            # e_{p+1,q}
    theta <- crossprod(Rp * w_h, ((xs - c) / h)^(p + 1L)) # X' A(h) S_{p+1}
    Aq_b <- t(Rq * w_b)                                   # X' A(b)
    Qmat <- t(Rp * w_h) - h^(p + 1L) * (theta %*% (t(e_p1) %*% invG_q %*% Aq_b))

    # coefficients
    beta_p <- invG_p %*% crossprod(Rp * w_h, ys)         # conventional
    beta_q <- invG_q %*% crossprod(Rq * w_b, ys)         # order-q (for bc residuals)
    beta_bc <- invG_p %*% (Qmat %*% ys)                  # bias-corrected

    # intercept influence rows (e_0' M), then g = influence * residual
    a_c  <- as.numeric(invG_p[1, ] %*% t(Rp * w_h))      # conventional intercept
    a_bc <- as.numeric(invG_p[1, ] %*% Qmat)             # bias-corrected intercept
    res_c <- sqrt(length(ys) / (length(ys) - (p + 1L))) * (ys - Rp %*% beta_p)
    res_b <- sqrt(length(ys) / (length(ys) - (q + 1L))) * (ys - Rq %*% beta_q)

    list(
      beta0    = beta_p[1L],
      beta0_bc = beta_bc[1L],
      slope    = as.numeric(beta_p[2L]),
      id       = ids,
      g        = a_c * as.numeric(res_c),
      g_bc     = a_bc * as.numeric(res_b),
      # influence on (conventional - bias-corrected): both are linear in the same
      # Y, so the difference has per-unit weight (a_c - a_bc); used to estimate
      # Var(B-hat) for bandwidth regularization.
      g_diff   = (a_c - a_bc) * as.numeric(res_c)
    )
  }

  R_side <- side_fit(x >= c)   # (+)
  L_side <- side_fit(x <  c)   # (-)

  D    <- R_side$beta0    - L_side$beta0
  D_bc <- R_side$beta0_bc - L_side$beta0_bc
  V_D    <- sum(R_side$g^2)    + sum(L_side$g^2)
  V_D_bc <- sum(R_side$g_bc^2) + sum(L_side$g_bc^2)

  # plug-in asymptotic constants (Appendix B.3, eq:per-period-orders):
  #   V(D_t(h)) = V_t / (n_t h)               ->  V-hat_t = n h V_D
  #   B_t(h)    = h^{p+1}/(p+1)! * B_t         ->  B-hat_t = (p+1)! (D - D_bc) / h^{p+1},
  # since D - D_bc is the estimated bias B-hat_t(h, b) of Appendix B.1 (with the
  # finite-sample B_{t,(side),p}(h) in place of the kernel constant nu_{(side),p}).
  v_const <- n * h * V_D
  b_const <- factorial(p + 1L) * (D - D_bc) / h^(p + 1L)

  structure(
    list(D = D, D_bc = D_bc, V_D = V_D, V_D_bc = V_D_bc,
         b_const = b_const, v_const = v_const,
         n = n, h = h, b = b, c = c, p = p, q = q, kernel = kernel,
         sides = list(`+` = R_side, `-` = L_side)),
    class = "rd_period")
}

#' @export
print.rd_period <- function(x, ...) {
  cat(sprintf("Single-period RD (p=%d, h=%.4g, b=%.4g, kernel=%s, n=%d)\n",
              x$p, x$h, x$b, x$kernel, x$n))
  cat(sprintf("  D (conventional)   = %+.5g  (se %.4g)\n", x$D, sqrt(x$V_D)))
  cat(sprintf("  D (bias-corrected) = %+.5g  (se %.4g)\n", x$D_bc, sqrt(x$V_D_bc)))
  invisible(x)
}
