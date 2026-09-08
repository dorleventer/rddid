# Kernel constants of Appendix B.3 ("Single-period asymptotics" and
# "Additional definitions") of the paper (Leventer & Nevo, RD-DID).
#
# Implements, for r_p(u) = [1, u, ..., u^p]' and a kernel K supported on
# [-1, 1], with side "+" integrating over [0, 1] and side "-" over [-1, 0]:
#
#   Gamma~_(side),p       = int K(u) r_p(u) r_p(u)' du
#   Psi~_(side),p         = int K(u)^2 r_p(u) r_p(u)' du
#   vartheta~_(side),p    = int K(u) u^(p+1) r_p(u) du
#   Omega~_(side),p(rho)  = int K(v) K(rho v) r_p(v) r_p(rho v)' dv
#
#   b_(side),p      = e0' Gamma~^{-1} vartheta~            (eq:per-period-orders)
#   v_(side),p      = e0' Gamma~^{-1} Psi~ Gamma~^{-1} e0  (eq:per-period-orders)
#   c_(side),p(rho) = e0' Gamma~^{-1} Omega~(rho) Gamma~^{-1} e0  (lem:cov-pc)
#
# Names: the paper writes these three constants as nu_(side),p (bias),
# omega_(side),p (variance) and omega_(side),p(rho) (cross-period, with
# omega(1) = omega); the package keeps the mnemonics .kc_b, .kc_v, .kc_c.
#
# Every one of these is a PURE function of (K, p, rho) -- no data, no
# randomness -- so each is computed once by quadrature and memoised (see
# .kc_cache below), never refit inside a bootstrap or coordinate-descent loop.
#
# Ported from the validated base-R module
# rd-did/code/simulations/appb_verify/R/kernel_constants.R (Module K of the
# Appendix B verification study), renamed to the package's dot-prefixed
# internal convention. K itself is obtained from this package's own
# .rd_kweight() (R/kernels.R) rather than a re-typed copy, so the three
# kernel shapes used here can never drift from the estimator's.
#
# Every object function below carries a quadrature `rule` ("simpson", the
# default composite Simpson rule at n = 4000 subintervals, or "gl" for
# 400-node Gauss-Legendre); the test suite certifies Simpson against GL
# wherever no closed form is convenient to hand-derive. The quadrature is
# vectorised (nodes x weights x basis matrices), so a fresh c(rho) costs well
# under a millisecond and can be evaluated inside optimisation loops.

## ---------------------------------------------------------------- kernel ----

## Canonical kernel name (case-insensitive, partial match) -- the same
## resolution .rd_kweight() itself performs -- used both to build cache keys
## and to fetch the kernel function below.
.kc_kernel_name <- function(kernel) {
  match.arg(tolower(kernel), c("triangular", "epanechnikov", "uniform"))
}

#' Kernel function K(u) for the kernel-constants module
#'
#' Returns the weight function `K(u)` used throughout this module, obtained by
#' currying `.rd_kweight()` (`R/kernels.R`) on `kernel` -- so the shape used
#' here is *exactly* the package's own kernel (triangular, Epanechnikov, or
#' uniform), never a second, independently-typed copy of those formulas.
#'
#' @param kernel one of `"triangular"` (default), `"epanechnikov"`, `"uniform"`.
#' @return a function `K(u)` returning the kernel weight, vectorised over `u`,
#'   zero outside `[-1, 1]`.
#' @keywords internal
#' @noRd
.kc_kernel <- function(kernel = "triangular") {
  kernel <- .kc_kernel_name(kernel)
  function(u) .rd_kweight(u, kernel)
}

## Polynomial basis r_p(u) = [1, u, ..., u^p]'. Returns length(u) x (p+1).
## p >= 0 is supported (p = 0 is used by no estimator in this package, which
## requires p >= 1, but the kernel constants are well defined there too).
.kc_r_p <- function(u, p) {
  stopifnot(length(p) == 1L, p >= 0, p == as.integer(p))
  outer(as.numeric(u), 0:as.integer(p), "^")
}

## side_idx("+") = 1, side_idx("-") = 2. Accepts "+"/"-" (vectorised) or
## numeric 1/2; a pasted Unicode minus (U+2212) maps to "-".
.kc_side_idx <- function(side) {
  if (is.numeric(side)) {
    if (!all(side %in% c(1, 2))) stop("numeric side must be 1 (\"+\") or 2 (\"-\")")
    return(as.integer(side))
  }
  s   <- as.character(side)
  out <- match(s, c("+", "-"))
  bad <- is.na(out)
  if (any(bad)) {
    mn <- as.raw(c(0xe2, 0x88, 0x92))
    out[bad] <- ifelse(vapply(s[bad], function(z) identical(charToRaw(z), mn),
                              logical(1)), 2L, NA_integer_)
  }
  if (anyNA(out)) stop("side must be \"+\" or \"-\", got: ", paste(s, collapse = ", "))
  out
}

## Integration limits for a side, assuming supp(K) is contained in [-1, 1].
.kc_side_limits <- function(side) if (.kc_side_idx(side) == 1L) c(0, 1) else c(-1, 0)

## e_{0,p}: first unit vector of length p+1.
.kc_e0 <- function(p) {
  e <- numeric(p + 1L)
  e[1L] <- 1
  e
}

## ------------------------------------------------------------ quadrature ----
## Both rules are used in VECTORISED form: a rule returns its nodes x and weights
## w on [a, b], and every integral below is a weighted matrix product over the
## nodes (sum_k w_k g(x_k) r_p(x_k) r_p(y_k)' = crossprod(R1 * (w g), R2)). This
## is what makes c(rho) cheap enough to sit inside the coordinate-descent loop of
## .bw_joint_iter(), where rho = h_t/h_s changes at every objective evaluation.

## Gauss-Legendre nodes/weights on [-1, 1] by the Golub-Welsch construction:
## eigen-decomposition of the symmetric tridiagonal Jacobi matrix of the
## Legendre recursion (off-diagonal beta_k = k / sqrt(4k^2 - 1)). Nodes are
## the eigenvalues, weights 2 * (first eigenvector component)^2. Cached by n.
.kc_gl_cache <- new.env(parent = emptyenv())

.kc_gl_rule <- function(n) {
  n <- as.integer(n)
  stopifnot(n >= 1L)
  key <- paste0("n", n)
  hit <- .kc_gl_cache[[key]]
  if (!is.null(hit)) return(hit)
  J <- matrix(0, n, n)
  if (n > 1L) {
    k <- seq_len(n - 1L)
    beta <- k / sqrt(4 * k^2 - 1)
    J[cbind(k, k + 1L)] <- beta
    J[cbind(k + 1L, k)] <- beta
  }
  eg <- eigen(J, symmetric = TRUE)
  nodes <- eg$values
  wts <- 2 * eg$vectors[1L, ]^2
  o <- order(nodes)
  out <- list(x = nodes[o], w = wts[o])
  if (abs(sum(out$w) - 2) > 1e-10)
    stop("Gauss-Legendre weights do not sum to 2 (n = ", n, ")")
  .kc_gl_cache[[key]] <- out
  out
}

#' Quadrature nodes and weights on `[a, b]`
#'
#' `rule = "simpson"`: composite Simpson with `n = 4000` subintervals (the
#' default everywhere). `rule = "gl"`: `n = 400` Gauss-Legendre nodes, an
#' independent rule the test suite certifies Simpson against.
#'
#' @param rule `"simpson"` or `"gl"`.
#' @param a,b integration limits (`b > a`).
#' @return list with numeric vectors `x` (nodes) and `w` (weights).
#' @keywords internal
#' @noRd
.kc_nodes <- function(rule, a, b) {
  stopifnot(is.finite(a), is.finite(b), b > a)
  if (rule == "simpson") {
    n <- 4000L
    x <- seq.int(a, b, length.out = n + 1L)
    w <- rep(2, n + 1L)
    w[seq.int(2L, n, by = 2L)] <- 4        # odd-indexed interior nodes
    w[c(1L, n + 1L)] <- 1
    list(x = x, w = w * (b - a) / (3 * n))
  } else {
    nw <- .kc_gl_rule(400L)
    list(x = (b - a) / 2 * nw$x + (a + b) / 2, w = (b - a) / 2 * nw$w)
  }
}

## sum_k w_k g_k r_p(x_k) r_p(y_k)'  with y = rho * x  (rho = 1: r_p(x) r_p(x)').
.kc_int_rr <- function(nodes, p, g, rho = 1) {
  R1 <- .kc_r_p(nodes$x, p)
  R2 <- if (rho == 1) R1 else .kc_r_p(rho * nodes$x, p)
  crossprod(R1 * (nodes$w * g), R2)
}

## ------------------------------------------------------------- memoisation --
## The kernel objects are constants; recomputing them inside a bootstrap or
## coordinate-descent loop is pure waste. Keyed by (object type, kernel name,
## p, canonical side index, ..., rule) -- `rule` is always part of the key, so
## Simpson and Gauss-Legendre results are cached under distinct entries and a
## Simpson-vs-GL comparison is a genuine agreement check between two
## independently computed numbers, never two reads of the same cached value.
## Unlike the source module (which cached only the built-in triangular kernel,
## keyed on function identity), every kernel is cached here, keyed on its
## canonical name string -- robust to `.kc_kernel()` returning a fresh closure
## on every call.
.kc_cache <- new.env(parent = emptyenv())

.kc_key <- function(...) paste(..., sep = "|")

.kc_memo <- function(key, compute) {
  hit <- .kc_cache[[key]]
  if (!is.null(hit)) return(hit)
  val <- compute()
  .kc_cache[[key]] <- val
  val
}

## -------------------------------------------------------- kernel objects ----

#' Gamma~ kernel matrix: int K(u) r_p(u) r_p(u)' du over one side
#'
#' `Gamma~_(side),p = \int K(u) r_p(u) r_p(u)' du`, `r_p(u) = [1, u, ..., u^p]'`,
#' integrated over `[0, 1]` (side `"+"`) or `[-1, 0]` (side `"-"`). Appendix
#' B.3, "Single-period asymptotics".
#'
#' @param p polynomial order (`p >= 0`; the package's own estimator requires
#'   `p >= 1`, but this module also supports `p = 0`).
#' @param side `"+"` (integrate over `[0, 1]`) or `"-"` (`[-1, 0]`).
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`
#'   (Gauss-Legendre); the two agree to high precision and exist so the test
#'   suite can certify one against the other.
#' @return a `(p+1) x (p+1)` numeric matrix.
#' @keywords internal
#' @noRd
.kc_gamma <- function(p, side, kernel = "triangular", rule = c("simpson", "gl")) {
  rule   <- match.arg(rule)
  kernel <- .kc_kernel_name(kernel)
  key <- .kc_key("gamma", kernel, p, .kc_side_idx(side), rule)
  .kc_memo(key, function() {
    K   <- .kc_kernel(kernel)
    lim <- .kc_side_limits(side)
    nd  <- .kc_nodes(rule, lim[1L], lim[2L])
    .kc_int_rr(nd, p, K(nd$x))                 # sum w K(u) r_p(u) r_p(u)'
  })
}

#' Psi~ kernel matrix: int K(u)^2 r_p(u) r_p(u)' du over one side
#'
#' `Psi~_(side),p = \int K(u)^2 r_p(u) r_p(u)' du`. Appendix B.3.
#'
#' @param p polynomial order (`p >= 0`).
#' @param side `"+"` or `"-"`.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`.
#' @return a `(p+1) x (p+1)` numeric matrix.
#' @keywords internal
#' @noRd
.kc_psi <- function(p, side, kernel = "triangular", rule = c("simpson", "gl")) {
  rule   <- match.arg(rule)
  kernel <- .kc_kernel_name(kernel)
  key <- .kc_key("psi", kernel, p, .kc_side_idx(side), rule)
  .kc_memo(key, function() {
    K   <- .kc_kernel(kernel)
    lim <- .kc_side_limits(side)
    nd  <- .kc_nodes(rule, lim[1L], lim[2L])
    .kc_int_rr(nd, p, K(nd$x)^2)               # sum w K(u)^2 r_p(u) r_p(u)'
  })
}

#' vartheta~ kernel vector: int K(u) u^order r_p(u) du over one side
#'
#' `vartheta~_(side),p = \int K(u) u^order r_p(u) du`; the paper's leading-bias
#' term uses `order = p + 1`. Appendix B.3.
#'
#' @param p polynomial order (`p >= 0`).
#' @param side `"+"` or `"-"`.
#' @param order power of `u` in the integrand (default `p + 1`).
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`.
#' @return a numeric vector of length `p + 1`.
#' @keywords internal
#' @noRd
.kc_theta <- function(p, side, order = p + 1, kernel = "triangular",
                      rule = c("simpson", "gl")) {
  rule   <- match.arg(rule)
  kernel <- .kc_kernel_name(kernel)
  key <- .kc_key("theta", kernel, p, .kc_side_idx(side), order, rule)
  .kc_memo(key, function() {
    K   <- .kc_kernel(kernel)
    lim <- .kc_side_limits(side)
    nd  <- .kc_nodes(rule, lim[1L], lim[2L])
    as.numeric(crossprod(.kc_r_p(nd$x, p), nd$w * K(nd$x) * nd$x^order))
  })
}

#' Omega~ kernel matrix: int K(v) K(rho v) r_p(v) r_p(rho v)' dv over one side
#'
#' `Omega~_(side),p(rho) = \int K(v) K(rho v) r_p(v) r_p(rho v)' dv`. Appendix
#' B.3; feeds the cross-period covariance constant `c(rho)` (`lem:cov-pc`).
#' For `rho > 1` the factor `K(rho v)` leaves the kernel support at
#' `|v| = 1/rho`, a kink in the integrand, so the integration interval is
#' split there -- both quadrature rules then see a smooth (in fact
#' polynomial) integrand on each piece. The integral is still taken over the
#' full side: the integrand is exactly zero beyond `1/rho`. The same split
#' applies for every kernel in this module (triangular, Epanechnikov,
#' uniform), since all have support exactly `[-1, 1]`.
#'
#' By construction, the break point is exactly where `|rho * v|` first
#' reaches 1, so on whichever side of it `K(rho v)` is analytically zero, it
#' is zero throughout that whole piece's *interior* -- that piece is skipped
#' rather than quadrature-integrated (each piece's classification is decided
#' once, from its midpoint). This matters only for the uniform kernel: it is
#' the sole kernel here with a *closed*, discontinuous support indicator
#' (`K(1) = 0.5`, not the continuous tapers of the triangular/Epanechnikov
#' shapes), so composite Simpson -- unlike Gauss-Legendre, which never
#' samples an interval's endpoints -- would otherwise evaluate the shared
#' kink node at its closed (in-support) value and assign it positive
#' quadrature weight inside the zero piece, a spurious O(1/n) contribution.
#' Skipping the zero piece sidesteps that without touching `.rd_kweight()`'s
#' own (correct) closed boundary convention, and is a no-op for the
#' triangular/Epanechnikov kernels, which already integrate to ~0 there.
#'
#' @param p polynomial order (`p >= 0`).
#' @param side `"+"` or `"-"`.
#' @param rho positive scalar frequency ratio; `Omega~(1) == Psi~`.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`.
#' @return a `(p+1) x (p+1)` numeric matrix.
#' @keywords internal
#' @noRd
.kc_omega <- function(p, side, rho, kernel = "triangular", rule = c("simpson", "gl")) {
  stopifnot(length(rho) == 1L, is.finite(rho), rho > 0)
  rule   <- match.arg(rule)
  kernel <- .kc_kernel_name(kernel)
  key <- .kc_key("omega", kernel, p, .kc_side_idx(side), format(rho, digits = 17), rule)
  .kc_memo(key, function() {
    K   <- .kc_kernel(kernel)
    brk <- if (.kc_side_idx(side) == 1L) {
      sort(unique(c(0, min(1, 1 / rho), 1)))
    } else {
      sort(unique(c(-1, max(-1, -1 / rho), 0)))
    }
    acc <- matrix(0, p + 1L, p + 1L)
    for (j in seq_len(length(brk) - 1L)) {
      a <- brk[j]; b <- brk[j + 1L]
      if (b <= a) next
      in_support <- abs(rho * ((a + b) / 2)) <= 1   # piece's own midpoint
      if (!in_support) next                        # K(rho v) == 0 on this piece
      nd  <- .kc_nodes(rule, a, b)
      acc <- acc + .kc_int_rr(nd, p, K(nd$x) * K(rho * nd$x), rho = rho)
    }
    acc
  })
}

## ------------------------------------------------------ kernel constants ----

#' Leading-bias kernel constant b_(side),p
#'
#' `b_(side),p = e0' Gamma~^{-1} vartheta~`. Appendix B.3, `eq:per-period-orders`.
#'
#' @param p polynomial order (`p >= 0`).
#' @param side `"+"` or `"-"`.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`.
#' @return a scalar.
#' @keywords internal
#' @noRd
.kc_b <- function(p, side, kernel = "triangular", rule = c("simpson", "gl")) {
  rule <- match.arg(rule)
  G  <- .kc_gamma(p, side, kernel = kernel, rule = rule)
  th <- .kc_theta(p, side, order = p + 1, kernel = kernel, rule = rule)
  as.numeric(solve(G, th)[1L])
}

#' Variance kernel constant v_(side),p
#'
#' `v_(side),p = e0' Gamma~^{-1} Psi~ Gamma~^{-1} e0`. Appendix B.3,
#' `eq:per-period-orders`.
#'
#' @param p polynomial order (`p >= 0`).
#' @param side `"+"` or `"-"`.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`.
#' @return a scalar.
#' @keywords internal
#' @noRd
.kc_v <- function(p, side, kernel = "triangular", rule = c("simpson", "gl")) {
  rule <- match.arg(rule)
  G  <- .kc_gamma(p, side, kernel = kernel, rule = rule)
  Ps <- .kc_psi(p, side, kernel = kernel, rule = rule)
  g <- solve(G, .kc_e0(p))
  as.numeric(crossprod(g, Ps %*% g))
}

#' Cross-period covariance kernel constant c_(side),p(rho)
#'
#' `c_(side),p(rho) = e0' Gamma~^{-1} Omega~(rho) Gamma~^{-1} e0`. Appendix
#' B.3, `lem:cov-pc`. Satisfies `c(1) = v` and `c(1/rho) = rho * c(rho)`.
#'
#' @param p polynomial order (`p >= 0`).
#' @param side `"+"` or `"-"`.
#' @param rho positive scalar frequency ratio.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param rule quadrature rule, `"simpson"` (default) or `"gl"`.
#' @return a scalar.
#' @keywords internal
#' @noRd
.kc_c <- function(p, side, rho, kernel = "triangular", rule = c("simpson", "gl")) {
  rule <- match.arg(rule)
  G  <- .kc_gamma(p, side, kernel = kernel, rule = rule)
  Om <- .kc_omega(p, side, rho, kernel = kernel, rule = rule)
  g <- solve(G, .kc_e0(p))
  as.numeric(crossprod(g, Om %*% g))
}
