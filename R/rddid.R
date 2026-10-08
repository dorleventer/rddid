#' Period coefficients from a trend / weighting scheme
#' @keywords internal
#' @noRd
.rddid_weights <- function(weights, comps, t_rd) {
  m <- length(comps)
  if (is.numeric(weights)) {
    if (length(weights) != m)
      stop("numeric `weights` must have one entry per comparison period.")
    w <- weights
    if (abs(sum(w) - 1) > 1e-8)
      warning("comparison weights do not sum to 1 (constant-trend admissibility).")
  } else {
    weights <- match.arg(weights, c("constant", "linear"))
    if (weights == "constant") {
      w <- rep(1 / m, m)                 # equal weights: constant confounding trend
    } else {
      if (m < 2) stop("linear weights need at least 2 comparison periods.")
      X <- cbind(1, comps)               # extrapolate the line through the D_t to t_rd
      w <- as.numeric(c(1, t_rd) %*% solve(crossprod(X), t(X)))
    }
  }
  stats::setNames(w, as.character(comps))
}

#' Classify the sampling scheme from a long (period, id, side) table
#'
#' No repeated unit across periods → `"cs"`; repeated units that switch side at
#' least once → `"pv"`; repeated units that never switch → `"pc"`.
#' @keywords internal
#' @noRd
.scheme_from_long <- function(long) {
  rep_ids <- names(which(table(unique(long[, c("period", "id")])$id) >= 2L))
  if (length(rep_ids) == 0L) return("cs")
  sub <- long[long$id %in% rep_ids, , drop = FALSE]
  switches <- tapply(sub$side, sub$id, function(s) length(unique(s)) > 1L)
  if (any(switches)) "pv" else "pc"
}

#' Detect the sampling scheme from the id / side structure
#' @keywords internal
#' @noRd
.detect_scheme <- function(plist, c = 0) {
  # side is 1 if x >= c (treated at the cutoff), 0 otherwise — no third "side"
  # at exactly x == c, so a unit sitting on the cutoff does not read as a switch.
  long <- do.call(rbind, lapply(names(plist), function(k)
    data.frame(period = k, id = plist[[k]]$id,
               side = as.integer(plist[[k]]$x >= c))))
  .scheme_from_long(long)
}

#' RD-DID estimation and inference
#'
#' Estimates a treatment effect in a regression-discontinuity
#' difference-in-discontinuities design. The period-\eqn{t_{\mathrm{RD}}}
#' discontinuity is contaminated by a confounding policy that switches at the
#' same cutoff; comparison periods, where the confounding is present but the
#' treatment of interest is uniform at the cutoff, identify and net out that
#' confounding. The estimator is
#' \eqn{\widehat{\mathrm{ATT}} = \widehat D_{t_{\mathrm{RD}}} - \sum_t w_t \widehat D_t},
#' with each \eqn{\widehat D_t} a standard local-linear RD. This estimates the
#' ATT, or the ATU when the comparison periods are uniformly treated
#' (`estimand = "atu"`).
#'
#' ## Targeting the ATU
#'
#' When the treatment of interest is uniformly present (equal to one) in the
#' comparison periods rather than uniformly absent, the same difference of
#' discontinuities identifies the ATU: Leventer and Nevo, Section 6, show that
#' the ATU design is the ATT design with the sides of the cutoff exchanged
#' (mirror \eqn{\tilde R = c - R} and apply the ATT procedure unchanged). The
#' point estimate, standard errors, and bandwidth rules are numerically the
#' same either way, so `estimand` only labels the output here; the argument is
#' also passed through to the validation tests (`?rd_typecont`, `?rd_homog`,
#' `?rd_trendcell`, `?rd_compstable`), where only [rd_compstable()] computes
#' differently.
#'
#' @param data a long data frame, one row per unit-period (repeated
#'   cross-section or panel; a panel need not be balanced).
#' @param y,x,time column names (strings) for the outcome, running variable, and
#'   period.
#' @param id column name for the unit id; `NULL` (default) treats every row as a
#'   distinct unit (repeated cross-section).
#' @param t_rd the value of `time` identifying the RD (treated-at-cutoff) period.
#' @param comparisons values of `time` to use as comparison periods: periods
#'   in which the treatment of interest does not switch at the cutoff. `NULL`
#'   (default) uses every other period present, so with more than one RD
#'   period pass `comparisons` explicitly.
#' @param trend how the confounding jump is assumed to move over time, which
#'   fixes the comparison-period weights: `"constant"` (default; the jump is
#'   the same in every period, so the comparison discontinuities get equal
#'   weights), `"linear"` (the jump moves linearly in time; a line through the
#'   comparison discontinuities is extrapolated to `t_rd`, which needs at least
#'   two comparison periods), or a numeric vector of weights, one per entry of
#'   `comparisons` in the order given.
#' @param weights the old name of `trend`; still accepted, with a message.
#' @param estimand `"att"` (default) when the treatment of interest is
#'   uniformly ZERO in the comparison periods (targets the ATT), `"atu"` when
#'   it is uniformly ONE (targets the ATU). See "Targeting the ATU" below.
#' @param bwselect `"joint"` (default; one common bandwidth for every period,
#'   chosen to minimise the asymptotic MSE of the aggregate estimator — the
#'   paper's common rule), `"cct"` (a separate CCT MSE-optimal bandwidth per
#'   period, [rd_bw_cct()]), or `"iter"` (period-specific bandwidths found by
#'   coordinate descent on the aggregate AMSE, started at the common
#'   bandwidth; not used in the paper, kept for the simulations). Ignored if
#'   `h` is supplied.
#'   Both joint rules estimate each period's bias and variance constants at
#'   that period's own CCT pilot, so neither depends on which period is
#'   labelled `t_rd`; under `"joint"` the pilot `b_t` keeps each period's CCT
#'   ratio `b_t^CCT / h_t^CCT` (Appendix B.4 of the paper).
#' @param start starting point of the iterative (`bwselect = "iter"`) coordinate
#'   descent. `"hstar"` (default) starts all periods at the common joint-optimal
#'   h*; `"cct"` starts each period at its own CCT h; or supply a named
#'   numeric vector/list with one entry per period. Ignored unless
#'   `bwselect = "iter"`.
#' @param regularize logical; if `TRUE` (default) the joint AMSE-optimal
#'   bandwidth adds an `rdrobust`-style regularization term to the squared bias
#'   so a near-zero estimated curvature cannot blow the bandwidth up. Ignored
#'   if `h` is supplied.
#' @param reg_const regularization constant for `regularize` (default 3,
#'   matching the CCT convention).
#' @param h,b optional common point / pilot bandwidths; if `h` is given it is
#'   used for every period (with `b` defaulting to `h`).
#' @param scheme sampling scheme: `"cs"` (repeated cross-section), `"pc"`
#'   (panel, time-constant running variable), `"pv"` (panel, time-varying
#'   running variable), or `"auto"` (default), which reads it off the data: no
#'   `id` gives `"cs"`; units observed in more than one period, each on the
#'   same side of the cutoff in every period, give `"pc"`; any unit on
#'   different sides in different periods gives `"pv"`. The scheme selects
#'   which standard error and CI are printed; all three are always returned.
#' @param c cutoff (default 0).
#' @param p,q point / bias-correction polynomial orders (default 1, 2).
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param level confidence level (default 0.95).
#'
#' @return An object of class `"rddid"`, a list with:
#'   \describe{
#'     \item{`estimates`}{data frame with rows `Conventional` (the local-linear
#'       estimate with its conventional standard error) and `Robust` (the
#'       bias-corrected estimate with its robust standard error) and columns
#'       `est`, `se`, `ci_l`, `ci_u`, `z`, `p` (at `scheme`), and `se_cs`,
#'       `se_pc`, `se_pv` (the robust standard error under each scheme).}
#'     \item{`scheme`}{the sampling scheme used; `scheme_requested` is the
#'       argument as passed.}
#'     \item{`weights`, `weights_type`}{the comparison-period weights and
#'       their kind.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`bandwidth`}{list with `method` (the `bwselect` value, or
#'       `"fixed"`), `h_by_period` and `b_by_period` (named numeric vectors: the
#'       main and pilot bandwidth used in each period, whatever the rule),
#'       `niter` for `"iter"`, and, when the rule produces one common value,
#'       `h` and `b` (`"fixed"`, `"joint"`).}
#'     \item{`fits`}{named list of [rd_period()] objects by period, the RD
#'       period first (index by name).}
#'     \item{`t_rd`, `comparisons`, `level`, `c`, `p`, `q`, `kernel`, `call`}{as passed.}
#'     \item{`coef`}{the coefficient of each period's discontinuity in the
#'       aggregate: +1 for `t_rd`, minus the weight for each comparison period.}
#'     \item{`n_by_period`}{rows used in each period.}
#'     \item{`scheme_detected`}{the scheme read off the data, whatever was requested.}
#'   }
#' @export
rddid <- function(data, y, x, time, id = NULL, t_rd,
                  comparisons = NULL, trend = "constant",
                  estimand = c("att", "atu"),
                  bwselect = c("joint", "iter", "cct"), h = NULL, b = NULL,
                  scheme = c("auto", "cs", "pc", "pv"),
                  c = 0, p = 1L, q = 2L, kernel = "triangular", level = 0.95,
                  start = "hstar", regularize = TRUE, reg_const = 3,
                  weights = NULL) {
  cl       <- match.call()
  bwselect <- match.arg(bwselect)
  scheme   <- match.arg(scheme)
  estimand <- match.arg(estimand)
  kernel   <- match.arg(kernel, c("triangular", "epanechnikov", "uniform"))
  if (!is.null(weights)) {                 # `weights` is the old name of `trend`
    if (!missing(trend))
      stop("supply either `trend` or its old name `weights`, not both.")
    message("rddid(): `weights` is now called `trend`; the old name still works.")
    trend <- weights
  }
  for (nm in c(y, x, time, id)) if (!nm %in% names(data))
    stop("column '", nm, "' not found in `data`.")

  tt <- data[[time]]
  if (!t_rd %in% tt) stop("t_rd = ", t_rd, " not present in `", time, "`.")
  if (is.null(comparisons)) comparisons <- sort(setdiff(unique(tt), t_rd))
  if (length(comparisons) < 1L) stop("need at least one comparison period.")
  if (t_rd %in% comparisons) stop("`t_rd` must not be one of the `comparisons`.")
  if (is.null(id) && scheme == "auto")
    message("rddid(): no `id` given, so every row is treated as a different unit ",
            "(repeated cross-section standard errors).")
  periods <- c(t_rd, comparisons)

  ii <- if (is.null(id)) seq_len(nrow(data)) else data[[id]]
  plist <- stats::setNames(lapply(periods, function(tv) {
    rows <- which(tt == tv)
    data.frame(y = data[[y]][rows], x = data[[x]][rows], id = ii[rows])
  }), as.character(periods))

  w    <- .rddid_weights(trend, comparisons, t_rd)
  coef <- c(stats::setNames(1, as.character(t_rd)), -w)

  detected <- .detect_scheme(plist, c = c)
  use_scheme <- if (scheme == "auto") detected else scheme
  if (scheme %in% c("pc", "pv") && detected == "cs")
    warning("scheme = \"", scheme, "\" requested but no unit id repeats across periods; ",
            "all cross-period covariances are zero, so the reported SE equals the CS one.")

  # ---- bandwidth ----
  if (!is.null(h)) {
    if (is.null(b)) b <- h
    bws <- stats::setNames(rep(list(c(h = h, b = b)), length(periods)),
                           as.character(periods))
    bw_info <- list(method = "fixed", h = h, b = b)
  } else if (bwselect == "cct") {
    bws <- .bw_cct(plist, c = c, p = p, kernel = kernel)
    bw_info <- list(method = "cct", bws = bws)
  } else if (bwselect == "iter") {
    ib <- .bw_joint_iter(plist, coef, as.character(t_rd), scheme = use_scheme,
                         start = start,
                         c = c, p = p, q = q, kernel = kernel,
                         regularize = regularize, reg_const = reg_const)
    bws <- ib$bws
    bw_info <- list(method = "iter", bws = bws, niter = ib$niter)
  } else {
    jb <- .bw_joint(plist, coef, as.character(t_rd), scheme = use_scheme,
                    c = c, p = p, q = q, kernel = kernel,
                    regularize = regularize, reg_const = reg_const)
    # one common h*; the pilot b keeps each period's own CCT ratio (B.4)
    bws <- stats::setNames(lapply(as.character(periods), function(k)
      c(h = jb$h, b = unname(jb$b[[k]]))), as.character(periods))
    bw_info <- list(method = "joint", h = jb$h, b = jb$b, B = jb$B,
                    Veff = jb$Veff, reg = jb$reg, pilot_bws = jb$pilot_bws)
  }
  # the bandwidth actually used in each period, whatever the rule
  bw_info$h_by_period <- vapply(bws, function(v) unname(v[["h"]]), numeric(1))
  bw_info$b_by_period <- vapply(bws, function(v) unname(v[["b"]]), numeric(1))

  # ---- per-period fits at chosen bandwidth(s) ----
  fits <- stats::setNames(lapply(as.character(periods), function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = bws[[k]]["h"], b = bws[[k]]["b"],
              id = plist[[k]]$id, c = c, p = p, q = q, kernel = kernel)),
    as.character(periods))

  ac  <- .aggregate_fits(fits, coef, bc = FALSE)
  abc <- .aggregate_fits(fits, coef, bc = TRUE)

  sefld <- c(cs = "V_cs", pc = "V_pc", pv = "V_pv")
  zc <- stats::qnorm(1 - (1 - level) / 2)
  mkrow <- function(agg) {
    se <- unname(sqrt(agg[sefld[use_scheme]]))
    z  <- agg[["est"]] / se
    c(est = agg[["est"]],
      se = se,
      se_cs = sqrt(agg[["V_cs"]]), se_pc = sqrt(agg[["V_pc"]]), se_pv = sqrt(agg[["V_pv"]]),
      ci_l = agg[["est"]] - zc * se, ci_u = agg[["est"]] + zc * se,
      z = z, p = 2 * stats::pnorm(-abs(z)))
  }
  est_tab <- rbind(Conventional = mkrow(ac), Robust = mkrow(abc))

  structure(list(
    estimates = as.data.frame(est_tab),
    coef = coef, weights = w, weights_type = if (is.numeric(trend)) "custom" else trend,
    estimand = estimand,
    t_rd = t_rd, comparisons = comparisons,
    scheme = use_scheme, scheme_detected = detected, scheme_requested = scheme,
    bandwidth = bw_info, fits = fits, level = level,
    p = p, q = q, kernel = kernel, c = c,
    n_by_period = vapply(fits, function(f) f$n, numeric(1)),
    call = cl
  ), class = "rddid")
}

#' @export
print.rddid <- function(x, digits = 4, ...) {
  est <- if (is.null(x$estimand)) "att" else x$estimand
  cat(sprintf("RD-DID estimate of the %s in period %s\n", toupper(est), x$t_rd))
  cat(sprintf("  Comparison periods: %s   (%s; weights %s)\n",
              paste(x$comparisons, collapse = ", "), .trend_label(x$weights_type),
              paste(trimws(formatC(x$weights, digits = 3, format = "g")), collapse = ", ")))
  cat(sprintf("  Sampling scheme: %s%s\n", .scheme_label(x$scheme),
              if (identical(x$scheme_requested, "auto")) " (detected from the data)" else ""))
  cat(sprintf("  Bandwidth: %s\n\n", .bandwidth_label(x$bandwidth)))
  .print_estimates(x, digits = digits)
  cat("\n  summary() shows the per-period fits and the s.e. under every sampling scheme.\n")
  invisible(x)
}

# ---- print helpers shared by print.rddid / summary.rddid ------------------------------------
.scheme_label <- function(s) {
  c(cs = "repeated cross-section",
    pc = "panel, running variable fixed over time",
    pv = "panel, running variable varies over time")[[s]]
}
.trend_label <- function(wt) {
  switch(wt, constant = "constant confounding trend", linear = "linear confounding trend",
         custom = "user-supplied weights", wt)
}
.bandwidth_label <- function(bw) {
  by_t <- function(v) paste(trimws(formatC(v, digits = 4, format = "g")), collapse = ", ")
  switch(bw$method,
    fixed = sprintf("h = %.4g in every period (fixed), pilot b = %.4g", bw$h, bw$b),
    joint = sprintf("common h = %.4g (rule \"joint\", AMSE-optimal for the aggregate)\n  Pilot bandwidth b by period: %s",
                    bw$h, by_t(bw$b_by_period)),
    cct   = sprintf("per-period CCT MSE-optimal (rule \"cct\"): h by period %s", by_t(bw$h_by_period)),
    iter  = sprintf("period-specific (rule \"iter\", %d iterations): h by period %s", bw$niter,
                    by_t(bw$h_by_period)))
}
.print_estimates <- function(x, digits = 4) {
  e <- x$estimates
  lab <- c(Conventional = "Conventional", Robust = "Robust (bias-corrected)")
  fmt <- function(v) formatC(v, digits = digits, format = "f")
  pfmt <- function(p) ifelse(p < 1e-3, "<0.001", formatC(p, digits = 3, format = "f"))
  cat(sprintf("  %-24s %10s %10s %7s %8s   %s%% CI\n", "", "Estimate", "Std. err.", "z",
              "p-value", format(100 * x$level)))
  for (r in rownames(e))
    cat(sprintf("  %-24s %10s %10s %7.2f %8s   [%s, %s]\n", lab[[r]], fmt(e[r, "est"]),
                fmt(e[r, "se"]), e[r, "z"], pfmt(e[r, "p"]), fmt(e[r, "ci_l"]), fmt(e[r, "ci_u"])))
}
