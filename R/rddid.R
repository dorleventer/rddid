# rddid.R -- the main estimator rddid() and its print method. The pieces it calls live in:
# trend_weights.R (weights), sampling_scheme.R (scheme), bandwidth_*.R (bandwidth rules),
# rd_period.R (per-period fits), cross_period_covariance.R (aggregation and variances).

#' Estimate the effect of a treatment at a cutoff shared with a confounding policy
#'
#' A treatment of interest switches on at a cutoff of a running variable in one
#' period, the **RD period**. A **confounding policy** switches at the same
#' cutoff, in every period, so the jump in the outcome at the cutoff in the RD
#' period mixes the treatment effect with the **confounding jump**. In the
#' **comparison periods** the treatment of interest is uniform at the cutoff
#' (nobody treated, or everybody treated), so, provided the
#' treatment of interest has no anticipation or carry-over effects there (which the
#' paper assumes), the jump there *is* the confounding jump. `rddid()` estimates the jump in every period by
#' local-linear RD and subtracts a weighted average of the comparison-period
#' jumps from the RD-period jump. How the weights are set is the
#' **confounding-trend assumption**: constant (equal weights) or linear in time.
#'
#' @details
#' ## The estimate
#'
#' \deqn{\mathrm{ATT}(t_{\mathrm{RD}}) = D_{t_{\mathrm{RD}}} - \sum_t w_t D_t,}{ATT(t_RD) = D_{t_RD} - sum_t w_t D_t,}
#' where \eqn{D_t} is the jump in the outcome at the cutoff in period \eqn{t},
#' estimated by a local-linear RD on that period's observations;
#' \eqn{t_{\mathrm{RD}}}{t_RD} is the RD period (`t_rd`); the sum runs over the
#' comparison periods; and \eqn{w_t} is the weight of comparison period
#' \eqn{t}. `summary()` lists every \eqn{D_t} with its coefficient in the sum.
#'
#' ## The confounding-trend assumption
#'
#' With `trend = "constant"` the confounding jump is the same in every period,
#' so the comparison periods get equal weights, \eqn{w_t = 1/m} with \eqn{m}
#' comparison periods. With `trend = "linear"` the confounding jump moves
#' linearly in time; the weights then extrapolate the least-squares line through the
#' comparison-period jumps to the RD period (with two comparison periods, the
#' line through them). They sum to one and can be
#' negative: with comparison periods 1 and 2 and RD period 3 they are -1 and 2.
#' A numeric `trend` supplies the weights directly.
#'
#' ## Sampling schemes
#'
#' The scheme sets the standard error and, under `"joint"` and `"iter"`, also
#' the bandwidth, and with it the estimate, because those rules balance bias
#' against the variance under the scheme in use; with a fixed `h` or `"cct"` it
#' changes only the standard error. In a repeated
#' cross-section (`"cs"`) the periods' samples are independent and the variance
#' is the weighted sum of the period variances. In a panel the same units
#' appear in several periods, so the period jumps are correlated and the
#' variance adds their covariances, computed by matching units on `id`: under
#' `"pc"` from units on the same side of the cutoff in both periods; under
#' `"pv"` also from the units that change side, which enter with the opposite
#' sign. The standard errors under all
#' three schemes, at the bandwidths actually used, are in `estimates` and in
#' `summary()`.
#'
#' ## Bandwidth rules
#'
#' * `"joint"` (default), the common bandwidth: one `h` in every period, chosen
#'   to minimize the asymptotic mean squared error of the RD-DID estimate, not
#'   of each period's jump, so the biases of the period jumps can partly cancel.
#' * `"cct"`, per-period CCT bandwidths: each period gets its own MSE-optimal
#'   bandwidth from [rdrobust::rdbwselect()], as if it were a stand-alone RD
#'   (Calonico, Cattaneo and Titiunik, 2014).
#' * `"iter"`: a separate bandwidth in each period, chosen together to minimize
#'   the asymptotic mean squared error of the RD-DID estimate; not used in the
#'   paper, kept for simulations.
#' * A numeric `h` is used in every period, with pilot bandwidth `b` (default
#'   `h`).
#'
#' The bandwidths used in each period are in `bandwidth$h_by_period` and
#' `bandwidth$b_by_period`, and in `summary()`.
#'
#' ## Targeting the ATU
#'
#' When the treatment of interest is uniformly present in the comparison periods
#' (everybody at the cutoff is treated), the same difference of jumps
#' identifies the ATU: the effect for the units just below the cutoff, which
#' are untreated in the RD period. The paper shows that this design is the ATT
#' design with the two sides of the cutoff exchanged (the running variable
#' mirrored around the cutoff). Mirroring changes neither the jump estimates
#' nor their standard errors nor the bandwidth rules, so `estimand = "atu"`
#' returns the same numbers as `"att"` and only labels the output. The four
#' tests of the assumptions take the same argument; only [rd_compstable()]
#' computes differently under `"atu"`.
#'
#' ## Technical details
#'
#' The `"joint"` rule is AMSE-optimal: it minimizes the asymptotic mean squared
#' error (AMSE) of the RD-DID estimate over a common `h`. Each period is first
#' fitted at its own CCT bandwidths to estimate that period's bias constant and
#' variance constant; these are combined, with the coefficients of the sum above
#' (and, under `"pc"`, the covariances between periods), into the AMSE, which
#' is then minimized in closed form. Each period's pilot bandwidth keeps that
#' period's CCT ratio \eqn{b/h}{b/h}. Each period's bias and variance
#' constants are estimated at that period's own CCT pilot bandwidths, not at the
#' RD period's, so an aggregate of several RD periods gets the same bandwidth
#' whichever of them carries the `t_rd` label.
#' The `"iter"` rule minimizes the same objective over one bandwidth per period
#' by coordinate descent, starting from `start`. With `regularize = TRUE` both
#' rules add `reg_const` times the estimated variance of the bias constants to
#' the squared bias, as \pkg{rdrobust} does, so a near-zero estimated bias
#' cannot make the bandwidth very large. Standard errors use the HC1 convention
#' of \pkg{rdrobust}. The derivations are in the paper's appendix on estimation
#' and bandwidth choice.
#'
#' @param data a data frame in long format, one row per unit and period: a
#'   repeated cross-section (different units in each period) or a panel (the
#'   same units in several periods; it need not be balanced).
#' @param y name of the outcome column (a string).
#' @param x name of the running-variable column (a string).
#' @param time name of the period column (a string). The column is usually
#'   numeric (a year); character or factor labels work as well, except under
#'   `trend = "linear"` (in [rddid()] or [rd_trendcell()]), where the line is
#'   fitted on the period values: these then need to be numeric and are the
#'   time scale of the line (so 2015, 2017, 2018 are unequally spaced);
#'   `as.numeric()` converts character labels such as `"2019"`. In [rddid()]
#'   the default comparison periods are every period other than `t_rd`, in
#'   sorted order of the values (alphabetical for character or factor labels),
#'   the order a numeric `trend` follows.
#' @param id name of the unit-identifier column (a string), needed for panel
#'   standard errors. With `NULL` (default) every row is treated as a different
#'   unit, which gives repeated cross-section standard errors (with a message).
#' @param t_rd the RD period: the value of `time` in which the treatment of
#'   interest switches on at the cutoff.
#' @param comparisons the comparison periods: values of `time` in which the
#'   treatment of interest is uniform at the cutoff (nobody treated, or, with
#'   `estimand = "atu"`, everybody treated). `NULL` (default) uses every period
#'   other than `t_rd`; pass the periods explicitly when the data contain
#'   periods that are neither (a second RD period, say).
#' @param trend the confounding-trend assumption, which sets the
#'   comparison-period weights: `"constant"` (default; the confounding jump is
#'   the same in every period, equal weights) or `"linear"` (the confounding
#'   jump moves linearly in time; needs at least two comparison periods). A
#'   numeric vector gives the weights directly, one per entry of `comparisons`
#'   in that order. Weights that sum to one cancel a constant confounding
#'   jump; `rddid()` warns when they do not.
#' @param estimand `"att"` (default) when the comparison periods are uniformly
#'   untreated: the ATT, the effect of the treatment on the units just above the
#'   cutoff that are treated in the RD period. `"atu"` when the comparison
#'   periods are uniformly treated: the ATU, the effect for the units just below
#'   the cutoff, which are untreated in the RD period. The numbers are the same
#'   either way; see "Targeting the ATU".
#' @param bwselect the bandwidth rule, used when `h` is not given: `"joint"`
#'   (default; the common bandwidth, one `h` for every period, chosen to minimize
#'   the asymptotic mean squared error of the RD-DID estimate), `"cct"`
#'   (per-period CCT bandwidths: each period's own MSE-optimal bandwidth of
#'   Calonico, Cattaneo and Titiunik, from [rd_bw_cct()]), or `"iter"` (a
#'   separate bandwidth in each period, chosen together for the RD-DID
#'   estimate; not used in the paper, kept for simulations). See "Bandwidth
#'   rules".
#' @param h the main bandwidth (point estimate). If given, it is used in every
#'   period and `bwselect` is ignored.
#' @param b the pilot bandwidth (bias correction), used with a given `h`;
#'   defaults to `h`.
#' @param scheme the sampling scheme, which sets the standard error: `"cs"`
#'   (repeated cross-section: different units in each period), `"pc"` (panel,
#'   running variable fixed over time: the same units, each on the same side of
#'   the cutoff in every period), `"pv"` (panel, running variable varies over
#'   time: some units change side between periods), or `"auto"` (default), which
#'   reads it off the data: no repeated `id` gives `"cs"`, repeated units that
#'   never change side give `"pc"`, any unit that changes side gives `"pv"`.
#' @param c the cutoff (default 0). A unit with `x >= c` is above the cutoff.
#' @param p order of the local polynomial for the point estimate (default 1,
#'   local linear).
#' @param q order of the local polynomial for the bias correction (default 2);
#'   must exceed `p`.
#' @param kernel the kernel: `"triangular"` (default), `"epanechnikov"` or
#'   `"uniform"`.
#' @param level confidence level of the reported intervals (default 0.95).
#' @param start where the `"iter"` rule starts: `"hstar"` (default; the common
#'   bandwidth in every period), `"cct"` (each period's CCT bandwidth), or a
#'   named numeric vector or list with one starting bandwidth per period. Used
#'   only with `bwselect = "iter"`.
#' @param regularize logical. If `TRUE` (default), the `"joint"` and `"iter"`
#'   rules add a regularization term to the estimated squared bias, as
#'   \pkg{rdrobust} does, so that a near-zero estimated bias cannot make the
#'   bandwidth very large. Not used with a given `h` or `bwselect = "cct"`.
#' @param reg_const the regularization constant: the multiple of the estimated
#'   variance of the bias constants added to the squared bias (default 3).
#' @param weights the old name of `trend`, still accepted with a message. It is
#'   not a vector of observation weights; `rddid()` has none.
#'
#' @return An object of class `"rddid"`, a list with:
#'   \describe{
#'     \item{`estimates`}{a data frame with two rows, `Conventional` (the
#'       local-linear estimate with its conventional standard error) and
#'       `Robust` (the bias-corrected estimate with its robust standard error,
#'       printed as "Robust (bias-corrected)"), and columns `est`, `se`, `ci_l`,
#'       `ci_u`, `z`, `p` (estimate, standard error, confidence limits, z
#'       statistic and p-value under `scheme`) and `se_cs`, `se_pc`, `se_pv`
#'       (that row's standard error under each sampling scheme, at the same
#'       bandwidths).}
#'     \item{`coef`}{named numeric vector, the coefficient of each period's jump
#'       in the estimate: 1 for the RD period, minus its weight for each
#'       comparison period. [coef()] returns the estimate itself.}
#'     \item{`weights`}{named numeric vector of the comparison-period weights.}
#'     \item{`weights_type`}{`"constant"`, `"linear"`, or `"custom"` for numeric
#'       weights.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`t_rd`, `comparisons`}{the RD period and the comparison periods
#'       used.}
#'     \item{`scheme`}{the sampling scheme behind `se`, `ci_l`, `ci_u`, `z` and
#'       `p`; `scheme_detected` is the scheme read off the data and
#'       `scheme_requested` the argument as passed.}
#'     \item{`bandwidth`}{a list: `method` (the `bwselect` value, or `"fixed"`
#'       when `h` is given); `h_by_period` and `b_by_period` (named numeric
#'       vectors, the main and pilot bandwidth used in each period, whatever the
#'       rule); for `"fixed"` and `"joint"`, the common main bandwidth `h` and
#'       the pilot `b` (one number for `"fixed"`, one per period for
#'       `"joint"`); for `"iter"`, the number of iterations `niter`; and
#'       intermediate quantities of the rule (`bws` for `"cct"` and `"iter"`;
#'       `B`, `Veff`, `reg`, `pilot_bws` for `"joint"`).}
#'     \item{`fits`}{named list of the per-period [rd_period()] fits, the RD
#'       period first (index by name).}
#'     \item{`n_by_period`}{the number of observations in each period.}
#'     \item{`level`, `c`, `p`, `q`, `kernel`}{as passed.}
#'     \item{`data`}{the per-period data used (a named list of data frames with
#'       columns `y`, `x`, `id`), for [plot.rddid()].}
#'     \item{`call`}{the matched call.}
#'   }
#'
#' @references
#' Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression Discontinuity
#' Designs Using Multiple Time-Period Data.* arXiv:2408.05847.
#' \url{https://arxiv.org/abs/2408.05847}
#'
#' Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust nonparametric
#' confidence intervals for regression-discontinuity designs. *Econometrica*
#' 82(6), 2295-2326.
#'
#' @seealso [summary.rddid()] and [rddid-methods] (`coef()`, `confint()`,
#'   `nobs()`) for the fitted object; the tests of the assumptions
#'   [rd_typecont()], [rd_compstable()], [rd_homog()] and [rd_trendcell()]; the
#'   example data [rddid_sim] and [rddid_sim_pv].
#' @family RD-DID estimation
#'
#' @examples
#' # rddid_sim: confounding jump 0.5 in every year, treatment effect 1 in year 3
#' fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' fit            # the estimate, its standard error and confidence interval
#' summary(fit)   # the jump in every period, and the s.e. under each scheme
#' coef(fit)
#' confint(fit, "Robust")
#' # the confounding jump moves linearly in time
#' rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
#'       trend = "linear")
#' # comparison periods uniformly treated: add estimand = "atu" (same numbers)
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
  kernel   <- match.arg(tolower(kernel), c("triangular", "epanechnikov", "uniform"))
  cutoff   <- c                            # `c` stays the argument name, as in rdrobust
  if (!is.null(weights)) {                 # `weights` is the old name of `trend`
    if (!missing(trend))
      stop("supply either `trend` or its old name `weights`, not both.")
    message("rddid(): `weights` is now called `trend`; the old name still works.")
    trend <- weights
  }

  # ---- periods and data -------------------------------------------------------------------
  .check_columns(data, c(y, x, time, id))
  time_col <- data[[time]]
  if (!t_rd %in% time_col) stop("t_rd = ", t_rd, " not present in `", time, "`.")
  if (is.null(comparisons)) comparisons <- sort(setdiff(unique(time_col), t_rd))
  if (length(comparisons) < 1L) stop("need at least one comparison period.")
  if (t_rd %in% comparisons) stop("`t_rd` must not be one of the `comparisons`.")
  if (anyDuplicated(comparisons)) stop("`comparisons` contains a period more than once.")
  if (!is.null(id)) {
    dup <- duplicated(data[, c(id, time)]) & !is.na(data[[id]]) & !is.na(time_col)
    if (any(dup))
      stop(sum(dup), " row(s) repeat a (", id, ", ", time, ") pair: each unit may appear once ",
           "per period (aggregate or deduplicate first).")
  }
  missing_comp <- setdiff(comparisons, time_col)
  if (length(missing_comp))
    stop("comparison period(s) not in `", time, "`: ", paste(missing_comp, collapse = ", "))
  if (!is.numeric(time_col) && !is.character(time_col) && !is.factor(time_col))
    stop("the period column `", time, "` must be numeric, character or factor (a Date can be ",
         "converted with as.numeric(format(x, \"%Y\")) or similar).")
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a single number strictly between 0 and 1 (e.g. 0.95).")
  if (!is.null(b) && is.null(h))
    stop("`b` is the pilot bandwidth that goes with a fixed `h`; supply `h` too, or neither.")
  if (!is.null(h) && (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0))
    stop("`h` must be a single positive number.")
  if (is.factor(t_rd)) t_rd <- as.character(t_rd)
  if (is.factor(comparisons)) comparisons <- as.character(comparisons)
  if (is.null(id) && scheme == "auto")
    message("rddid(): no `id` given, so every row is treated as a different unit ",
            "(repeated cross-section standard errors).")
  periods <- c(t_rd, comparisons)          # the RD period first, then the comparisons
  plist   <- .rddid_period_frames(data, y, x, id, time_col, periods)

  # ---- the estimator's coefficients: +1 on the RD period, minus the weight on each comparison
  comp_weights <- .rddid_weights(trend, comparisons, t_rd)
  period_coef  <- c(stats::setNames(1, as.character(t_rd)), -comp_weights)

  # ---- sampling scheme (sets the standard error and, for the joint rules, the bandwidth) --
  detected   <- .detect_scheme(plist, c = cutoff)
  use_scheme <- if (scheme == "auto") detected else scheme
  if (scheme %in% c("pc", "pv") && detected == "cs")
    warning("scheme = \"", scheme, "\" requested but no unit id repeats across periods; ",
            "all cross-period covariances are zero, so the reported SE equals the CS one.")
  else if (scheme != "auto" && scheme != detected)
    message("rddid(): scheme = \"", scheme, "\" as requested; the data look like \"", detected,
            "\" (", .scheme_label(detected), ").")

  # ---- bandwidths, then one local-linear fit per period at those bandwidths ---------------
  bw <- .rddid_bandwidths(plist, period_coef, t_rd, periods, use_scheme, h, b, bwselect, start,
                          cutoff, p, q, kernel, regularize, reg_const)
  fits <- stats::setNames(lapply(as.character(periods), function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = bw$bws[[k]]["h"], b = bw$bws[[k]]["b"],
              id = plist[[k]]$id, c = cutoff, p = p, q = q, kernel = kernel)),
    as.character(periods))

  # ---- aggregate: estimate and variance under every scheme, conventional and bias-corrected
  agg_conv <- .aggregate_fits(fits, period_coef, bc = FALSE)
  agg_bc   <- .aggregate_fits(fits, period_coef, bc = TRUE)

  structure(list(
    estimates = .rddid_estimate_table(agg_conv, agg_bc, use_scheme, level),
    coef = period_coef, weights = comp_weights,
    weights_type = if (is.numeric(trend)) "custom" else match.arg(trend, c("constant", "linear")),
    estimand = estimand,
    t_rd = t_rd, comparisons = comparisons,
    scheme = use_scheme, scheme_detected = detected, scheme_requested = scheme,
    bandwidth = bw$info, fits = fits, level = level,
    p = p, q = q, kernel = kernel, c = cutoff,
    n_by_period = vapply(fits, function(f) f$n, numeric(1)),
    data = plist,
    call = cl
  ), class = "rddid")
}

#' One (y, x, id) data frame per period, the RD period first
#' @keywords internal
#' @noRd
.rddid_period_frames <- function(data, y, x, id, time_col, periods) {
  unit_id <- if (is.null(id)) seq_len(nrow(data)) else data[[id]]
  frames <- stats::setNames(lapply(periods, function(tv) {
    rows <- which(time_col == tv)
    data.frame(y = data[[y]][rows], x = data[[x]][rows], id = unit_id[rows])
  }), as.character(periods))
  # only complete (y, x, id) rows enter: the fits drop them anyway, and the sampling scheme,
  # n_by_period and nobs() are read off the rows actually used (an NA running variable used
  # to read as a side switch and turn a "pc" panel into "pv")
  dropped <- 0L
  frames <- lapply(frames, function(d) {
    ok <- stats::complete.cases(d)
    dropped <<- dropped + sum(!ok)
    d[ok, , drop = FALSE]
  })
  if (dropped > 0L)
    message("rddid(): ", dropped, " row(s) with a missing outcome, running variable or id dropped.")
  frames
}

#' Bandwidth pair (h, b) for every period under the requested rule
#'
#' @return list with `bws` (named list of `c(h, b)` per period, used for the fits) and `info`
#'   (the `bandwidth` element of the returned object: method, rule-specific details, and the
#'   per-period `h_by_period` / `b_by_period`).
#' @keywords internal
#' @noRd
.rddid_bandwidths <- function(plist, period_coef, t_rd, periods, use_scheme, h, b, bwselect,
                              start, cutoff, p, q, kernel, regularize, reg_const) {
  labels <- as.character(periods)
  if (!is.null(h)) {
    if (is.null(b)) b <- h
    bws  <- stats::setNames(rep(list(c(h = h, b = b)), length(periods)), labels)
    info <- list(method = "fixed", h = h, b = b)
  } else if (bwselect == "cct") {
    bws  <- .bw_cct(plist, c = cutoff, p = p, q = q, kernel = kernel)
    info <- list(method = "cct", bws = bws)
  } else if (bwselect == "iter") {
    ib <- .bw_joint_iter(plist, period_coef, as.character(t_rd), scheme = use_scheme,
                         start = start,
                         c = cutoff, p = p, q = q, kernel = kernel,
                         regularize = regularize, reg_const = reg_const)
    bws  <- ib$bws
    info <- list(method = "iter", bws = bws, niter = ib$niter)
  } else {
    jb <- .bw_joint(plist, period_coef, as.character(t_rd), scheme = use_scheme,
                    c = cutoff, p = p, q = q, kernel = kernel,
                    regularize = regularize, reg_const = reg_const)
    # one common h; each period's pilot b keeps its own CCT ratio b/h
    bws  <- stats::setNames(lapply(labels, function(k)
      c(h = jb$h, b = unname(jb$b[[k]]))), labels)
    info <- list(method = "joint", h = jb$h, b = jb$b, B = jb$B,
                 Veff = jb$Veff, reg = jb$reg, pilot_bws = jb$pilot_bws)
  }
  # the bandwidths actually used in each period, whatever the rule
  info$h_by_period <- vapply(bws, function(v) unname(v[["h"]]), numeric(1))
  info$b_by_period <- vapply(bws, function(v) unname(v[["b"]]), numeric(1))
  list(bws = bws, info = info)
}

#' The two-row estimate table (Conventional, Robust) at the scheme's standard error
#' @keywords internal
#' @noRd
.rddid_estimate_table <- function(agg_conv, agg_bc, use_scheme, level) {
  var_field <- c(cs = "V_cs", pc = "V_pc", pv = "V_pv")
  z_crit    <- stats::qnorm(1 - (1 - level) / 2)
  make_row  <- function(agg) {
    se <- unname(sqrt(agg[var_field[use_scheme]]))
    z  <- agg[["est"]] / se
    c(est = agg[["est"]],
      se = se,
      se_cs = sqrt(agg[["V_cs"]]), se_pc = sqrt(agg[["V_pc"]]), se_pv = sqrt(agg[["V_pv"]]),
      ci_l = agg[["est"]] - z_crit * se, ci_u = agg[["est"]] + z_crit * se,
      z = z, p = 2 * stats::pnorm(-abs(z)))
  }
  as.data.frame(rbind(Conventional = make_row(agg_conv), Robust = make_row(agg_bc)))
}

#' @export
print.rddid <- function(x, digits = 4, hint = TRUE, ...) {
  est <- if (is.null(x$estimand)) "att" else x$estimand
  cat(sprintf("RD-DID estimate of the %s in period %s\n", toupper(est), x$t_rd))
  cat(sprintf("  Comparison periods: %s   (%s; weights %s)\n",
              paste(x$comparisons, collapse = ", "), .trend_label(x$weights_type),
              paste(trimws(formatC(x$weights, digits = 3, format = "g")), collapse = ", ")))
  cat(sprintf("  Sampling scheme: %s%s\n", .scheme_label(x$scheme),
              if (identical(x$scheme_requested, "auto")) " (detected)" else ""))
  cat(sprintf("  Bandwidth: %s\n\n", .bandwidth_label(x$bandwidth)))
  .print_estimates(x, digits = digits)
  if (hint)
    cat("\n  summary() shows the per-period fits and the s.e. under every sampling scheme.\n")
  invisible(x)
}

# ---- print helpers shared by print.rddid / summary.rddid ------------------------------------
#' Sampling scheme in words, for printouts
#' @keywords internal
#' @noRd
.scheme_label <- function(s) {
  labels <- c(cs = "repeated cross-section",
              pc = "panel, no unit changes side of the cutoff",
              pv = "panel, some units change side of the cutoff")
  if (s %in% names(labels)) labels[[s]] else s   # e.g. "mixed" (rd_compstable pairs differ)
}
#' Confounding-trend assumption in words, for printouts
#' @keywords internal
#' @noRd
.trend_label <- function(wt) {
  switch(wt, constant = "constant confounding trend", linear = "linear confounding trend",
         custom = "user-supplied weights", wt)
}
#' Bandwidth rule and the values used, for print.rddid
#' @keywords internal
#' @noRd
.bandwidth_label <- function(bw) {
  by_t <- function(v) {
    v <- v[.period_order(names(v))]            # time order, not the RD-period-first storage
    paste(sprintf("%s = %s", names(v), trimws(formatC(v, digits = 4, format = "g"))),
          collapse = ", ")
  }
  switch(bw$method,
    fixed = sprintf("h = %.4g in every period (fixed), pilot b = %.4g", bw$h, bw$b),
    joint = sprintf("common h = %.4g (rule \"joint\": one bandwidth, chosen for the RD-DID estimate)\n  Pilot bandwidth b (period = value): %s",
                    bw$h, by_t(bw$b_by_period)),
    cct   = sprintf("per-period CCT MSE-optimal (rule \"cct\"): h (period = value) %s", by_t(bw$h_by_period)),
    iter  = sprintf("period-specific (rule \"iter\", %d iterations): h (period = value) %s", bw$niter,
                    by_t(bw$h_by_period)))
}
#' The two-row estimate table of print.rddid
#' @keywords internal
#' @noRd
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
