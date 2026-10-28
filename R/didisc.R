# didisc.R -- the main estimator didisc() and its print method. The pieces it calls live in:
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
#' paper assumes), the jump there *is* the confounding jump. `didisc()` estimates the jump in every period by
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
#' \eqn{t}. `summary()` lists every \eqn{D_t} with its signed weight in the sum.
#'
#' ## Several RD periods
#'
#' With a vector `t_rd`, `didisc()` estimates the weighted average of the RD
#' periods' effects, \eqn{\mathrm{ATT}^\phi = \sum_t \phi_t \mathrm{ATT}(t)}{ATT^phi = sum_t phi_t ATT(t)},
#' with weights \eqn{\phi_t}{phi_t} from `rd_weights` (equal by default; any
#' finite weights are allowed, as they define the estimand). Each period's jump
#' is estimated once: each RD period \eqn{t} enters with \eqn{\phi_t}{phi_t},
#' and each comparison period \eqn{s} with the aggregated weight
#' \eqn{w^\phi_s = \sum_t \phi_t w_s(t)}{w^phi_s = sum_t phi_t w_s(t)}, where
#' \eqn{w_s(t)}{w_s(t)} is the weight that RD period \eqn{t} alone gives
#' comparison period \eqn{s} under the trend assumption. The standard errors
#' include the covariances between RD periods, and `"joint"` chooses one
#' bandwidth for the average, not for any single RD period. `weighting =
#' "min_variance"` is available for a single RD period only.
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
#' A numeric `trend` supplies the weights directly. These are the weights under
#' the default `weighting = "ols"`; `weighting = "min_variance"` chooses others
#' the same assumption allows (see below).
#'
#' ## Minimum-variance weights
#'
#' All weights the confounding-trend assumption allows estimate the same
#' effect; they differ in precision. `weighting = "min_variance"` picks, for each
#' row of `estimates`, the ones with the smallest variance of that estimate. For
#' the conventional estimate it uses the estimated
#' variances of the per-period jumps and, in a panel (`scheme` `"pc"` or
#' `"pv"`), their covariances (the quantities behind the conventional standard
#' error). In a
#' repeated cross-section with `trend = "constant"` these are inverse-variance
#' weights: a comparison period whose jump is estimated more precisely gets more
#' weight. Under `bwselect = "joint"` the weights and the bandwidth are chosen in
#' three steps: the common bandwidth with the `"ols"` weights, the
#' minimum-variance weights at that bandwidth, and the bandwidths (common `h`
#' and pilot `b`) again with those weights. With `"cct"` or a fixed `h` the
#' bandwidths do not depend on the weights, so the weights are computed once.
#' The bias-corrected estimate gets its own weights, computed in the same way
#' from the variances and covariances of the bias-corrected jumps (those behind
#' the robust standard error), at the same bandwidths; they are in
#' `weights_detail$weights_bc`, and the common bandwidth is chosen with the
#' conventional ones. The standard errors treat the weights as known; in large samples,
#' estimating them does not change the distribution of the estimate, because
#' every allowed set of weights cancels the confounding jump exactly, so an
#' error in the weights multiplies only the estimation error of the jumps. With
#' one comparison period under `"constant"`, or two under `"linear"`, the
#' assumption allows only one set of weights; `didisc()` returns it and sets
#' `weights_detail$pinned` to `TRUE`.
#'
#' ## Sampling schemes
#'
#' The scheme sets the standard error and, under `"joint"`, also
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
#'   to minimize the asymptotic mean squared error of the diff-in-disc estimate, not
#'   of each period's jump, so the biases of the period jumps can partly cancel.
#' * `"cct"`, per-period CCT bandwidths: each period gets its own MSE-optimal
#'   bandwidth from [rdrobust::rdbwselect()], as if it were a stand-alone RD
#'   (Calonico, Cattaneo and Titiunik, 2014).
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
#' error (AMSE) of the diff-in-disc estimate over a common `h`, in the way
#' \pkg{rdrobust} chooses a bandwidth for one period. Each period's variance
#' constant, under `"pc"` the covariances between periods, and the design
#' moment in each period's bias constant are estimated at one rule-of-thumb
#' bandwidth shared by every period (the preliminary bandwidth of
#' \pkg{rdrobust}, computed from the running variables of all periods pooled;
#' `bandwidth$h_rot` in the result). The other factor of the bias constant, the
#' curvature of the conditional mean, is estimated in every period from a fit at
#' one pilot bandwidth `b`, which is itself chosen for the diff-in-disc estimate as
#' \pkg{rdrobust} chooses its pilot for one period: `b` is AMSE-optimal for the
#' same combination of the periods' \eqn{(p+1)}{(p+1)}-th derivatives, with its
#' own curvature from fits at a bandwidth `d` chosen the same way one order up
#' (`bandwidth$d`). These are combined, with the signed weights of the sum above,
#' into the AMSE, which is then minimized in closed form. Every period's bias
#' correction uses the same `b`. With one period this is
#' [rdrobust::rdbwselect()] (`bwselect = "mserd"`, `vce = "hc1"`); in a panel
#' with a running variable fixed over time, \pkg{rdrobust}'s selector applied
#' to the combined outcome. Nothing in the rule depends on which period is the
#' RD period, so an aggregate of several RD periods gets the same bandwidth
#' whatever the order of `t_rd`.
#' With `regularize = TRUE` the rule adds `reg_const` times the estimated
#' variance of the aggregate bias constant, including its covariances between
#' periods, to the squared bias of `h` and of `b` (not of `d`), as
#' \pkg{rdrobust} does, so a near-zero estimated bias cannot make a bandwidth
#' very large. Standard errors use the HC1 convention
#' of \pkg{rdrobust}. The derivations are in the paper's appendix on estimation
#' and bandwidth choice.
#'
#' @param data a data frame in long format, one row per unit and period: a
#'   repeated cross-section (different units in each period) or a panel (the
#'   same units in several periods; it need not be balanced). Only the rows of
#'   `t_rd` and `comparisons` are used; among them, rows with a missing outcome,
#'   running variable, period or id are dropped, with a message.
#' @param y name of the outcome column (a string).
#' @param x name of the running-variable column (a string).
#' @param time name of the period column (a string). The column is usually
#'   numeric (a year); character or factor labels work as well, except under
#'   `trend = "linear"` (in [didisc()] or [rd_trendcell()]), where the line is
#'   fitted on the period values: these then need to be numeric and are the
#'   time scale of the line (so 2015, 2017, 2018 are unequally spaced);
#'   `as.numeric()` converts character labels such as `"2019"`. In [didisc()]
#'   the default comparison periods are every period other than `t_rd`, in
#'   sorted order of the values (alphabetical for character or factor labels),
#'   the order a numeric `trend` follows.
#' @param id name of the unit-identifier column (a string), needed for panel
#'   standard errors. With `NULL` (default) every row is treated as a different
#'   unit, which gives repeated cross-section standard errors (with a message).
#' @param t_rd the RD period: the value of `time` in which the treatment of
#'   interest switches on at the cutoff. A vector of several RD periods estimates
#'   the weighted average of their effects (see "Several RD periods").
#' @param comparisons the comparison periods: values of `time` in which the
#'   treatment of interest is uniform at the cutoff (nobody treated, or, with
#'   `estimand = "atu"`, everybody treated). `NULL` (default) uses every period
#'   other than `t_rd`; pass the periods explicitly when the data contain
#'   periods that belong in neither.
#' @param trend the confounding-trend assumption, which restricts the
#'   comparison-period weights (`weighting` picks them within it): `"constant"`
#'   (default; the confounding jump is the same in every period; equal weights
#'   under the default `weighting`) or `"linear"` (the confounding
#'   jump moves linearly in time; needs at least two comparison periods). A
#'   numeric vector gives the weights directly, one per entry of `comparisons`
#'   in that order. Weights that sum to one cancel a constant confounding
#'   jump; `didisc()` warns when they do not. With several RD periods a numeric
#'   vector applies to every RD period, and a numeric matrix gives one row of
#'   weights per RD period (rows in `t_rd` order or named by period; columns in
#'   `comparisons` order or named by period).
#' @param rd_weights the weights \eqn{\phi_t}{phi_t} of the RD periods in the
#'   average, one per entry of `t_rd` (in that order, or named by period).
#'   `NULL` (default) gives equal weights. Any finite weights are allowed;
#'   weights that sum to one make the estimand a weighted average. Only for
#'   several RD periods (a single RD period has weight 1).
#' @param weighting how the comparison-period weights are chosen within the
#'   `trend` assumption: `"ols"` (default) uses the weights `trend` implies
#'   (equal weights for `"constant"`, the least-squares line for `"linear"`);
#'   `"min_variance"` uses the minimum-variance weights, the comparison-period
#'   weights the confounding-trend assumption allows that make the variance of
#'   the estimate smallest (see Details). For a single RD period and a
#'   character `trend`.
#' @param estimand `"att"` (default) when the comparison periods are uniformly
#'   untreated: the ATT, the effect of the treatment on the units just above the
#'   cutoff that are treated in the RD period. `"atu"` when the comparison
#'   periods are uniformly treated: the ATU, the effect for the units just below
#'   the cutoff, which are untreated in the RD period. The numbers are the same
#'   either way; see "Targeting the ATU".
#' @param bwselect the bandwidth rule, used when `h` is not given: `"joint"`
#'   (default; the common bandwidth, one `h` for every period, chosen to minimize
#'   the asymptotic mean squared error of the diff-in-disc estimate), `"cct"`
#'   (per-period CCT bandwidths: each period's own MSE-optimal bandwidth of
#'   Calonico, Cattaneo and Titiunik, from [rd_bw_cct()]). See "Bandwidth
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
#' @param regularize logical. If `TRUE` (default), the `"joint"` rule adds a
#'   regularization term to the estimated squared bias, as
#'   \pkg{rdrobust} does, so that a near-zero estimated bias cannot make the
#'   bandwidth very large. Not used with a given `h` or `bwselect = "cct"`.
#' @param reg_const the regularization constant: the multiple of the estimated
#'   variance of the bias constants added to the squared bias (default 3).
#'
#' @return An object of class `"didisc"`, a list with:
#'   \describe{
#'     \item{`estimates`}{a data frame with two rows, `Conventional` (the
#'       local-linear estimate with its conventional standard error) and
#'       `Robust` (the bias-corrected estimate with its robust standard error,
#'       printed as "Robust (bias-corrected)"), and columns `est`, `se`, `ci_l`,
#'       `ci_u`, `z`, `p` (estimate, standard error, confidence limits, z
#'       statistic and p-value under `scheme`) and `se_cs`, `se_pc`, `se_pv`
#'       (that row's standard error under each sampling scheme, at the same
#'       bandwidths).}
#'     \item{`period_weights`}{named numeric vector, the signed weight of each period's
#'       jump in the estimate: \eqn{\phi_t}{phi_t} for each RD period (1 for a single
#'       RD period), minus its (aggregated) weight for each comparison period.
#'       [coef()] returns the estimate itself.}
#'     \item{`period_weights_bc`}{the same for the bias-corrected estimate; it
#'       differs from `period_weights` only under `weighting = "min_variance"`.}
#'     \item{`weights`}{named numeric vector of the comparison-period weights;
#'       with several RD periods, the aggregated weights
#'       \eqn{w^\phi_s}{w^phi_s}; under `weighting = "min_variance"`, those of
#'       the conventional estimate (the bias-corrected estimate's are in
#'       `weights_detail$weights_bc`).}
#'     \item{`rd_weights`}{named numeric vector of the RD-period weights
#'       \eqn{\phi_t}{phi_t} (1 for a single RD period).}
#'     \item{`weights_by_target`}{numeric matrix of the comparison-period
#'       weights for each RD period, one row per RD period.}
#'     \item{`weights_type`}{`"constant"`, `"linear"`, or `"custom"` for numeric
#'       weights.}
#'     \item{`estimand`}{`"att"` or `"atu"`, as passed.}
#'     \item{`t_rd`, `comparisons`}{the RD period(s) and the comparison
#'       periods used.}
#'     \item{`scheme`}{the sampling scheme behind `se`, `ci_l`, `ci_u`, `z` and
#'       `p`; `scheme_detected` is the scheme read off the data and
#'       `scheme_requested` the argument as passed.}
#'     \item{`bandwidth`}{a list: `method` (the `bwselect` value, or `"fixed"`
#'       when `h` is given); `h_by_period` and `b_by_period` (named numeric
#'       vectors, the main and pilot bandwidth used in each period, whatever the
#'       rule); for `"fixed"` and `"joint"`, the common main bandwidth `h` and
#'       the pilot `b`, one number each; and intermediate quantities of the
#'       rule (`bws` for `"cct"`; `h_rot`, `d`, `B`, `Veff`, `reg` for
#'       `"joint"`).}
#'     \item{`fits`}{named list of the per-period [rd_period()] fits, the RD
#'       period(s) first (index by name).}
#'     \item{`n_by_period`}{the number of observations in each period.}
#'     \item{`level`, `c`, `p`, `q`, `kernel`}{as passed.}
#'     \item{`data`}{the per-period data used (a named list of data frames with
#'       columns `y`, `x`, `id`), for [plot.didisc()].}
#'     \item{`call`}{the matched call.}
#'     \item{`weighting`}{`"ols"` or `"min_variance"`, as passed.}
#'     \item{`weights_detail`}{`NULL` for `"ols"`; for `"min_variance"` a list
#'       with `pilot` (the `"ols"` weights), `Xi` (the estimated covariance
#'       matrix of the comparison-period jumps), `xi` (their estimated
#'       covariances with the RD-period jump), `h` (the main bandwidths at which
#'       these were estimated; under `"joint"`, the first-step bandwidths, not
#'       the final ones in `bandwidth`) and `pinned` (`TRUE` when the assumption
#'       allows only one set of weights), and `weights_bc`, `Xi_bc`, `xi_bc`, the
#'       same for the bias-corrected estimate.}
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
#' @seealso [summary.didisc()] and [didisc-methods] (`coef()`, `confint()`,
#'   `nobs()`) for the fitted object; the tests of the assumptions
#'   [rd_typecont()], [rd_compstable()], [rd_homog()] and [rd_trendcell()]; the
#'   example data [didisc_sim] and [didisc_sim_pv].
#' @family diff-in-disc estimation
#'
#' @examples
#' # didisc_sim: confounding jump 0.5 in every year, treatment effect 1 in year 3
#' fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' fit            # the estimate, its standard error and confidence interval
#' summary(fit)   # the jump in every period, and the s.e. under each scheme
#' coef(fit)
#' confint(fit, "Robust")
#' # the confounding jump moves linearly in time
#' didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
#'       trend = "linear")
#' # minimum-variance comparison-period weights
#' didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
#'       weighting = "min_variance")
#' # several RD periods: the average of their effects (year 3 copied as a
#' # second RD year, 4, for illustration)
#' sim4 <- rbind(didisc_sim, transform(didisc_sim[didisc_sim$year == 3, ], year = 4))
#' didisc(sim4, y = "Y", x = "R", time = "year", id = "id", t_rd = c(3, 4))
#' # comparison periods uniformly treated: add estimand = "atu" (same numbers)
#' @export
didisc <- function(data, y, x, time, id = NULL, t_rd,
                  comparisons = NULL, trend = "constant", rd_weights = NULL,
                  weighting = c("ols", "min_variance"),
                  estimand = c("att", "atu"),
                  bwselect = c("joint", "cct"), h = NULL, b = NULL,
                  scheme = c("auto", "cs", "pc", "pv"),
                  c = 0, p = 1L, q = 2L, kernel = "triangular", level = 0.95,
                  regularize = TRUE, reg_const = 3) {
  cl       <- match.call()
  bwselect <- .match_arg(bwselect)
  scheme   <- .match_arg(scheme)
  estimand <- .match_arg(estimand)
  weighting <- .match_arg(weighting)
  kernel   <- .match_arg(tolower(kernel), c("triangular", "epanechnikov", "uniform"), "kernel")
  cutoff   <- c                            # `c` stays the argument name, as in rdrobust
  if (is.character(trend)) {
    if (length(trend) == 1L && trend %in% c("ols", "min_variance"))
      .abort("\"", trend, "\" is a value of `weighting`, not of `trend` (\"constant\" or ",
             "\"linear\"); use weighting = \"", trend, "\".")
    trend <- .match_arg(trend, c("constant", "linear"))
  } else if (!is.numeric(trend)) {
    .abort("`trend` must be \"constant\", \"linear\" or numeric weights.")
  }

  # ---- periods and data -------------------------------------------------------------------
  # one shared prepare step (input_checks.R): checks, factor periods, duplicates, missing rows.
  # Without an id every row is its own unit, numbered by its row in the data as passed.
  if (is.null(id)) data$.didisc_row_id <- seq_len(nrow(data))
  if (missing(t_rd)) t_rd <- NULL
  prep <- .prepare_input(data, "didisc()", x = x, time = time, id = id, y = y, t_rd = t_rd,
                         comparisons = comparisons, t_rd_scalar = FALSE)
  data <- prep$data; t_rd <- prep$t_rd; comparisons <- prep$comparisons
  time_col <- data[[time]]
  if (identical(trend, "linear") && !is.numeric(time_col))
    .abort("trend = \"linear\" fits a line in the period values, so the period column `", time,
           "` must be numeric (it is ", class(time_col)[1L], ").")
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    .abort("`level` must be a single number strictly between 0 and 1 (e.g. 0.95).")
  if (!is.null(b) && is.null(h))
    .abort("`b` is the pilot bandwidth that goes with a fixed `h`; supply `h` too, or neither.")
  if (!is.null(h) && (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0))
    .abort("`h` must be a single positive number.")
  if (is.null(id) && scheme == "auto")
    message("didisc(): no `id` given, so every row is treated as a different unit ",
            "(repeated cross-section standard errors).")
  periods <- c(t_rd, comparisons)          # the RD period first, then the comparisons
  period_order <- as.character(if (is.factor(time_col)) intersect(levels(time_col), periods)
                               else sort(periods))   # time order, for printouts and plots
  plist   <- .didisc_period_frames(data, y, x, id, time_col, periods)

  # ---- the signed period weights: phi_t on each RD period (1 for a single one), minus the
  # aggregated weight w^phi_t0 = sum_t phi_t w_t0(t) on each comparison period (Cor. cor:agg)
  rd_w <- .didisc_rd_weights(rd_weights, t_rd)
  weights_by_target <- .didisc_weights_by_target(trend, comparisons, t_rd)
  comp_weights <- colSums(rd_w * weights_by_target)
  period_weights  <- c(rd_w, -comp_weights)
  if (weighting == "min_variance") {
    if (is.numeric(trend))
      .abort("`weighting = \"min_variance\"` chooses the comparison weights and a numeric `trend` ",
             "fixes them; use one or the other.")
    if (length(t_rd) > 1L)
      .abort("`weighting = \"min_variance\"` is available for a single RD period; with several ",
             "RD periods use the default `weighting = \"ols\"`.")
  }

  # ---- sampling scheme (sets the standard error and, for the joint rules, the bandwidth) --
  detected   <- .detect_scheme(plist, c = cutoff)
  use_scheme <- if (scheme == "auto") detected else scheme
  if (scheme %in% c("pc", "pv") && detected == "cs")
    .warn("scheme = \"", scheme, "\" requested but no unit id repeats across periods; ",
          "all cross-period covariances are zero, so the reported SE equals the CS one.")
  else if (scheme != "auto" && scheme != detected)
    message("didisc(): scheme = \"", scheme, "\" as requested; the data look like \"", detected,
            "\" (", .scheme_label(detected), ").")

  # ---- bandwidths, then one local-linear fit per period at those bandwidths ---------------
  fit_all <- function(bw) stats::setNames(lapply(as.character(periods), function(k)
    rd_period(plist[[k]]$y, plist[[k]]$x, h = bw$bws[[k]]["h"], b = bw$bws[[k]]["b"],
              id = plist[[k]]$id, c = cutoff, p = p, q = q, kernel = kernel)),
    as.character(periods))
  bw <- .didisc_bandwidths(plist, period_weights, t_rd, periods, use_scheme, h, b, bwselect,
                          cutoff, p, q, kernel, regularize, reg_const)
  fits <- fit_all(bw)

  # ---- minimum-variance weights: computed from these fits (the "ols" weights are the pilot);
  # under "joint" the common bandwidth is then chosen again with them -- three steps, no
  # iteration, as the weights do not depend on a common bandwidth to first order
  weights_detail <- NULL
  period_weights_bc <- NULL                       # the bias-corrected row's weights, if different
  if (weighting == "min_variance") {
    # the unscaled weights of the single RD period (setNames: [1, ] drops the names of a 1 x 1 matrix)
    w_pilot <- stats::setNames(as.numeric(weights_by_target[1, ]), colnames(weights_by_target))
    trend_type <- trend
    mv    <- .mv_weights(fits, comparisons, t_rd, trend_type, use_scheme, w_pilot)
    # the bias-corrected estimate minimizes its own (robust) variance, from the same fits
    mv_bc <- .mv_weights(fits, comparisons, t_rd, trend_type, use_scheme, w_pilot, bc = TRUE)
    weights_detail <- list(pilot = w_pilot, Xi = mv$Xi, xi = mv$xi,
                           h = bw$info$h_by_period, pinned = mv$pinned,
                           weights_bc = rd_w[[1]] * mv_bc$w, Xi_bc = mv_bc$Xi, xi_bc = mv_bc$xi)
    comp_weights <- rd_w[[1]] * mv$w
    weights_by_target[1, ] <- mv$w
    period_weights    <- c(rd_w, -comp_weights)
    period_weights_bc <- c(rd_w, -rd_w[[1]] * mv_bc$w)
    if (is.null(h) && bwselect == "joint" && !mv$pinned) {
      bw <- .didisc_bandwidths(plist, period_weights, t_rd, periods, use_scheme, h, b, bwselect,
                              cutoff, p, q, kernel, regularize, reg_const)
      fits <- fit_all(bw)
    }
  }

  # ---- aggregate: estimate and variance under every scheme, conventional and bias-corrected
  agg_conv <- .aggregate_fits(fits, period_weights, bc = FALSE)
  if (is.null(period_weights_bc)) period_weights_bc <- period_weights
  agg_bc   <- .aggregate_fits(fits, period_weights_bc, bc = TRUE)

  structure(list(
    estimates = .didisc_estimate_table(agg_conv, agg_bc, use_scheme, level),
    period_weights = period_weights, weights = comp_weights,
    rd_weights = rd_w, weights_by_target = weights_by_target, period_weights_bc = period_weights_bc,
    weights_type = if (is.numeric(trend)) "custom" else trend,
    estimand = estimand,
    t_rd = t_rd, comparisons = comparisons, period_order = period_order,
    scheme = use_scheme, scheme_detected = detected, scheme_requested = scheme,
    bandwidth = bw$info, fits = fits, level = level,
    p = p, q = q, kernel = kernel, c = cutoff,
    n_by_period = vapply(fits, function(f) f$n, numeric(1)),
    data = plist,
    call = cl,
    weighting = weighting, weights_detail = weights_detail
  ), class = "didisc")
}

#' One (y, x, id) data frame per period, the RD period first
#' @keywords internal
#' @noRd
.didisc_period_frames <- function(data, y, x, id, time_col, periods) {
  unit_id <- if (is.null(id)) data$.didisc_row_id else data[[id]]
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
    message("didisc(): ", dropped, " row(s) with a missing outcome, running variable or id dropped.")
  frames
}

#' Bandwidth pair (h, b) for every period under the requested rule
#'
#' @return list with `bws` (named list of `c(h, b)` per period, used for the fits) and `info`
#'   (the `bandwidth` element of the returned object: method, rule-specific details, and the
#'   per-period `h_by_period` / `b_by_period`).
#' @keywords internal
#' @noRd
.didisc_bandwidths <- function(plist, period_weights, t_rd, periods, use_scheme, h, b, bwselect,
                              cutoff, p, q, kernel, regularize, reg_const) {
  labels <- as.character(periods)
  if (!is.null(h)) {
    if (is.null(b)) b <- h
    bws  <- stats::setNames(rep(list(c(h = h, b = b)), length(periods)), labels)
    info <- list(method = "fixed", h = h, b = b)
  } else if (bwselect == "cct") {
    bws  <- .bw_cct(plist, c = cutoff, p = p, q = q, kernel = kernel)
    info <- list(method = "cct", bws = bws)
  } else {
    jb <- .bw_joint(plist, period_weights, as.character(t_rd), scheme = use_scheme,
                    c = cutoff, p = p, q = q, kernel = kernel,
                    regularize = regularize, reg_const = reg_const)
    # one common h and one common pilot b, both chosen for the aggregate
    bws  <- stats::setNames(rep(list(c(h = jb$h, b = jb$b)), length(labels)), labels)
    info <- list(method = "joint", h = jb$h, b = jb$b, d = jb$d, B = jb$B,
                 Veff = jb$Veff, reg = jb$reg, h_rot = jb$h_rot)
  }
  # the bandwidths actually used in each period, whatever the rule
  info$h_by_period <- vapply(bws, function(v) unname(v[["h"]]), numeric(1))
  info$b_by_period <- vapply(bws, function(v) unname(v[["b"]]), numeric(1))
  list(bws = bws, info = info)
}

#' The two-row estimate table (Conventional, Robust) at the scheme's standard error
#' @keywords internal
#' @noRd
.didisc_estimate_table <- function(agg_conv, agg_bc, use_scheme, level) {
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
print.didisc <- function(x, digits = 4, hint = TRUE, ...) {
  est <- if (is.null(x$estimand)) "att" else x$estimand
  # weights as "period: weight", three decimals
  fmt_w <- function(w) {
    v <- formatC(unname(w), digits = 3, format = "f")
    paste(if (is.null(names(w))) v else paste0(names(w), ": ", v), collapse = ", ")
  }
  several <- length(x$t_rd) > 1L
  if (several)
    cat(sprintf("Diff-in-disc estimate of the %s %s RD periods %s\n", toupper(est),
                if (abs(sum(x$rd_weights) - 1) < 1e-8) "averaged over" else "weighted sum over",
                paste(x$t_rd, collapse = ", ")),
        sprintf("  RD-period weights: %s\n", fmt_w(x$rd_weights)), sep = "")
  else
    cat(sprintf("Diff-in-disc estimate of the %s in period %s\n", toupper(est), x$t_rd))
  cat(sprintf("  Comparison periods: %s   (%s)\n", paste(x$comparisons, collapse = ", "),
              .trend_label(x$weights_type)))
  w_bc <- x$weights_detail$weights_bc
  own_bc <- !is.null(w_bc) && !isTRUE(all.equal(unname(w_bc), unname(x$weights)))
  cat(sprintf("  Comparison weights%s: %s%s\n",
              if (identical(x$weighting, "min_variance")) " (minimum variance)"
              else if (several) " (summed over the RD periods)" else "",
              fmt_w(x$weights),
              if (own_bc) paste0("; Robust row: ", fmt_w(stats::setNames(w_bc, names(x$weights))))
              else ""))
  cat(sprintf("  Sampling scheme: %s%s\n", .scheme_label(x$scheme),
              if (identical(x$scheme_requested, "auto")) " (detected)" else ""))
  cat(sprintf("  Bandwidth: %s\n\n", .bandwidth_label(x$bandwidth)))
  .print_estimates(x, digits = digits)
  if (hint)
    cat("\n  summary() shows the per-period fits and the s.e. under every sampling scheme.\n")
  invisible(x)
}

# ---- print helpers shared by print.didisc / summary.didisc ------------------------------------
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
#' Bandwidth rule and the values used, for print.didisc
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
    joint = sprintf("common h = %.4g (rule \"joint\": one bandwidth, chosen for the diff-in-disc estimate)\n  Pilot bandwidth b = %.4g (also chosen for the diff-in-disc estimate)",
                    bw$h, bw$b),
    cct   = sprintf("per-period CCT MSE-optimal (rule \"cct\"): h (period = value) %s", by_t(bw$h_by_period)))
}
#' The two-row estimate table of print.didisc
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
