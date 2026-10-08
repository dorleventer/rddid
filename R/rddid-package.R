#' rddid: treatment effects at a cutoff shared with a confounding policy
#'
#' A treatment of interest switches on at a cutoff of a running variable in one
#' period, the **RD period**. A **confounding policy** switches at the same
#' cutoff, in every period, so the jump in the outcome at the cutoff in the RD
#' period mixes the treatment effect with the **confounding jump**. In the
#' **comparison periods** the treatment of interest is uniform at the cutoff
#' (nobody treated, or everybody treated), so the jump there *is* the
#' confounding jump. [rddid()] estimates the jump in every period by
#' local-linear RD and subtracts a weighted average of the comparison-period
#' jumps from the RD-period jump. How the weights are set is the
#' **confounding-trend assumption**: constant (equal weights) or linear in time.
#'
#' @section Workflow:
#' 1. **Estimate.** `fit <- rddid(data, y = , x = , time = , id = , t_rd = )`;
#'    printing `fit` shows the estimate, its standard error and confidence
#'    interval, the comparison periods with their weights, the sampling scheme
#'    and the bandwidth. See [rddid()].
#' 2. **Look inside.** `summary(fit)` shows the jump in every period and the
#'    standard error under each sampling scheme ([summary.rddid()]); `coef()`,
#'    `confint()` and `nobs()` work as usual ([rddid-methods]), and `tidy()`
#'    and `glance()` feed table makers ([rddid-tidiers]).
#' 3. **Check the assumptions.** When the running variable moves over time,
#'    units can change side of the cutoff between periods, and four tests check
#'    the assumptions this adds: type continuity ([rd_typecont()]), composition
#'    stability ([rd_compstable()]), homogeneous confounding ([rd_homog()]) and
#'    constant within-type confounding ([rd_trendcell()]).
#'
#' The building blocks [rd_period()] (the local-linear RD in one period) and
#' [rd_bw_cct()] (its CCT bandwidths) are exported for users who want to work
#' period by period. [rddid_sim] and [rddid_sim_pv] are simulated panels for
#' trying things out. `citation("rddid")` gives the reference below.
#'
#' @references
#' Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression Discontinuity
#' Designs Using Multiple Time-Period Data.* arXiv:2408.05847.
#' \url{https://arxiv.org/abs/2408.05847}
#'
#' @seealso The package website, with a getting-started guide and articles on
#'   the computation, the options and the tests:
#'   \url{https://dorleventer.github.io/rddid/}
#' @keywords internal
"_PACKAGE"
