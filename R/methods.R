# methods.R -- standard methods for the objects the package returns: summary(), coef(), confint(),
# nobs(), broom-style tidy()/glance() (registered for the generics package), and the print helpers
# shared by the four assumption tests. Nothing here re-estimates anything.

# Standard methods for the objects the package returns: summary(), coef(), confint(), nobs(),
# and broom-style tidy()/glance() (registered for the `generics` package when it is loaded).
# Everything here is computed from fields that rddid()/the tests already store; nothing is
# re-estimated.

#' Summary of an RD-DID fit
#'
#' `summary()` of an [rddid()] object adds, to what `print()` shows, a per-period table (the
#' local-linear jump in every period with its bandwidths and sample size, and its coefficient
#' in the estimate) and the standard error of the estimate under each of the three sampling
#' schemes.
#'
#' @param object an object of class `"rddid"`.
#' @param digits number of decimals in the printed tables.
#' @param x a `"summary.rddid"` object.
#' @param ... unused.
#' @return `summary()` returns an object of class `"summary.rddid"`: a list with `fit` (the
#'   object) and `per_period`, a data frame with one row per period and columns `period`,
#'   `role` (`"RD"` or `"comparison"`), `coef` (its coefficient in the estimate, 1 for the RD
#'   period and minus its weight for a comparison period), `n`, `h`, `b`, `jump` and `se` (the
#'   conventional local-linear jump and its standard error) and `jump_bc`, `se_rb`
#'   (bias-corrected jump, robust standard error).
#' @examples
#' fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' summary(fit)
#' summary(fit)$per_period
#' @export
summary.rddid <- function(object, ...) {
  structure(list(fit = object, per_period = .rddid_period_table(object)),
            class = "summary.rddid")
}

#' @rdname summary.rddid
#' @export
print.summary.rddid <- function(x, digits = 4, ...) {
  fit <- x$fit
  print(fit, digits = digits)
  cat("\n  Per-period local-linear fits (estimate = sum of coef x jump):\n")
  tab <- x$per_period
  fmt <- function(v, d = digits) formatC(v, digits = d, format = "f")
  cat(sprintf("  %-8s %-11s %6s %6s %8s %8s %10s %9s %10s %9s\n", "period", "role", "coef",
              "n", "h", "b", "jump", "s.e.", "jump (bc)", "s.e. (rb)"))
  for (i in seq_len(nrow(tab)))
    cat(sprintf("  %-8s %-11s %6s %6d %8s %8s %10s %9s %10s %9s\n", tab$period[i], tab$role[i],
                formatC(tab$coef[i], digits = 3, format = "g"), tab$n[i], fmt(tab$h[i]),
                fmt(tab$b[i]), fmt(tab$jump[i]), fmt(tab$se[i]), fmt(tab$jump_bc[i]),
                fmt(tab$se_rb[i])))
  e <- fit$estimates
  cat(sprintf("\n  Robust s.e. under each sampling scheme:  cross-section %s   panel, fixed R %s   panel, varying R %s\n",
              fmt(e["Robust", "se_cs"]), fmt(e["Robust", "se_pc"]), fmt(e["Robust", "se_pv"])))
  cat(sprintf("  (the printed s.e. is the one for scheme \"%s\"; the others are shown for comparison)\n",
              fit$scheme))
  invisible(x)
}

#' Per-period table behind summary.rddid(): role, coefficient, n, bandwidths, jumps
#' @keywords internal
#' @noRd
.rddid_period_table <- function(x) {
  per <- names(x$fits)
  f <- x$fits
  data.frame(
    period  = per,
    role    = ifelse(per == as.character(x$t_rd), "RD", "comparison"),
    coef    = unname(x$coef[per]),
    n       = vapply(f, function(z) as.integer(z$n), integer(1)),
    h       = vapply(f, function(z) unname(z$h), numeric(1)),
    b       = vapply(f, function(z) unname(z$b), numeric(1)),
    jump    = vapply(f, function(z) z$D, numeric(1)),
    se      = vapply(f, function(z) sqrt(z$V_D), numeric(1)),
    jump_bc = vapply(f, function(z) z$D_bc, numeric(1)),
    se_rb   = vapply(f, function(z) sqrt(z$V_D_bc), numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Coefficients, confidence intervals and sample size of an RD-DID fit
#'
#' `coef()` returns the two estimates of an [rddid()] fit, `Conventional` (local-linear) and
#' `Robust` (bias-corrected); `confint()` their confidence intervals under the fit's sampling
#' scheme; `nobs()` the number of observations used.
#'
#' @param object an object of class `"rddid"`.
#' @param parm which rows of the estimate table: `"Conventional"` (local-linear estimate,
#'   conventional standard error), `"Robust"` (bias-corrected estimate, robust standard
#'   error), or both (default).
#' @param level confidence level; `NULL` (default) returns the interval stored in the object
#'   (at the `level` given to [rddid()]).
#' @param ... unused.
#' @return `coef()` a named numeric vector, the `Conventional` and `Robust` estimates;
#'   `confint()` a matrix with one row per `parm` and
#'   columns giving the lower and upper limits; `nobs()` the number of unit-period rows used
#'   across all periods.
#' @examples
#' fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' coef(fit)
#' confint(fit)
#' confint(fit, "Robust", level = 0.9)
#' nobs(fit)
#' @name rddid-methods
#' @importFrom stats coef confint nobs
NULL

#' @rdname rddid-methods
#' @export
coef.rddid <- function(object, ...) {
  e <- object$estimates
  stats::setNames(e$est, rownames(e))
}

#' @rdname rddid-methods
#' @export
confint.rddid <- function(object, parm = c("Conventional", "Robust"), level = NULL, ...) {
  parm <- match.arg(parm, several.ok = TRUE)
  e <- object$estimates[parm, , drop = FALSE]
  if (is.null(level)) {
    level <- object$level
    ci <- cbind(e$ci_l, e$ci_u)
  } else {
    zc <- stats::qnorm(1 - (1 - level) / 2)
    ci <- cbind(e$est - zc * e$se, e$est + zc * e$se)
  }
  a <- (1 - level) / 2
  dimnames(ci) <- list(parm, paste(format(100 * c(a, 1 - a), trim = TRUE), "%"))
  ci
}

#' @rdname rddid-methods
#' @export
nobs.rddid <- function(object, ...) as.integer(sum(object$n_by_period))

# ---- tidy / glance (broom conventions) -----------------------------------------------------

#' Tidy output for RD-DID fits and validation tests
#'
#' `tidy()` and `glance()` methods following the `broom` conventions, so that fits and tests
#' can be passed to table makers such as `modelsummary`. The generics come from the
#' `generics` package; the methods are registered when that package (or `broom`) is loaded.
#'
#' @param x an object returned by [rddid()], [rd_typecont()], [rd_compstable()],
#'   [rd_homog()] or [rd_trendcell()].
#' @param ... unused.
#' @return `tidy()` returns a data frame with one row per estimate (`term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`) for a fit, and one row per
#'   test (`test`, `statistic`, `df`, `p.value`) for a validation test. `glance()` returns a
#'   one-row data frame describing the fit (`nobs`, `t_rd`, `comparisons`, `trend`, `weights`,
#'   `bwselect`, `h` (the common bandwidth under `"joint"`/fixed `h`, `NA` otherwise),
#'   `scheme`, `level`).
#' @examples
#' fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   generics::tidy(fit)
#'   generics::glance(fit)
#'   generics::tidy(rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id"))
#' }
#' @name rddid-tidiers
NULL

#' @rdname rddid-tidiers
#' @exportS3Method generics::tidy
tidy.rddid <- function(x, ...) {
  e <- x$estimates
  data.frame(term = rownames(e), estimate = e$est, std.error = e$se, statistic = e$z,
             p.value = e$p, conf.low = e$ci_l, conf.high = e$ci_u,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @rdname rddid-tidiers
#' @exportS3Method generics::glance
glance.rddid <- function(x, ...) {
  data.frame(nobs = nobs.rddid(x), t_rd = x$t_rd,
             comparisons = paste(x$comparisons, collapse = ", "),
             trend = x$weights_type, weights = paste(signif(x$weights, 3), collapse = ", "),
             bwselect = x$bandwidth$method,
             h = if (!is.null(x$bandwidth$h)) unname(x$bandwidth$h) else NA_real_,
             scheme = x$scheme, level = x$level, row.names = NULL, stringsAsFactors = FALSE)
}

#' One-row tidy() data frame for a validation test
#' @keywords internal
#' @noRd
.tidy_test <- function(x, name) {
  data.frame(test = name, statistic = x$statistic, df = x$df, p.value = x$p_value,
             row.names = NULL, stringsAsFactors = FALSE)
}
#' @rdname rddid-tidiers
#' @exportS3Method generics::tidy
tidy.rd_typecont <- function(x, ...) .tidy_test(x, "type continuity")
#' @rdname rddid-tidiers
#' @exportS3Method generics::tidy
tidy.rd_compstable <- function(x, ...) .tidy_test(x, "composition stability")
#' @rdname rddid-tidiers
#' @exportS3Method generics::tidy
tidy.rd_homog <- function(x, ...) .tidy_test(x, "homogeneous confounding")
#' @rdname rddid-tidiers
#' @exportS3Method generics::tidy
tidy.rd_trendcell <- function(x, ...) .tidy_test(x, "constant within-type confounding")

# ---- shared print helpers for the validation tests ------------------------------------------
#' Header shared by the four tests' print methods: title, H0, scheme, estimand
#' @keywords internal
#' @noRd
.print_test_header <- function(title, fn, h0, scheme, scheme_detected = TRUE, estimand = "att",
                               atu_note = "test unchanged") {
  cat(sprintf("Test of %s  [%s()]\n", title, fn))
  cat(sprintf("  H0: %s\n", h0))
  if (!is.null(scheme))
    cat(sprintf("  Sampling scheme: %s%s\n", .scheme_label(scheme),
                if (scheme_detected) " (detected from the data)" else ""))
  if (identical(estimand, "atu")) cat(sprintf("  Estimand: ATU (%s)\n", atu_note))
}
#' One line: Wald chi-squared(df) = stat, p = p (or 'not testable')
#' @keywords internal
#' @noRd
.print_wald <- function(stat, df, p, label = "Wald", indent = "  ") {
  if (is.na(stat)) {
    cat(sprintf("%s%s: not testable (df = %d)\n", indent, label, df))
  } else {
    cat(sprintf("%s%s chi-squared(%d) = %.3f,  p = %s\n", indent, label, as.integer(df), stat,
                if (p < 1e-3) "<0.001" else formatC(p, digits = 3, format = "f")))
  }
}
#' Bandwidth description for the tests' printouts
#' @keywords internal
#' @noRd
.bw_label_test <- function(h, bwselect) {
  if (!is.null(h) && !is.na(h)) sprintf("h = %.4g in every cell", h)
  else if (bwselect == "cct") "CCT MSE-optimal, chosen per cell"
  else "rule of thumb, per cell"
}
#' Per-cell jump table shared by print.rd_homog / print.rd_trendcell
#' @keywords internal
#' @noRd
.print_jump_table <- function(df, cols, labels) {
  if (is.null(df) || nrow(df) == 0L) return(invisible())
  cat(sprintf("\n  Per-cell local-linear jumps (%s):\n", "comparison periods"))
  cat(sprintf("    %-10s %-12s %10s %10s %7s\n", labels[1], labels[2], "jump", "s.e.", "n"))
  for (k in seq_len(nrow(df)))
    cat(sprintf("    %-10s %-12s %10.4f %10.4f %7d%s\n", df[[cols[1]]][k], df[[cols[2]]][k],
                df$jump[k], df$se[k], df$n[k], if (isTRUE(df$reference[k])) "  (reference)" else ""))
}
#' NULL-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
