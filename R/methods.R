# methods.R -- standard methods for the objects the package returns: summary(), coef(), confint(),
# nobs(), broom-style tidy()/glance() (registered for the generics package), and the print helpers
# shared by the four assumption tests. Nothing here re-estimates anything.

# Standard methods for the objects the package returns: summary(), coef(), confint(), nobs(),
# and broom-style tidy()/glance() (registered for the `generics` package when it is loaded).
# Everything here is computed from fields that didisc()/the tests already store; nothing is
# re-estimated.

#' Summary of a diff-in-disc fit
#'
#' `summary()` of a [didisc()] object adds, to what `print()` shows, a per-period table (the
#' local-linear jump in every period with its standard error, conventional and bias-corrected,
#' its bandwidths and sample size) and a table of the standard error of the estimate,
#' conventional and robust, under each of the three sampling schemes (RCS, PC, PV).
#'
#' @param object an object of class `"didisc"`.
#' @param digits number of decimals in the printed tables.
#' @param x a `"summary.didisc"` object.
#' @param ... unused.
#' @return `summary()` returns an object of class `"summary.didisc"`: a list with `fit` (the
#'   object) and `per_period`, a data frame with one row per period in time order and columns
#'   `period`, `role` (`"RD"` or `"comparison"`), `weight` (the period's signed weight in the
#'   conventional estimate: the RD period's weight in `rd_weights`, 1 for a single RD period,
#'   and minus its weight for a comparison period), `n`, `h`, `b`, `jump` and `se` (the
#'   conventional local-linear jump and its standard error), `jump_bc`, `se_rb`
#'   (bias-corrected jump, robust standard error) and, under `weighting = "min_variance"`,
#'   `weight_bc` (the signed weights of the bias-corrected estimate, which has its own weights).
#' @examples
#' fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' summary(fit)
#' summary(fit)$per_period
#' @export
summary.didisc <- function(object, ...) {
  structure(list(fit = object, per_period = .didisc_period_table(object)),
            class = "summary.didisc")
}

#' @rdname summary.didisc
#' @export
print.summary.didisc <- function(x, digits = 4, ...) {
  fit <- x$fit
  print(fit, digits = digits, hint = FALSE)
  tab <- x$per_period
  fmt    <- function(v) formatC(v, digits = digits, format = "f")
  center <- function(s, w) { l <- (w - nchar(s)) %/% 2; paste0(strrep(" ", l), s, strrep(" ", w - nchar(s) - l)) }
  conv <- "Conventional"; rob <- "Robust (bias-corrected)"
  # per-period fits: the jump and its standard error, conventional and bias-corrected (the
  # weights that combine them are printed above)
  cat("\n  Per-period local-linear fits:\n")
  cat(sprintf("  %45s   %s   %s\n", "", center(conv, 23), center(rob, 23)))
  cat(sprintf("  %-8s %-11s %6s %8s %8s   %11s %11s   %11s %11s\n", "Period", "Role", "n", "h", "b",
              "Jump", "Std. err.", "Jump", "Std. err."))
  for (i in seq_len(nrow(tab)))
    cat(sprintf("  %-8s %-11s %6d %8s %8s   %11s %11s   %11s %11s\n", tab$period[i], tab$role[i],
                tab$n[i], fmt(tab$h[i]), fmt(tab$b[i]), fmt(tab$jump[i]), fmt(tab$se[i]),
                fmt(tab$jump_bc[i]), fmt(tab$se_rb[i])))
  # standard error of the estimate under each sampling scheme, at the same bandwidths
  e <- fit$estimates
  code <- c(cs = "RCS", pc = "PC", pv = "PV")
  cat(sprintf("\n  Std. err. of the estimate under each sampling scheme (the estimate table at the top uses %s):\n",
              if (fit$scheme %in% names(code)) code[[fit$scheme]] else fit$scheme))
  cat(sprintf("  %-48s %12s   %s\n", "", conv, rob))
  for (s in names(code))
    cat(sprintf("  %-4s %-43s %12s   %23s\n", code[[s]], .scheme_label(s),
                fmt(e["Conventional", paste0("se_", s)]), fmt(e["Robust", paste0("se_", s)])))
  invisible(x)
}

#' Per-period table behind summary.didisc(): role, weight, n, bandwidths, jumps
#' @keywords internal
#' @noRd
.didisc_period_table <- function(x) {
  per <- .fit_period_order(x)              # rows in time order ...
  f   <- x$fits[per]                       # ... and every column taken in that same order
  out <- data.frame(
    period  = per,
    role    = ifelse(per %in% as.character(x$t_rd), "RD", "comparison"),
    weight  = unname(x$period_weights[per]),
    n       = vapply(f, function(z) as.integer(z$n), integer(1)),
    h       = vapply(f, function(z) unname(z$h), numeric(1)),
    b       = vapply(f, function(z) unname(z$b), numeric(1)),
    jump    = vapply(f, function(z) z$D, numeric(1)),
    se      = vapply(f, function(z) sqrt(z$V_D), numeric(1)),
    jump_bc = vapply(f, function(z) z$D_bc, numeric(1)),
    se_rb   = vapply(f, function(z) sqrt(z$V_D_bc), numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(x$period_weights_bc) && !isTRUE(all.equal(x$period_weights_bc, x$period_weights)))
    out$weight_bc <- unname(x$period_weights_bc[per])
  out
}

#' Coefficients, confidence intervals and sample size of a diff-in-disc fit
#'
#' `coef()` returns the two estimates of a [didisc()] fit, `Conventional` (local-linear) and
#' `Robust` (bias-corrected); `confint()` their confidence intervals under the fit's sampling
#' scheme; `nobs()` the number of observations used.
#'
#' @param object an object of class `"didisc"`.
#' @param parm which rows of the estimate table: `"Conventional"` (local-linear estimate,
#'   conventional standard error), `"Robust"` (bias-corrected estimate, robust standard
#'   error), or both (default).
#' @param level confidence level; `NULL` (default) returns the interval stored in the object
#'   (at the `level` given to [didisc()]).
#' @param ... unused.
#' @return `coef()` a named numeric vector, the `Conventional` and `Robust` estimates;
#'   `confint()` a matrix with one row per `parm` and
#'   columns giving the lower and upper limits; `nobs()` the number of unit-period rows used
#'   across all periods.
#' @examples
#' fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' coef(fit)
#' confint(fit)
#' confint(fit, "Robust", level = 0.9)
#' nobs(fit)
#' @name didisc-methods
#' @importFrom stats coef confint nobs
NULL

#' @rdname didisc-methods
#' @export
coef.didisc <- function(object, ...) {
  e <- object$estimates
  stats::setNames(e$est, rownames(e))
}

#' @rdname didisc-methods
#' @export
confint.didisc <- function(object, parm = c("Conventional", "Robust"), level = NULL, ...) {
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

#' @rdname didisc-methods
#' @export
nobs.didisc <- function(object, ...) as.integer(sum(object$n_by_period))

# ---- tidy / glance (broom conventions) -----------------------------------------------------

#' Tidy output for diff-in-disc fits and validation tests
#'
#' `tidy()` and `glance()` methods following the `broom` conventions, so that fits and tests
#' can be passed to table makers such as `modelsummary`. The generics come from the
#' `generics` package; the methods are registered when that package (or `broom`) is loaded.
#'
#' @param x an object returned by [didisc()], [rd_typecont()], [rd_compstable()],
#'   [rd_homog()] or [rd_trendcell()].
#' @param conf.int,conf.level for `tidy()` of a fit: include the confidence interval (default
#'   `TRUE`) and at which level (default the fit's `level`; another value recomputes the
#'   interval from the estimate and its standard error), as `modelsummary` and other table
#'   makers pass them.
#' @param ... unused.
#' @return `tidy()` returns a data frame with one row per estimate (`term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`) for a fit, and one row per
#'   test for a validation test: `test` (the function, as in [rd_checks()]), `assumption`,
#'   `statistic` (the largest |t| under `joint = "supt"` for [rd_typecont()] and
#'   [rd_compstable()], whose row is the joint test over pairs; the chi-squared statistic for
#'   [rd_homog()] and [rd_trendcell()]), `df` (the chi-squared degrees of freedom of a Wald
#'   statistic; `NA` for sup-t and sWald), `p.value` and `method`.
#'   `glance()` returns a
#'   one-row data frame describing the fit (`nobs`, `t_rd`, `rd_weights`, `comparisons`,
#'   `trend`, `weights`,
#'   `bwselect`, `h` (the common bandwidth under `"joint"`/fixed `h`, `NA` otherwise),
#'   `scheme`, `level`, `weighting`).
#' @examples
#' fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' tidy(fit)
#' glance(fit)
#' tidy(rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3))
#' @name didisc-tidiers
NULL

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @rdname didisc-tidiers
#' @export
tidy.didisc <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  e <- x$estimates
  out <- data.frame(term = rownames(e), estimate = e$est, std.error = e$se, statistic = e$z,
                    p.value = e$p, row.names = NULL, stringsAsFactors = FALSE)
  if (isTRUE(conf.int)) {
    if (is.null(conf.level) || isTRUE(all.equal(conf.level, x$level))) {
      out$conf.low  <- e$ci_l
      out$conf.high <- e$ci_u
    } else {
      zc <- stats::qnorm(1 - (1 - conf.level) / 2)
      out$conf.low  <- e$est - zc * e$se
      out$conf.high <- e$est + zc * e$se
    }
  }
  out
}

#' @rdname didisc-tidiers
#' @export
glance.didisc <- function(x, ...) {
  data.frame(nobs = nobs.didisc(x),
             t_rd = if (length(x$t_rd) == 1L) x$t_rd else paste(x$t_rd, collapse = ", "),
             comparisons = paste(x$comparisons, collapse = ", "),
             rd_weights = paste(signif(x$rd_weights, 3), collapse = ", "),
             trend = x$weights_type, weights = paste(signif(x$weights, 3), collapse = ", "),
             bwselect = x$bandwidth$method,
             h = if (!is.null(x$bandwidth[["h"]])) unname(x$bandwidth[["h"]]) else NA_real_,
             scheme = x$scheme, level = x$level,
             weighting = if (is.null(x$weighting)) "ols" else x$weighting,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' One-row tidy() data frame for a validation test
#' @keywords internal
#' @noRd
.tidy_test <- function(x, fn, assumption) {
  # df is a chi-squared degrees of freedom: defined for the Wald statistics only, not for sup-t
  # or sWald (whose critical values are simulated)
  wald <- identical(x$joint_method %||% "wald", "wald")
  data.frame(test = fn, assumption = assumption, statistic = x$statistic,
             df = if (wald) x$df else NA_integer_, p.value = x$p_value,
             method = .checks_method(x), row.names = NULL, stringsAsFactors = FALSE)
}
#' @rdname didisc-tidiers
#' @export
tidy.rd_typecont <- function(x, ...) .tidy_test(x, "rd_typecont", "type continuity")
#' @rdname didisc-tidiers
#' @export
tidy.rd_compstable <- function(x, ...) .tidy_test(x, "rd_compstable", "composition stability")
#' @rdname didisc-tidiers
#' @export
tidy.rd_homog <- function(x, ...) .tidy_test(x, "rd_homog", "homogeneous confounding")
#' @rdname didisc-tidiers
#' @export
tidy.rd_trendcell <- function(x, ...)
  .tidy_test(x, "rd_trendcell", paste(x$trend, "within-type confounding"))

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
                if (scheme_detected) " (detected)" else ""))
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
#' One line for a joint test under `method` ("supt", "swald" or "wald"). `label` is "Joint" or
#' a prefix such as "Period 1:"; the Wald line is .print_wald()'s ("Joint Wald chi-squared(...)").
#' @keywords internal
#' @noRd
.print_joint <- function(test, method, label = "Joint", indent = "  ", cell_label = NULL,
                         studentize = FALSE) {
  if (method == "wald")
    return(.print_wald(test$stat, test$df, test$p,
                       label = if (identical(label, "Joint")) "Joint Wald" else label,
                       indent = indent))
  name <- sprintf("%s (%s critical values)", if (method == "supt") "sup-t" else "sWald",
                  if (studentize) "t" else "normal")
  if (is.null(test$stat) || is.na(test$stat)) {
    cat(sprintf("%s%s %s: not testable\n", indent, label, name))
    return(invisible())
  }
  pf    <- if (test$p < 1e-3) "<0.001" else formatC(test$p, digits = 3, format = "f")
  jumps <- sprintf("%d jump%s", as.integer(test$df), if (test$df == 1L) "" else "s")
  if (method == "supt") {
    cat(sprintf("%s%s %s max |t| = %.3f over %s,  p = %s%s\n", indent, label, name, test$stat,
                jumps, pf, if (!is.null(cell_label) && test$df > 1L)
                             sprintf("  (largest: %s)", cell_label) else ""))
  } else {
    cat(sprintf("%s%s %s sum of t^2 = %.3f over %s,  p = %s\n", indent, label, name, test$stat,
                jumps, pf))
  }
}
#' Bandwidth description for the tests' printouts
#' @keywords internal
#' @noRd
.bw_label_test <- function(h, bwselect) {
  if (identical(bwselect, "pooled") && !is.null(h) && !is.na(h))
    sprintf("h = %.4g in every cell (CCT on the pooled comparison periods)", h)
  else if (!is.null(h) && !is.na(h)) sprintf("h = %.4g in every cell", h)
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
#' @rdname didisc-tidiers
#' @export
tidy.rd_checks <- function(x, ...) {
  out <- x$table
  names(out)[names(out) == "p_value"] <- "p.value"
  out
}
