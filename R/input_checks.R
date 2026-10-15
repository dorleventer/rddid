# input_checks.R -- argument checks shared by every exported function.

#' Stop if a named column is missing from `data`
#' @keywords internal
#' @noRd
.check_columns <- function(data, cols) {
  for (col in cols)
    if (!col %in% names(data))
      stop("column '", col, "' not found in `data`.", call. = FALSE)
  invisible(TRUE)
}

#' Validate and prepare the input data of an exported function
#'
#' One place for the checks every entry point needs: `data` is a data frame with the named
#' columns, `x` (and `y`) numeric, the period column numeric, character or factor; `t_rd` and
#' `comparisons` exist, do not overlap and do not repeat; each unit appears once per period;
#' rows with a missing value in a used column are dropped, with a message. Errors carry the
#' call of the exported function (`Error in rd_homog(...)`), so they name it once; the
#' missing-value message starts with `fn` (e.g. `"rddid()"`, `"rd_homog"`).
#'
#' @param fn the calling function's label for the message.
#' @param complete_cols the columns a row must have to be kept (default: every used column).
#' @param y,x,time,id column names; `y` and `id` may be `NULL`.
#' @param t_rd,comparisons the periods, as passed (factors become character).
#' @param t_rd_required stop if `t_rd` is `NULL`.
#' @param t_rd_scalar stop if `t_rd` has more than one period (only `rddid()` takes several).
#' @param default_comparisons if `comparisons` is `NULL`: `TRUE` uses every other period
#'   (every period when `t_rd` is `NULL`); `FALSE` leaves it `NULL`.
#' @param restrict keep only the rows of `t_rd` and `comparisons` (before the missing-value
#'   drop, so the message counts only those rows).
#' @return list: `data` (checked, possibly restricted, complete rows), `t_rd`, `comparisons`.
#' @keywords internal
#' @noRd
.prepare_input <- function(data, fn, x, time, id = NULL, y = NULL, t_rd = NULL,
                           comparisons = NULL, t_rd_required = TRUE, t_rd_scalar = TRUE,
                           default_comparisons = TRUE, restrict = FALSE,
                           complete_cols = base::c(y, x, time, id)) {
  caller <- sys.call(-1L)   # the exported function's call, shown in the error header
  err <- function(...) stop(simpleError(paste0(...), caller))
  if (!is.data.frame(data)) err("`data` must be a data frame.")
  .check_columns(data, base::c(y, x, time, id))
  if (!is.null(y) && !is.numeric(data[[y]]) && !is.logical(data[[y]]))
    err("the outcome column '", y, "' must be numeric (or logical).")
  if (!is.numeric(data[[x]])) err("the running variable column '", x, "' must be numeric.")
  time_col <- data[[time]]
  if (!is.numeric(time_col) && !is.character(time_col) && !is.factor(time_col))
    err("the period column `", time, "` must be numeric, character or factor (a Date can be ",
        "converted with as.numeric(format(x, \"%Y\")) or similar).")

  # ---- the RD period(s) ----
  if (is.factor(t_rd)) t_rd <- as.character(t_rd)
  if (is.null(t_rd)) {
    if (t_rd_required) err("`t_rd` (the RD period) is required.")
  } else {
    if (length(t_rd) < 1L) err("`t_rd` must give at least one RD period.")
    if (t_rd_scalar && length(t_rd) > 1L)
      err("`t_rd` must be a single period (several RD periods are supported by rddid() only).")
    if (anyDuplicated(t_rd)) err("`t_rd` contains a period more than once.")
    missing_rd <- t_rd[!t_rd %in% time_col]
    if (length(missing_rd))
      err("`t_rd` = ", paste(missing_rd, collapse = ", "), " is not a period in `", time, "`.")
  }

  # ---- the comparison periods ----
  if (is.factor(comparisons)) comparisons <- as.character(comparisons)
  if (is.null(comparisons) && default_comparisons)
    comparisons <- sort(setdiff(unique(time_col), t_rd))
  if (!is.null(comparisons) || default_comparisons) {
    if (length(comparisons) < 1L) err("need at least one comparison period.")
    if (!is.null(t_rd) && any(t_rd %in% comparisons))
      err("`t_rd` must not be one of the `comparisons`.")
    if (anyDuplicated(comparisons)) err("`comparisons` contains a period more than once.")
    missing_comp <- setdiff(comparisons, time_col)
    if (length(missing_comp))
      err("comparison period(s) not in `", time, "`: ", paste(missing_comp, collapse = ", "))
  }

  # ---- rows used: the requested periods (if restricted), complete in the used columns ----
  if (restrict) data <- data[!is.na(time_col) & time_col %in% base::c(t_rd, comparisons), , drop = FALSE]
  cols <- complete_cols
  keep <- stats::complete.cases(data[, cols, drop = FALSE])
  if (any(!keep)) {
    message(fn, ": ", sum(!keep), " row(s) with a missing value in ",
            paste0("`", cols, "`", collapse = ", "), " dropped.")
    data <- data[keep, , drop = FALSE]
  }
  empty <- setdiff(base::c(t_rd, comparisons), data[[time]])
  if (length(empty))
    err("no complete row in period(s) ", paste(empty, collapse = ", "), ": every row there has a ",
        "missing value in ", paste0("`", cols, "`", collapse = ", "), ".")

  # ---- each unit once per period, among the rows used ----
  if (!is.null(id)) {
    dup <- duplicated(data[, base::c(id, time)])
    if (any(dup))
      err(sum(dup), " row(s) repeat a (", id, ", ", time, ") pair: each unit may appear once ",
          "per period (aggregate or deduplicate first).")
  }
  list(data = data, t_rd = t_rd, comparisons = comparisons)
}
