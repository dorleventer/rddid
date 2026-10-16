# input_checks.R -- argument checks shared by every exported function.

#' Stop if a named column is missing from `data`
#' @keywords internal
#' @noRd
.check_columns <- function(data, cols) {
  for (col in cols)
    if (!col %in% names(data))
      .abort("column '", col, "' not found in `data`.")
  invisible(TRUE)
}

#' Validate and prepare the input data of an exported function
#'
#' One place for the checks every entry point needs: `data` is a data frame with the named
#' columns, `x` (and `y`) numeric, the period column numeric, character or factor; `t_rd` and
#' `comparisons` exist, do not overlap and do not repeat; each unit appears once per period.
#' Only the rows of `t_rd` and `comparisons` are kept, and of those the rows complete in
#' `complete_cols` (dropped with a message). Periods given as text or factor are read in the
#' period column's type. Errors name the user's call (.abort()); the message starts with `fn`.
#'
#' @param fn the calling function's label for the message, e.g. `"rd_homog()"`.
#' @param y,x,time,id column names; `y` and `id` may be `NULL`.
#' @param t_rd,comparisons the periods, as passed; `comparisons = NULL` is every other period
#'   (in time order; in level order for a factor period column).
#' @param t_rd_required stop if `t_rd` is `NULL` (only [rd_typecont()] runs without it, on
#'   `comparisons`, by default every period).
#' @param t_rd_scalar stop if `t_rd` has more than one period (only `rddid()` takes several).
#' @param complete_cols the columns a row must have to be kept (default: every used column).
#' @return list: `data` (the rows used), `t_rd`, `comparisons`.
#' @keywords internal
#' @noRd
.prepare_input <- function(data, fn, x, time, id = NULL, y = NULL, t_rd = NULL,
                           comparisons = NULL, t_rd_required = TRUE, t_rd_scalar = TRUE,
                           complete_cols = base::c(y, x, time, id)) {
  if (!is.data.frame(data)) .abort("`data` must be a data frame.")

  # ---- columns: present, and of a usable type ----
  given  <- base::c(y = y, x = x, time = time, id = id)
  absent <- given[!given %in% names(data)]
  if (length(absent))
    .abort("not a column of `data`: ", paste0("`", names(absent), "` = \"", absent, "\"",
                                             collapse = ", "),
           ". The columns are ", .list_values(names(data)), ".")
  if (!is.null(y) && !is.numeric(data[[y]]) && !is.logical(data[[y]]))
    .abort("the outcome column `y` = \"", y, "\" must be numeric", .type_hint(data[[y]], y), ".")
  if (!is.numeric(data[[x]]))
    .abort("the running variable column `x` = \"", x, "\" must be numeric",
           .type_hint(data[[x]], x), ".")
  time_col <- data[[time]]
  if (!is.numeric(time_col) && !is.character(time_col) && !is.factor(time_col))
    .abort("the period column `time` = \"", time, "\" must be numeric, character or factor (a ",
           "Date converts with as.numeric(format(x, \"%Y\")) or similar).")
  periods_in_data <- if (is.factor(time_col)) levels(droplevels(time_col)) else
    sort(unique(time_col[!is.na(time_col)]))
  # periods given as text or factor are read in the period column's type
  as_period <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.numeric(time_col) && is.character(v)) {
      num <- suppressWarnings(as.numeric(v))
      if (!anyNA(num)) v <- num
    }
    v
  }

  # ---- the RD period(s) ----
  t_rd <- as_period(t_rd)
  if (is.null(t_rd) && t_rd_required) .abort("`t_rd` (the RD period) is required.")
  if (!is.null(t_rd)) {
    if (length(t_rd) < 1L) .abort("`t_rd` must give at least one RD period.")
    if (t_rd_scalar && length(t_rd) > 1L)
      .abort("`t_rd` must be a single period here (several RD periods are for rddid()).")
    if (anyDuplicated(t_rd)) .abort("`t_rd` contains a period more than once.")
    missing_rd <- t_rd[!t_rd %in% time_col]
    if (length(missing_rd))
      .abort("`t_rd` = ", paste(missing_rd, collapse = ", "), " is not a period in `", time,
             "` (its periods: ", .list_values(periods_in_data), ").")
  }

  # ---- the comparison periods: by default every other period ----
  comparisons <- if (is.null(comparisons)) setdiff(periods_in_data, t_rd) else as_period(comparisons)
  if (length(comparisons) < 1L) .abort("need at least one comparison period.")
  overlap <- intersect(t_rd, comparisons)
  if (length(overlap))
    .abort("period ", paste(overlap, collapse = ", "), " is in both `t_rd` and `comparisons`.")
  if (anyDuplicated(comparisons)) .abort("`comparisons` contains a period more than once.")
  missing_comp <- setdiff(comparisons, time_col)
  if (length(missing_comp))
    .abort("comparison period(s) ", paste(missing_comp, collapse = ", "), " not in `", time,
           "` (its periods: ", .list_values(periods_in_data), ").")

  # ---- rows used: those of the periods used, complete in the columns used ----
  data <- data[!is.na(time_col) & time_col %in% base::c(t_rd, comparisons), , drop = FALSE]
  keep <- stats::complete.cases(data[, complete_cols, drop = FALSE])
  if (any(!keep)) {
    na_cols <- complete_cols[vapply(complete_cols, function(cc) anyNA(data[[cc]]), logical(1))]
    message(fn, ": ", sum(!keep), " row(s) with a missing ",
            paste0("`", na_cols, "`", collapse = " or "), " dropped.")
    data <- data[keep, , drop = FALSE]
  }
  empty <- setdiff(base::c(t_rd, comparisons), data[[time]])
  if (length(empty))
    .abort("no complete row in period(s) ", paste(empty, collapse = ", "), ": every row there has ",
           "a missing value in ", paste0("`", complete_cols, "`", collapse = ", "), ".")

  # ---- each unit once per period, among the rows used ----
  if (!is.null(id)) {
    dup <- duplicated(data[, base::c(id, time)])
    if (any(dup))
      .abort(sum(dup), " row(s) repeat a (", id, ", ", time, ") pair: each unit may appear once ",
             "per period (aggregate or deduplicate first).")
  }
  list(data = data, t_rd = t_rd, comparisons = comparisons)
}

#' Up to `max` values, comma-separated, with "..." when there are more
#' @keywords internal
#' @noRd
.list_values <- function(v, max = 12L) {
  v <- as.character(v)
  if (length(v) > max) paste0(paste(v[seq_len(max)], collapse = ", "), ", ...")
  else paste(v, collapse = ", ")
}

#' "; it is <type> (<how to convert>)" for a column that should be numeric
#' @keywords internal
#' @noRd
.type_hint <- function(v, col) {
  if (is.factor(v))
    paste0("; it is a factor (as.numeric(as.character(data$", col, ")) converts its labels, ",
           "while as.numeric() alone gives the level codes)")
  else if (is.character(v)) paste0("; it is character (as.numeric(data$", col, ") converts it)")
  else paste0("; it is ", class(v)[1L])
}
