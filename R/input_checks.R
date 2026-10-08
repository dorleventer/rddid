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
