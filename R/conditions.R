# conditions.R -- errors and warnings that name the function the user called.
# Every user-facing error goes through .abort(): its header is the user's own call
# ("Error in didisc(...)"), wherever in the package the check sits, so no internal helper's
# name reaches the user and no message needs a "fn: " prefix.

#' The outermost call to an exported didisc function on the call stack (or NULL)
#' @keywords internal
#' @noRd
.user_call <- function() {
  exported <- getNamespaceExports("didisc")
  for (cl in sys.calls()) {
    fn <- cl[[1L]]
    if (is.call(fn) && identical(fn[[1L]], as.name("::"))) fn <- fn[[3L]]
    if (is.name(fn) && as.character(fn) %in% exported) return(cl)
    if (is.function(fn)) {                 # do.call(didisc, ...): name the function, not its data
      ns  <- asNamespace("didisc")
      hit <- Filter(function(e) identical(fn, get(e, envir = ns)), exported)
      if (length(hit)) return(as.call(list(as.name(hit[[1L]]))))
    }
  }
  NULL
}

#' Stop with an error whose header is the user's call
#' @keywords internal
#' @noRd
.abort <- function(...) stop(simpleError(paste0(...), .user_call()))

#' Warn with a warning whose header is the user's call
#' @keywords internal
#' @noRd
.warn <- function(...) warning(simpleWarning(paste0(...), .user_call()))

#' match.arg() with an error that names the argument, its choices and the value given
#'
#' `choices` defaults to the calling function's default for the argument, as in match.arg();
#' `name` is the argument's name in messages (pass it when `arg` is an expression).
#' @keywords internal
#' @noRd
.match_arg <- function(arg, choices, name = deparse(substitute(arg))) {
  if (missing(choices)) choices <- eval(formals(sys.function(sys.parent()))[[name]])
  if (identical(arg, choices)) return(arg[[1L]])
  i <- if (is.character(arg) && length(arg) == 1L && !is.na(arg)) pmatch(arg, choices) else NA
  if (is.na(i))
    .abort("`", name, "` must be one of ", paste0("\"", choices, "\"", collapse = ", "),
           " (got ", paste(deparse(arg), collapse = ""), ").")
  choices[[i]]
}
