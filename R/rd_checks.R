# rd_checks.R -- rd_checks(): the four tests of the identification assumptions in one call, with one
# set of options, and one table; print and tidy methods.

#' Run the four tests of the identification assumptions
#'
#' When the running variable moves over time, the estimate of [rddid()] is
#' identified under type continuity, constant within-type confounding, and
#' either composition stability or homogeneous confounding (Leventer and Nevo,
#' 2024). `rd_checks()` runs the four tests, [rd_typecont()],
#' [rd_compstable()], [rd_homog()] and [rd_trendcell()], with one set of
#' options, and prints one table. Each test's own result is kept, so every
#' detail stays available.
#'
#' @details
#' The arguments shared by the four tests (`t_rd`, `comparisons`, `estimand`,
#' `c`, `kernel`, `scheme`, `bc`) are set once. `joint` and `studentize` apply
#' to the two tests on type shares, `outcome_bwselect` to the two tests on the
#' outcome; see the help of each test for what they do. Anything else, or a
#' different value for one test only, goes in that test's list, for example
#' `typecont = list(joint = "wald")` or `homog = list(bwselect = "cct", min_n = 20)`.
#' A test that cannot run on the data (for example, no unit changes side of
#' the cutoff) gives a row with its error message and `NA` results.
#'
#' The printout ends with the reading at the level `level`: which assumptions
#' a test rejects. Not rejecting an assumption does not establish it; constant
#' within-type confounding concerns the RD period and is checked in the
#' comparison periods only (a pre-trends check).
#'
#' @param data,y,x,time,id as in [rddid()]; `id` is required.
#' @param t_rd the RD period (required).
#' @param comparisons,estimand,c,kernel,scheme,bc shared by the four tests;
#'   defaults as in each test.
#' @param joint the joint statistic of [rd_typecont()] and [rd_compstable()]:
#'   `"supt"` (default), `"swald"` or `"wald"` (the chi-squared Wald).
#' @param studentize for `joint = "supt"` or `"swald"`: `TRUE` (default; t
#'   critical values) or `FALSE` (normal critical values).
#' @param outcome_bwselect the bandwidth rule of [rd_homog()] and
#'   [rd_trendcell()]: `"common"` (default; one bandwidth for every cell, the
#'   CCT bandwidth of the outcome on the pooled comparison periods), `"cct"`
#'   (one CCT bandwidth per cell) or `"rot"` (a rule of thumb per cell). The
#'   share tests' own `bwselect` (per-fit CCT by default) can be changed in
#'   `typecont` / `compstable`; their rule of thumb is a different one
#'   (`0.5 * IQR(x)` in every regression).
#' @param typecont,compstable,homog,trendcell lists of further arguments for
#'   one test, which override the shared ones.
#' @param level the level of the reading in the printout (default 0.05).
#' @return An object of class `"rd_checks"`: a list with `table` (one row per
#'   test: `test`, `assumption`, `statistic`, `p_value`, `method`, `bandwidth`,
#'   `error`), the four test results `typecont`, `compstable`, `homog`,
#'   `trendcell` (`NULL` where a test failed), `level` and `call`.
#' @seealso [rd_typecont()], [rd_compstable()], [rd_homog()], [rd_trendcell()];
#'   `tidy()` in [rddid-tidiers].
#' @family tests of the assumptions
#' @examples
#' ch <- rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' ch
#' # switch options: the chi-squared Wald and per-cell bandwidths everywhere
#' rd_checks(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
#'           joint = "wald", outcome_bwselect = "cct")
#' @export
rd_checks <- function(data, y, x, time, id, t_rd,
                      comparisons = NULL,
                      estimand = c("att", "atu"),
                      c = 0,
                      kernel = "triangular",
                      scheme = c("auto", "cs", "pc", "pv"),
                      bc = TRUE,
                      joint = c("supt", "swald", "wald"),
                      studentize = TRUE,
                      outcome_bwselect = c("common", "cct", "rot"),
                      typecont = list(), compstable = list(), homog = list(), trendcell = list(),
                      level = 0.05) {
  cl               <- match.call()
  estimand         <- match.arg(estimand)
  scheme           <- match.arg(scheme)
  joint            <- match.arg(joint)
  outcome_bwselect <- match.arg(outcome_bwselect)
  if (missing(t_rd)) stop("rd_checks: `t_rd` (the RD period) is required.", call. = FALSE)

  shared <- list(data = data, x = x, time = time, id = id, t_rd = t_rd,
                 comparisons = comparisons, estimand = estimand, c = c, kernel = kernel,
                 scheme = scheme, bc = bc)
  share  <- list(joint = joint, studentize = studentize)
  outc   <- list(y = y, bwselect = outcome_bwselect)
  .check_studentize(studentize, "rd_checks")
  run <- function(fun, args) tryCatch(do.call(fun, args), error = function(e) e)
  # a test's own list replaces the shared values (not merged: `data` may itself be a list)
  over <- function(base, own) { base[names(own)] <- own; base }
  res <- list(
    typecont   = run(rd_typecont,   over(c(shared, share), typecont)),
    compstable = run(rd_compstable, over(c(shared, share), compstable)),
    homog      = run(rd_homog,      over(c(shared, outc), homog)),
    trendcell  = run(rd_trendcell,  over(c(shared, outc), trendcell)))

  labels <- c(typecont = "type continuity", compstable = "composition stability",
              homog = "homogeneous confounding", trendcell = "constant within-type confounding")
  rows <- lapply(names(res), function(k) {
    r <- res[[k]]
    if (inherits(r, "error"))
      return(data.frame(test = paste0("rd_", k), assumption = labels[[k]], statistic = NA_real_,
                        p_value = NA_real_, method = NA_character_, bandwidth = NA_character_,
                        error = sub("^(.+?\\.)(\\s.*)?$", "\\1", conditionMessage(r), perl = TRUE),
                        stringsAsFactors = FALSE))
    data.frame(test = paste0("rd_", k), assumption = labels[[k]], statistic = r$statistic,
               p_value = r$p_value, method = .checks_method(r), bandwidth = .checks_bw(r),
               error = NA_character_, stringsAsFactors = FALSE)
  })
  structure(list(table = do.call(rbind, rows),
                 typecont   = if (inherits(res$typecont, "error")) NULL else res$typecont,
                 compstable = if (inherits(res$compstable, "error")) NULL else res$compstable,
                 homog      = if (inherits(res$homog, "error")) NULL else res$homog,
                 trendcell  = if (inherits(res$trendcell, "error")) NULL else res$trendcell,
                 level = level, call = cl),
            class = "rd_checks")
}

#' The statistic of a test result, in words (for the rd_checks table)
#' @keywords internal
#' @noRd
.checks_method <- function(r) {
  m <- r$joint_method %||% "wald"
  if (m == "wald") return(sprintf("Wald chi-squared(%d)", as.integer(r$df)))
  sprintf("%s, %s critical values", if (m == "supt") "sup-t" else "sWald",
          if (isTRUE(r$studentize)) "t" else "normal")
}

#' The bandwidth rule of a test result, in words (for the rd_checks table)
#' @keywords internal
#' @noRd
.checks_bw <- function(r) {
  if (!is.null(r$bandwidth)) {
    return(switch(r$bandwidth$method,
                  common = sprintf("common (pooled CCT), h = %.3g", r$bandwidth$h),
                  fixed  = sprintf("fixed, h = %.3g", r$bandwidth$h),
                  cct    = "CCT per cell",
                  rot    = "rule of thumb per cell"))
  }
  if (identical(r$meta$bwselect, "rot")) return(sprintf("rule of thumb, h = %.3g", r$meta$h))
  if (!is.null(r$meta$h) && !is.na(r$meta$h)) return(sprintf("fixed, h = %.3g", r$meta$h))
  "CCT per cell"
}

#' @export
print.rd_checks <- function(x, ...) {
  cat("Tests of the identification assumptions  [rd_checks()]\n\n")
  tab <- x$table
  fmt_p <- function(p) ifelse(is.na(p), "     NA", ifelse(p < 1e-3, " <0.001", formatC(p, digits = 3, format = "f", width = 7)))
  cat(sprintf("  %-34s %10s %8s   %-26s %s\n", "Assumption", "Statistic", "p", "Statistic type", "Bandwidth"))
  for (i in seq_len(nrow(tab))) {
    if (!is.na(tab$error[i])) {
      cat(sprintf("  %-34s not run: %s\n", tab$assumption[i], tab$error[i]))
      next
    }
    cat(sprintf("  %-34s %10s %8s   %-26s %s\n", tab$assumption[i],
                if (is.na(tab$statistic[i])) "NA" else formatC(tab$statistic[i], digits = 3, format = "f"),
                fmt_p(tab$p_value[i]), tab$method[i], tab$bandwidth[i]))
  }
  rej <- function(k) { p <- tab$p_value[tab$test == paste0("rd_", k)]
    if (!length(p) || is.na(p)) NA else p < x$level }
  word <- function(r) if (is.na(r)) "not tested" else if (r) "rejected" else "not rejected"
  cat(sprintf("\n  Reading at the %g%% level (not rejecting does not establish an assumption):\n",
              100 * x$level))
  cat(sprintf("    %-34s %s\n", "type continuity", word(rej("typecont"))))
  cat(sprintf("    %-34s %-13s} identification needs\n", "composition stability", word(rej("compstable"))))
  cat(sprintf("    %-34s %-13s} only one of these two\n", "homogeneous confounding", word(rej("homog"))))
  cat(sprintf("    %-34s %s\n", "constant within-type confounding",
              paste(word(rej("trendcell")), "(checked in the comparison periods)")))
  cat("  Each test's details: x$typecont, x$compstable, x$homog, x$trendcell.\n")
  invisible(x)
}
