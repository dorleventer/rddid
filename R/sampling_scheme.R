# sampling_scheme.R -- detection of the sampling scheme from the (period, id, side) structure:
# "cs" (no unit repeats), "pc" (units repeat, never change side), "pv" (some unit changes side).
# Used by rddid(), rd_homog(), rd_trendcell() and rd_typecont(); rd_compstable() has its own
# pair-level rule (a unit above the cutoff in both periods of a pair).

#' Classify the sampling scheme from a long (period, id, side) table
#'
#' No repeated unit across periods → `"cs"`; repeated units that switch side at
#' least once → `"pv"`; repeated units that never switch → `"pc"`.
#' @keywords internal
#' @noRd
.scheme_from_long <- function(long) {
  rep_ids <- names(which(table(unique(long[, c("period", "id")])$id) >= 2L))
  if (length(rep_ids) == 0L) return("cs")
  sub <- long[long$id %in% rep_ids, , drop = FALSE]
  switches <- tapply(sub$side, sub$id, function(s) length(unique(s)) > 1L)
  if (any(switches)) "pv" else "pc"
}

#' Detect the sampling scheme from the id / side structure
#' @keywords internal
#' @noRd
.detect_scheme <- function(plist, c = 0) {
  # side is 1 if x >= c (treated at the cutoff), 0 otherwise — no third "side"
  # at exactly x == c, so a unit sitting on the cutoff does not read as a switch.
  long <- do.call(rbind, lapply(names(plist), function(k)
    data.frame(period = k, id = plist[[k]]$id,
               side = as.integer(plist[[k]]$x >= c))))
  .scheme_from_long(long)
}
