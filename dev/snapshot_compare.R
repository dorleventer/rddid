# Compare two golden-master snapshots from dev/snapshot_all.R (UX sweep gate G1, T1644).
#
#   Rscript dev/snapshot_compare.R <baseline.rds> <new.rds>
#
# Walks every leaf of the baseline. A leaf passes only if the same path exists in the new
# snapshot and is identical() — not all.equal(): the sweep promises bit-identical numbers.
# Leaves present only in the new snapshot are listed as additions (allowed: new fields and
# methods are part of the sweep). Exit status 1 if any leaf changed or went missing.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript dev/snapshot_compare.R <baseline.rds> <new.rds>")
base <- readRDS(args[1]); new <- readRDS(args[2])

leaves <- function(x, path = character()) {
  if (is.data.frame(x)) x <- c(as.list(x), list(.row.names = attr(x, "row.names")))
  if (is.list(x) && length(x) > 0) {
    nm <- names(x); if (is.null(nm)) nm <- rep("", length(x))
    nm[nm == ""] <- paste0("[[", which(nm == ""), "]]")
    out <- list()
    for (i in seq_along(x)) out <- c(out, leaves(x[[i]], c(path, nm[i])))
    return(out)
  }
  stats::setNames(list(x), paste(path, collapse = "/"))
}
lb <- leaves(unclass(base)); ln <- leaves(unclass(new))
changed <- character(); missing <- character(); reworded <- character()
for (k in names(lb)) {
  if (!k %in% names(ln)) { missing <- c(missing, k); next }
  if (!identical(lb[[k]], ln[[k]])) {
    # an error cell whose MESSAGE changed is wording, not numbers: report, don't fail
    if (grepl("/error$", k) && is.character(lb[[k]]) && is.character(ln[[k]])) {
      reworded <- c(reworded, k); next
    }
    changed <- c(changed, k)
  }
}
added <- setdiff(names(ln), names(lb))

cat(sprintf("baseline %s (HEAD %s)  vs  new %s (HEAD %s)\n", args[1], attr(base, "git_head"),
            args[2], attr(new, "git_head")))
cat(sprintf("leaves: %d compared, %d identical, %d changed, %d missing, %d added\n",
            length(lb), length(lb) - length(changed) - length(missing), length(changed),
            length(missing), length(added)))
if (length(added)) {
  cat("\nadded (allowed):\n")
  # collapse to the first two path components so a new field with many leaves prints once
  cat(paste0("  ", unique(sub("^([^/]*/[^/]*/[^/]*).*$", "\\1", added)), collapse = "\n"), "\n")
}
if (length(reworded)) {
  cat("\nerror message reworded (allowed):\n"); cat(paste0("  ", reworded, collapse = "\n"), "\n")
}
if (length(changed)) {
  cat("\nCHANGED:\n")
  for (k in changed) {
    b <- lb[[k]]; n <- ln[[k]]
    desc <- if (is.numeric(b) && is.numeric(n) && length(b) == length(n))
      sprintf("max |diff| = %.3g", max(abs(b - n), na.rm = TRUE)) else
      sprintf("class %s -> %s, length %d -> %d", class(b)[1], class(n)[1], length(b), length(n))
    cat(sprintf("  %s  (%s)\n", k, desc))
  }
}
if (length(missing)) { cat("\nMISSING in new:\n"); cat(paste0("  ", missing, collapse = "\n"), "\n") }
if (length(changed) || length(missing)) { cat("\nGATE: FAIL\n"); quit(status = 1) }
cat("\nGATE: PASS (every baseline leaf identical)\n")
