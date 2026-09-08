# Check that every paper label cited in dev/appB_map.md and in R/*.R still exists
# in the paper's main.tex. Labels get renamed during line-passes; a stale citation
# in the map or in a code comment is silent drift.
#
# Usage (from the package root):
#   Rscript dev/check_appB_labels.R [path/to/main.tex]
# Default path: ../rd-did/writing/main_revamp/main.tex (sibling repo).
# Exit status 1 if any cited label is missing.

args <- commandArgs(trailingOnly = TRUE)
tex_path <- if (length(args) >= 1L) args[1L] else
  file.path("..", "rd-did", "writing", "main_revamp", "main.tex")
if (!file.exists(tex_path)) stop("main.tex not found at: ", tex_path)

tex <- readLines(tex_path, warn = FALSE)
defined <- unique(regmatches(tex, gregexpr("\\\\label\\{[^}]+\\}", tex)) |> unlist())
defined <- sub("^\\\\label\\{(.*)\\}$", "\\1", defined)

# Cited labels: prefix:name tokens, prefixes used by the paper.
lab_re <- "\\b(eq|lem|prop|thm|cor|ass|sec|app|alg|tab|fig):[A-Za-z0-9_-]+"
cite_in <- function(files) {
  out <- character(0)
  for (f in files) {
    x <- readLines(f, warn = FALSE)
    m <- regmatches(x, gregexpr(lab_re, x, perl = TRUE))
    hits <- unlist(m)
    if (length(hits)) out <- rbind(out, cbind(file = basename(f), label = hits))
  }
  if (length(out)) as.data.frame(out, stringsAsFactors = FALSE) else
    data.frame(file = character(0), label = character(0))
}
cited <- cite_in(c(list.files("dev", full.names = TRUE, pattern = "\\.md$"),
                   list.files("R", full.names = TRUE, pattern = "\\.R$")))
cited <- unique(cited)

# Labels the map lists as deliberately removed are exempt when cited only in the map's
# "must not be cited" sentence; flag them anywhere else.
removed <- c("lem:coercive", "eq:amse-ps")
missing <- cited[!cited$label %in% defined, , drop = FALSE]
missing <- missing[!(missing$label %in% removed & grepl("_map\\.md$", missing$file)), , drop = FALSE]

cat(sprintf("main.tex: %d labels defined; %d distinct labels cited in dev/*.md + R/\n",
            length(defined), length(unique(cited$label))))
if (nrow(missing)) {
  cat("MISSING (cited but not defined in main.tex):\n")
  print(missing, row.names = FALSE)
  quit(status = 1L)
}
cat("OK: every cited label exists in main.tex\n")
