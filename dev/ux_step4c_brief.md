# Step 4c brief — readability refactor of ONE assumption-test file (UX sweep, 2026-10-08)

You are refactoring exactly one file, `R/rd_<name>.R` (named in your task), for human readability. **Numerical behaviour must not change at all**, and the golden-master gate decides: your work is accepted only when it passes.

## What you may touch
- Only `R/rd_<name>.R`. Not `R/assumption_tests_helpers.R` (three other agents share it; the lead dedupes afterwards), not any other file, not tests, not docs.
- Inside the file: the function body, the print method, plain `#` comments, and roxygen only if a `@param` name must follow a local rename (there should be none — argument names are frozen).

## What must stay identical
- The signature: argument names, order, defaults, `match.arg` vectors (default first).
- The returned object: class, element names, element ORDER, types; every number, string, NA.
- Every `stop()`, `warning()`, `message()` text and the condition under which it fires.
- The arithmetic and its order: the sequence of statements that touch numbers, loop order (`for (tp ...) for (vt ...)`), `sort(..., method = "radix")`, `decreasing =`, tolerances, `tryCatch(..., error = function(e) NULL)` cell skipping, the reference-type rule, how `Sigma` is assembled (which index is row, which is column), `.wald_eigen` vs `.joint_wald` choice.
- Calls into shared helpers (`.build_types`, `.cell_bandwidth`, `.cross_cov`, `.cov_scheme`, `.match_sum`, `.detect_scheme`, `.scheme_from_long`, `.joint_wald`, `.wald_eigen`, `rd_period`, `rd_bw_cct`) with the same arguments.

## What to do (in this order, saving after each)
1. **Alias the cutoff.** First line of the body after the `match.arg` block: `cutoff <- c   # the cutoff; \`c\` stays the argument name for rdrobust users`. Replace every *numeric* use of `c` in the body by `cutoff` (comparisons `x >= c`, `c = c` passed to helpers becomes `c = cutoff`, `c_orig <- c` etc.). Leave the function-call uses `c(...)` and `base::c(...)` exactly as they are.
2. **Split the body into named steps.** The function is a pipeline: validate inputs → select periods → build types/cells → detect the scheme → fit every (period, type) cell → assemble the covariance → Wald → assemble the output. Extract each non-trivial step into a helper defined **in this file**, named `.<name>_<step>()` (e.g. `.homog_fit_cells()`), placed above the exported function, with a one-line `#' ... @noRd` roxygen saying what it returns. The exported function should read as a short sequence of those calls (target: under 80 lines). Each helper contains the **same statements in the same order** as the code it replaces; pass what it needs as arguments and return what the caller needs (a list if several things). Do not merge two loops, do not vectorise a loop, do not reorder independent statements "for clarity".
3. **One statement per line.** Break `a <- 1; b <- 2` lines; break expressions longer than ~100 characters at natural points without changing grouping (parentheses stay).
4. **Names.** Rename local variables that collide with scheme codes or are cryptic, consistently within the file: a p-value held in `pv` → `p_value`; `vt`/`ck` loop variables → `type_v` / `cell_k` (or keep if already clear); `bt` (the `.build_types()` result) → `types`; `b_t`, `b_v`, `b_idx` second-loop indices → `idx_t2`, `idx_v2`, `type_idx2`; `ew` → `wald`; `df` (degrees of freedom) → `wald_df` where it shadows `stats::df`; `sub`, `labels` where they shadow base functions. Do not rename anything that is part of the returned object or printed.
5. **Comments.** Delete comments that restate the code. Keep and, where thin, add one-line WHY comments at: the reference-type drop (why one type per period is dropped — the indicators sum to one), the `method = "radix"` sorts (locale independence), the `tryCatch` cell skipping (which failures it hides), the scheme's effect on `Sigma` (which covariance terms enter under cs/pc/pv), any `1e-8`/`eps` tolerance, and any place the code deliberately differs from the obvious (the file's existing good comments mark several). Dated history ("this pass", "was ...", "before 2026-..", "current behaviour", "original behavior") goes; it lives in NEWS.
6. **Print method.** Leave its output byte-identical; you may tidy its code.

## Acceptance (run all; paste the results in your report)
From the package root:
```
Rscript dev/snapshot_all.R /tmp/step4c_<name>.rds > /dev/null 2>&1 && Rscript dev/snapshot_compare.R dev/snapshots/baseline_step1.rds /tmp/step4c_<name>.rds | tail -3
Rscript -e 'devtools::test(filter = "<name>|estimand|mirror-invariance|s44-conformance|bwselect-paths", reporter = "summary")' 2>&1 | tail -5
Rscript -e 'devtools::document()' 2>&1 | grep -iv "resolve link|Writing|Loading|^i|^v" | head
```
The compare must print `GATE: PASS`; the tests must show no failures. If the gate fails, find the statement you changed, restore it, and re-run — never "fix" a number. Report: line count of the exported function before and after, the helpers created with one line each, every rename, and the acceptance output. Do not commit. End with DONE / DONE_WITH_CONCERNS / BLOCKED.
