# dev/ — developer notes for rddid

Not part of the package build (`.Rbuildignore`). Contents: `appB_map.md` and `tests_map.md` (code ↔ paper maps), `check_appB_labels.R`, `snapshot_rddid.R`, `site_plan.md` and `site_dgp_check.R` (pkgdown site rebuild, 2026-09).

## Pre-commit hook

After cloning, enable the doc-sync pre-commit hook (regenerates `man/*.Rd` from roxygen and blocks commits where the generated docs are stale — the mismatch that otherwise fails `R CMD check`):

``` sh
git config core.hooksPath .githooks
```
