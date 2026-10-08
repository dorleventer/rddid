# dev/ — developer notes for rddid

**2026-10-08:** the composition-adjusted family (`rd_att`, `rd_sadjust`, `rd_c`, `rd_adjust`) was removed from the package (UX sweep, T1644); references to it in `appB_map.md`, `site_plan.md`, `atu_estimand_plan.md` are historical. Recover at tag `v0.4.0.9000-composition`.

Not part of the package build (`.Rbuildignore`).

**Numerical gate (UX sweep, 2026-10-08).** `dev/snapshot_all.R` records the full return object of every exported function on seeded data under every option (134 cells) to an `.rds`; `dev/snapshot_compare.R <baseline> <new>` requires every baseline leaf to be `identical()` (new leaves allowed, reworded error messages allowed). Current baseline: `dev/snapshots/baseline_ux2.rds` (2026-10-08, after the review round: the no-switcher panel cells of the four tests now error). Earlier ones kept for history: `baseline_bugfix.rds` (after the five bug fixes; differs from `baseline_step1.rds` only in the `rd_homog(p = 2)` cells). Run both before and after any change to `R/`.

Contents: `appB_map.md` and `tests_map.md` (code ↔ paper maps), `check_appB_labels.R`, `snapshot_rddid.R`, `site_plan.md` and `site_dgp_check.R` (pkgdown site rebuild, 2026-09), `atu_estimand_plan.md` (the `estimand` argument, 2026-09-10).

## Pre-commit hook

After cloning, enable the doc-sync pre-commit hook (regenerates `man/*.Rd` from roxygen and blocks commits where the generated docs are stale — the mismatch that otherwise fails `R CMD check`):

``` sh
git config core.hooksPath .githooks
```
