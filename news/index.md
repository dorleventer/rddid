# Changelog

## rddid 0.4.0.9000 (development)

### 2026-10-08: `tidy()`/`glance()` always available; a no-Suggests check

- `generics` moved from Suggests to Imports and its
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)
  are re-exported: `tidy(fit)` works after
  [`library(rddid)`](https://github.com/dorleventer/rddid) without
  `generics::` and without installing anything else. No numerical
  change.
- New workflow `check-no-suggests.yaml`: R CMD check with only the hard
  dependencies installed (plus testthat, knitr, rmarkdown), so every use
  of ggplot2, modelsummary or broom stays guarded. `broom` added to
  Suggests (the Get-started modelsummary chunk needs it).

### 2026-10-08: second bug-hunt round (input handling; no numerical change on valid input)

- [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  now works on complete (outcome, running variable, id) rows only and
  says how many it dropped: an NA running variable used to read as a
  side switch, turning a “no unit changes side” panel into “some units
  change side” and moving the common bandwidth and the estimate;
  `n_by_period` and [`nobs()`](https://rdrr.io/r/stats/nobs.html) now
  count the rows actually used.
- Errors instead of silent misuse: a unit appearing twice in a period
  ([`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  and the four tests), a period repeated in `comparisons`, `comparisons`
  that include `t_rd` or periods not in the data (all four tests), `b`
  without `h`, an `h` that is not a single positive number, a `t_rd` not
  in the data
  ([`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)),
  `type_by = "pattern"` in
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  (a cell fixed across the comparison periods cannot use a comparison
  period’s own side).
- Named numeric `trend` weights are matched to the comparison periods by
  name; `q` now reaches
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  (the pilot bandwidth used to be chosen for `q = p + 1` whatever `q`
  was); `kernel` accepts any case; a factor `t_rd` works.
- Untestable results (no testable contrast, or a type share whose jump
  has a zero standard error) are reported as `statistic = NA`, `df = 0`,
  `p = NA` instead of a chi-squared of 0 with p = 1 or a rounding-noise
  rejection; the joint composition-stability test skips untestable
  pairs.
- The no-switcher guard looks only at the periods a test uses, and names
  a repeated cross-section for what it is.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  `rd_homog`/`rd_trendcell` works when `t_rd` was passed as a variable
  (the objects now carry `t_rd`); `plot(cs, pair = )` checks its index;
  [`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md)
  errors on a cross-section; more than five types get distinct colours.
- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)/[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  choose a cell’s bandwidth only after the `min_n` check (no fallback
  noise from skipped cells);
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  prints the skipped cells before it stops for lack of contrasts; its
  `reference` flag is set only in periods that contributed a contrast;
  comparison periods print in time order.
- Tests: console snapshots (`tests/testthat/_snaps/`) and one regression
  test per fix.
- Also: [`tidy()`](https://generics.r-lib.org/reference/tidy.html) of a
  fit honours `conf.int` and `conf.level`; DESCRIPTION cites the paper
  and has a shorter Title; `CITATION.cff`; lifecycle badge; the
  How-it-works article is guarded on ggplot2; rdrobust skips removed
  from the tests.

### 2026-10-08: review round on the plots, the console and the site

- **Figures** (after a review against the paper’s own figures): no
  titles, subtitles or statistics in any figure; the paper’s palette and
  `theme_bw`; bins never straddle the cutoff; share plots clip at 0–1
  without dropping line segments; the composition-stability axis shows
  the running variable on both sides (folded ticks) and marks the
  mirrored period correctly under `estimand = "atu"`;
  [`plot.rd_trendcell()`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md)
  draws the least-squares line under `trend = "linear"`; switcher
  colours follow the side in the second period; periods in time order
  everywhere; legends keep no internal names when a theme is replaced.
- **Console**: [`summary()`](https://rdrr.io/r/base/summary.html)’s
  per-period table is aligned with its (time-ordered) labels — the
  previous commit had scrambled the numbers;
  [`glance()`](https://generics.r-lib.org/reference/glance.html) is one
  row under every bandwidth rule (`$h` used to partial-match
  `h_by_period`); shorter scheme labels;
  [`summary()`](https://rdrr.io/r/base/summary.html) no longer repeats
  the print hint; `(detected)` only when the scheme was auto-detected;
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  checks `comparisons`, `level` and the period column’s type and says
  so.
- **Degenerate panels**: when no unit changes side of the cutoff between
  the periods (every unit has the same type in every period,
  e.g. `rddid_sim`), the four tests now stop with one clear message
  instead of returning a mechanical rejection, a `chi-squared(0)` or
  internal errors. Every other numerical result is unchanged.
- **Site and help pages**: *Get started* reads the four test printouts
  line by line, says what the period column may be (`trend = "linear"`
  needs numeric periods), estimates on `rddid_sim_pv` before testing its
  assumptions, defines a unit’s type once in a table, and shows how to
  export a table; the ATU section moved after the tests, and “What to
  report” became “What the paper reports”. Advice to the researcher is
  rewritten as description across the articles and help pages, with “can
  be biased” throughout; the joint test over pairs of
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  is flagged as approximate wherever it is tabulated. Articles menu: the
  tests, the plots (now *Plots of the estimate and the checks*, with a
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  example), the options, then *How rddid() computes the estimate* (for
  referees). The README links the tests and plots articles and
  `citation("rddid")`.

### 2026-10-08: plots

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods,
  mirroring the validation figures of the paper’s application:
  `plot(fit)` draws the per-period RD plots behind an
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  estimate (binned outcome, local-linear fits on each side, the jump at
  the cutoff); [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  of an `rd_typecont` or `rd_compstable` result draws the binned share
  of each type with its fitted lines on both sides of the (artificial)
  cutoff; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an
  `rd_homog` or `rd_trendcell` result draws the within-type confounding
  jumps with 95% intervals by period.
  [`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md)
  draws the running variable in one period against another and counts
  the units that change side. All return ggplot objects (ggplot2 in
  Suggests). New article *Plots of the estimate and the checks*.
- To feed the plots, the objects now also carry what the methods read:
  `rddid()$data` (the per-period data used),
  `rd_typecont()$fits`/`$data`, each
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  pair’s `fits` and `sample`, and `rd_typecont()$meta$c`. Nothing else
  changed.

### 2026-10-08: five bug fixes found during the UX sweep

- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  /
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  with `p > 1`: the per-cell CCT bandwidth is now chosen for the
  order-`p` fit (it used to be chosen at `p = 1` whatever `p` the fit
  used). Default calls (`p = 1`) are unchanged.
- `rd_trendcell(trend = "linear")`: the second differences are now taken
  in time, using the period values, so a jump that is linear in time is
  annihilated whatever the spacing of the comparison periods (the rows
  used to be index-based and could reject a true linear trend with
  unequally spaced periods). Equally spaced periods give exactly the
  same contrasts as before. Non-numeric period values now error under
  `"linear"`.
- [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md):
  a type whose local-linear fit fails is skipped as documented; when the
  last type (in radix order) failed, the function used to stop with
  “subscript out of bounds”.
- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  /
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md):
  when the contrasts’ covariance is zero to working precision (no
  residual variation within the cells, e.g. a constant outcome) the
  functions now stop with a clear message instead of reporting a
  chi-squared statistic made of rounding noise.
- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md):
  the `reference` column of `period_type_jumps` marks the type actually
  used as the reference in each period (the all-below type when its cell
  was fitted, else the next type in order); it used to flag the
  all-below type only, leaving some periods with no reference row.
- Numerical results are unchanged except `rd_homog(p = 2)` (first fix).

### 2026-10-08: source files renamed by what they hold

- **Files are named after what they hold.** The four assumption tests
  live in `rd_typecont.R`, `rd_compstable.R`, `rd_homog.R`,
  `rd_trendcell.R` (they used to be `test_*.R`, which read as unit
  tests); their shared internals in `assumption_tests_helpers.R`; the
  covariance and aggregation in `cross_period_covariance.R`; the
  bandwidth rules in `bandwidth_cct.R` and `bandwidth_joint.R`; the
  trend weights and scheme detection in `trend_weights.R` and
  `sampling_scheme.R`. Every file opens with a header saying what it
  holds and who calls it.
- **Functions read as pipelines.**
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  and the four tests are short sequences of named steps (validate,
  periods, types, scheme, fit cells, covariance, Wald, output), each a
  small helper with a one-line note; the side-of-the-cutoff fit of
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  is `.rd_side_fit()`; the coordinate-descent bandwidth rule is split
  into objective, start and descent.
- **Names and comments.** The cutoff is `cutoff` inside every function
  (`c` stays the argument name, as in rdrobust); locals that collided
  with the scheme codes or base functions were renamed; comments say
  *why* at the numerically sensitive lines (the active-set floor, the
  bias-correction matrix grouping, the HC1 factors, the radix sorts, the
  reference-type drop, the summation order of the pair loops) and dated
  history moved here.
- **Nothing numerical changed.** Every statement that touches a number
  is the same, in the same order; a snapshot of 134 calls (7,763 values)
  is identical before and after every commit of the sweep. Column-check
  errors no longer carry an `Error in <fn>` prefix (the message is
  unchanged).

### 2026-10-08: UX sweep, step 3 — documentation and site

- Every help page rewritten in one vocabulary with a runnable example on
  the shipped data; package help page `?rddid-package`;
  `citation("rddid")`.
- Vignettes: *Get started* reaches
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) in
  four lines of code and reads the printout line by line; *How rddid()
  computes the estimate* holds the hand-built reconstructions;
  *Bandwidth rules and sampling schemes* and *Checking the
  identification assumptions* rewritten on the shipped data. Reference
  index grouped Estimate / Check the assumptions / Building blocks /
  Data.
- Example data `rddid_sim` (running variable fixed) and `rddid_sim_pv`
  (moving; composition stability fails there by design while the other
  assumptions hold).

### 2026-10-08: UX sweep, step 2 — the API

- **Methods.** [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) for `rddid` objects, and
  broom-style
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)
  for fits and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) for the
  four validation tests (registered for the `generics` package, so
  `modelsummary` tables work).
- **Printing.** [`print()`](https://rdrr.io/r/base/print.html) of a fit
  now shows the estimand and period in words, the sampling scheme in
  words, the bandwidth actually used, and z and p-values; the per-period
  fits and the standard error under every sampling scheme moved to
  [`summary()`](https://rdrr.io/r/base/summary.html). The four tests
  print the same way: the assumption, its null in words, then the Wald
  statistic.
- **Fields added (nothing removed or renamed):**
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  returns `call`, `bandwidth$h_by_period`/`$b_by_period` (the bandwidths
  used in every period, whatever the rule) and `estimates$z`/`$p`;
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  and
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  return top-level `statistic`, `df`, `p_value`, `scheme`, `estimand`,
  `call` like the other two tests.
- **Arguments.** `rddid(weights=)` is now `trend=` (the old name still
  works, with a message): it names the assumption on the confounding
  jump, and `weights` means observation weights elsewhere in R.
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  gains `t_rd`/`comparisons` so the same call works for every function;
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)/[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  take `p`, `q` explicitly instead of `...` (passing `b` through `...`
  used to fail with “argument matches multiple formal arguments”), and
  list `t_rd` before `comparisons` like
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)/[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  no longer silently accept unknown arguments. `kernel` is validated on
  entry everywhere.
- **Messages.**
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  says so when no `id` is given (rows are then treated as separate
  units);
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)/[`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  list the cells they skip for `min_n`.
- Numerical output is unchanged: every value of the step-1 snapshot is
  identical.

### 2026-10-08: UX sweep, step 1

- **The default bandwidth rule is now the common bandwidth,
  `bwselect = "joint"`** (one `h` for every period, minimising the
  asymptotic MSE of the aggregate estimator). The iterative
  period-specific rule `"iter"` remains available but is no longer the
  default and is not used in the paper. Calls that pass `bwselect`
  explicitly are unaffected.
- **`rdrobust` moved from Suggests to Imports.**
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  no longer falls back to `0.5 * IQR(x)` when rdrobust is not installed
  (the fallback on an `rdbwselect()` error is unchanged), so the default
  path gives the same numbers on every machine.
- **The composition-adjusted family is removed**: `rd_att()`,
  `rd_sadjust()`, `rd_c()`, `rd_adjust()` and their tests. They belong
  to a companion paper, not to the RD-DID paper, and are recoverable at
  the git tag `v0.4.0.9000-composition`.

### Earlier in 0.4.0.9000 (before the 2026-10-08 sweep; the default rule was `"iter"` then)

- **`bwselect = "joint"` no longer depends on which period is labelled
  `t_rd`.** The common AMSE-optimal bandwidth used to fit *every* period
  at the RD period’s CCT pilot pair to estimate the aggregate bias and
  variance constants, so an estimator that aggregates several RD periods
  (the same linear combination of per-period discontinuities, whichever
  RD period carries the `t_rd` label) got a different `h*` for each
  labelling (Grembi, aggregate ATU over 2001-2004: 482 vs 320). Both
  joint rules now estimate each period’s constants at that period’s
  **own** CCT pilot through one shared helper (`.bw_constants()`), and
  the common `h*` is the exact scalar minimizer of the iterative rule’s
  objective (pinned in `test-appB-bandwidth.R`, together with a
  `t_rd`-relabelling invariance test). Under `"joint"` the pilot `b_t`
  now keeps each period’s CCT ratio, `b_t = h* b_t^CCT / h_t^CCT` (was
  the RD period’s ratio for all periods); `rddid()$bandwidth$b` is
  therefore a named per-period vector. `"joint"` numbers move (a few
  percent in the snapshot DGPs; Grembi Table 1 common-`h*` rows change);
  `"iter"` moves only through its seed (below 1e-4 relative); `"cct"` is
  unchanged.

- **New argument `estimand = c("att", "atu")`** on
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md),
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  and
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  (default `"att"`, nothing existing moves). Set `"atu"` when the
  comparison periods are uniformly *treated* (the paper’s Section 6: the
  ATU design is the ATT design with the sides of the cutoff exchanged).
  For
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md),
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  and
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  the estimates and tests are numerically identical, so `estimand` only
  labels the output (pinned by `test-mirror-invariance.R`).
  `rd_compstable(estimand = "atu")` mirrors the running variable
  (`x -> c - x`) and tests composition stability on the **below**-cutoff
  shares (type indicator unchanged, “above the cutoff in the other
  period”, so the jump is `pi_{t_RD,(-)}(1) - pi_{t_0,(-)}(1)`);
  observations at `x == c` error under `"atu"` (place the cutoff between
  support points).

- **Canay-Kamat permutation test and McCrary tests removed** from
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  and
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  (arguments `q`, `S`; outputs `ck_perm`, `mccrary_within`,
  `mccrary_pooled`; internal helpers `.q_rot()`, `.mccrary()`). The
  paper reports the local-linear Wald tests only. The remaining outputs
  are numerically unchanged.

- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  default `type_by` is now `"rd_side"` (the unit’s side of the cutoff in
  the RD period, the partition of the paper’s Section 4.4), matching
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md);
  was `"pattern"`.

- Print methods echo option values:
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  prints `bwselect: iter (5 iterations)` and
  `scheme: pc (auto-detected)` (were
  `bandwidth: period-specific joint AMSE ...` and
  `sampling scheme: PC`), and `Robust SE by scheme: cs= pc= pv=`; the
  test print methods drop the LaTeX assumption labels and the
  necessary/sufficient notes.

- Documentation pass after a user-focused review of the site
  (`rd-did/docs/reviews/ 2026-09-09_oren-persona_rddid-site.md`): what a
  comparison period is and the scheme detection rule stated on the home
  page / in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  the tests vignette opens with its scope (a time-varying running
  variable);
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  leads with its `rdrobust` equivalence;
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)’s
  Value lists the element names; assumption references by name, not
  label; the normative wording (“suggestive”, “neither necessary nor
  sufficient”, “preferred”) removed throughout.

- **[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  default bandwidth rule is now `bwselect = "iter"`** (period-specific
  bandwidths by coordinate descent on the aggregate AMSE, started at the
  common h\*), the rule Section 5.3 of the paper states as preferred;
  previously `"joint"`. Calls that pass `bwselect` explicitly (all the
  paper’s scripts do) are unaffected.

- **Site rebuilt around three vignettes** (`vignettes/`, plan in
  `dev/site_plan.md`): *Estimating the ATT with rddid()* (per-period
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md) +
  [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
  then
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  under constant, linear and custom weights), *Bandwidth rules and
  sampling schemes* (`bwselect = "cct" / "joint" / "iter"`, `start`,
  fixed `h`; `scheme` detection and the three standard errors) and
  *Composition validation tests*
  ([`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  on simulated scenarios that satisfy or violate each assumption). Every
  vignette hand-codes its DGP and prints the truth next to each
  estimate; `dev/site_dgp_check.R` verifies the scenario claims.
  `_pkgdown.yml` groups the reference into *Estimation* and *Validation
  tests*; README gains a quick start (now generated from `README.Rmd`).

- The composition-adjusted estimators (`rd_adjust()`, `rd_sadjust()`,
  `rd_c()`, `rd_att()`) are tagged `@keywords internal`: still exported
  and tested, but no longer listed on the site. They belong to a
  companion paper and are not part of the current manuscript.

- `Suggests` gains `knitr`, `rmarkdown`, `ggplot2`;
  `VignetteBuilder: knitr`.

## rddid 0.3.5.9000 (development)

- **Synced with Appendix B of the paper (rewritten 2026-09-08).** A code
  \<-\> equation map now lives in `dev/appB_map.md`: one row per object
  the package computes, with the paper’s exact expression, the
  implementing symbol, conventions, audit status, and the test that pins
  it. `dev/check_appB_labels.R` verifies that every cited paper label
  still exists in `main.tex`; `dev/snapshot_rddid.R` is a numerical
  regression snapshot.

- Bandwidth selectors follow App. B.4 at general polynomial order `p`
  (previously the exponents and constants were hard-coded for `p = 1`):
  `rd_period()$b_const` is `(p+1)! (D - D_bc) / h^{p+1}`; the common
  `h*` (`eq:common_h_opt`) uses the constant `((p+1)!)^2 / (2(p+1))` and
  exponent `1/(2p+3)`; the period-specific objective and its
  regularization use `h^{p+1}/(p+1)!`. Numerically identical at `p = 1`.

- `bwselect = "iter"` under `scheme = "pc"`: the same-side cross-period
  covariance term of the aggregate AMSE now scales as
  `omega(h_t/h_s)/h_s` per side, with `omega(rho)` the kernel constant
  of Lemma `cov-pc` (new internal module `R/kernel_constants.R`, ported
  from the paper’s simulation-verification code), instead of the
  previous `1/max(h_t, h_s)` approximation (exact only for the uniform
  kernel at `p = 0`). `.bw_joint_iter()` additionally returns the
  objective value and the objective function for diagnostics.

- [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md):
  the active set is the union of the pilot and main windows (was the
  pilot window alone, which silently truncated the main fit when
  `b < h`). Identical whenever `b >= h`.

- New tests: `test-appB-conformance.R` (from-the-equations reference
  implementation of the estimator, bias correction and the three
  sampling-scheme variances, 1e-10), `test-appB-bandwidth.R` (B.4
  objectives and selectors), `test-kernel-constants.R`.

- Regularization of the bandwidth selectors: the curvature-variance
  estimate `Var(B-hat_t)` now pairs the pilot-window influence weights
  with the residuals of the order-`q` pilot fit at `b` (the same
  convention as the bias-corrected variance) instead of the order-`p`
  residuals at `h`. Moves `bwselect = "joint"`/`"iter"` bandwidths by a
  fraction of a percent at the CCT pilot ratio; at wide pilots
  (`b/h >= 3`) the old estimate over-stated the variance.

- Validation tests synced with Section 4.4 of the paper
  (`dev/tests_map.md`): documentation cites the assumptions by label
  (`ass:type-cont`, `ass:comp-stable`, `ass:homog`, `ass:trend-cell`)
  instead of numbers that changed in the paper;
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  drops one reference type instead of pseudo-inverting the structurally
  singular covariance (same statistic for binary types, exact df); a
  comparison-period unit exactly at the cutoff now stays on the
  reflected side; the joint-over-pairs result is documented as
  approximate. New from-the-text conformance tests
  `test-s44-conformance.R`.

- Guards: `bwselect = "cct"` and the joint pilot go through the guarded
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  (fallback + finiteness checks) instead of calling
  [`rdrobust::rdbwselect()`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
  directly; `.bw_joint()` warns when `h*` exceeds the running-variable
  radius; the coordinate descent warns when a bandwidth ends on the
  search boundary; `rddid(scheme = "pc"/"pv")` warns when no unit id
  repeats across periods (all cross-period covariances are then zero);
  the search cap is NA-safe.

- Stale references in code comments to the removed “Appendix C”,
  `lem:coercive` and `eq:amse-ps` replaced by the current labels
  (`app:est-bw`, `eq:amse-att`, `eq:update`, `alg:coorddesc`).

## rddid 0.3.0.9000 (development)

- `rddid(..., bwselect = "iter")` gains a `start` argument controlling
  the coordinate-descent seed: `"hstar"` (default) seeds all periods at
  the common joint-optimal bandwidth h*; `"cct"` seeds each period at
  its own CCT/IK pilot h; or supply a named numeric/list of per-period
  bandwidths for a manual seed. Each run weakly improves the joint AMSE
  over its own start; seeding from h* therefore weakly dominates the
  common-h rule. (Seed from `"cct"` when the per-period biases nearly
  cancel, where h\* is inflated.)

- Bug fix (reproducibility):
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  called `set.seed(NULL)` immediately before its Canay-Kamat
  permutation, which **re-initialised** the RNG from system entropy and
  discarded any seed the caller had set — so the CK p-value changed on
  every call. Removed; the permutation now inherits the caller’s RNG
  state, so [`set.seed()`](https://rdrr.io/r/base/Random.html) before
  the call makes the CK p reproducible.
  ([`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  never had this and was already reproducible.)

- [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md):
  the Canay-Kamat permutation test now chooses the number of nearest
  observations per side, `q`, by the Canay & Kamat (2018) rule of thumb
  **by default** (`q = NULL`), per period. A fixed `q` over-rejects in
  finite samples when the type distribution varies steeply in the
  running variable at the cutoff; the rule of thumb shrinks `q` as that
  association strengthens. Pass an integer `q` to force a fixed value;
  the per-period `q` used is returned in `meta$q_used`.

- [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md):
  same change — its Canay-Kamat permutation shares the identical
  fixed-`q` exposure, so `q` now defaults to the rule of thumb
  (`q = NULL`), chosen per `(t_RD, t_0)` pair on the pooled reflected
  sample. The per-pair `q` used is returned in `meta$q_used` (and echoed
  in each `pairs[[...]]$q`). Pass an integer `q` to force a fixed value.

- Bug fix (numerical robustness): the joint Wald pseudo-inverse
  (`.joint_wald()`, used by
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  and
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md))
  now uses the [`MASS::ginv`](https://rdrr.io/pkg/MASS/man/ginv.html)
  relative tolerance `sqrt(eps)*max(sv)`. The previous, tighter
  tolerance could leave the structural-zero singular value (the
  per-period type indicators sum to 1) just above the cut on some LAPACK
  builds, inflating the statistic into a platform-dependent false
  rejection. Caught by the new CI on Ubuntu.

- Continuous integration: added an `R-CMD-check` GitHub Action (standard
  multi-OS matrix) and R-CMD-check / MIT-license badges to the README.
  Removed a placeholder ORCID from `DESCRIPTION`.

- `R CMD check` is now clean (0 errors / 0 warnings). Fixes: replaced
  the `\insertCite{}` macros in
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  (Rdpack was not a dependency) with the plain-text citations already in
  the References; documented the `regularize`/`reg_const` arguments of
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  and dropped the `VignetteBuilder: knitr` declaration (and the
  `knitr`/`rmarkdown` Suggests) since the package ships no vignettes.

- Internal: consolidated duplicated logic into shared helpers. A new
  `.cov_scheme()` (the `cs`/`pc`/`pv` scheme-combine of `.cross_cov()`)
  now backs the cross-period covariance in
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md),
  and
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
  replacing three near-identical inline copies (including the former
  `.cross_cov_homog()`/`.match_sum_homog()`). A new
  `.scheme_from_long()` primitive backs both `.detect_scheme()` and
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)’s
  scheme detection. No change in results (verified to machine
  precision).

## rddid 0.2.1

- Internal: removed a duplicate `.build_types()` (a second copy lived in
  `test_homog.R` and shadowed the canonical one in `test_helpers.R` at
  load time). There is now a single shared implementation used by
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  and
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md).
- Convention: units exactly at the cutoff are now treated as above it
  (`V_i = 1{R_i >= c}`) everywhere, including sampling-scheme detection.
  A unit sitting on the cutoff no longer registers as a separate “side”
  and so cannot be misread as a side-switch.
- Internal: dropped the unused single-cell `.ck_perm()` helper (and its
  test).
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)/[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  use an inlined *joint* Canay–Kamat permutation with one shared
  per-period shuffle; the standalone helper was a dead parallel path.
- Documentation: corrected the assumption numbering in the Section 3.4
  test functions to match the manuscript’s `\begin{assumption}` ordering
  —
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  is Assumption A7 (was mislabelled A6),
  [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  is A8 (was A7), and
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  is A9 (was A8). Affects titles, `print` output, and cross-references
  only; no behaviour change.
- Documentation: fixed an unmatched apostrophe in the
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  `@examples` comment that was silently dropping the entire example from
  the rendered help page.

## rddid 0.2.0

- Added tests for the Section 3.4 identifying assumptions of the
  time-varying-running-variable design:
  - [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
    — continuity of the type distribution (LL-Wald and Canay–Kamat
    permutation, both necessary & sufficient; McCrary within-type
    sufficient-not-necessary; McCrary pooled neither).
  - [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
    — composition stability across periods, via the reflection
    construction (LL-Wald and Canay–Kamat permutation, both necessary &
    sufficient). The permutation uses the partially-overlapping-samples
    scheme for units above the cutoff in more than one period.
  - [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
    — type-homogeneous confounding, tested in comparison periods (only
    suggestive of the assumption at the RD period: neither necessary nor
    sufficient there).

## rddid 0.1.0

- From-scratch rewrite: per-period local-linear RD engine
  ([`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
  validated against `rdrobust` to machine precision), the
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  aggregate estimator with constant/linear/custom weights, CS/PC/PV
  sampling-scheme variances, and joint / CCT / period-specific bandwidth
  selection.
