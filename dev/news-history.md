<!-- The development log of rddid up to its first public version (kept for the authors; the
package's user-facing NEWS is inst/NEWS.md). -->

# rddid 0.4.0.9000 (development)

## 2026-10-29: no "pilot" in user-facing text (0.4.0.9010)

* T1793 part 2, Dor's rule for the paper ("no pilot"): exported docs, vignettes, README, NEWS,
  printed output and error messages say "bias-correction bandwidth" for `b`; `weights_detail$pilot`
  -> `weights_detail$initial` (the paper's "initial weights"). Internal names (`.bw_pilots()`,
  `pilot_fits`, `w_pilot`) and @noRd developer docs keep "pilot".
* Golden master against master 34a05c4 (0.4.0.9009): every leaf identical except
  `weights_detail/pilot` -> `weights_detail/initial` in the min_variance cells.

## 2026-10-28: the regularization's cross-period covariances under PC only (0.4.0.9009)

* App. B.4 Regularization (`eq:vhat-bbar`, Overleaf 88120e3; Dor's decision (a), T1793): the
  variance of the aggregate bias constant adds the same-side cross-period covariances only under
  PC, the same rule as the variance V-hat; `b` is regularized the same way, `d` not at all.
  `.bw_constants()` `var_bias` (for `h`) and `.bw_stage()`'s `R` (for `b`) now use the same-side
  covariances under `"pc"` and none otherwise; before, `"pv"` added same minus opposite side.
  `.bw_agg_var()` loses its `"same_opp"` rule.
* Golden master against master aec94cb (0.4.0.9008): exactly 11 cells change, every one with
  `scheme = "pv"` and `bwselect = "joint"` (`pc|joint|scheme-pv`, the pv constant / linear / p2 /
  atu / numeric-w / scheme-pv / min_variance cells, `pv4|joint|linear-3comp`, the two `pv4|agg`
  cells); `pv|joint|noreg`, PV data under `scheme = "cs"`/`"pc"`, and CS data under `"pv"` (no
  unit in two periods) are identical, as are all other cells. Test: PV equals CS on the same
  data for `d`, `b`, `h`, `reg`.

## 2026-10-27: the pilot bandwidths chosen for the aggregate (PR #12: rddid 0.4.0.9004 on its branch, lands as didisc 0.4.0.9008)

* `bwselect = "joint"`: the pilot bandwidth `b` is now one bandwidth for every period, chosen for
  the RD-DID estimate as rdrobust's `rdbwselect` (mserd) chooses its pilot for one period: `b` is
  AMSE-optimal for the aggregate `sum_t coef_t` of the periods' `(p+1)`-th derivative jumps (order-q
  fits), with its own curvature from order-(q+1) fits at a bandwidth `d`, itself AMSE-optimal for
  the aggregate `(q+1)`-th derivative jump (order-(q+2) fits at each side's full range; not
  regularized). Every variance and design moment is at the shared `h_rot`; `b` is regularized
  like `h`; `d` and `b` are capped at the largest distance to the cutoff (rdrobust's
  `bwrestrict`). The curvature of `h`'s bias is then taken at this `b` in every period, and every
  period's bias correction uses it. Before, each period's curvature came from its own CCT pilot
  `b_t^CCT`, and the bias correction used `b_t = h* b_t^CCT / h_t^CCT`.
* New internals `.bw_stage()` (one rdrobust_bw stage for the aggregate, in rdrobust's units) and
  `.bw_pilots()`; `.bw_constants()` / `.bw_joint()` take one `b` instead of `pilot_bws`;
  `bandwidth` in the result has `d` and one `b` (no `pilot_bws`). `.bw_rot()` takes the IQR from
  type-2 quantiles, as rdrobust (was `IQR()`, type 7; matters only when IQR / 1.349 < sd).
* rdrobust's handling of the data, ported so that one period is exactly `rdbwselect`:
  `masspoints = "adjust"` (new `.bw_mass()`: `h_rot` uses the average number of distinct
  running-variable values per period instead of rows; with mass points, a side with
  `1 - M/N >= 0.2`, `h_rot` and `d` are raised to at least the 10th distinct distance to the cutoff,
  rdrobust's `bwcheck = 10`) and `bwrestrict` (`h_rot`, `d`, `b` and now also `h` capped at the
  largest distance to the cutoff; `h` still warns when capped). Without this the order-(q+1) fit
  of the `d` stage at `h_rot` aborted on discrete running variables that master handled (review
  finding).
* Why: paper App. B.4 "Bandwidth estimation" (Overleaf `c0b890c`, 2026-10-27) and the T1787
  simulations (rd-did `code/simulations/bw_routes`, `docs/bw_routes_note.md` section 8): joint
  pilots give the most accurate common bandwidth of the routes compared.
* Checks: one period, `h` and `b` equal `rdbwselect(mserd, hc1)` to 1e-8 (p = 1, 2; both branches
  of the rule of thumb); balanced PC panel, `b` and `h` equal `rdbwselect` on the composite
  outcome at the same `h_rot`; `.bw_stage()`'s third stage equals `.bw_joint()`'s `h` (two
  implementations of eq:common_h_opt); `d` (one period and CS) and `b` (CS) equal a from-scratch
  `lm()` + HC1 reference; a discrete running variable runs, and one period equals `rdbwselect`
  with its mass-point adjustment.
* rdrobust 4.1 (CRAN 2026-10-07; CI installs it): the selector follows 4.1.1, whose d stage takes
  each side's range times 1 + sqrt(eps) (the farthest unit keeps a positive weight) and whose
  `bwcheck` floor is the 10th distinct distance times 1 + sqrt(eps). Against rdrobust 2.2 the
  one-period equality then holds to ~1e-3 on heavy-tailed data (tests use 1e-8 for rdrobust >= 4.1,
  1e-3 below).
* Golden master, all recorded under rdrobust 4.1.1: `baseline_bmomrot_rd41.rds` (master
  d25dd5e) against the new `baseline_jointpilots.rds`: exactly the 31 joint cells change (h by a
  factor 0.86-1.52); the other 98 are identical. Upgrading rdrobust 2.2 -> 4.1.1 alone moves 80
  cells of master by at most 6e-7 (relative; h and estimates), through `rdbwselect`. On `rddid_sim` the print snapshot's common h moves 0.2694 -> 0.3986, with
  b = 0.5936 in every period.

## 2026-10-27: the bias design moment at the shared rule-of-thumb bandwidth (0.4.0.9003)

* `bwselect = "joint"`: on each side, the bias constant is the design moment
  `bmom = e_0'(X'A(h)X)^{-1} X'A(h) ((x - c)/h)^{p+1}` times the curvature (the order-q
  coefficient on `(x - c)^{p+1}` at `b`). The design moment is now evaluated at `h_rot`, like the
  variance constants; the curvature still comes from each period's CCT pilot at `b_t`. Before,
  the design moment was taken at each period's CCT `h_t` (via `b_const`, the gap between the
  conventional and bias-corrected jumps at the CCT pair). This is what rdrobust's `rdbwselect`
  does (`BConst` at the preliminary bandwidth `c_bw`, `beta_B` at `b_bw`), so every quantity
  that depends only on the running variable and the residuals is estimated at the one
  bandwidth `h_rot`, and only the curvature needs a pilot. The regularization's `Var(B-hat_t)`
  uses the same factorization.
* Why: paper App. B.4 "Bandwidth estimation" (2026-10-27): one rule for the variance and design
  constants, as CCT; the limit is the same (both design moments tend to the kernel constant).
* `rd_period()` side fields: new `bmom`, `curv`, `g_curv` (with `g_diff = h^{p+1} bmom g_curv`).
* Golden master (new baseline `dev/snapshots/baseline_bmomrot.rds`, against `baseline_bwrot.rds`):
  exactly the 31 joint cells change; the other 98 are identical apart from the new fields. On
  `rddid_sim` the print snapshot's common h moves 0.2734 -> 0.2694.

## 2026-10-24: the common bandwidth estimates its constants as CCT does (0.4.0.9002)

* `bwselect = "joint"`: the variance constants, and under `"pc"` the same-side cross-period
  covariances, are now estimated at one rule-of-thumb bandwidth `h_rot` shared by every period
  (rdrobust's preliminary bandwidth, `2.576 min(sd, IQR/1.349) N^{-1/5}` for the triangular
  kernel, on the pooled running variables; `bandwidth$h_rot` in the result). Before, each period
  was fitted at its own CCT pair and the PC covariance was rescaled to a bandwidth ratio of one by
  the kernel constant. The bias constants still come from each period's CCT pilot.
* The regularization is now `reg_const` times the variance of the aggregate bias estimate
  `sum_t coef_t B-hat_t`, with the cross-period covariances of the `B-hat_t` (same side under
  `"pc"`, same minus opposite side under `"pv"`). Before, it summed the per-period variances,
  which overstates it 1.7 to 7.7 times when the periods' errors are correlated.
* Why: paper App. B.4 "Bandwidth estimation" (Overleaf 3762b95..ba9a6c4); simulation evidence in
  rd-did `docs/bw_routes_note.md` (T1755): the variance side was already accurate, the
  regularization was the defect; the new rule is closer to the AMSE-optimal common bandwidth under
  PC and PV with little movement, unchanged under CS.
* Golden master (new baseline `dev/snapshots/baseline_bwrot.rds`, against a fresh snapshot of
  8300e37): exactly the 31 `rddid(bwselect = "joint")` cells change (bandwidths and everything
  computed at them); every `cct`, fixed-`h` and test cell is identical. The print snapshot moves
  with the joint bandwidth (h 0.2672 -> 0.2734 on `rddid_sim`).
* Internal: `.bw_constants()` returns `var_bias` and `var_cross` (K x K) and `h_rot`, no longer
  `var_bias_const`, `kappa_plus`, `kappa_minus`; `.bw_rot()` is new. The kernel constants in
  `kernel_constants.R` are no longer used by the selector (only by the conformance tests).

## 2026-10-16: the tests use only the periods they are given

* `rd_compstable()`, `rd_homog()` and `rd_trendcell()` now use only the rows of `t_rd` and
  `comparisons`, as `rddid()` does. Before, they built each unit's type from every period in
  `data`: a unit missing from a period outside the test was dropped, and under
  `rd_compstable()` the unused periods split the units into finer types. With `comparisons =
  NULL` (every period) nothing changes. `rd_typecont()` already did this.
* Golden master: one cell changes, `rd_compstable(t_rd = 3, comparisons = 2)` on three
  periods: 2 types and 1 degree of freedom (before: 4 and 3, period 1 splitting the types). The
  Grembi application calls (A7 to A10) are identical.
* With a single contrast, the sup-t / sWald p-value no longer carries the contrast's name.

## 2026-10-15: one input step for `rddid()` and the four tests

* `rddid()`, `rd_typecont()`, `rd_compstable()`, `rd_homog()` and `rd_trendcell()` now check and
  prepare their input in one shared internal step, so the same mistake gives the same message
  everywhere: `data` must be a data frame, the running variable numeric (the outcome numeric or
  logical), the period column numeric, character or factor; `t_rd` and `comparisons` must be
  periods in the data, not overlap, and not repeat; each unit may appear once per period; a
  factor `t_rd` or `comparisons` is read as character. The four tests say that they take one RD
  period, and `kernel` is not case sensitive in any of them. Errors name the function in their
  header (`Error in rd_homog(...)`), not in the message (before: an `rd_homog: ` prefix in some).
* Missing values are handled the same way everywhere: rows with a missing value in a used column
  are dropped, with a message. Before, `rd_typecont()` and `rd_compstable()` dropped them without
  a message, and `rd_homog()` and `rd_trendcell()` kept them (a missing outcome counted toward
  `min_n` and a cell's `n`). In these two tests a row with a missing outcome still gives its
  unit's side of the cutoff and is left out of the fits only.
* `rddid()` checks for repeated (id, period) rows among the rows it uses (before: in all periods
  of `data`), and names a period in which every row has a missing value.
* Results on data without missing values are unchanged (golden master: every result identical).

## 2026-10-15: exact linear-trend weights with calendar years

* `trend = "linear"`: the comparison weights are now computed with time centred at the RD period.
  The weights are the same in exact arithmetic, but with calendar years (e.g. 1999, 2000) the old
  computation lost about 1e-9 to rounding. Linear-trend results move only at that level (largest
  relative change 7.6e-10 across the reference calls); constant-trend and numeric weights are
  unchanged.

## 2026-10-15: minimum-variance weights for the bias-corrected estimate

* With `weighting = "min_variance"`, the bias-corrected (robust) row of `estimates` now uses its
  own minimum-variance weights, computed from the variances and covariances of the bias-corrected
  jumps (the paper, Section 5.5). Before, it reused the weights that minimize the conventional
  variance. The conventional row, the bandwidths and `weighting = "ols"` are unchanged.
* New fields `coef_bc` and `weights_detail$weights_bc` (with `Xi_bc`, `xi_bc`); `print()` shows
  the bias-corrected weights when they differ.

## 2026-10-15: several RD periods in `rddid()`

* `t_rd` may be a vector of RD periods. `rddid()` then estimates the weighted average of their
  effects, ATT^phi = sum_t phi_t ATT(t) (the paper's Corollary 2 and Section 5.6), with weights
  from the new argument `rd_weights` (equal weights by default; any finite weights). Each
  period's jump is estimated once; each comparison period enters with the phi-weighted sum of
  the weights each RD period gives it. Standard errors include the covariances between RD periods,
  and `bwselect = "joint"` chooses one bandwidth for the average.
* New return fields `rd_weights` and `weights_by_target`; with several RD periods `weights` holds
  the aggregated comparison weights. A numeric matrix `trend` sets the comparison weights per RD
  period.
* `weighting = "min_variance"` is available for a single RD period only.
* Results for a single RD period are unchanged.

## 2026-10-15: `bwselect = "iter"` removed

* The period-specific bandwidth rule `bwselect = "iter"` (one bandwidth per period, by coordinate
  descent on the RD-DID estimate's AMSE) and its argument `start` are removed; the paper no longer
  uses the rule. `bwselect` is now `"joint"` (default) or `"cct"`. Estimates under `"joint"`,
  `"cct"` and a fixed `h` are unchanged. The last version with the rule is tagged
  `v0.4.0.9000-iter`.

## 2026-10-15: size corrections for the validation tests (changes the tests' default results)

A size study (clean null, n = 1,000, three periods, 5% level) found the four tests rejecting too
often: type continuity 9%, composition stability 9%, homogeneous confounding 8%, constant
within-type confounding 8%. The defaults bring this to 6-7% (5.7% for composition stability):
sup-t avoids inverting the estimated covariance, but it still uses estimated standard errors, so
a variance deficit or noisy standard errors still inflate it, by less than the Wald. Estimates of
`rddid()` are unchanged.

* `rd_typecont()` and `rd_compstable()`: new arguments `joint = c("supt", "swald", "wald")` and
  `studentize = TRUE`.
  * The default `"supt"` compares the largest |t| of the tested jumps with a sup-t critical value
    simulated from their estimated correlation.
  * `studentize = TRUE` (default) makes that critical value t-based. Each standard error is an
    estimate, with the degrees of freedom of its robust variance: Satterthwaite / Bell-McCaffrey,
    from the jump's regression weights. With normal critical values (`studentize = FALSE`) sup-t
    rejected a true null 7% of the time at three periods; studentized, 6%, at the same
    size-adjusted power.
  * `"swald"` uses the sum of squared z's (Fusejima, Ishihara and Sawada, 2025).
  * `"wald"` is the chi-squared Wald statistic used so far, which inverts the estimated covariance
    and over-rejects as the number of jumps grows (7-9% at 9 jumps, 15% at 28).
  * The critical values are simulated with 100,000 draws under a fixed seed (reproducible; the
    caller's random-number state is restored).
  * With one jump all three are the same two-sided normal test, so two-period calls give the
    Wald's p-value.
  * New fields `joint_method`, `studentize`, `joint_tests` (every statistic, with t, normal or
    chi-squared critical values) and, for `rd_typecont()`, `t_stats` and `nu`. The printouts name the
    statistic and its critical values, and point to `x$joint_tests`.
    `ll_wald` and the per-period/per-pair `ll_wald` keep the Wald; the per-period/per-pair `test`
    holds the test under `joint`.
* `rd_compstable()`: the joint test over pairs no longer adds up the pair statistics as if the
  pairs were independent; they share the RD-period units.
  * With `joint = "wald"` each pair's test is unchanged, but the joint test over pairs differs from
    the old summed statistic. It is now the Wald over all pairs' jumps with their full covariance.
  * It now tests all pairs' jumps together, with their full covariance.
  * The covariance across pairs comes from the shared units, matched on id, following `scheme`.
  * `joint$ll_wald` is now the Wald over the stacked jumps, and `joint$theta` / `joint$Sigma`
    hold them.
  * Each pair's `kept` and `cov_matrix` are returned.
* New `rd_checks()` runs the four tests with one set of options and prints one table, ending with
  which assumptions are rejected. Options for one test only go in that test's list, for example
  `homog = list(bwselect = "cct")`. `tidy()` returns the table.
* `rd_period()`: the side lists keep the intercept weights `a` and `a_bc`. This is append-only;
  no number changes.
* `rd_homog()` and `rd_trendcell()`: new default `bwselect = "common"`.
  * It uses one bandwidth in every cell: the CCT bandwidth of the outcome on the running variable
    with all comparison periods pooled.
  * Per-cell CCT bandwidths (`"cct"`, the former default) made every cell variance 10-15% too
    small.
  * The common bandwidth gives about 6% instead of 8%, at about 8 points of power against a
    violation confined to one cell.
  * New field `bandwidth` (`method`, `h`); the printout reports it.

## 2026-10-09: minimum-variance comparison-period weights

* New argument `weighting` in `rddid()`. `"ols"` (default) keeps the weights used so far (equal under
  `trend = "constant"`, the least-squares line under `"linear"`). `"min_variance"` uses the
  minimum-variance weights: the comparison-period weights the confounding-trend assumption allows
  that make the variance of the conventional estimate smallest.
  * They are computed from the per-period variances and, in a panel, covariances behind the
    standard error.
  * In a repeated cross-section with a constant trend they are inverse-variance weights.
  * Under `bwselect = "joint"` the weights and the common bandwidth are chosen in three steps; with
    `"cct"` or a fixed `h` the bandwidths do not depend on the weights, so the weights are computed
    once.
  * Not available with a numeric `trend`, or with `bwselect = "iter"` unless `h` is given.
  * New fields `weighting` and `weights_detail` in the returned object; `print()` and `glance()`
    report the weighting.
  * Results with the default `weighting = "ols"` are unchanged.

## 2026-10-09: verification fixes (no change on the paper's numbers)

* `rd_period()` requires a single main bandwidth `h` and a single pilot `b`. A vector used to be
  recycled across the observations by the kernel weights, silently giving each unit its own
  bandwidth (the joint rule's per-period pilots passed as one `b` did exactly that).
* `rd_typecont()` and `rd_compstable()`: a type share that is constant near the cutoff on both
  sides with no jump (the type is absent, or universal, in the window) carries no information
  and is now left out of the Wald test, which is run on the remaining shares; `ll_wald$dropped`
  counts the shares left out. A share that jumps deterministically (zero standard error, jump
  not zero) still makes the test undefined (`NA`). The 2026-10-08 guard made both cases `NA`,
  which turned informative tests into `NA` whenever one type was absent near the cutoff.

## 2026-10-08: `tidy()`/`glance()` always available; a no-Suggests check

* `generics` moved from Suggests to Imports and its `tidy()`/`glance()` are re-exported:
  `tidy(fit)` works after `library(rddid)` without `generics::` and without installing anything
  else. No numerical change.
* New workflow `check-no-suggests.yaml`: R CMD check with only the hard dependencies installed
  (plus testthat, knitr, rmarkdown), so every use of ggplot2, modelsummary or broom stays
  guarded. `broom` added to Suggests (the Get-started modelsummary chunk needs it).

## 2026-10-08: second bug-hunt round (input handling; no numerical change on valid input)

* `rddid()` now works on complete (outcome, running variable, id) rows only and says how many
  it dropped: an NA running variable used to read as a side switch, turning a "no unit changes
  side" panel into "some units change side" and moving the common bandwidth and the estimate;
  `n_by_period` and `nobs()` now count the rows actually used.
* Errors instead of silent misuse: a unit appearing twice in a period (`rddid()` and the four
  tests), a period repeated in `comparisons`, `comparisons` that include `t_rd` or periods not in
  the data (all four tests), `b` without `h`, an `h` that is not a single positive number, a
  `t_rd` not in the data (`rd_typecont()`), `type_by = "pattern"` in `rd_trendcell()` (a cell
  fixed across the comparison periods cannot use a comparison period's own side).
* Named numeric `trend` weights are matched to the comparison periods by name; `q` now reaches
  `rd_bw_cct()` (the pilot bandwidth used to be chosen for `q = p + 1` whatever `q` was);
  `kernel` accepts any case; a factor `t_rd` works.
* Untestable results (no testable contrast, or a type share whose jump has a zero standard
  error) are reported as `statistic = NA`, `df = 0`, `p = NA` instead of a chi-squared of 0
  with p = 1 or a rounding-noise rejection; the joint composition-stability test skips
  untestable pairs.
* The no-switcher guard looks only at the periods a test uses, and names a repeated
  cross-section for what it is.
* `plot()` of `rd_homog`/`rd_trendcell` works when `t_rd` was passed as a variable (the objects
  now carry `t_rd`); `plot(cs, pair = )` checks its index; `plot_switchers()` errors on a
  cross-section; more than five types get distinct colours.
* `rd_homog()`/`rd_trendcell()` choose a cell's bandwidth only after the `min_n` check (no
  fallback noise from skipped cells); `rd_homog()` prints the skipped cells before it stops
  for lack of contrasts; its `reference` flag is set only in periods that contributed a
  contrast; comparison periods print in time order.
* Tests: console snapshots (`tests/testthat/_snaps/`) and one regression test per fix.
* Also: `tidy()` of a fit honours `conf.int` and `conf.level`; DESCRIPTION cites the paper and has
  a shorter Title; `CITATION.cff`; lifecycle badge; the How-it-works article is guarded on
  ggplot2; rdrobust skips removed from the tests.

## 2026-10-08: review round on the plots, the console and the site

* **Figures** (after a review against the paper's own figures): no titles, subtitles or
  statistics in any figure; the paper's palette and `theme_bw`; bins never straddle the cutoff;
  share plots clip at 0–1 without dropping line segments; the composition-stability axis shows
  the running variable on both sides (folded ticks) and marks the mirrored period correctly
  under `estimand = "atu"`; `plot.rd_trendcell()` draws the least-squares line under
  `trend = "linear"`; switcher colours follow the side in the second period; periods in time
  order everywhere; legends keep no internal names when a theme is replaced.
* **Console**: `summary()`'s per-period table is aligned with its (time-ordered) labels — the
  previous commit had scrambled the numbers; `glance()` is one row under every bandwidth rule
  (`$h` used to partial-match `h_by_period`); shorter scheme labels; `summary()` no longer
  repeats the print hint; `(detected)` only when the scheme was auto-detected; `rddid()`
  checks `comparisons`, `level` and the period column's type and says so.
* **Degenerate panels**: when no unit changes side of the cutoff between the periods (every
  unit has the same type in every period, e.g. `rddid_sim`), the four tests now stop with one
  clear message instead of returning a mechanical rejection, a `chi-squared(0)` or internal
  errors. Every other numerical result is unchanged.
* **Site and help pages**: *Get started* reads the four test printouts line by line, says what
  the period column may be (`trend = "linear"` needs numeric periods), estimates on
  `rddid_sim_pv` before testing its assumptions, defines a unit's type once in a table, and
  shows how to export a table; the ATU section moved after the tests, and "What to report"
  became "What the paper reports". Advice to the researcher is rewritten as description
  across the articles and help pages, with "can be biased" throughout; the joint test over
  pairs of `rd_compstable()` is flagged as approximate wherever it is tabulated. Articles menu:
  the tests, the plots (now *Plots of the estimate and the checks*, with a `ggsave()`
  example), the options, then *How rddid() computes the estimate* (for referees). The README
  links the tests and plots articles and `citation("rddid")`.

## 2026-10-08: plots

* `plot()` methods, mirroring the validation figures of the paper's application: `plot(fit)`
  draws the per-period RD plots behind an `rddid()` estimate (binned outcome, local-linear
  fits on each side, the jump at the cutoff); `plot()` of an `rd_typecont` or `rd_compstable`
  result draws the binned share of each type with its fitted lines on both sides of the
  (artificial) cutoff; `plot()` of an `rd_homog` or `rd_trendcell` result draws the within-type
  confounding jumps with 95% intervals by period. `plot_switchers()` draws the running
  variable in one period against another and counts the units that change side. All return
  ggplot objects (ggplot2 in Suggests). New article *Plots of the estimate and the checks*.
* To feed the plots, the objects now also carry what the methods read: `rddid()$data` (the
  per-period data used), `rd_typecont()$fits`/`$data`, each `rd_compstable()` pair's `fits`
  and `sample`, and `rd_typecont()$meta$c`. Nothing else changed.

## 2026-10-08: five bug fixes found during the UX sweep

* `rd_homog()` / `rd_trendcell()` with `p > 1`: the per-cell CCT bandwidth is now chosen for
  the order-`p` fit (it used to be chosen at `p = 1` whatever `p` the fit used). Default calls
  (`p = 1`) are unchanged.
* `rd_trendcell(trend = "linear")`: the second differences are now taken in time, using the
  period values, so a jump that is linear in time is annihilated whatever the spacing of the
  comparison periods (the rows used to be index-based and could reject a true linear trend with
  unequally spaced periods). Equally spaced periods give exactly the same contrasts as before.
  Non-numeric period values now error under `"linear"`.
* `rd_compstable()`: a type whose local-linear fit fails is skipped as documented; when the
  last type (in radix order) failed, the function used to stop with "subscript out of bounds".
* `rd_homog()` / `rd_trendcell()`: when the contrasts' covariance is zero to working precision
  (no residual variation within the cells, e.g. a constant outcome) the functions now stop with
  a clear message instead of reporting a chi-squared statistic made of rounding noise.
* `rd_homog()`: the `reference` column of `period_type_jumps` marks the type actually used as
  the reference in each period (the all-below type when its cell was fitted, else the next
  type in order); it used to flag the all-below type only, leaving some periods with no
  reference row.
* Numerical results are unchanged except `rd_homog(p = 2)` (first fix).

## 2026-10-08: source files renamed by what they hold

* **Files are named after what they hold.** The four assumption tests live in `rd_typecont.R`,
  `rd_compstable.R`, `rd_homog.R`, `rd_trendcell.R` (they used to be `test_*.R`, which read as
  unit tests); their shared internals in `assumption_tests_helpers.R`; the covariance and
  aggregation in `cross_period_covariance.R`; the bandwidth rules in `bandwidth_cct.R` and
  `bandwidth_joint.R`; the trend weights and scheme detection in `trend_weights.R` and
  `sampling_scheme.R`. Every file opens with a header saying what it holds and who calls it.
* **Functions read as pipelines.** `rddid()` and the four tests are short sequences of named
  steps (validate, periods, types, scheme, fit cells, covariance, Wald, output), each a small
  helper with a one-line note; the side-of-the-cutoff fit of `rd_period()` is `.rd_side_fit()`;
  the coordinate-descent bandwidth rule is split into objective, start and descent.
* **Names and comments.** The cutoff is `cutoff` inside every function (`c` stays the argument
  name, as in rdrobust); locals that collided with the scheme codes or base functions were
  renamed; comments say *why* at the numerically sensitive lines (the active-set floor, the
  bias-correction matrix grouping, the HC1 factors, the radix sorts, the reference-type drop,
  the summation order of the pair loops) and dated history moved here.
* **Nothing numerical changed.** Every statement that touches a number is the same, in the same
  order; a snapshot of 134 calls (7,763 values) is identical before and after every commit of
  the sweep. Column-check errors no longer carry an `Error in <fn>`
  prefix (the message is unchanged).

## 2026-10-08: UX sweep, step 3 — documentation and site

* Every help page rewritten in one vocabulary with a runnable example on the shipped data;
  package help page `?rddid-package`; `citation("rddid")`.
* Vignettes: *Get started* reaches `rddid()` in four lines of code and reads the printout line
  by line; *How rddid() computes the estimate* holds the hand-built reconstructions; *Bandwidth
  rules and sampling schemes* and *Checking the identification assumptions* rewritten on the
  shipped data. Reference index grouped Estimate / Check the assumptions / Building blocks / Data.
* Example data `rddid_sim` (running variable fixed) and `rddid_sim_pv` (moving; composition
  stability fails there by design while the other assumptions hold).

## 2026-10-08: UX sweep, step 2 — the API

* **Methods.** `summary()`, `coef()`, `confint()`, `nobs()` for `rddid` objects, and broom-style
  `tidy()`/`glance()` for fits and `tidy()` for the four validation tests (registered for the
  `generics` package, so `modelsummary` tables work).
* **Printing.** `print()` of a fit now shows the estimand and period in words, the sampling scheme
  in words, the bandwidth actually used, and z and p-values; the per-period fits and the standard
  error under every sampling scheme moved to `summary()`. The four tests print the same way:
  the assumption, its null in words, then the Wald statistic.
* **Fields added (nothing removed or renamed):** `rddid()` returns `call`,
  `bandwidth$h_by_period`/`$b_by_period` (the bandwidths used in every period, whatever the
  rule) and `estimates$z`/`$p`; `rd_typecont()` and `rd_compstable()` return top-level
  `statistic`, `df`, `p_value`, `scheme`, `estimand`, `call` like the other two tests.
* **Arguments.** `rddid(weights=)` is now `trend=` (the old name still works, with a message):
  it names the assumption on the confounding jump, and `weights` means observation weights
  elsewhere in R. `rd_typecont()` gains `t_rd`/`comparisons` so the same call works for every
  function; `rd_homog()`/`rd_trendcell()` take `p`, `q` explicitly instead of `...` (passing
  `b` through `...` used to fail with "argument matches multiple formal arguments"), and list
  `t_rd` before `comparisons` like `rddid()`. `rd_typecont()`/`rd_compstable()` no longer
  silently accept unknown arguments. `kernel` is validated on entry everywhere.
* **Messages.** `rddid()` says so when no `id` is given (rows are then treated as separate
  units); `rd_homog()`/`rd_trendcell()` list the cells they skip for `min_n`.
* Numerical output is unchanged: every value of the step-1 snapshot is identical.

## 2026-10-08: UX sweep, step 1

* **The default bandwidth rule is now the common bandwidth, `bwselect = "joint"`** (one `h` for
  every period, minimising the asymptotic MSE of the aggregate estimator). The iterative
  period-specific rule `"iter"` remains available but is no longer the default and is not used in
  the paper. Calls that pass `bwselect` explicitly are unaffected.
* **`rdrobust` moved from Suggests to Imports.** `rd_bw_cct()` no longer falls back to
  `0.5 * IQR(x)` when rdrobust is not installed (the fallback on an `rdbwselect()` error is
  unchanged), so the default path gives the same numbers on every machine.
* **The composition-adjusted family is removed**: `rd_att()`, `rd_sadjust()`, `rd_c()`,
  `rd_adjust()` and their tests. They belong to a companion paper, not to the RD-DID paper, and are
  recoverable at the git tag `v0.4.0.9000-composition`.

## Earlier in 0.4.0.9000 (before the 2026-10-08 sweep; the default rule was `"iter"` then)

* **`bwselect = "joint"` no longer depends on which period is labelled `t_rd`.** The common
  AMSE-optimal bandwidth used to fit *every* period at the RD period's CCT pilot pair to estimate
  the aggregate bias and variance constants, so an estimator that aggregates several RD periods
  (the same linear combination of per-period discontinuities, whichever RD period carries the
  `t_rd` label) got a different `h*` for each labelling (Grembi, aggregate ATU over 2001-2004:
  482 vs 320). Both joint rules now estimate each period's constants at that period's **own** CCT
  pilot through one shared helper (`.bw_constants()`), and the common `h*` is the exact scalar
  minimizer of the iterative rule's objective (pinned in `test-appB-bandwidth.R`, together with
  a `t_rd`-relabelling invariance test). Under `"joint"` the pilot `b_t` now keeps each period's
  CCT ratio, `b_t = h* b_t^CCT / h_t^CCT` (was the RD period's ratio for all periods);
  `rddid()$bandwidth$b` is therefore a named per-period vector. `"joint"` numbers move (a few
  percent in the snapshot DGPs; Grembi Table 1 common-`h*` rows change); `"iter"` moves only
  through its seed (below 1e-4 relative); `"cct"` is unchanged.

* **New argument `estimand = c("att", "atu")`** on `rddid()`, `rd_typecont()`, `rd_compstable()`,
  `rd_homog()` and `rd_trendcell()` (default `"att"`, nothing existing moves). Set `"atu"` when the
  comparison periods are uniformly *treated* (the paper's Section 6: the ATU design is the ATT
  design with the sides of the cutoff exchanged). For `rddid()`, `rd_typecont()`, `rd_homog()` and
  `rd_trendcell()` the estimates and tests are numerically identical, so `estimand` only labels the
  output (pinned by `test-mirror-invariance.R`). `rd_compstable(estimand = "atu")` mirrors the
  running variable (`x -> c - x`) and tests composition stability on the **below**-cutoff shares
  (type indicator unchanged, "above the cutoff in the other period", so the jump is
  `pi_{t_RD,(-)}(1) - pi_{t_0,(-)}(1)`);
  observations at `x == c` error under `"atu"` (place the cutoff between support points).

* **Canay-Kamat permutation test and McCrary tests removed** from `rd_typecont()` and
  `rd_compstable()` (arguments `q`, `S`; outputs `ck_perm`, `mccrary_within`, `mccrary_pooled`;
  internal helpers `.q_rot()`, `.mccrary()`). The paper reports the local-linear Wald tests only.
  The remaining outputs are numerically unchanged.

* `rd_homog()` default `type_by` is now `"rd_side"` (the unit's side of the cutoff in the RD
  period, the partition of the paper's Section 4.4), matching `rd_trendcell()`; was `"pattern"`.

* Print methods echo option values: `rddid()` prints `bwselect: iter (5 iterations)` and
  `scheme: pc (auto-detected)` (were `bandwidth: period-specific joint AMSE ...` and
  `sampling scheme: PC`), and `Robust SE by scheme: cs= pc= pv=`; the test print methods drop
  the LaTeX assumption labels and the necessary/sufficient notes.

* Documentation pass after a user-focused review of the site (`rd-did/docs/reviews/
  2026-09-09_oren-persona_rddid-site.md`): what a comparison period is and the scheme detection
  rule stated on the home page / in `rddid()`; the tests vignette opens with its scope (a
  time-varying running variable); `rd_period()` leads with its `rdrobust` equivalence;
  `rddid()`'s Value lists the element names; assumption references by name, not label; the
  normative wording ("suggestive", "neither necessary nor sufficient", "preferred") removed
  throughout.

* **`rddid()` default bandwidth rule is now `bwselect = "iter"`** (period-specific bandwidths by
  coordinate descent on the aggregate AMSE, started at the common h*), the rule Section 5.3 of the
  paper states as preferred; previously `"joint"`. Calls that pass `bwselect` explicitly (all the
  paper's scripts do) are unaffected.

* **Site rebuilt around three vignettes** (`vignettes/`, plan in `dev/site_plan.md`): *Estimating
  the ATT with rddid()* (per-period `rd_bw_cct()` + `rd_period()`, then `rddid()` under constant,
  linear and custom weights), *Bandwidth rules and sampling schemes* (`bwselect = "cct" / "joint" /
  "iter"`, `start`, fixed `h`; `scheme` detection and the three standard errors) and *Composition
  validation tests* (`rd_typecont()`, `rd_compstable()`, `rd_homog()`, `rd_trendcell()` on
  simulated scenarios that satisfy or violate each assumption). Every vignette hand-codes its DGP
  and prints the truth next to each estimate; `dev/site_dgp_check.R` verifies the scenario claims.
  `_pkgdown.yml` groups the reference into *Estimation* and *Validation tests*; README gains a
  quick start (now generated from `README.Rmd`).

* The composition-adjusted estimators (`rd_adjust()`, `rd_sadjust()`, `rd_c()`, `rd_att()`) are
  tagged `@keywords internal`: still exported and tested, but no longer listed on the site. They
  belong to a companion paper and are not part of the current manuscript.

* `Suggests` gains `knitr`, `rmarkdown`, `ggplot2`; `VignetteBuilder: knitr`.

# rddid 0.3.5.9000 (development)

* **Synced with Appendix B of the paper (rewritten 2026-09-08).** A code <-> equation map now
  lives in `dev/appB_map.md`: one row per object the package computes, with the paper's exact
  expression, the implementing symbol, conventions, audit status, and the test that pins it.
  `dev/check_appB_labels.R` verifies that every cited paper label still exists in `main.tex`;
  `dev/snapshot_rddid.R` is a numerical regression snapshot.

* Bandwidth selectors follow App. B.4 at general polynomial order `p` (previously the
  exponents and constants were hard-coded for `p = 1`): `rd_period()$b_const` is
  `(p+1)! (D - D_bc) / h^{p+1}`; the common `h*` (`eq:common_h_opt`) uses the constant
  `((p+1)!)^2 / (2(p+1))` and exponent `1/(2p+3)`; the period-specific objective and its
  regularization use `h^{p+1}/(p+1)!`. Numerically identical at `p = 1`.

* `bwselect = "iter"` under `scheme = "pc"`: the same-side cross-period covariance term of the
  aggregate AMSE now scales as `omega(h_t/h_s)/h_s` per side, with `omega(rho)` the kernel
  constant of Lemma `cov-pc` (new internal module `R/kernel_constants.R`, ported from the paper's
  simulation-verification code), instead of the previous `1/max(h_t, h_s)` approximation (exact
  only for the uniform kernel at `p = 0`). `.bw_joint_iter()` additionally returns the objective
  value and the objective function for diagnostics.

* `rd_period()`: the active set is the union of the pilot and main windows (was the pilot window
  alone, which silently truncated the main fit when `b < h`). Identical whenever `b >= h`.

* New tests: `test-appB-conformance.R` (from-the-equations reference implementation of the
  estimator, bias correction and the three sampling-scheme variances, 1e-10),
  `test-appB-bandwidth.R` (B.4 objectives and selectors), `test-kernel-constants.R`.

* Regularization of the bandwidth selectors: the curvature-variance estimate `Var(B-hat_t)`
  now pairs the pilot-window influence weights with the residuals of the order-`q` pilot fit at
  `b` (the same convention as the bias-corrected variance) instead of the order-`p` residuals at
  `h`. Moves `bwselect = "joint"`/`"iter"` bandwidths by a fraction of a percent at the CCT pilot
  ratio; at wide pilots (`b/h >= 3`) the old estimate over-stated the variance.

* Validation tests synced with Section 4.4 of the paper (`dev/tests_map.md`): documentation
  cites the assumptions by label (`ass:type-cont`, `ass:comp-stable`, `ass:homog`,
  `ass:trend-cell`) instead of numbers that changed in the paper; `rd_compstable()` drops one
  reference type instead of pseudo-inverting the structurally singular covariance (same
  statistic for binary types, exact df); a comparison-period unit exactly at the cutoff now
  stays on the reflected side; the joint-over-pairs result is documented as approximate. New
  from-the-text conformance tests `test-s44-conformance.R`.

* Guards: `bwselect = "cct"` and the joint pilot go through the guarded `rd_bw_cct()`
  (fallback + finiteness checks) instead of calling `rdrobust::rdbwselect()` directly;
  `.bw_joint()` warns when `h*` exceeds the running-variable radius; the coordinate descent
  warns when a bandwidth ends on the search boundary; `rddid(scheme = "pc"/"pv")` warns when
  no unit id repeats across periods (all cross-period covariances are then zero); the search
  cap is NA-safe.

* Stale references in code comments to the removed "Appendix C", `lem:coercive` and
  `eq:amse-ps` replaced by the current labels (`app:est-bw`, `eq:amse-att`, `eq:update`,
  `alg:coorddesc`).

# rddid 0.3.0.9000 (development)

* `rddid(..., bwselect = "iter")` gains a `start` argument controlling the
  coordinate-descent seed: `"hstar"` (default) seeds all periods at the common
  joint-optimal bandwidth h*; `"cct"` seeds each period at its own CCT/IK pilot
  h; or supply a named numeric/list of per-period bandwidths for a manual seed.
  Each run weakly improves the joint AMSE over its own start; seeding from h*
  therefore weakly dominates the common-h rule. (Seed from `"cct"` when the
  per-period biases nearly cancel, where h* is inflated.)

* Bug fix (reproducibility): `rd_typecont()` called `set.seed(NULL)` immediately
  before its Canay-Kamat permutation, which **re-initialised** the RNG from
  system entropy and discarded any seed the caller had set — so the CK p-value
  changed on every call. Removed; the permutation now inherits the caller's RNG
  state, so `set.seed()` before the call makes the CK p reproducible.
  (`rd_compstable()` never had this and was already reproducible.)

* `rd_typecont()`: the Canay-Kamat permutation test now chooses the number of
  nearest observations per side, `q`, by the Canay & Kamat (2018) rule of thumb
  **by default** (`q = NULL`), per period. A fixed `q` over-rejects in finite
  samples when the type distribution varies steeply in the running variable at
  the cutoff; the rule of thumb shrinks `q` as that association strengthens.
  Pass an integer `q` to force a fixed value; the per-period `q` used is
  returned in `meta$q_used`.

* `rd_compstable()`: same change — its Canay-Kamat permutation shares the
  identical fixed-`q` exposure, so `q` now defaults to the rule of thumb
  (`q = NULL`), chosen per `(t_RD, t_0)` pair on the pooled reflected sample.
  The per-pair `q` used is returned in `meta$q_used` (and echoed in each
  `pairs[[...]]$q`). Pass an integer `q` to force a fixed value.

* Bug fix (numerical robustness): the joint Wald pseudo-inverse (`.joint_wald()`,
  used by `rd_typecont()` and `rd_compstable()`) now uses the `MASS::ginv`
  relative tolerance `sqrt(eps)*max(sv)`. The previous, tighter tolerance could
  leave the structural-zero singular value (the per-period type indicators sum
  to 1) just above the cut on some LAPACK builds, inflating the statistic into a
  platform-dependent false rejection. Caught by the new CI on Ubuntu.

* Continuous integration: added an `R-CMD-check` GitHub Action (standard
  multi-OS matrix) and R-CMD-check / MIT-license badges to the README. Removed
  a placeholder ORCID from `DESCRIPTION`.
* `R CMD check` is now clean (0 errors / 0 warnings). Fixes: replaced the
  `\insertCite{}` macros in `rd_compstable()` (Rdpack was not a dependency)
  with the plain-text citations already in the References; documented the
  `regularize`/`reg_const` arguments of `rddid()`; and dropped the
  `VignetteBuilder: knitr` declaration (and the `knitr`/`rmarkdown` Suggests)
  since the package ships no vignettes.

* Internal: consolidated duplicated logic into shared helpers. A new
  `.cov_scheme()` (the `cs`/`pc`/`pv` scheme-combine of `.cross_cov()`) now
  backs the cross-period covariance in `rd_typecont()`, `rd_compstable()`, and
  `rd_homog()`, replacing three near-identical inline copies (including the
  former `.cross_cov_homog()`/`.match_sum_homog()`). A new `.scheme_from_long()`
  primitive backs both `.detect_scheme()` and `rd_typecont()`'s scheme
  detection. No change in results (verified to machine precision).

# rddid 0.2.1

* Internal: removed a duplicate `.build_types()` (a second copy lived in
  `test_homog.R` and shadowed the canonical one in `test_helpers.R` at load
  time). There is now a single shared implementation used by `rd_typecont()`
  and `rd_homog()`.
* Convention: units exactly at the cutoff are now treated as above it
  (`V_i = 1{R_i >= c}`) everywhere, including sampling-scheme detection. A unit
  sitting on the cutoff no longer registers as a separate "side" and so cannot
  be misread as a side-switch.
* Internal: dropped the unused single-cell `.ck_perm()` helper (and its test).
  `rd_typecont()`/`rd_compstable()` use an inlined *joint* Canay–Kamat
  permutation with one shared per-period shuffle; the standalone helper was a
  dead parallel path.
* Documentation: corrected the assumption numbering in the Section 3.4 test
  functions to match the manuscript's `\begin{assumption}` ordering —
  `rd_typecont()` is Assumption A7 (was mislabelled A6), `rd_compstable()` is
  A8 (was A7), and `rd_homog()` is A9 (was A8). Affects titles, `print` output,
  and cross-references only; no behaviour change.
* Documentation: fixed an unmatched apostrophe in the `rd_homog()` `@examples`
  comment that was silently dropping the entire example from the rendered help
  page.

# rddid 0.2.0

* Added tests for the Section 3.4 identifying assumptions of the
  time-varying-running-variable design:
  * `rd_typecont()` — continuity of the type distribution (LL-Wald and
    Canay–Kamat permutation, both necessary & sufficient; McCrary within-type
    sufficient-not-necessary; McCrary pooled neither).
  * `rd_compstable()` — composition stability across periods, via the reflection
    construction (LL-Wald and Canay–Kamat permutation, both necessary &
    sufficient). The permutation uses the partially-overlapping-samples scheme
    for units above the cutoff in more than one period.
  * `rd_homog()` — type-homogeneous confounding, tested in comparison periods
    (only suggestive of the assumption at the RD period: neither necessary nor
    sufficient there).

# rddid 0.1.0

* From-scratch rewrite: per-period local-linear RD engine (`rd_period()`,
  validated against `rdrobust` to machine precision), the `rddid()` aggregate
  estimator with constant/linear/custom weights, CS/PC/PV sampling-scheme
  variances, and joint / CCT / period-specific bandwidth selection.
