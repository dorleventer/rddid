# `estimand = c("att", "atu")` — plan

*Written 2026-09-10 (RA-project task T1507). Status: APPROVED 2026-09-10 (Dor: plan as written; ATU on the site "just this once, to show how to use this feature"). §7 decisions = the first option in each.*

Goal: let a user say that the comparison periods are uniformly **treated** (the paper's new Section 6, "Identification and estimation of the ATU"), so that the package (a) labels the estimand correctly and (b) runs the composition-stability test on the side the ATU requires. The paper's statement: mirror the design ($\tilde R=c-R$, $\tilde W=1-W$, $\tilde V=1-V$), apply the ATT procedure as it stands; the estimator returns $-\widehat{\atu}$ with identical SEs and bandwidths, three tests are unchanged, the composition-stability test uses the below-cutoff units.

## 1. What the package actually has to do

The package has no `W`/`V` columns; sides are derived from `x >= c`. So "mirroring" inside the package is `x <- c - x; c <- 0`, and nothing else. Working through each exported function on mirrored data (paper §6; RA-project session 2026-09-10-02 derivation):

| Function | On mirrored data | Reported in original orientation? | Numerically different from `"att"`? |
|---|---|---|---|
| `rddid()` | $\widehat{\tilde D}_t=-\widehat D_t$ exactly (LL intercepts, residuals, pilot second derivative, and the bias design constant $B$ are invariant to reflecting $x$ within a side); est $=-\widehat{\atu}$; SEs, cross-period covariances, all three bandwidth rules identical | after negation: identical to the `"att"` run | **No** |
| `rd_typecont()` | outcome $1-\mathbf 1\{V_s=1\}$ on $-R$: jump $(1-\hat\beta_-)-(1-\hat\beta_+)=\hat D$; same SE | identical | **No** |
| `rd_homog()` | $\tilde D_{t_0}(\tilde v)=-D_{t_0}(1-\tilde v)$; types relabel; equality contrasts and Wald identical | after negation + relabel: identical | **No** |
| `rd_trendcell()` | as `rd_homog()` | identical | **No** |
| `rd_compstable()` | the recipe's "above-cutoff units" are the original **below**-cutoff units; jump $=\tilde\pi_{t_\mathrm{RD},(+)}(1)-\tilde\pi_{t_0,(+)}(1)=\pi_{t_\mathrm{RD},(-)}(0)-\pi_{t_0,(-)}(0)$ | reported as is (share of below-cutoff units that are below in the other period) | **Yes** |
| `rd_period()`, `rd_bw_cct()` | single-period primitives, no estimand | — | not touched |
| `rd_adjust`, `rd_sadjust`, `rd_c`, `rd_att` | composition-adjusted (internal, companion paper) | — | not touched; documented ATT-only |

So the *only* computation that changes is `rd_compstable()`. Everywhere else `estimand` is a label, and the invariance is a theorem that the test suite pins numerically.

## 2. Design (recommended)

**Argument.** `estimand = c("att", "atu")`, default `"att"` (no existing number moves), added to `rddid()`, `rd_typecont()`, `rd_compstable()`, `rd_homog()`, `rd_trendcell()`. Placed after `comparisons` (tests) / after `weights` (`rddid()`). Stored in the returned object (`$estimand` for `rddid`, `$meta$estimand` for the tests).

**`rddid()`, `rd_typecont()`, `rd_homog()`, `rd_trendcell()`: label only.** No mirroring in code. Print methods show `ATU(t_RD)` / a one-line "estimand: ATU (comparison periods uniformly treated)". Rd `@details` state the identity and point to the paper's Section 6 and to the invariance tests. Rationale: mirroring and translating back (`fits` per period carry `sides[["+"]]/[["-"]]`, slopes, per-side residual vectors) adds surface with zero numerical content; the claim is instead *verified* by tests that mirror the data by hand and compare (§4).

**`rd_compstable(estimand = "atu")`: real mirror.** Internally `x <- c - x; c <- 0` before the existing construction; everything downstream unchanged. Output fields keep their names; the Rd documents both readings of `jumps`, `n_trd`, `n_t0`, `n_both` ("above-cutoff" ↔ "below-cutoff" under `"atu"`). Print says which shares were tested.

**Ties at the cutoff under `"atu"`.** Units with `x == c` are treated in the original design (`>= c`) but land at $\tilde x=0$ and would be read as mirrored-above (= original untreated). Exact handling needs a strict side rule threaded through `rd_period()`'s split; instead: `stop()` with guidance — *"estimand = 'atu': units at x == c are treated in the original design; place the cutoff between support points (e.g. c = 4999.5 for integer populations) so no unit sits on it."* Same advice as for any discrete running variable; Grembi has no ties (checked: 0 rows at 5,000). Note for later: the reflected construction already has a tie edge case under `"att"` (`test_compstable.R:235` comment); unchanged here.

**Not exported:** no `rd_mirror()` helper (site rule: omit rather than add). A one-line mirror helper lives in the test helpers only.

Alternative considered: mirror inside every function and translate outputs back. Rejected for the reason above; revisit only if a user-facing per-period object needs mirrored orientation.

## 3. Files touched

- `R/rddid.R`, `R/test_typecont.R`, `R/test_homog.R`, `R/test_trendcell.R`: argument, `match.arg`, stored field, print label, roxygen `@param estimand` + `@details`.
- `R/test_compstable.R`: argument, tie check, mirror lines, roxygen (both readings), print.
- `man/*.Rd` regenerated (pre-commit hook).
- `tests/testthat/test-estimand.R` (new; §4), `helper-mirror.R` (new, 3 lines).
- `NEWS.md` bullet under the current dev version (bump is Dor's call, §7).
- `dev/appB_map.md`: one row (paper §6 mirror ↔ `rd_compstable(estimand = "atu")`, "label only" rows for the others); `dev/tests_map.md`: the new test file; `dev/README.md`: list this plan.
- Vignettes and `README.Rmd` (§5), `_pkgdown.yml` article title.

## 4. Verification (before anything ships)

**Automated — `tests/testthat/test-estimand.R`.** Data: the vignettes' DGP-A (PC, 3 periods) and DGP-B S1 (PV, drift) at a fixed seed; the mirror is applied by hand in the test (`x <- -x`, `c = 0`), independent of the package's `estimand` code.

1. `rddid()`: for `bwselect` in `cct`, `joint`, `iter` and `scheme` in `cs`, `pc`, `pv`: `est(mirrored) == -est(original)`, `se_cs/se_pc/se_pv` equal, per-period `h_t`, `b_t` equal, Conventional and Robust rows. Tolerance `1e-10` on estimates/SEs; `1e-8` on bandwidths (rdrobust's `rdbwselect` on reflected data — if it is not symmetric to `1e-8` that is a finding to record, not to paper over).
2. `rd_typecont()`: `ll_wald$stat`, per-period jumps and SEs equal (conv and BC).
3. `rd_homog()`, `rd_trendcell()`: statistics equal; `period_type_jumps` negated with type labels swapped.
4. `rd_compstable(estimand = "atu")` on original data `==` `rd_compstable(estimand = "att")` on hand-mirrored data (exact); `!=` `estimand = "att"` on original data for S1 (they test different units); equals a from-the-equations below-side reflected fit built in the test with `rd_bw_cct` + `rd_period` (point estimate; the package SE carries the id-matched cross-side term, so compare $z^2$ to the package `ll_wald$stat` as the existing S5 cross-check does).
5. Ties: `estimand = "atu"` with any `x == c` errors with the guidance message; `"att"` unaffected.
6. `estimand = "att"` default: `dev/snapshot_rddid.R` byte-identical; full suite (currently 2252) passes; `R CMD check --as-cran` clean.

**Manual, on Grembi (in `rd-did/code/application/`, not part of the package).**
7. `s6_estimates.R` with `estimand = "atu"` on every `rddid()` call: `outputs/s6_main.tex`, `s6_numbers.tex` byte-identical.
8. `s5_application.R`: A7/A9/A10 macros byte-identical with `estimand = "atu"`; `rd_compstable(estimand = "atu")` jump/SE/$\chi^2$ reproduced by an independent below-side reflected fit written in the script (mirror of the existing above-side cross-check at `s5_application.R:556-575`); the old above-side numbers kept in the script log for the record.
9. Read the printed output of all five functions with `estimand = "atu"` on Grembi and check every label by eye.

## 5. Site and README

The site plan's rule "ATT framing throughout; **no ATU mention**" (2026-09-09) predates the paper's Section 6 and needs Dor's explicit lift (§7). Proposed, minimal:

- **`rddid-estimation`** (Get started): one closing section *"Comparison periods that are uniformly treated"* — two sentences (the paper's §6 in one line: same estimator, label ATU; the composition-stability test changes side), one `rddid(..., estimand = "atu")` call on a DGP-A variant with $W\equiv1$ in the comparison periods and a $W\times V$ interaction so ATT $\ne$ ATU at the cutoff, truth printed next to the estimate.
- **`rddid-validation-tests`**: one section after `rd_compstable()` — *"ATU designs"*: `rd_compstable(S1, ..., estimand = "atu")` next to the `"att"` call, showing the two jumps and that they are estimated on different units; one sentence that the other three tests are unchanged (shown by calling one with `estimand = "atu"`).
- **`rddid-options`**: nothing (bandwidths and schemes are estimand-invariant; one sentence saying so in the Sampling schemes section, optional).
- **README quick start**: one line: `estimand = "atu"` when the comparison periods are uniformly treated.
- `_pkgdown.yml`: article title "Estimating the ATT with rddid()" → "Estimating the ATT or ATU with rddid()".
- `dev/site_plan.md` §0: amend the no-ATU rule with a dated note pointing here.

## 6. Order of work (one chunk per agent, commit per chunk, lead verifies each)

1. Package code + Rd + NEWS (implementer; contract = the five signatures and the tie error message; self-test = existing suite passes).
2. `test-estimand.R` (separate implementer; contract = the six checks of §4; PASS lines).
3. Lead: run suite + `R CMD check`; manual Grembi checks 7–9; record results here (§8).
4. Vignettes + README + pkgdown (implementer; after Dor's decision on §5); lead renders the site locally and reads the new sections.
5. Push; CI green; live site checked.
6. Then the paper side (separate task step): rerun `s5_application.R`, regenerate the A8 macros/table/figure panel (b), rewrite the §7.3 composition-stability paragraph and the summary.

## 7. Decisions for Dor

1. Argument name: `estimand` (proposed) vs `target`.
2. Label-only for `rddid()` and three tests + real mirror in `rd_compstable()` (proposed), vs mirror everywhere.
3. Ties under `"atu"`: error with the between-support-points guidance (proposed), vs exact strict-side handling.
4. Lift the site's no-ATU rule for the two sections in §5.
5. Version: add the NEWS bullet under 0.4.0.9000 (proposed) vs bump to 0.5.0.9000 for the new argument.

## 8. Results log

*(filled in as steps of §6 complete)*
