# rddid site rebuild — plan

*Written 2026-09-09 (RA-project task T1503); revised the same day after Dor's decisions (§7). Status: APPROVED, building.*

Goal: the pkgdown site (`https://dorleventer.github.io/rddid/`) shows how to use each function of the package on hand-coded simulated data with known truth, following the arc of the paper's Section 6 (per-period discontinuities → RD-DID estimate → bandwidth rules → standard errors by sampling scheme → composition validation tests).

## 0. Audience and style (Dor, 2026-09-09)

- The reader **has read the paper and understands it**; the vignettes show off the package functions, they do not re-teach the method. No derivations, no equations for the bandwidth rules, no lemma-level exposition.
- **Simple as possible.** Short sections; one dataset per vignette where possible; one idea per section.
- **Notation as in the paper, with a one-line reminder at first use** (e.g. "`D_t` is the observed outcome discontinuity in period `t`", "`alpha_{t,0}` the confounding discontinuity", "`w_t` the comparison-period weights", "CS / PC / PV = repeated cross-section / panel with time-constant / time-varying running variable", "a unit's type in period `t` is its side of the cutoff in the other period").
- **Descriptive, not normative.** "For standard errors under the PV sampling scheme, set `scheme = "pv"`", not "you should use PV". Each bandwidth rule: what it is, how to call it, how to read the output — no ranking.
- ATT framing throughout; **no ATU mention**. Generic periods `t = 1, 2, 3`.
- **No composition-adjusted estimators anywhere on the site** (`rd_adjust`, `rd_sadjust`, `rd_c`, `rd_att`): tagged `@keywords internal` so they stay exported and tested but leave the reference index; not mentioned in any vignette or the README.

## 1. Where the site stands

- The pkgdown workflow is green on every push to `master` and deploys to `gh-pages`; the live site is current with the code but has **no usage content**: home = README (install + dev hook only), a flat reference index, no articles.
- The two Feb-2025 articles were deleted in the June 2026 rewrite (`1c16b2a`); they used an API that no longer exists. `vignettes/` holds only a `.gitignore`; no `_pkgdown.yml`; no `VignetteBuilder`; `Suggests` = `rdrobust, testthat`.

## 2. Package exports → where they appear

| Export | Paper object | Vignette |
|---|---|---|
| `rd_bw_cct()` | single-period CCT bandwidths `h`, `b` | V1 |
| `rd_period()` | per-period local-linear RD, conventional + robust bias-corrected `D_t`; per-side intercept/slope for plots | V1 |
| `rddid()` — `weights` | `eq:est_agg_bc`, `cor:trend` (constant / linear / custom `w_t`) | V1 |
| `rddid()` — `bwselect`, `start`, `regularize`, `h`, `b` | §5.3 rules: per-period CCT, common `h*`, iterative period-specific | V2 |
| `rddid()` — `scheme`, `id` | `eq:var-cs/pc/pv`, auto-detection | V2 |
| `rd_typecont()` | `ass:type-cont` test (side-share RD per period + joint Wald) | V3 |
| `rd_compstable()` | `ass:comp-stable` test (reflected-cutoff RD) | V3 |
| `rd_homog()` | `ass:homog` test (within-type comparison-period jumps equal across types) | V3 |
| `rd_trendcell()` | `ass:trend-cell` test (within-type jumps equal across comparison periods) | V3 |
| `rd_adjust`, `rd_sadjust`, `rd_c`, `rd_att` | companion paper | not on the site |

## 3. The paper's illustration, as the template

§6.2: per-period RD table (conv, BC, SE, `h`, `n`) → RD-DID under constant and linear trend with the arithmetic written out → the three bandwidth rules listed with their bandwidths and estimates. §6.3: one paragraph per assumption, always *estimates with SEs → test statistic → reading*, plus the type-share RD plots, the reflected RD, and per-type jumps ± CI. All of it is `rd_bw_cct` + `rd_period` + `rddid()` + the four test functions (`rd-did/code/application/s6_estimates.R`, `s5_application.R`, `s5_validation_figure.R`); the vignettes use the same call recipe.

## 4. DGPs (hand-coded, shown in the vignette that uses them)

Truths are checked by the lead against a large oracle draw *before* parameters are frozen; the vignettes print the true ATT next to every estimate but do not derive anything.

**DGP-A (V1, V2): time-invariant running variable.** `n = 2000` units, `t ∈ {1,2,3}`, comparison periods `{1,2}`, RD period `3`, `c = 0`.
```
R_i ~ U(-1,1);  V_it = 1{R_i >= 0};  W_it = V_it 1{t = 3}
Y_it = m_t(R_i) + alpha_t V_it + tau W_it + u_i + eps_it,   m_t(r) = r + (theta_t/2) r^2 1{r >= 0}
u_i ~ N(0, 0.5^2),  eps_it ~ N(0, 0.5^2),  tau = 0.5
```
`D_t = alpha_t + tau 1{t=3}` exactly. Constant confounding: `alpha_t = 1`. Linear: `alpha_t = 0.5 + 0.5 t`. Curvature `theta_t = (1, 6, 3)` in V2 so the per-period CCT bandwidths differ. Sampling variants in V2: CS (fresh units each period), PC (as written), PV (`R_it = R_i + nu_it`, `nu ~ N(0, 0.1^2)`; composition-free because `alpha` is homogeneous and `u_i ⟂ R`).

**DGP-B (V3): time-varying running variable.** `n = 4000`, `t0 = 1`, `t_RD = 2` (three periods `{1,2}` → `3` for the trend-cell section).
```
eta_i ~ N(0,1);  R_i1 = eta_i + e_i1;  R_i2 = eta_i + d + e_i2;  e ~ N(0, 0.5^2)
V_it = 1{R_it >= 0};  type in period t = side in the other period
Y_it = m(R_it) + alpha(V_is) V_it + tau V_it 1{t = 2} + kappa eta_i + eps_it,   m(r) = r + r^2 1{r >= 0}
```
Scenarios, one dial each: **S0** all assumptions hold (`d = 0`, `alpha(0) = alpha(1) = 1`, `kappa = 0`); **S1** composition stability fails, confounding homogeneous (`d = 0.5`) — estimate still recovers `tau`; **S2** composition stability and homogeneity fail (`d = 0.5`, `alpha(0) = 0.5`, `alpha(1) = 1.5`) — estimate biased; **S3** type continuity fails in period 2 (`kappa = 1`, units just below the cutoff with `V_1 = 1` reflect above with probability `0.5` inside a window of `0.5`; `dgp_s3.R` device); **S4** within-type trend (three periods, `alpha_t(1) = 1 + 0.3 t`, `alpha_t(0) = 1`).

## 5. Vignettes

Conventions: `set.seed` at the top; every number in prose is inline R; true ATT printed next to every estimate; ggplot2 figures in the paper's style (binned means, local-linear lines ending at the bandwidth); runtime ≤ 60 s each; a hidden `stopifnot` gate on truth-vs-estimate (deterministic at the seed).

### V1 `rddid-estimation.Rmd` — "Estimating the ATT with rddid()"  *(Get started)*
1. DGP-A code (constant confounding), five lines of context, notation reminder (`D_t`, `alpha_{t,0}`, `tau = ATT(t_RD)`).
2. Per period: `rd_bw_cct()` → `rd_period()`; table `D`, SE, `D_bc`, SE, `h`, `n`; the three RD plots. Reading: `D_1 = D_2 = alpha`, `D_3 = tau + alpha`.
3. The arithmetic (average of the comparison discontinuities, subtract), then `rddid(weights = "constant")` reproduces it; walk through the print (Conventional / Robust rows, SE, CI, the three scheme SEs).
4. Linear confounding: regenerate; `rddid(weights = "linear")`; show `r$weights = (-1, 2)`; the slope-and-extrapolate arithmetic; numeric `weights` in one line. Descriptive note: with two comparison periods the linear trend is just-identified.

### V2 `rddid-options.Rmd` — "Bandwidth rules and sampling schemes"
1. **Bandwidths.** DGP-A with `theta_t = (1, 6, 3)`. For each of `bwselect = "cct"`, `"joint"`, `"iter"`: one sentence on what it is (per-period CCT/IK; one common `h*` for the aggregate estimator; period-specific bandwidths from coordinate descent started at `h*`), the call, the resulting bandwidths (`r$bandwidth`, `r$fits[[t]]$h`) and estimate/SE in one table. Also `start`, `regularize`, and fixing `h`/`b` by hand.
2. **Sampling schemes.** Reminder of CS / PC / PV. Same DGP sampled three ways; `scheme = "auto"` and what it detected; the three SE columns on each dataset; how to set `scheme` explicitly and what `id = NULL` means.

### V3 `rddid-validation-tests.Rmd` — "Composition validation tests"
DGP-B code; the type; scatter of `(R_1, R_2)` coloured by type. Then one section per test, each shown on S0 (holds) and the scenario that breaks it, in the paper's order *estimates with SEs → test → reading*:
1. `rd_typecont()` (S0 vs S3): per-period side-share jumps, joint Wald (conventional and `bc = TRUE`); the two type-share RD plots.
2. `rd_compstable()` (S0 vs S1): reflected RD jump, SE, `p`; the reflected-RD plot; the shared-unit dependence is inside the SE.
3. `rd_homog(type_by = "rd_side")` (S1 vs S2): the within-type comparison-period jumps with SEs (`period_type_jumps`), the equality test; per-type jumps ± CI.
4. `rd_trendcell(trend = "constant")` (S0-3p vs S4): same jumps compared across comparison periods; runs only with ≥ 2 comparison periods.
5. **Summary table**: scenario × p-values × `rddid()` bias, read as in the paper's summary paragraph (which route the estimate rests on).

### README
12-line quick start (DGP-A, `rddid()` call, print) above the install block; dev-hook section kept.

## 6. Site mechanics
- `DESCRIPTION`: `Suggests` += `knitr, rmarkdown, ggplot2`; `VignetteBuilder: knitr`; version `0.4.0.9000`; `NEWS.md` entry.
- `@keywords internal` on the four CA functions (still exported; tests untouched).
- `_pkgdown.yml`: Bootstrap 5; navbar Get started (V1) · Articles (V2, V3) · Reference · News; reference groups *Estimation* (`rddid`, `rd_period`, `rd_bw_cct`) and *Validation tests* (four functions).
- Both workflows build vignettes; all new Suggests are CRAN packages.

## 7. Decisions (Dor, 2026-09-09)
D1 no composition adjustment on the site at all · D2 Monte Carlo in V2: dropped (simplicity) · D3 three vignettes grouped by function family (estimator / estimator options / validation tests) · D4 no ATU mention · D5 generic periods.

## 8. Build sequence
1. Lead: DGP-A/B generators + oracle checks of the scenario claims (S0 all hold; S1 `C = 0`; S2, S3, S4 biases visible at `n`; tests reject/hold at the seed and across 20 seeds); freeze parameters.
2. Lead: scaffolding (DESCRIPTION, `_pkgdown.yml`, keywords, README), local `pkgdown::build_site()`.
3. One vignette per write → render → inspect → fix → commit cycle; build and check as separate dispatches.
4. `R CMD check --as-cran`, `pkgdown::build_site()`, push, both workflows green, live articles opened, brief + task board.
