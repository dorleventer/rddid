# rddid ↔ Section 4.4 map (validation tests)

**Purpose.** Same contract as `dev/appB_map.md`, for the four validation tests of Section 4.4
"Testing the identification assumptions" (`sec:tvrv-assess`) of Leventer & Nevo. Section 4.4
describes each test in prose, without displayed formulas; the rows below quote the operative
sentence, state the estimand and null it pins down, and give the statistic the package computes.
The variance machinery is Section 5.2 / Appendix B.2 (mapped in `dev/appB_map.md` §2.3).

**Last synced.** Paper: `rd-did` commit `c1a9fbf` (2026-09-08; Section 4.4 lean rewrite of
2026-08-28 `445ee8f`, scope sentence and A7 type clause of 2026-09-08). Package: this file's commit.

**Assumption numbers.** The paper renumbered on 2026-09-03; the package's documentation still
carried the old numbers until this pass. Cite by label, not number.

| Label | Name | Paper number now | Old package number |
|---|---|---|---|
| `ass:trend-cell` | Constant confounding discontinuity within types | A7 | A10 |
| `ass:type-cont` | Continuous type distribution | A8 | A7 |
| `ass:comp-stable` | Composition stability | A9 | A8 |
| `ass:homog` | Homogeneous confounding | A10 | A9 |

**Update 2026-09-09:** the permutation and McCrary rows below are REMOVED from the package (see NEWS 0.4.0.9000); `rd_homog()` now defaults to `type_by = "rd_side"`. The restructure note that follows is otherwise unchanged.

**Paper restructure in progress (2026-09-08, 18:06–18:13, uncommitted while this was written).** Session 02 deleted Appendix C (simulation evidence) and the general-$P$ Appendix A, and moved the body proofs to a new Appendix A "Proofs". Once that lands: the McCrary rows below revert to NOT-IN-PAPER (the simulation appendix that ran them is gone), the general-$P$ sign-pattern types and $(\mathbf u,b)$ types are no longer licensed by a paper statement (`ass:type-cont-gen`, `ass:comp-stable-gen` deleted), and the labels `app:sim-val`, `app:sim-val-typecont`, `app:est-adjust` cited in the two maps must be re-pointed. Re-run `dev/check_appB_labels.R` against the new HEAD and update the sync stamps then; this file is pinned to `c1a9fbf`.

**Status legend.** As in `dev/appB_map.md`: `MATCH` (verified by `tests/testthat/test-s44-conformance.R`
at 1e-10 up to the stated conventions), `MATCH*`, `MISMATCH`, `NOT-IN-PAPER` (package computes it, the
paper does not state it), `DECISION`.

---

## 1. The four tests

### 1.1 `ass:type-cont` — `rd_typecont()`

| Item | Paper (Section 4.4) | Package |
|---|---|---|
| Text | "in each period $t$, to run a local-linear RD of the side indicator $\mathbf 1\{V_{i,s}=1\}$ on $R_{i,t}$, which estimates the side-share jump, $\widehat\pi_{t,(+)}(1)-\widehat\pi_{t,(-)}(1)$. One can evaluate whether the jump is equal to zero in each period separately, or in both periods via a Wald statistic. [...] the off-diagonal element in the variance-covariance matrix of the estimators is non-zero in panel data, as the same units enter both periods' regressions." | `.build_types()`: the type of unit $i$ in period $t$ is the sign pattern of its sides in the OTHER periods (`"+"`/`"-"` string); at $P=2$ this is $\mathbf 1\{V_{i,s}=1\}$. For each (period, type value) `rd_period(y = 1{type == v}, x = R_t)` gives the share jump and its influence vectors. |
| Estimand / null | $\pi_{t,(+)}(1)-\pi_{t,(-)}(1)$ per period; $H_0$: all zero | `theta` = stacked jumps (conventional `D` or BC `D_bc`) |
| Statistic | per period: the jump's $t$/$\chi^2(1)$; jointly: Wald over both periods | `per_period[[t]]$ll_wald` ($\chi^2(1)$ with binary types); `ll_wald` = joint Wald `theta' Sigma^{-1} theta`, $\chi^2$(number of kept contrasts); one reference type dropped per period (the indicators sum to 1), so the covariance is full rank — equivalent to the paper's single jump per period at $P=2$ |
| Covariance | "non-zero off-diagonal … same units enter both regressions" | within period, across types: same-side term `.cross_cov()$pc` (same units, same cutoff); across periods: `.cov_scheme()` = `eq:cross-decomp` by sampling scheme (`auto` detection on every typed unit of each period — this pass; was on a rule-of-thumb window unrelated to the per-cell CCT bandwidths, which could under-detect switching; PV keeps all four side pairs) |
| Conventions | none stated | per-cell CCT bandwidth (`bwselect = "cct"`, via `rd_bw_cct()`), HC1 factor, `bc = TRUE` default (BC jumps + BC variance) |
| "Necessary and sufficient" | stated | holds for binary types (one jump per period is the whole assumption); for general $P$ the package tests all sign-pattern cells (`ass:type-cont-gen`, Appendix A) |
| Status | MATCH | pinned by `ass:type-cont —` tests in `test-s44-conformance.R` |

### 1.2 `ass:comp-stable` — `rd_compstable()`

| Item | Paper (Section 4.4) | Package |
|---|---|---|
| Text | "take the units above the cutoff in period $t_0$, center their running variable and flip its sign [...] $-(R_{i,t_0}-c)$. Then, append the units above the cutoff in period $t_{\mathrm{RD}}$, and center their running variable, i.e., $R_{i,t_{\mathrm{RD}}}-c$. [...] for the period-$t_0$ units create the indicator $\mathbf 1\{V_{i,t_{\mathrm{RD}}}=1\}$, and for the period-$t_{\mathrm{RD}}$ units the indicator $\mathbf 1\{V_{i,t_0}=1\}$. Finally, run a local-linear RD with the indicator as the outcome on the constructed data, which estimates $\widehat\pi_{t_{\mathrm{RD}},(+)}(1)-\widehat\pi_{t_0,(+)}(1)$ [...] Since a unit above the cutoff in both periods appears in both groups, the variance of the estimated jump must account for this dependence." | For each pair $(t_{\mathrm{RD}},t_0)$: `xref_trd = R_trd - c` (above-cutoff units of $t_{\mathrm{RD}}$, right of 0), `xref_t0 = -(R_t0 - c)` (above-cutoff units of $t_0$, left of 0); type = partner side ($b$), prefixed by the shared other-period sides $\mathbf u$ when $P\ge3$; `rd_period(y = 1{type == v}, x = xref, c = 0)` on the stacked data, so its "+" side is the $t_{\mathrm{RD}}$ group and its "−" side the $t_0$ group |
| Estimand / null | $\pi_{t_{\mathrm{RD}},(+)}(1)-\pi_{t_0,(+)}(1)=0$ | jump of the indicator at the artificial cutoff (`D` / `D_bc`) |
| Statistic | test the discontinuity $=0$ | `pairs[[pair]]$ll_wald`: Wald on the kept type jumps — the reference type (partner side 0, first in radix order) is dropped when every type fitted (this pass; was a Moore–Penrose pseudo-inverse of the structurally singular full set, T1086); with binary types the kept jump is $\widehat\pi_{t_{\mathrm{RD}},(+)}(1)-\widehat\pi_{t_0,(+)}(1)$ itself, $\chi^2(1)$ = the paper's test. At a common bandwidth drop-one and the pseudo-inverse coincide exactly; with per-type CCT bandwidths the jumps do not sum exactly to zero, so at $P\ge3$ drop-one is a different (still valid) full-rank test |
| Covariance | "must account for this dependence" | `scheme = "auto"` → `"pv"` when any unit is above in both periods: $\mathrm{Var}=V_{(+)}+V_{(-)}-2\,\mathrm{Cov}_{(+,-)}$ with the shared-unit cross term `.match_sum()` of the two sides' influence vectors (the `eq:cross-decomp` opposite-side term); `"cs"` treats the sides as independent. Off-diagonal between two type indicators ($P\ge3$ only): same-side term always (both fits use the same reflected sample), opposite-side term under `"pv"` (this pass; was zeroed under `"cs"`) |
| Conventions | none stated | per-type CCT bandwidth in the reflected space, HC1, `bc = TRUE` default; a $t_0$ unit exactly at the cutoff is kept on the reflected (left) side (this pass; was assigned to the right side by the `x >= 0` split) |
| Status | MATCH | pinned by `ass:comp-stable —` tests |

### 1.3 `ass:homog` — `rd_homog()`

| Item | Paper (Section 4.4) | Package |
|---|---|---|
| Text | "estimate both $\widehat D_{t_0}(0)$ within the subpopulation $V_{i,t_{\mathrm{RD}}}=0$ and $\widehat D_{t_0}(1)$ within the subpopulation $V_{i,t_{\mathrm{RD}}}=1$ via local-linear RD, and test whether the within-type jumps are equal." | `type_by = "rd_side"`: type = side in `t_rd` (labels `"+"`/`"-"`); for each comparison period and type, `rd_period()` on that subsample; contrasts = type minus reference type |
| Estimand / null | $\alpha_{t_0,0}(0)=\alpha_{t_0,0}(1)$, tested as $D_{t_0}(0)=D_{t_0}(1)$ | `contrasts` $=\widehat D_{t_0}(1)-\widehat D_{t_0}(0)$ per comparison period (reference = the all-below type `"-"`, i.e. $V_{i,t_{\mathrm{RD}}}=0$, in locale-independent radix order; before this pass the reference was `sort()`'s first element, whose order for `"+"`/`"-"` depends on the locale, so the exported sign could flip between machines — the Wald never did) |
| Statistic | equality of the two jumps | `statistic` = Wald on the stacked contrasts, `.wald_eigen()` (eigen pseudo-inverse dropping non-positive directions; `df` = retained directions); one comparison period: $(\widehat D(1)-\widehat D(0))^2/(V(1)+V(0))$, $\chi^2(1)$ |
| Covariance | (subsamples are disjoint) | within a period: 0 (disjoint); across comparison periods: `.cov_scheme()` (shared units), so the joint over periods is $\chi^2(P_0)$ — Section 6 reports this joint |
| Conventions | none stated | per-cell CCT, HC1, `bc = TRUE` default, `min_n = 10` per side per cell; default `type_by = "pattern"` uses the full other-period sign pattern (identical at $P=2$; Section 6 uses `"rd_side"`) |
| Status | MATCH | pinned by `ass:homog —` tests |

### 1.4 `ass:trend-cell` — `rd_trendcell()`

| Item | Paper (Section 4.4) | Package |
|---|---|---|
| Text | "a suggestive test can be conducted when multiple comparison time periods are available. For example, if $\mathcal T_0=\{t_1,t_2\}$, and since $D_t(v)=\alpha_{t,0}(v)$ for $t\in\mathcal T_0$, we can test $\alpha_{t_1,0}(v)=\alpha_{t_2,0}(v)$ via $D_{t_1}(v)=D_{t_2}(v)$, with $v$ the period-$t_{\mathrm{RD}}$ side." | cells fixed across comparison periods: `type_by = "rd_side"` (default) = side in `t_rd`; per (cell, comparison period) `rd_period()`; `trend = "constant"`: contrasts = later period minus first period within each cell |
| Estimand / null | $D_{t_1}(v)=D_{t_2}(v)$ for each $v$ | `contrasts` $=\widehat D_{t_2}(v)-\widehat D_{t_1}(v)$, $v\in\{0,1\}$ |
| Statistic | (not specified beyond the equality) | `statistic` = Wald over all cells' contrasts via `.wald_eigen()`, $\chi^2(\#\text{cells}\times(P_0-1))$; Section 6: $\chi^2(2)$ |
| Covariance | not mentioned (Section 5 owns it) | within a cell across periods: `.cov_scheme()` (shared units); across cells: 0 (cells partition units and are fixed over time) |
| Extras | — | `trend = "linear"`: second differences (needs $P_0\ge3$; the paper notes only the constant form is testable with two comparison periods) |
| Status | MATCH | pinned by `ass:trend-cell —` tests |

---

## 2. What the package computes that Section 4.4 does not state

| Package object | Where | Paper | Status / note |
|---|---|---|---|
| Canay–Kamat permutation test | **REMOVED 2026-09-09** (args `q`, `S`; outputs `ck_perm`, `per_period$ck_p`) | dropped from the paper's validation table on 2026-09-09 (`27cf985`); the paper no longer mentions a permutation test | REMOVED |
| McCrary density tests | **REMOVED 2026-09-09** (outputs `mccrary_within`, `mccrary_pooled`, helper `.mccrary`) | not in 4.4 (which argues McCrary "asks a different question"); the simulation appendix that ran them was deleted on 2026-09-08; Section 6's own McCrary $p$-values come from `rddensity` in `s5_application.R` | REMOVED |
| Bias-corrected variants (`bc = TRUE`, the default) | all four | 4.4 silent; Section 6.3 reports the bias-corrected column with the caution paragraph | NOT-IN-PAPER as a definition; the construction is the App. B.1–B.2 BC jump/variance applied to the same regressions (MATCH by composition) |
| General-$P$ types: sign pattern of all other periods (`"pattern"`); $(\mathbf u,b)$ types in the reflection | `.build_types()`, `rd_compstable` | Appendix A (`ass:type-cont-gen`, `ass:comp-stable-gen`) | consistent with Appendix A; identical to 4.4 at $P=2$ |
| Joint over several $(t_{\mathrm{RD}},t_0)$ pairs (`joint$ll_wald` sums the per-pair $\chi^2$; `joint$ck_perm` Fisher-combines) | `rd_compstable` | 4.4 is one pair; Section 6 reports the 2004::2000 pair | NOT-IN-PAPER; the sum assumes independent pairs, which fails when the pairs share the RD-period above-cutoff group (they always do) — documented as approximate this pass; the paper's test is per pair |
| Joint homogeneity over comparison periods | `rd_homog` | 4.4 is one $t_0$; Section 6.3 reports the joint $\chi^2(2)$ with the shared-unit correlation | consistent (Section 6 explains it) |
| `.wald_eigen()` drops non-positive eigen-directions and reduces `df` silently | `rd_homog`, `rd_trendcell` | — | design choice for contrast covariances assembled from estimated pieces; cannot bite at $df=1$; documented |
| Rule-of-thumb bandwidth when `bwselect = "rot"`: $0.5\,\mathrm{IQR}$ (`rd_typecont`, `rd_compstable`) vs $0.2\times$range per cell (`rd_homog`, `rd_trendcell`) | `.cell_bandwidth()` | — | inconsistent fallbacks, non-default paths; documented |

---

## 3. Findings of this pass (2026-09-08)

1. **Documentation used the pre-2026-09-03 assumption numbers throughout** (titles, `print` methods, comments in the four test files, `adjust.R`, `sadjust.R`, `test_helpers.R`). Replaced by names and labels; numbers are not repeated in the package.
2. **T1086 closed.** `rd_compstable` now drops one reference type (as `rd_typecont` does) instead of pseudo-inverting the structurally singular full covariance. Numerically identical for binary types (the two jumps are exact negatives); exact `df` on every platform.
3. **Cutoff tie in the reflection.** A $t_0$-above unit with $R_{i,t_0}=c$ got `xref = 0` and was assigned to the right ($t_{\mathrm{RD}}$) group by `rd_period`'s `x >= c` split. Now kept on the left. No unit sits at the cutoff in the Grembi data, so Section 6 is unaffected.
4. **Joint over pairs** in `rd_compstable` assumes independent pairs; documented as approximate (the paper's test is per pair).
5. **Paper-side (recorded, not changed):** D5 — the permutation test is defined only in the notes to Table `tab:app-validation` (with the Canay–Kamat citation), not in 4.4. Also stale after the 2026-09-03 renumbering: the same table notes (`main.tex` ≈ line 587) still say "A7 (type-share continuity) … A8 (composition stability) … A9 (type-homogeneous confounding) … A10", now pointing at the wrong assumptions; and the simulation appendix (≈ line 1559) says "the ranking of the three tests in Section 4.4", which 4.4 no longer contains. Scope sentence of 4.4 now reads "balanced panel data with time varying running variables" (Dor, 2026-09-08); the A7 type clause is in.
7. **Adversarial review (2026-09-08) fixes:** (a) `rd_compstable` off-diagonal between type indicators under `"cs"` was zero; now the same-side term (reachable at $P\ge3$ only); (b) type/cell orders are locale-independent (radix), with the reference type documented as the all-below pattern — the exported `rd_homog` contrast sign no longer depends on the machine; (c) drop-one applies only when every type fitted, so a failed fit no longer yields $p=1$ silently; (d) units unobserved in a shared period get no type in the reflection (were pasted as `"NA"` types, inflating df); (e) `rd_typecont`'s scheme detection uses every typed unit, not a rule-of-thumb window.
6. `docs/s34-assumption-tests.md` (untracked, Dropbox-only) was the June spec for these functions; it is superseded by this file (header added).

---

## 4. Drift control

- `tests/testthat/test-s44-conformance.R` + `helper-s44-reference.R`: from-the-text references (reusing the App. B reference sandwiches) vs the four functions at 1e-10, fixed bandwidths, PV and census-like PC panels.
- `dev/check_appB_labels.R` now scans every `dev/*.md`, so the labels cited here are checked against `main.tex` too.
