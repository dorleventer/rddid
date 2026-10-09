# rddid documentation style and vocabulary (UX sweep, 2026-10-08)

Shared brief for everyone writing user-facing text (roxygen, vignettes, README, pkgdown). The reader is an applied economist who knows `rdrobust` and wants an RD-DID estimate plus the validation tests. Numbers are never changed by documentation work; **in `R/*.R` files only lines starting with `#'` may be edited** (plus the removal of stale `# Manuscript ref:` header comments).

## The one-paragraph story (use it, in these words)
A treatment of interest switches on at a cutoff of a running variable in one period, the **RD period**. A **confounding policy** switches at the same cutoff, in every period, so the jump in the outcome at the cutoff in the RD period mixes the treatment effect with the **confounding jump**. In the **comparison periods** the treatment of interest is uniform at the cutoff (nobody treated, or everybody treated), so, provided the treatment of interest has no anticipation or carry-over effects there (which the paper assumes), the jump there *is* the confounding jump. `rddid()` estimates the jump in every period by local-linear RD and subtracts a weighted average of the comparison-period jumps from the RD-period jump. How the weights are set is the **confounding-trend assumption**: constant (equal weights) or linear in time.

Display equation (use once per document, define every symbol in words right after it):
`ATT(t_RD) = D_{t_RD} - sum_t w_t D_t`, where `D_t` is the jump in the outcome at the cutoff in period `t` and `w_t` are the comparison-period weights.

## Vocabulary (one name per thing; never the paper's LaTeX labels)
| Thing | Say | Do not say |
|---|---|---|
| period in which the treatment of interest switches at the cutoff | the RD period (`t_rd`) | treated period, t_RD in prose |
| periods where it is uniform at the cutoff | comparison periods (`comparisons`) | control periods, pre-periods, T_0 |
| the discontinuity in the mean *untreated / treated* outcome at the cutoff | the confounding jump | alpha, bias parameter |
| `trend = "constant"` / `"linear"` | the confounding jump is the same in every period / moves linearly in time | g_0, discontinuity trend |
| estimand when the comparison periods are uniformly **untreated** | the ATT: the effect of the treatment on the units just above the cutoff that are treated in the RD period | |
| estimand when the comparison periods are uniformly **treated** (`estimand = "atu"`) | the ATU: the effect for the units just below the cutoff, which are untreated in the RD period | "mirror" only in a Details paragraph |
| `scheme = "cs"` | repeated cross-section: different units in each period | |
| `scheme = "pc"` | panel, running variable fixed over time: the same units, each on the same side of the cutoff in every period | time-constant |
| `scheme = "pv"` | panel, running variable varies over time: some units change side between periods ("switchers") | time-varying |
| what the scheme does | sets which standard error is reported; under `"joint"`/`"iter"` it also enters the common bandwidth (so the estimate); with a fixed `h` or `"cct"` the estimate is unchanged | "only changes the s.e." |
| `bwselect = "joint"` (default) | the common bandwidth: one `h` for every period, chosen to minimise the asymptotic mean squared error of the RD-DID estimate | AMSE-optimal (only in the Details) |
| `bwselect = "cct"` | per-period CCT bandwidths: each period's own MSE-optimal bandwidth from `rdrobust::rdbwselect` (Calonico, Cattaneo and Titiunik, 2014) | |
| `bwselect = "iter"` | period-specific bandwidths by coordinate descent on the aggregate's asymptotic MSE; not used in the paper, kept for simulations | preferred rule |
| `h`, `b` | the main bandwidth (point estimate) and the pilot bandwidth (bias correction) | |
| `weighting = "min_variance"` | the minimum-variance weights: the comparison-period weights the confounding-trend assumption allows that make the variance of the estimate smallest | GLS weights, MV, optimal weights |
| estimates row `Conventional` | the local-linear estimate with its conventional standard error | |
| estimates row `Robust` | the bias-corrected estimate with its robust standard error (Calonico, Cattaneo and Titiunik, 2014); printed as "Robust (bias-corrected)" | |
| the four tests | type continuity (`rd_typecont`), composition stability (`rd_compstable`), homogeneous confounding (`rd_homog`), constant within-type confounding (`rd_trendcell`) | assumption numbers (A7–A10 — they change) |
| a unit's **type** | the side of the cutoff it is on in the *other* period(s); with two periods, "above in the other period" or "below in the other period" | sign pattern (only in Details) |

Nulls of the four tests, in words (the print methods use exactly these):
- type continuity: the share of each type jumps by zero at the cutoff, in every period;
- composition stability: the share of each type among the units just above the cutoff is the same in the RD period and in each comparison period;
- homogeneous confounding: in each comparison period the confounding jump is the same for every type;
- constant within-type confounding: within each type, the confounding jump is the same in every comparison period (`trend = "linear"`: moves linearly).

When the running variable does not move (`rddid_sim`), the types are degenerate and the tests are not informative; show them on `rddid_sim_pv`.

## Citation (use this everywhere)
Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression Discontinuity Designs Using Multiple Time-Period Data.* arXiv:2408.05847. Refer to parts of the paper by topic ("the paper's treatment of a time-varying running variable"), not by section number.

Also: Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica* 82(6), 2295–2326.

## Rules for examples and vignettes
- Every example uses `rddid_sim` (fixed running variable) or `rddid_sim_pv` (moving; for the tests), is at most 10 lines, runs in under 5 seconds, no `\dontrun`.
- Name every argument in examples (`y = "Y"`, never positional beyond `data`).
- Defaults are the headline: the first `rddid()` call in any document uses defaults except `t_rd`.
- Show output and read it: after a printout, say in one or two sentences what each line means.
- Get started reaches the first estimate within its first 10 lines of code.
- No LaTeX labels (`eq:`, `thm:`, `ass:`, `app:`, `prop:`), no `dev/` paths, no theorem numbers in titles.
- Jargon is defined at first use, in the sentence that uses it.
- Shared `@param` text: `data`, `y`, `x`, `time`, `id`, `t_rd`, `comparisons`, `estimand`, `c`, `h`, `b`, `kernel`, `scheme`, `bwselect`, `p`, `q`, `level` are documented once on `rddid()` and pulled in with `@inheritParams rddid` where the meaning is the same; where a function's option set differs (e.g. `bwselect` in the tests: `"cct"`/`"rot"`), document it locally.
- British/American: American spelling ("minimize" is fine either way; be consistent within a file).

## Site structure (pkgdown)
Navbar: Get started · Articles (How rddid() computes the estimate · Bandwidth rules and sampling schemes · Checking the identification assumptions) · Reference · Changelog · GitHub.
Reference groups: **Estimate** (`rddid`, `summary.rddid`, `rddid-methods`, `rddid-tidiers`) · **Check the assumptions** (`rd_typecont`, `rd_compstable`, `rd_homog`, `rd_trendcell`) · **Building blocks** (`rd_period`, `rd_bw_cct`) · **Data** (`rddid_sim`, `rddid_sim_pv`).
