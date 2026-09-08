# rddid ↔ Appendix B map

**Purpose.** One row per object the package computes that has a definition in the paper
(Leventer & Nevo, RD-DID; `rd-did/writing/main_revamp/main.tex`, Section 5 and Appendix B).
Each row gives the paper's exact expression, the package symbol that implements it, the
scaling/convention difference (if any), the audit status, and the test that pins it.
The paper is the source of truth; the package follows.

**Last synced.** Paper: `rd-did` commit `6eae87d` (2026-09-08; App. B.3–B.4 rewrite `795ca4a`..`e187613`, then the App. B review fixes `b3ae549`, pilot-bandwidth rates `bdd08e8`, §5.3 sync `6eae87d`). Symbols: the paper's kernel constants are $\nu_{(\blacktriangle),p}$ (bias), $\omega_{(\blacktriangle),p}$ (variance), $\omega_{(\blacktriangle),p}(\rho)$ (cross-period, $\omega(1)=\omega$); the package's mnemonics are `.kc_b`, `.kc_v`, `.kc_c`.
Package: this file's commit. Update both stamps whenever either side changes a mapped object.

**How to keep it honest.**
- Every `MATCH` row names a test in `tests/testthat/` whose title starts with the paper label.
  A row without a test is a claim, not a fact.
- `Rscript dev/check_appB_labels.R` checks that every paper label cited here and in `R/*.R`
  still exists in `main.tex` (labels get renamed; this catches silent drift).
- When a formula changes in the paper, change the package function in the same sitting, then
  re-run `tests/testthat/test-appB-conformance.R` and update the row here.

**Status legend.** `MATCH` = same object, verified by test to 1e-10 (up to the stated
convention). `MATCH*` = same object up to a documented finite-sample convention. `MISMATCH` =
package computed something else (fixed in this pass unless marked OPEN). `PAPER-GAP` = the
package needs a choice the paper does not state. `DECISION` = needs Dor's call. `DELEGATED` =
the package calls `rdrobust` for this object. `OUT-OF-SCOPE` = mapped to a section outside
Appendix B, not audited here.

---

## 0. Conventions that differ from the paper's notation (all cancel or are documented)

| Topic | Paper | Package | Consequence |
|---|---|---|---|
| Kernel weight | $A_{t,(\blacktriangle)}(h_t)=\mathrm{diag}(\mathbf 1\{i\in\mathcal N_t\}\mathbf 1\{\text{side}\}K_{h_t}(R_{i,t}-c))$, $K_h(u)=K(u/h)/h$ | `.rd_kweight((x-c)/h, kernel)`: $K(u)$ **not** divided by $h$; only side-active units enter | $1/h$ and $1/n$ cancel in $\widehat\beta$, in $Q$-form BC, and in every sandwich; verified to 1e-10 (`eq:wls_rd` test) |
| Kernel shapes | generic $K$ (Assumption `ass:reg`(b)) | triangular $(1-\lvert u\rvert)_+$, Epanechnikov $0.75(1-u^2)_+$, uniform $0.5\cdot\mathbf 1\{\lvert u\rvert\le1\}$ | matches `rdrobust`; the kernel constants of §2.4 are computed for the same shapes |
| Sample | $i=1,\dots,n$ units, $\mathcal N_t$, $n_t$, $\mathcal N_{t,s}$ | one `rd_period` fit per period on that period's rows; cross-period objects matched on `id` | $n_t$ = `fit$n`; $\mathcal N_{t,s}$ = shared `id`s; CS = distinct ids |
| Active set within a period | all $i\in\mathcal N_t$ (zero weight outside the window) | units with $K((R-c)/b_t)>0$ **or** $K((R-c)/h_t)>0$ (union of pilot and main windows; before this pass the pilot window alone, which silently truncated the main fit when $b_t<h_t$) | identical when $b_t\ge h_t$ (the case for every selector in the package) |
| Residual for $\widehat{\boldsymbol\Sigma}$ | $\widehat\varepsilon_{i,t}$ = own-side $p$-fit at $h_t$ (preamble) for the conventional variance; own-side $q$-fit at $b_t$ for the BC variance and the curvature-variance (B.2 "Estimation"; sentence to be added per D1/D4, settled 2026-09-08) | conventional: own-side $p$-fit at $h_t$; BC and `g_diff`: own-side $q$-fit at $b_t$ (`rdrobust` convention) | D1/D4 settled (Dor, 2026-09-08) |
| Finite-sample factor | none yet (HC1 clause to be added to B.2 "Estimation" per D2, settled 2026-09-08) | HC1: residual $\times\sqrt{n_{\blacktriangle}/(n_{\blacktriangle}-(p+1))}$ (conv) or $\sqrt{n_{\blacktriangle}/(n_{\blacktriangle}-(q+1))}$ (BC), $n_{\blacktriangle}$ = side-$\blacktriangle$ units in the pilot window (`rdrobust` `vce="hc1"`, matched to 1e-10 in `test-rd_period.R`) | variances $\approx+0.5$–$0.7\%$ (Study 5c of the App. B simulation verification); D2 settled (Dor, 2026-09-08) |
| Sides | $\blacktriangle\in\{+,-\}$, $(+)$: $R\ge c$ | `sides[["+"]]`, `sides[["-"]]`; `x >= c` is $(+)$ | same |
| Signed weights | $\tilde w_{t_{\mathrm{RD}}}=1$, $\tilde w_{t_0}=-w_{t_0}$ | `coef` (named vector; RD period $+1$, comparisons $-w$) | same object |

---

## 1. Label registry

Paper labels cited in this file or in `R/*.R` comments (checked by `dev/check_appB_labels.R`):
`sec:estimation`, `sec:est-point`, `eq:wls_rd`, `eq:est_agg_bc`, `sec:est-inf`, `eq:agg-var`,
`eq:var-cs`, `eq:var-pc`, `eq:var-pv`, `sec:est-bw`, `eq:joint_amse`, `alg:coorddesc`,
`sec:est-adjust`, `app:est`, `ass:reg`, `app:est-bias`, `eq:bias_single`, `eq:bc_intercept`,
`eq:bc_Q`, `eq:bias-agg`, `app:est-var`, `eq:cross-decomp`, `app:est-asy`,
`eq:per-period-orders`, `eq:cov-rescaled`, `lem:aux-gamma`, `lem:aux-psi-pc`, `lem:aux-psi-pv`,
`lem:cov-cs`, `lem:cov-pc`, `lem:cov-pv`, `lem:agg-var`, `app:est-bw`, `eq:amse-att`,
`eq:common_h_opt`, `eq:update`, `app:est-adjust`.

Labels that were removed from the paper on 2026-09-08 and must not be cited: `lem:coercive`,
`eq:amse-ps`, "Appendix C" (the old bandwidth appendix; its content is now `app:est-bw`).

---

## 2. Object-by-object map

### 2.1 Preamble of Appendix B (`app:est`) and the estimator (`sec:est-point`)

| Paper object (label) | Exact expression in the paper | Package | Convention | Status | Test |
|---|---|---|---|---|---|
| $\boldsymbol X_p(r)$ | $[1,(r-c),\dots,(r-c)^p]'$ | `Rq <- outer(xs - c, 0:q, "^")`, `Rp <- Rq[, 1:(p+1)]` (`rd_period`) | — | MATCH | `eq:wls_rd` |
| $\boldsymbol\Gamma_{t,(\blacktriangle),p}(h_t)$ | $\tfrac1n\boldsymbol X_{t,p}'\boldsymbol A_{t,(\blacktriangle)}(h_t)\boldsymbol X_{t,p}$ | `invG_p <- .qrXXinv(sqrt(w_h) * Rp)` $=(\boldsymbol X'\boldsymbol W_h\boldsymbol X)^{-1}$ | $\boldsymbol\Gamma^{-1}=nh_t\cdot$`invG_p` | MATCH | `eq:wls_rd` |
| $\widehat\beta_{t,(\blacktriangle),p}(h_t)$ (`eq:wls_rd`) | $\tfrac1n\boldsymbol\Gamma^{-1}\boldsymbol X_{t,p}'\boldsymbol A_{t,(\blacktriangle)}(h_t)\boldsymbol Y_t$ | `beta_p <- invG_p %*% crossprod(Rp * w_h, ys)` | scale-free | MATCH | `eq:wls_rd` |
| $\widehat\beta^{(v)}_{t,(\blacktriangle),p}$ | $v!\,e_{v,p}'\widehat\beta$ | intercept `beta_p[1]`; slope `beta_p[2]` (only $v=0,1$ exposed; $v=p+1$ from the $q$-fit enters $Q$ implicitly) | — | MATCH | `eq:bc_intercept == eq:bc_Q` |
| $\widehat D_t(h_t)$ | $\widehat\beta^{(0)}_{t,(+),p}(h_t)-\widehat\beta^{(0)}_{t,(-),p}(h_t)$ | `D <- R_side$beta0 - L_side$beta0` | — | MATCH | `eq:wls_rd` |
| $\widehat{\att}(t_{\mathrm{RD}}\mid\{h_t\})$ | $\widehat D_{t_{\mathrm{RD}}}-\sum_{t_0\in\mathcal T_0}w_{t_0}\widehat D_{t_0}$ | `.aggregate_fits()$est` $=\sum_\tau$`coef`$_\tau\widehat D_\tau$ | — | MATCH | `eq:var-cs, eq:var-pc, eq:var-pv` |
| $\widehat\varepsilon_{i,t}$ | $Y_{i,t}-\boldsymbol X_p(R_{i,t})'\widehat\beta_{t,(\pm),p}(h_t)$, own side | `res_c <- sqrt(n_s/(n_s-(p+1))) * (ys - Rp %*% beta_p)` | HC1 factor (§0) | MATCH* | `B.2 single-period sandwich` |
| $n_t$ | $\lvert\mathcal N_t\rvert$ | `fit$n` (complete cases of the period's rows) | — | MATCH | — |

### 2.2 B.1 First-order bias (`app:est-bias`)

| Paper object (label) | Exact expression in the paper | Package | Convention | Status | Test |
|---|---|---|---|---|---|
| $\boldsymbol Z_t(h_t)$, $\boldsymbol\vartheta_{t,(\blacktriangle),p}(h_t)$ | $\boldsymbol Z_t=(z_{i,t}^{p+1})_i$, $z_{i,t}=(R_{i,t}-c)/h_t$; $\boldsymbol\vartheta=\tfrac1n\boldsymbol X_{t,p}'\boldsymbol A_{t,(\blacktriangle)}(h_t)\boldsymbol Z_t(h_t)$ | `theta <- crossprod(Rp * w_h, ((xs - c)/h)^(p+1))` | $\boldsymbol\vartheta=$`theta`$/(nh_t)$ | MATCH | `eq:bc_intercept == eq:bc_Q` |
| $B_{t,(\blacktriangle),p}(h_t)$ | $e_{0,p}'\boldsymbol\Gamma_{t,(\blacktriangle),p}(h_t)^{-1}\boldsymbol\vartheta_{t,(\blacktriangle),p}(h_t)$ | not formed explicitly; enters through $Q$ | — | MATCH (implicit) | `B.1 bias constant` |
| BC intercept (`eq:bc_intercept`) | $\widehat\beta^{(0)}_{t,(\blacktriangle),p}(h_t)-\tfrac{h_t^{p+1}}{(p+1)!}\widehat\beta^{(p+1)}_{t,(\blacktriangle),q}(b_t)B_{t,(\blacktriangle),p}(h_t)$ | computed in the equivalent $Q$-form below | $(p+1)!$ cancels against $\widehat\beta^{(p+1)}=(p+1)!e_{p+1,q}'\widehat\beta_q$ | MATCH | `eq:bc_intercept == eq:bc_Q` |
| $\boldsymbol Q_{t,(\blacktriangle),p,q}(h_t,b_t)$ (`eq:bc_Q`) | $\boldsymbol X_{t,p}'\boldsymbol A_{t,(\blacktriangle)}(h_t)-h_t^{p+1}\boldsymbol\vartheta_{t,(\blacktriangle),p}(h_t)e_{p+1,q}'\boldsymbol\Gamma_{t,(\blacktriangle),q}(b_t)^{-1}\boldsymbol X_{t,q}'\boldsymbol A_{t,(\blacktriangle)}(b_t)$ | `Qmat <- t(Rp * w_h) - h^(p+1) * (theta %*% (t(e_p1) %*% invG_q %*% Aq_b))`, `e_p1` $=e_{p+1,q}$ (index `p+2`) | `Qmat` $=h_t\boldsymbol Q$ | MATCH | `eq:bc_intercept == eq:bc_Q` |
| $\widehat\beta^{(0),\mathrm{BC}}$ via `eq:bc_Q` | $e_{0,p}'\boldsymbol\Gamma_{t,(\blacktriangle),p}(h_t)^{-1}\tfrac1n\boldsymbol Q\boldsymbol Y_t$ | `beta_bc <- invG_p %*% (Qmat %*% ys)`; `beta0_bc <- beta_bc[1]` | scale-free | MATCH (Study 5c: 1e-14) | `eq:bc_intercept == eq:bc_Q` |
| $\widehat D^{\mathrm{BC}}_t(h_t,b_t)$ | $\widehat\beta^{(0),\mathrm{BC}}_{t,(+)}-\widehat\beta^{(0),\mathrm{BC}}_{t,(-)}$ | `D_bc` | — | MATCH | same |
| $\widehat{\att}^{\mathrm{BC}}$ (`eq:est_agg_bc`) | $\widehat D^{\mathrm{BC}}_{t_{\mathrm{RD}}}-\sum_{t_0}w_{t_0}\widehat D^{\mathrm{BC}}_{t_0}$ | `.aggregate_fits(fits, coef, bc = TRUE)$est` | — | MATCH | `eq:var-cs, eq:var-pc, eq:var-pv` (bc branch) |
| $\widehat{\mathtt B}_t(h_t,b_t)$ | $\tfrac{h_t^{p+1}}{(p+1)!}\bigl[\widehat\beta^{(p+1)}_{t,(+),q}(b_t)B_{t,(+),p}(h_t)-\widehat\beta^{(p+1)}_{t,(-),q}(b_t)B_{t,(-),p}(h_t)\bigr]$ | $=\;$`D - D_bc` (identity, both linear in $\boldsymbol Y_t$) | — | MATCH | `B.1 bias constant` |
| $\mathtt B_t(h_t)$ leading form (`eq:bias_single`, `eq:per-period-orders`) | $\tfrac{h_t^{p+1}}{(p+1)!}\mathcal B_t+o_p(h_t^{p+1})$, $\mathcal B_t=m^{(p+1)}_{t,(+)}\nu_{(+),p}-m^{(p+1)}_{t,(-)}\nu_{(-),p}$ | `b_const <- factorial(p+1) * (D - D_bc) / h^(p+1)` $=\widehat{\mathcal B}_t$ with $B_{t,(\blacktriangle),p}(h_t)$ in place of $\nu_{(\blacktriangle),p}$ | finite-sample $B(h)\to_p\nu$ (`lem:aux-gamma`); **was hard-coded to $p=1$ (`2*(D-D_bc)/h^2`)** | MISMATCH → fixed (general $p$) | `B.4 constants general p` |

### 2.3 B.2 Variance (`app:est-var`) and Section 5.2 (`sec:est-inf`)

| Paper object (label) | Exact expression in the paper | Package | Convention | Status | Test |
|---|---|---|---|---|---|
| Sandwich $V(\widehat\beta_{t,(\blacktriangle),p}(h_t)\mid\boldsymbol R_t)$ | $\tfrac1n\boldsymbol\Gamma^{-1}\boldsymbol\Psi_{t,t,(\blacktriangle,\blacktriangle),p}\boldsymbol\Gamma^{-1}$, $\boldsymbol\Psi_{t,t}=\tfrac1n\boldsymbol X'\boldsymbol A\boldsymbol\Sigma_t\boldsymbol A\boldsymbol X$ | intercept row only: `a_c <- invG_p[1,] %*% t(Rp * w_h)`; `g <- a_c * res_c`; $V(\widehat\beta^{(0)})=$`sum(g^2)` | $=e_{0,p}'\boldsymbol\Gamma^{-1}\boldsymbol\Psi\boldsymbol\Gamma^{-1}e_{0,p}$ with $\widehat{\boldsymbol\Sigma}_t=\mathrm{diag}(\widehat\varepsilon^2)$, HC1 | MATCH* | `B.2 single-period sandwich` |
| $V(\widehat D_t\mid\boldsymbol R_t)$ | $V(\widehat\beta^{(0)}_{t,(+)})+V(\widehat\beta^{(0)}_{t,(-)})$ | `V_D <- sum(R$g^2) + sum(L$g^2)` | — | MATCH* | same |
| $\boldsymbol\Psi^{\mathrm{BC}}$, BC sandwich | $\tfrac1n\boldsymbol Q_t\boldsymbol\Sigma_{t,s}\boldsymbol Q_s'$ in place of $\boldsymbol\Psi$ | `a_bc <- invG_p[1,] %*% Qmat`; `g_bc <- a_bc * res_b`; `V_D_bc <- sum(g_bc^2)` over sides | residual: $q$-fit at $b_t$ (D1); HC1 with $q+1$ | MATCH* | `B.2 single-period sandwich` (bc) |
| $\widehat{\boldsymbol\Sigma}_t$, $\widehat{\boldsymbol\Sigma}_{t,s}$ ("Estimation") | $\widehat\varepsilon_{i,t}^2$; $\widehat\varepsilon_{i,t}\widehat\varepsilon_{i,s}$ on $\mathcal N_{t,s}$ | implicit in `g` products: `.match_sum(idA, gA, idB, gB)` $=\sum_{i\in\mathcal N_{t,s}}g_{t,i}g_{s,i}$ | — | MATCH* | `eq:cross-decomp` |
| Cross-period sandwich, $\mathrm{Cov}(\widehat\beta^{(0)}_{t,(\blacktriangle)},\widehat\beta^{(0)}_{s,(\blacktriangledown)})$ | $e_{0,p}'\tfrac1n\boldsymbol\Gamma_t^{-1}\boldsymbol\Psi_{t,s,(\blacktriangle,\blacktriangledown),p}(h_t,h_s)\boldsymbol\Gamma_s^{-1}e_{0,p}$ | `.match_sum` of the two sides' `g` vectors (`aggregate.R`) | — | MATCH* | `eq:cross-decomp` |
| $C^{\mathrm{same}}_{t,s}$ | $\mathrm{Cov}(\widehat\beta^{(0)}_{t,(+)},\widehat\beta^{(0)}_{s,(+)})+\mathrm{Cov}(\widehat\beta^{(0)}_{t,(-)},\widehat\beta^{(0)}_{s,(-)})$ | `.cross_cov(ft, fs)$pc` (per-side parts `pc_p`, `pc_m` added this pass for §2.5) | — | MATCH* | `eq:cross-decomp` |
| $C^{\mathrm{opp}}_{t,s}$ | $\mathrm{Cov}(\widehat\beta^{(0)}_{t,(+)},\widehat\beta^{(0)}_{s,(-)})+\mathrm{Cov}(\widehat\beta^{(0)}_{t,(-)},\widehat\beta^{(0)}_{s,(+)})$ | `.cross_cov(ft, fs)$pv` | — | MATCH* | `eq:cross-decomp` |
| $\mathrm{Cov}(\widehat D_t,\widehat D_s)$ (`eq:cross-decomp`) | $C^{\mathrm{same}}_{t,s}-C^{\mathrm{opp}}_{t,s}$ | `.cov_scheme(ft, fs, "pv")` $=$`pc - pv` | — | MATCH* | `eq:cross-decomp` |
| $V^{\mathrm{CS}}$ (`eq:var-cs`) | $V(\widehat D_{t_{\mathrm{RD}}})+\sum_{t_0}w_{t_0}^2V(\widehat D_{t_0})$ | `.aggregate_fits()$V_cs` $=\sum_\tau$`coef`$_\tau^2V_D^\tau$ | — | MATCH* | `eq:var-cs, eq:var-pc, eq:var-pv` |
| $V^{\mathrm{PC}}$ (`eq:var-pc`) | $V^{\mathrm{CS}}-2\sum_{t_0}w_{t_0}C^{\mathrm{same}}_{t_{\mathrm{RD}},t_0}+\sum_{t_0}\sum_{s_0\ne t_0}w_{t_0}w_{s_0}C^{\mathrm{same}}_{t_0,s_0}$ | `V_pc = V_cs + Σ_{i<j} 2 coef_i coef_j pc_ij` | ordered double sum = $2\times$ unordered | MATCH* | same |
| $V^{\mathrm{PV}}$ (`eq:var-pv`) | $V^{\mathrm{PC}}+2\sum_{t_0}w_{t_0}C^{\mathrm{opp}}_{t_{\mathrm{RD}},t_0}-\sum_{t_0}\sum_{s_0\ne t_0}w_{t_0}w_{s_0}C^{\mathrm{opp}}_{t_0,s_0}$ | `V_pv = V_cs + Σ_{i<j} 2 coef_i coef_j (pc_ij - pv_ij)` | — | MATCH* | same |
| Sampling scheme ("Sampling scheme" paragraph) | CS: $\boldsymbol\Sigma_{t,s}=0$; PC: opposite-side covariances zero; PV: all four | `.detect_scheme()`: no repeated ids → `cs`; repeated, never switching side → `pc`; else `pv`. All three variances always returned; `scheme` picks the headline SE | the paper defines the scheme by design, the package infers it from the data | MATCH (structural) | `PC data: C^opp = 0 ...` |
| Inference (`sec:est-inf`) | $\widehat{\att}\pm z_{1-\alpha/2}\widehat V^{1/2}$ | `rddid()$estimates`: `se`, `ci_l`, `ci_u` at `scheme`; rows Conventional / Robust (BC point + BC variance) | — | MATCH | `test-rddid.R` |

### 2.4 B.3 Asymptotics (`app:est-asy`)

| Paper object (label) | Exact expression in the paper | Package | Convention | Status | Test |
|---|---|---|---|---|---|
| $\widetilde{\boldsymbol\Gamma}_{(\blacktriangle),p}$, $\widetilde{\boldsymbol\Psi}_{(\blacktriangle),p}$, $\widetilde{\boldsymbol\vartheta}_{(\blacktriangle),p}$ | $\int K(u)\boldsymbol r_p\boldsymbol r_p'du$, $\int K(u)^2\boldsymbol r_p\boldsymbol r_p'du$, $\int K(u)u^{p+1}\boldsymbol r_pdu$ over $[0,\infty)$ / $(-\infty,0]$ | `.kc_gamma(p, side, kernel)`, `.kc_psi`, `.kc_theta` (`R/kernel_constants.R`, composite Simpson, memoized) | support $[-1,1]$ so the integrals run over $[0,1]$ / $[-1,0]$ | MATCH (ported from the App. B simulation verification, Study 3) | `test-kernel-constants.R` |
| $\nu_{(\blacktriangle),p}$, $\omega_{(\blacktriangle),p}$ | $e_{0,p}'\widetilde{\boldsymbol\Gamma}^{-1}\widetilde{\boldsymbol\vartheta}$; $e_{0,p}'\widetilde{\boldsymbol\Gamma}^{-1}\widetilde{\boldsymbol\Psi}\widetilde{\boldsymbol\Gamma}^{-1}e_{0,p}$ | `.kc_b`, `.kc_v` | triangular $p=1$: $\nu_+=-0.1$, $\omega_+=4.8$ | MATCH | same |
| $\widetilde{\boldsymbol\Omega}_{(\blacktriangle),p}(\rho)$, $\omega_{(\blacktriangle),p}(\rho)$ | $\int K(v)K(\rho v)\boldsymbol r_p(v)\boldsymbol r_p(\rho v)'dv$; $e_{0,p}'\widetilde{\boldsymbol\Gamma}^{-1}\widetilde{\boldsymbol\Omega}(\rho)\widetilde{\boldsymbol\Gamma}^{-1}e_{0,p}$; $\omega(1)=\omega$, $\omega(1/\rho)=\rho\,\omega(\rho)$ | `.kc_omega(p, side, rho, kernel)`, `.kc_c` | kink split at $\lvert v\rvert=1/\rho$ | MATCH | same |
| $\mathcal B_t$, $\mathcal V_t$ | $m^{(p+1)}_{t,(+)}\nu_{(+),p}-m^{(p+1)}_{t,(-)}\nu_{(-),p}$; $\bigl(\sigma^2_{t,(+)}\omega_{(+),p}+\sigma^2_{t,(-)}\omega_{(-),p}\bigr)/f_t(c)$ | plug-ins: `b_const` (§2.2) and `v_const <- n * h * V_D` (`eq:per-period-orders`: $V(\widehat D_t)=\mathcal V_t/(n_th_t)$) | plug-in through the variance, not through $\sigma^2$, $f_t(c)$ separately | MATCH (plug-in) | `B.4 constants general p` |
| $\mathcal C_{t,s}$ (`lem:cov-pc`) | $\bigl(\sigma_{t,s,(+,+)}\omega_{(+),p}(\rho_{t,s})+\sigma_{t,s,(-,-)}\omega_{(-),p}(\rho_{t,s})\bigr)/f(c)$, $\rho_{t,s}=h_t/h_s$, with $nh_s\mathrm{Cov}(\widehat D_t,\widehat D_s)\to_p\mathcal C_{t,s}$ | per-side plug-in in `.bw_joint_iter`: $\widehat\kappa^{\blacktriangle}_{t,s}=$ `pc_side * h0_s / .kc_c(p, side, h0_t/h0_s)` $\approx\sigma_{t,s,(\blacktriangle,\blacktriangle)}/(f(c)\,n)$ from the pilot fits | the $h$-dependence $\omega(h_t/h_s)/h_s$ is applied exactly; **before this pass the package used $\kappa/\max(h_t,h_s)$** (exact only for the uniform kernel at $p=0$) | MISMATCH → fixed | `B.4 PC cross term` |
| `lem:cov-pv` | cross-period covariance $O_p(1/n)=o_p(1/(nh))$ under PV | dropped from bandwidth objectives (`Vfld pv → V_cs`; no `cov_term`); **kept** in the reported $V^{\mathrm{PV}}$ (finite-sample) | consistent with `lem:agg-var` ($\mathbf 1\{\mathrm S=\mathrm{PC}\}$) | MATCH | `B.4 scheme switch` |
| `lem:agg-var` | $V^{\mathrm S}=\sum_\tau\tilde w_\tau^2\mathcal V_\tau/(\pi_\tau nh_\tau)+\mathbf 1\{\mathrm S=\mathrm{PC}\}\sum_\tau\sum_{s\ne\tau}\tilde w_\tau\tilde w_s\mathcal C_{\tau,s}/(nh_s)+o_p(r_n)$ | `amse()` in `.bw_joint_iter`: `var_term = Σ coef² v_const/(n_t h_t)` ($\mathcal V_\tau/(\pi_\tau n h_\tau)=\mathcal V_\tau/(n_\tau h_\tau)$) `+ cov_term` (PC only) | uses $n_\tau$ directly instead of $\pi_\tau n$ | MATCH | `B.4 objective` |

### 2.5 B.4 MSE-optimal bandwidth (`app:est-bw`) and Section 5.3 (`sec:est-bw`)

| Paper object (label) | Exact expression in the paper | Package | Convention | Status | Test |
|---|---|---|---|---|---|
| $h_t^{\mathrm{CCT}}$ (P1) | $\bigl(\tfrac{((p+1)!)^2}{2(p+1)}\tfrac{\mathcal V_t}{\mathcal B_t^2}\bigr)^{1/(2p+3)}n_t^{-1/(2p+3)}$ | `.bw_cct()` / `rd_bw_cct()` → `rdrobust::rdbwselect(bwselect = "mserd")` (CCT's feasible version incl. their regularization; also returns $b_t$) | delegated | DELEGATED | `test-bw_cct.R` |
| $\bar{\mathtt B}(\{h_t\})$ (P2) | $\tfrac{1}{(p+1)!}\sum_\tau\tilde w_\tau h_\tau^{p+1}\mathcal B_\tau$ | `Bbar <- sum(cf * hv^(p+1) * bt) / factorial(p+1)`; **was `0.5 * sum(cf * hv^2 * bt)`** | plug-in $\widehat{\mathcal B}_\tau$ = `b_const` | MISMATCH → fixed | `B.4 objective` |
| $\mathrm{AMSE}^{\mathrm S}(\{h_t\})$ (`eq:amse-att`) | $\bar{\mathtt B}^2+\bar V^{\mathrm S}$, $\bar V^{\mathrm S}$ = `lem:agg-var` without remainder | `amse(hv) = Bbar^2 + pen + var_term + cov_term` | `pen` = regularization (below) | MATCH (after fix) | `B.4 objective` |
| $\mathcal B(t_{\mathrm{RD}})$, $\mathcal V^{\mathrm S}(t_{\mathrm{RD}})$ (P3) | $\sum_\tau\tilde w_\tau\mathcal B_\tau$; $\sum_\tau\tilde w_\tau^2\mathcal V_\tau/\pi_\tau+\mathbf 1\{\mathrm S=\mathrm{PC}\}\sum_\tau\sum_{s\ne\tau}\tilde w_\tau\tilde w_s\mathcal C_{\tau,s}$ | `.bw_joint()`: `B <- factorial(p+1) * (att_conv - att_bc) / h0^(p+1)` $=\sum_\tau$`coef`$_\tau\widehat{\mathcal B}_\tau$ at the common pilot $h_0$; `Veff <- h0 * V^S(h0)` $\approx\mathcal V^{\mathrm S}/n$ with `V^S` = `V_cs` (CS, PV) or `V_pc` (PC) at the common pilot ($\rho=1$, so $\mathcal C_{\tau,s}$ enters through the plug-in covariance) | **`B` was hard-coded to $p=1$** | MISMATCH → fixed | `B.4 common h` |
| $h^{\star,\mathrm S}$ (`eq:common_h_opt`) | $\bigl(\tfrac{((p+1)!)^2}{2(p+1)}\tfrac{\mathcal V^{\mathrm S}(t_{\mathrm{RD}})}{\mathcal B(t_{\mathrm{RD}})^2}\bigr)^{1/(2p+3)}n^{-1/(2p+3)}$ | `h_star <- (factorial(p+1)^2 / (2*(p+1)) * Veff / (B^2 + reg))^(1/(2*p+3))`; **was `(Veff/denom)^(1/5)`** | $n$ is inside `Veff`; `reg` = regularization below | MISMATCH → fixed (identical at $p=1$) | `B.4 common h` |
| Pilot bandwidth $b_t$ for the aggregate | not stated | `.bw_joint`: $b^\star=h^\star\cdot b_0/h_0$ (pilot ratio of the RD period's CCT pair); `.bw_joint_iter`: $b_t=h_t\cdot b^{\mathrm{CCT}}_t/h^{\mathrm{CCT}}_t$ per period | — | PAPER-GAP (G1) | — |
| Update step (`eq:update`, `alg:coorddesc`, P4) | $h_\tau\leftarrow\arg\min_{h_\tau}\mathrm{AMSE}^{\mathrm S}(h_\tau\mid\{h_t\}_{t\ne\tau})$, cycled until no bandwidth moves more than a tolerance or 50 sweeps | `.bw_joint_iter()`: `optimize(obj, c(lo, hmax))` per period, sweeps until `max|Δh| < tol*hmax` or `maxit = 50`; `lo = 0.03*hmax`, `hmax` = running-variable radius; a boundary solution now raises a warning | the paper's argmin is unconstrained; the search box is the package's (with `regularize = TRUE` no boundary hit was found in any design tried) | MATCH* | `B.4 objective` |
| Start of the descent (P4, `sec:est-bw`) | common $h^{\star,\mathrm S}$ or per-period $h^{\mathrm{CCT}}_t$ | `start = "hstar"` (default) / `"cct"` / numeric | — | MATCH | `test-rddid.R` |
| PC cross term in the update (P4) | $\mathbf 1\{\mathrm S=\mathrm{PC}\}\,2\sum_{s\ne\tau}\tilde w_\tau\tilde w_s\mathcal C_{\tau,s}/(nh_s)$ with $\mathcal C_{\tau,s}=\mathcal C_{\tau,s}(\rho_{\tau,s})$, $\rho_{\tau,s}=h_\tau/h_s$ | `cov_term = 2 Σ_{i<j} cf_i cf_j Σ_side κ̂_side_ij · .kc_c(p, side, hv_i/hv_j) / hv_j` (see §2.4) | symmetric in $(i,j)$ by $\omega(1/\rho)=\rho\omega(\rho)$; **was `P_ij * max(h0_i,h0_j) / max(hv_i,hv_j)`** | MISMATCH → fixed | `B.4 PC cross term` |
| Regularization | add $\lambda\sum_\tau\tilde w_\tau^2\bigl(h_\tau^{p+1}/(p+1)!\bigr)^2\operatorname{Var}(\widehat{\mathcal B}_\tau)$ to `eq:amse-att`, default $\lambda=3$; $\operatorname{Var}(\widehat{\mathcal B}_\tau)$ from the per-unit influence of $\widehat{\mathcal B}_\tau$ | `pen <- reg_const * sum(cf^2 * (hv^(p+1)/factorial(p+1))^2 * var_b)`, `var_b <- (factorial(p+1)/h0^(p+1))^2 * Σ g_diff^2`, `g_diff = (a_c - a_bc) * res_c` (influence of $\widehat D-\widehat D^{\mathrm{BC}}$); common-$h$ case enters as `denom <- B^2 + reg` | `reg_const` = $\lambda$; `g_diff` uses the $q$-fit-at-$b$ residuals (D4); **exponents were $p=1$** | MISMATCH → fixed | `B.4 common h` (reg term), `B.4 objective` (pen) |
| $\widehat{\mathcal B}_\tau$ (Regularization paragraph) | $\widehat\beta^{(p+1)}_{\tau,(+),q}(b_\tau)\nu_{(+),p}-\widehat\beta^{(p+1)}_{\tau,(-),q}(b_\tau)\nu_{(-),p}$ | `b_const` $=(p+1)!(\widehat D-\widehat D^{\mathrm{BC}})/h^{p+1}=\widehat\beta^{(p+1)}_{+}B_{+}(h)-\widehat\beta^{(p+1)}_{-}B_{-}(h)$ | $B_{(\blacktriangle),p}(h)$ (finite-sample) in place of $\nu_{(\blacktriangle),p}$; $B(h)\to_p\nu$ | MATCH* | `B.4 constants general p` |
| Body `eq:joint_amse` (`sec:est-bw`) | $\mathrm{AMSE}^{\mathrm S}(\{h_t\})=\tfrac14(\sum\tilde w_th_t^2\mathcal B_t)^2+\sum\tfrac{\tilde w_t^2\mathcal V_t}{n_th_t}+\mathbf 1\{\mathrm S=\mathrm{PC}\}\sum\sum_{s\ne t}\tfrac{\tilde w_t\tilde w_s\mathcal C_{t,s}}{nh_s}$ ($p=1$ form of `eq:amse-att`) | `amse()` in `.bw_joint_iter` at $p=1$ | synced to B.4 on 2026-09-08 (`6eae87d`); body is $p=1$ by design | MATCH | `B.4 objective` |

### 2.6 B.5 Composition-adjusted estimators (`app:est-adjust`, `sec:est-adjust`)

Function-level map only. B.5 has not changed in substance since the 2026-06-30 code↔math audit
(SOUND); only notation (`bc`→`BC`, $C^{\mathrm{PC}}/C^{\mathrm{PV}}$→$C^{\mathrm{same}}/C^{\mathrm{opp}}$).

| Paper object | Package | Status |
|---|---|---|
| $\widehat S_t=\sum_{\mathbf v_{-t}}\widehat{\Delta\pi}_t(\mathbf v_{-t})\widehat\mu_{t,(0,0),(-)}(\mathbf v_{-t})$; $P=2$ form $\widehat{\Delta\pi}_t(\mathrm{above})(\widehat\mu^{(-)}_{\mathrm{above}}-\widehat\mu^{(-)}_{\mathrm{below}})$ | `rd_sadjust()` (`R/sadjust.R`), blocks via `rd_period` | audited 2026-06-30 |
| $\widehat C_{t_0,t_{\mathrm{RD}}}=\sum_a(\widehat\pi_{t_{\mathrm{RD}},(+)}(v_{t_0}{=}a)-\widehat\pi_{t_0,(+)}(v_{t_{\mathrm{RD}}}{=}a))\widehat D_{t_0}(v_{t_{\mathrm{RD}}}{=}a)$ | `rd_c()` (`R/cterm.R`) | audited 2026-06-30 |
| family $\widehat{\att}_s,\widehat{\att}_c,\widehat{\att}_{sc}$; BC by blocks; one shared unit cluster bootstrap; per-block CCT bandwidths fixed across replications | `rd_att()` (`R/att.R`) | audited 2026-06-30 |
| within-type (reweighting) form | `rd_adjust()` (`R/adjust.R`) — legacy; the paper headlines the correction form | not in the paper |

### 2.7 Outside Appendix B (pointers only)

| Package | Paper | Note |
|---|---|---|
| `.rddid_weights()` (`constant` = equal $1/m$; `linear` = OLS line through the comparison periods extrapolated to $t_{\mathrm{RD}}$) | Corollary `cor:trend`, Section 3 | admissibility ($\sum w=1$ for constant, line weights for linear) |
| `rd_typecont`, `rd_compstable`, `rd_homog`, `rd_trendcell` | Section 4.4 (`sec:tvrv-assess`); simulations `app:sim-val` | OUT-OF-SCOPE; separate map (`docs/s34-assumption-tests.md`, untracked) |
| `.detect_scheme()` | Section 5.2 prose | data-driven; the paper defines schemes by design |

---

## 3. Audit findings of this pass (2026-09-08)

**Fixed in the package (numerically identical at $p=1$ except the PC cross term in `bwselect = "iter"`; adversarial verify-only review 2026-09-08: APPROVED, every load-bearing formula confirmed, several non-circularly against Monte Carlo and closed-form theory):**
1. `rd_period$b_const`, `.bw_joint` (`B`, `reg`, exponent and constant of $h^\star$), `.bw_joint_iter`
   (`Bbar`, `pen`, `var_b`) generalized from hard-coded $p=1$ to the paper's general-$p$ forms
   (T1492b). Regression: the frozen $p=1$ snapshot reproduces to machine precision.
2. PC cross term of the period-specific objective now scales as $\omega_{(\blacktriangle),p}(h_\tau/h_s)/h_s$
   per side (`lem:cov-pc`, B.4 P4) instead of $1/\max(h_\tau,h_s)$ (T1492a). This changes
   `bwselect = "iter"` under `scheme = "pc"` only (by construction the two coincide at $\rho=1$).
3. Stale references in `R/bandwidth.R` (Appendix C, `lem:coercive`, `eq:amse-ps`, "Section 4.3",
   the dropped simultaneous-degeneracy discussion) replaced by the current labels (T1492c).
4. `rd_period` active set = union of pilot and main windows (was pilot only; wrong when $b_t<h_t$).
5. Guards (behavior-preserving on the default path): the selectors call the guarded `rd_bw_cct()` instead of
   `rdrobust::rdbwselect` directly; `.bw_joint` warns when $h^\star$ exceeds the data radius; the descent warns on a
   boundary solution; `scheme = "pc"/"pv"` with no repeated ids warns (the cross-period terms are then identically 0).

**Decisions for the paper (not changed here):**
- **D1 (T1474) — SETTLED (Dor, 2026-09-08): paper adopts the `rdrobust` convention.** B.2 "Estimation" is to
  state that the BC variance uses the residuals of the order-$q$ pilot fit at $b_t$ (sentence handed to Dor; the tex
  was being edited by another session at the time). Both variants are consistent
  (Study 5 of the App. B verification: `BC/p_h` and `BC/q_b` ratios → 1; gap 0.1–0.9% in the conformance tests).
- **D2 — SETTLED (Dor, 2026-09-08): HC1 clause to be added to B.2 "Estimation" (handed to Dor).** The package applies `rdrobust`'s
  HC1 factor with the pilot-window side count.
- **D3 (T1475) — RESOLVED on the paper side (rd-did `6eae87d`, 2026-09-08).** Body `eq:joint_amse` now carries
  the scheme superscript, the $\mathbf 1\{\mathrm S=\mathrm{PC}\}$ indicator and the $\mathcal C_{t,s}/(nh_s)$
  normalization; the common-$h$ line uses $\mathcal B(t_{\mathrm{RD}})$, $\mathcal V^{\mathrm S}(t_{\mathrm{RD}})$ and
  cites `eq:common_h_opt`. The body stays at $p=1$ by design; App. B.4 is general $p$.
- **D4 — RESOLVED (Dor, 2026-09-08; rddid this commit).** The curvature-variance estimate
  $\widehat{\operatorname{Var}}(\widehat{\mathcal B}_\tau)=\bigl((p+1)!/h_0^{p+1}\bigr)^2\sum_i g^{\mathrm{diff}}_i{}^2$ now uses
  `g_diff = (a_c - a_bc) * res_b`: influence weights supported on the pilot window paired with the $q$-fit-at-$b$
  residuals, the same convention as the BC variance (D1). Before: `res_c` ($p$-fit at $h$). At the CCT pilot ratio
  ($b/h\approx1.5$) the two agree with Monte Carlo to within a few percent ($h^\star$ shift ~0.2%); at $b/h\ge3$ the
  $p$-fit version over-estimated (1.26× at 3, 2.2× at 6) while the $q$-fit version stays within 3% (adversarial review).
  Consequence: every `bwselect = "joint"`/`"iter"` bandwidth moves slightly; S5 macros regenerated in the same pass.
- **Plug-in choice disclosed (not a paper mismatch).** `.bw_joint` estimates every period's $\widehat{\mathcal B}_\tau$ at the
  RD period's CCT pilot $(h_0,b_0)$ rather than at period-specific pilots $b_\tau$; sensible for a common-$h$ target, but the
  Regularization paragraph writes $b_\tau$.
- **G1 (paper gap; partly closed).** Assumption `ass:reg`(f) now imposes $b_t\to0$, $n_tb_t^{2p+3}\to\infty$ and
  $h_t/b_t\to\bar\rho_t\in[0,\infty)$ (rd-did `bdd08e8`), which licenses the package's fixed pilot ratio. The paper
  still does not say *which* $b_t$ the aggregate selectors use in practice (package: the CCT ratio, $b^\star=h^\star b_0/h_0$
  in `joint`, $b_t=h_t\,b^{\mathrm{CCT}}_t/h^{\mathrm{CCT}}_t$ in `iter`). One sentence in §5.3 or B.4 would close it.

---

## 4. Drift control

- `tests/testthat/test-appB-conformance.R`: from-the-equations reference (matrix form, literal
  scaling, no HC1) vs `rd_period` / `.aggregate_fits`, 1e-10, three sampling designs, $p=1,2$.
- `tests/testthat/test-kernel-constants.R`: kernel constants vs Gauss–Legendre and closed forms.
- `tests/testthat/test-appB-bandwidth.R`: B.4 objectives and selectors vs hand-coded formulas.
- `dev/check_appB_labels.R`: label existence in `main.tex`.
- `dev/snapshot_rddid.R`: numerical regression snapshot (run before/after any change to
  `rd_period`, `aggregate.R`, `bandwidth.R`; diff the `.txt`).
