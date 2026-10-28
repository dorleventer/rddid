# Package index

## Package overview

- [`didisc-package`](https://dorleventer.github.io/didisc/reference/didisc-package.md)
  : didisc: treatment effects at a cutoff shared with a confounding
  policy

## Estimate

The diff-in-disc estimator and what to do with a fit.

- [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
  : Estimate the effect of a treatment at a cutoff shared with a
  confounding policy
- [`summary(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/summary.didisc.md)
  [`print(`*`<summary.didisc>`*`)`](https://dorleventer.github.io/didisc/reference/summary.didisc.md)
  : Summary of a diff-in-disc fit
- [`plot(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/plot.didisc.md)
  : Plot the per-period RD fits behind a diff-in-disc estimate
- [`coef(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-methods.md)
  [`confint(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-methods.md)
  [`nobs(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-methods.md)
  : Coefficients, confidence intervals and sample size of a diff-in-disc
  fit
- [`tidy(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  [`glance(`*`<didisc>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  [`tidy(`*`<rd_typecont>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  [`tidy(`*`<rd_compstable>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  [`tidy(`*`<rd_homog>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  [`tidy(`*`<rd_trendcell>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  [`tidy(`*`<rd_checks>`*`)`](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)
  : Tidy output for diff-in-disc fits and validation tests
- [`reexports`](https://dorleventer.github.io/didisc/reference/reexports.md)
  [`tidy`](https://dorleventer.github.io/didisc/reference/reexports.md)
  [`glance`](https://dorleventer.github.io/didisc/reference/reexports.md)
  : Objects exported from other packages

## Check the assumptions

Tests of the identification assumptions when the running variable moves
over time, so that units can change side of the cutoff between periods.

- [`rd_checks()`](https://dorleventer.github.io/didisc/reference/rd_checks.md)
  : Run the four tests of the identification assumptions
- [`rd_typecont()`](https://dorleventer.github.io/didisc/reference/rd_typecont.md)
  : Test of type continuity
- [`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md)
  : Test of composition stability
- [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md)
  : Test of homogeneous confounding
- [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)
  : Test of constant within-type confounding
- [`plot(`*`<rd_typecont>`*`)`](https://dorleventer.github.io/didisc/reference/plot.rd_typecont.md)
  : Plot a type-continuity test: one comparison period against one RD
  period
- [`plot(`*`<rd_compstable>`*`)`](https://dorleventer.github.io/didisc/reference/plot.rd_compstable.md)
  : Plot a composition-stability test: the reflected sample of one pair
  of periods
- [`plot(`*`<rd_homog>`*`)`](https://dorleventer.github.io/didisc/reference/plot.rd_homog.md)
  : Plot a homogeneous-confounding test: the confounding jump of each
  type, by comparison period
- [`plot(`*`<rd_trendcell>`*`)`](https://dorleventer.github.io/didisc/reference/plot.rd_trendcell.md)
  : Plot a constant-within-type-confounding test: each type's
  confounding jump over time
- [`plot_switchers()`](https://dorleventer.github.io/didisc/reference/plot_switchers.md)
  : Plot the switchers: the running variable in one period against
  another

## Building blocks

The local-linear RD in one period and its CCT bandwidths, for working
period by period.

- [`rd_period()`](https://dorleventer.github.io/didisc/reference/rd_period.md)
  : Local-linear RD in one period (building block)
- [`rd_bw_cct()`](https://dorleventer.github.io/didisc/reference/rd_bw_cct.md)
  : CCT bandwidths for one period (building block)

## Data

Simulated panels for trying things out.

- [`didisc_sim`](https://dorleventer.github.io/didisc/reference/didisc_sim.md)
  : Simulated diff-in-disc panel (running variable fixed over time)
- [`didisc_sim_pv`](https://dorleventer.github.io/didisc/reference/didisc_sim_pv.md)
  : Simulated diff-in-disc panel (running variable moves over time)
