# didisc: treatment effects at a cutoff shared with a confounding policy

A treatment of interest switches on at a cutoff of a running variable in
one period, the **RD period**. A **confounding policy** switches at the
same cutoff, in every period, so the jump in the outcome at the cutoff
in the RD period mixes the treatment effect with the **confounding
jump**. In the **comparison periods** the treatment of interest is
uniform at the cutoff (nobody treated, or everybody treated), so,
provided the treatment of interest has no anticipation or carry-over
effects there (which the paper assumes), the jump there *is* the
confounding jump.
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
estimates the jump in every period by local-linear RD and subtracts a
weighted average of the comparison-period jumps from the RD-period jump.
How the weights are set is the **confounding-trend assumption**:
constant (equal weights) or linear in time.

## Workflow

1.  **Estimate.**
    `fit <- didisc(data, y = , x = , time = , id = , t_rd = )`; printing
    `fit` shows the estimate, its standard error and confidence
    interval, the comparison periods with their weights, the sampling
    scheme and the bandwidth. See
    [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md).

2.  **Look inside.** `summary(fit)` shows the jump in every period and
    the standard error under each sampling scheme
    ([`summary.didisc()`](https://dorleventer.github.io/didisc/reference/summary.didisc.md));
    [`coef()`](https://rdrr.io/r/stats/coef.html),
    [`confint()`](https://rdrr.io/r/stats/confint.html) and
    [`nobs()`](https://rdrr.io/r/stats/nobs.html) work as usual
    ([didisc-methods](https://dorleventer.github.io/didisc/reference/didisc-methods.md)),
    and [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
    [`glance()`](https://generics.r-lib.org/reference/glance.html) feed
    table makers
    ([didisc-tidiers](https://dorleventer.github.io/didisc/reference/didisc-tidiers.md)).

3.  **Check the assumptions.** When the running variable moves over
    time, units can change side of the cutoff between periods, and four
    tests check the assumptions this adds: type continuity
    ([`rd_typecont()`](https://dorleventer.github.io/didisc/reference/rd_typecont.md)),
    composition stability
    ([`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md)),
    homogeneous confounding
    ([`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md))
    and constant within-type confounding
    ([`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)).

The building blocks
[`rd_period()`](https://dorleventer.github.io/didisc/reference/rd_period.md)
(the local-linear RD in one period) and
[`rd_bw_cct()`](https://dorleventer.github.io/didisc/reference/rd_bw_cct.md)
(its CCT bandwidths) are exported for users who want to work period by
period.
[didisc_sim](https://dorleventer.github.io/didisc/reference/didisc_sim.md)
and
[didisc_sim_pv](https://dorleventer.github.io/didisc/reference/didisc_sim_pv.md)
are simulated panels for trying things out. `citation("didisc")` gives
the reference below.

## References

Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data.*
arXiv:2408.05847. <https://arxiv.org/abs/2408.05847>

## See also

The package website, with a getting-started guide and articles on the
computation, the options and the tests:
<https://dorleventer.github.io/didisc/>

## Author

**Maintainer**: Dor Leventer <leventerdor@gmail.com> \[copyright
holder\]

Authors:

- Daniel Nevo
