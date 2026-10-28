# Checking the identification assumptions

[Get
started](https://dorleventer.github.io/didisc/articles/didisc-estimation.md)
runs the four tests and reads their printouts. This article adds how
each statistic is built, the tests under the ATU design, their options,
and how the four results read together.

## When the tests apply

The tests are for a panel in which the running variable moves over time,
so that some units, the **switchers**, are above the cutoff in one
period and below it in another. The units at the cutoff in the RD period
can then be a different mix from the units at the cutoff in a comparison
period, and the difference of the two jumps can mix the treatment effect
with that difference in who is at the cutoff. A unit’s **type** records
its side of the cutoff in other periods; the table in [Get
started](https://dorleventer.github.io/didisc/articles/didisc-estimation.html#when-the-running-variable-moves-over-time)
gives the two versions the tests use. The diff-in-disc estimate
identifies the ATT when

1.  **type continuity** holds,
2.  **constant within-type confounding** holds, and
3.  **composition stability** or **homogeneous confounding** holds (one
    of the two is enough).

This is the paper’s result for a time-varying running variable (Leventer
and Nevo, 2024). There is one test per assumption. When the running
variable does not move, as in `didisc_sim`, every unit has the same type
in every period and the tests are not informative; in a repeated
cross-section the types are not observed.

The examples use `didisc_sim_pv`: the `didisc_sim` design (a confounding
jump of 0.5 in every year, a treatment effect of 1 in year 3, the RD
period) with a running variable that drifts between years. Most units
stay on one side; some switch:

``` r

library(didisc)
above <- didisc_sim_pv$R >= 0
switcher <- tapply(above, didisc_sim_pv$id, function(a) length(unique(a)) > 1)
sum(switcher)
#> [1] 176
```

Every test takes the same long data frame and column names as
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md).
By default each test uses the bias-corrected jumps with their robust
variances (`bc = TRUE`), like the `Robust (bias-corrected)` row of
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md).
The two tests on type shares pick a separate CCT bandwidth for every
local-linear fit (`bwselect = "cct"`) and combine their jumps with a
studentized sup-t statistic (`joint = "supt"`, `studentize = TRUE`); the
two tests on the outcome use one bandwidth for every cell
(`bwselect = "pooled"`). The last section says why.

## Type continuity: `rd_typecont()`

**H0: the share of each type jumps by zero at the cutoff, in every
period.**

``` r

tc <- rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
tc
#> Test of type continuity  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#>   Type in a period: a unit's sides of the cutoff in the other periods, in time order
#>     ("+" at or above, "-" below).
#> 
#>   Joint sup-t (t critical values) max |t| = 1.232 over 9 jumps,  p = 0.869  (largest: period:type 1:+-)
#>     Period 1: sup-t (t critical values) max |t| = 1.232 over 3 jumps,  p = 0.499
#>     Period 2: sup-t (t critical values) max |t| = 0.582 over 3 jumps,  p = 0.902
#>     Period 3: sup-t (t critical values) max |t| = 1.048 over 3 jumps,  p = 0.631
#>   All joint statistics: x$joint_tests
```

- `Types: ++, +-, -+, --`: a unit’s sides of the cutoff in the other two
  years, in time order (the first row of the type table in Get started).
- The test runs, in each year and for each type, a local-linear RD of
  the indicator “the unit is of this type” on the running variable; the
  jump is the difference between the type’s share just above and just
  below the cutoff.
- `Joint sup-t (studentized) ... over 9 jumps`: in each year the four
  shares sum to one, so three jumps are free; three years give 9. Each
  jump is divided by its standard error, and the largest \|t\| is
  compared with the sup-t critical value: the quantile of the largest
  \|t\| when every share is continuous, from the estimated correlation
  of the nine jumps and, since each standard error is itself an
  estimate, with t rather than normal tails. The printout names the year
  and type with the largest \|t\|. The `Period` lines test each year
  alone.

Here p = 0.869. Type continuity holds in `didisc_sim_pv` by
construction, so a p-value below 0.05 would be a false rejection.

**If it rejects:** units sort around the cutoff according to their side
in other periods, so the jump in the outcome partly reflects a jump in
who the units are. The diff-in-disc estimate can then be biased,
whatever the other tests say.

## Composition stability: `rd_compstable()`

**H0: the share of each type among the units just above the cutoff is
the same in the RD period and in each comparison period.**

``` r

cs <- rd_compstable(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
cs
#> Test of composition stability  [rd_compstable()]
#>   H0: the share of each type among the units just above the cutoff is the same in the RD period and in each comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Pair 3::1: sup-t (t critical values) max |t| = 3.086 over 3 jumps,  p = 0.009
#>     n above the cutoff: 490 (RD period), 508 (comparison), 430 in both
#>   Pair 3::2: sup-t (t critical values) max |t| = 2.059 over 3 jumps,  p = 0.120
#>     n above the cutoff: 490 (RD period), 497 (comparison), 418 in both
#> 
#>   Joint over pairs: sup-t (t critical values) max |t| = 3.086 over 6 jumps,  p = 0.016  (largest: pair:type 3::1:-+)
#>   All joint statistics: x$joint_tests
```

- One block per pair of the RD period and a comparison period:
  `Pair 3::1` compares year 3 with year 1. The types are those of type
  continuity (a unit’s sides in the other two years), so there are four
  types and three free shares (`over 3 jumps`).
- The test stacks the units above the cutoff in the two years, with the
  comparison year’s running variable reflected below an artificial
  cutoff, and tests the jump in each type share there.
- `n above the cutoff`: the units above the cutoff in the RD year, in
  the comparison year, and in both. Units above the cutoff in both years
  appear on both sides of the stacked regression; the standard error
  accounts for them.
- `Joint over pairs`: one sup-t test over the jumps of both pairs (six
  here), with their full covariance. The pairs share the RD-period
  units, and in a panel also the units above the cutoff in both
  comparison years, so their jumps are correlated; the covariance across
  pairs comes from those shared units.

Here p = 0.016. Composition stability is false in `didisc_sim_pv` by
construction: the units at the year-3 cutoff were mostly on the same
side in years 1 and 2, while the units at the year-1 cutoff are spread
evenly over the four types. A rejection here is a correct detection;
with 1,000 units the test does not always detect it.

**If it rejects:** the units at the cutoff in the RD period are a
different mix of types from those in a comparison period. If the
confounding jump also differs across types, the comparison-period jump
measures the confounding jump of a different mix, and the estimate can
be biased. If the confounding jump is the same for every type (the next
test), the mix does not matter and the estimate is still valid.

## Homogeneous confounding: `rd_homog()`

**H0: in each comparison period the confounding jump is the same for
every type.**

Here a unit’s type is its side of the cutoff in the RD period (the
second row of the type table in Get started; `type_by = "rd_side"`, the
default).

``` r

hg <- rd_homog(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
hg
#> Test of homogeneous confounding  [rd_homog()]
#>   H0: in each comparison period the confounding jump is the same for every type
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Comparison periods: 1, 2   Bandwidth: h = 0.2104 in every cell (CCT on the pooled comparison periods)
#> 
#>   Type: a unit's side of the cutoff in RD period 3 ("+" at or above, "-" below).
#>   Wald chi-squared(2) = 0.247,  p = 0.884
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Period     Type               jump       s.e.       n
#>     1          -                0.5347     0.2640     510  (reference)
#>     1          +                0.7070     0.4047     490
#>     2          -                1.0987     0.4183     510  (reference)
#>     2          +                0.8699     0.5244     490
```

- The table gives, in each comparison year, the jump in the outcome at
  the cutoff among the units below (`-`) and above (`+`) the cutoff in
  year 3. In a comparison period the jump is the confounding jump, so
  these are the confounding jumps by type.
- `Wald chi-squared(2)`: one contrast (type `+` against the reference
  type `-`) per comparison year.
- `Bandwidth`: one bandwidth in every cell, the CCT bandwidth of the
  outcome with both comparison years pooled.

Here p = 0.884. The null is true in `didisc_sim_pv` (the confounding
jump is 0.5 for every unit), so a p-value below 0.05 would be a false
rejection. The per-type jumps are imprecise (standard errors 0.26 to
0.52).

**If it rejects:** the confounding jump differs across types. If
composition stability holds, a confounding jump that differs across
types does not bias the estimate. If both nulls are false, the estimate
can be biased.

## Constant within-type confounding: `rd_trendcell()`

**H0: within each type, the confounding jump is the same in every
comparison period.**

``` r

tr <- rd_trendcell(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
tr
#> Test of constant within-type confounding  [rd_trendcell()]
#>   H0: within each type, the confounding jump is the same in every comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Comparison periods: 1, 2   Trend: constant
#>   Bandwidth: h = 0.2104 in every cell (CCT on the pooled comparison periods)
#> 
#>   Type: a unit's side of the cutoff in RD period 3 ("+" at or above, "-" below).
#>   Wald chi-squared(2) = 1.445,  p = 0.485
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Type       Period             jump       s.e.       n
#>     +          1                0.7070     0.4047     490  (reference)
#>     +          2                0.8699     0.5244     490
#>     -          1                0.5347     0.2640     510  (reference)
#>     -          2                1.0987     0.4183     510
```

- The per-cell jumps are those of
  [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md),
  compared the other way:
  [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md)
  compares types within a comparison period,
  [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)
  compares comparison periods within a type.
- `Wald chi-squared(2)`: one contrast (year 2 against year 1) per type.

This is a pre-trends check. The assumption concerns the RD period, where
the confounding jump cannot be separated from the effect, so the test
asks whether the within-type confounding jump is stable across the
comparison periods. Not rejecting does not establish the assumption in
the RD period. Here p = 0.485.

With `trend = "linear"` the null is that the within-type confounding
jump moves linearly across the comparison periods. That needs at least
three comparison periods; with two, the function says so and returns no
statistic.

**If it rejects:** the within-type confounding jump changes over time,
and constant weights do not cancel it. `didisc(trend = "linear")` allows
a within-type confounding jump that moves linearly in time, and
`rd_trendcell(trend = "linear")` tests that version (three or more
comparison periods).

## ATU designs

When everybody is treated in the comparison periods
(`estimand = "atu"`), the estimate concerns the units just below the
cutoff, and composition stability concerns the below-cutoff shares.
`rd_compstable(estimand = "atu")` mirrors the running variable and runs
the same test on the units below the cutoff; the other three tests give
the same numbers with either estimand, so `estimand` only labels them.
`didisc_sim_pv` is an ATT design, so the call below only shows the
mechanics:

``` r

rd_compstable(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu")
#> Test of composition stability  [rd_compstable()]
#>   H0: the share of each type among the units just below the cutoff is the same in the RD period and in each comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Estimand: ATU (the units below the cutoff are the ones untreated in the RD period, so the test is on their shares (mirrored design))
#>   RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Pair 3::1: sup-t (t critical values) max |t| = 1.901 over 3 jumps,  p = 0.166
#>     n below the cutoff: 510 (RD period), 492 (comparison), 432 in both
#>   Pair 3::2: sup-t (t critical values) max |t| = 1.705 over 3 jumps,  p = 0.240
#>     n below the cutoff: 510 (RD period), 503 (comparison), 431 in both
#> 
#>   Joint over pairs: sup-t (t critical values) max |t| = 1.901 over 6 jumps,  p = 0.293  (largest: pair:type 3::1:++)
#>   All joint statistics: x$joint_tests
```

In this printout the H0 line and the counts refer to the units below the
cutoff: under the ATU design those are the units untreated in the RD
period, whose composition the estimate relies on.

## Run all four at once: `rd_checks()`

``` r

ch <- rd_checks(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
ch
#> Tests of the identification assumptions  [rd_checks()]
#> 
#>   Assumption                          Statistic        p   Statistic type             Bandwidth
#>   type continuity                         1.232    0.869   sup-t, t critical values   CCT per cell
#>   composition stability                   3.086    0.016   sup-t, t critical values   CCT per cell
#>   homogeneous confounding                 0.247    0.884   Wald chi-squared(2)        pooled CCT, h = 0.21
#>   constant within-type confounding        1.445    0.485   Wald chi-squared(2)        pooled CCT, h = 0.21
#> 
#>   Reading at level 0.95, rejecting when p < 0.05 (not rejecting does not establish an assumption):
#>     type continuity                    not rejected
#>     composition stability              rejected     } identification needs
#>     homogeneous confounding            not rejected } only one of these two
#>     constant within-type confounding   not rejected (checked in the comparison periods)
#>   Each test's details: x$typecont, x$compstable, x$homog, x$trendcell.
```

[`rd_checks()`](https://dorleventer.github.io/didisc/reference/rd_checks.md)
runs the four tests with one set of options and prints one row per test,
with the statistic, its type and the bandwidth rule. The four results
are inside (`ch$typecont`, …), and `tidy(ch)` returns the table for a
paper or a table maker such as modelsummary. The statistic column holds
the largest \|t\| for the two share tests (sup-t) and the chi-squared
statistic for the two outcome tests.

Together, the four results map onto the list at the top. A rejection of
type continuity, or of constant within-type confounding, would call the
estimate into question; so would a rejection of both composition
stability and homogeneous confounding. A rejection of only one of those
two would not on its own, because either one is enough for
identification. In `didisc_sim_pv` the truth is known: composition
stability fails and the other three assumptions hold, so the estimate is
valid:

``` r

didisc(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Bandwidth: common h = 0.242 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b = 0.3679 (also chosen for the diff-in-disc estimate)
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9689     0.1909    5.08   <0.001   [0.5948, 1.3430]
#>   Robust (bias-corrected)      0.9735     0.2332    4.17   <0.001   [0.5164, 1.4305]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

The true effect is 1.

## Why sup-t, and a pooled bandwidth for the outcome tests

In simulations with a clean null (every assumption true, 1,000 units,
three periods, 5% level), the chi-squared Wald statistics with per-cell
bandwidths rejected a true null too often: 9% for the share tests and 8%
for the outcome tests. The excess comes from the estimated covariance of
the jumps, not from the jumps themselves, in two ways.

- A bandwidth chosen on each cell’s own data is noisy, and the variance
  formula treats it as fixed, so each variance comes out 10-15% too
  small. For the outcome tests, whose cells are thin, one bandwidth
  chosen on the pooled comparison periods removes this: 6% instead of
  8%, at about 8 points of power against a violation confined to one
  cell.
- A chi-squared Wald inverts a 6-by-6 or 9-by-9 covariance estimated
  from a few hundred units near the cutoff, and the error in that
  inverse grows with the number of jumps (15% at 28 jumps, four
  periods). sup-t and sWald use the covariance only through its
  correlations, so they avoid that inversion: 6-7% at three periods,
  4-5% at four. They still divide each jump by its estimated standard
  error, so a variance that is too small, or noisy, inflates them as
  well, by less than the Wald; that is why the share tests stay slightly
  above 5% at three periods. For the share tests the per-fit bandwidth
  is kept, because it adapts to a sorting window and every pooled
  bandwidth tried halved the power against local sorting.

The default sup-t is studentized: each standard error is an estimate,
with the degrees of freedom of its robust variance (about 50-150 here,
computed from each jump’s regression weights), and the critical value
accounts for that noise. With the same simulations the studentized sup-t
rejected a true null about 6% of the time at three periods and 4-5% at
four, against 7% and 5-6% for sup-t with normal critical values, with
the same power once the sizes are matched.

With two periods (one jump per test and period) the statistics are the
same two-sided test, and the share tests did not over-reject.

## Switching options

The defaults are the settings recommended above; every alternative is
one argument away, in each test and in
[`rd_checks()`](https://dorleventer.github.io/didisc/reference/rd_checks.md).

``` r

# share tests: the chi-squared Wald, or sup-t with normal critical values
rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, joint = "wald")
rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, studentize = FALSE)
# outcome tests: one bandwidth per cell instead of one pooled bandwidth
rd_homog(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct")
# all four at once, and one test differently from the others
rd_checks(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
          joint = "wald", outcome_bwselect = "cct", homog = list(min_n = 20))
```

Every statistic of a share test is also in `$joint_tests`, so the
alternatives can be compared without running the test again:

``` r

ch$typecont$joint_tests
#>   method   reference statistic df   p_value
#> 1   supt           t  1.231775  9 0.8688813
#> 2   supt      normal  1.231775  9 0.8609214
#> 3  swald           t  5.110313  9 0.7953520
#> 4  swald      normal  5.110313  9 0.7861821
#> 5   wald chi-squared  4.781633  9 0.8529133
```

## Options

- `bwselect`: for
  [`rd_typecont()`](https://dorleventer.github.io/didisc/reference/rd_typecont.md)
  and
  [`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md),
  `"cct"` (default) chooses a CCT bandwidth for every local-linear fit
  and `"rot"` uses a rule of thumb; for
  [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md)
  and
  [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md),
  `"pooled"` (default) uses the CCT bandwidth of the outcome on the
  pooled comparison periods in every cell, `"cct"` one per cell, `"rot"`
  a rule of thumb. `h` fixes the bandwidth.
- `joint`
  ([`rd_typecont()`](https://dorleventer.github.io/didisc/reference/rd_typecont.md),
  [`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md)):
  `"supt"` (default), `"swald"` (the sum of squared t’s) or `"wald"`
  (the chi-squared Wald); `studentize`: `TRUE` (default, t-based
  critical values for sup-t and sWald) or `FALSE` (normal). Every
  combination is in `$joint_tests`.
- `bc`: `TRUE` (default) tests the bias-corrected jumps with their
  robust variances; `FALSE` uses the conventional ones.
- `scheme`: detected from the data as in
  [`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md);
  set it to override.
- [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md):
  `type_by = "rd_side"` (default, the side in the RD period) or
  `"pattern"` (the sides in the other periods used).
  [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md)
  always types units by their side in the RD period.
- [`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md)
  and
  [`rd_trendcell()`](https://dorleventer.github.io/didisc/reference/rd_trendcell.md):
  `min_n` drops cells with too few observations on either side of the
  cutoff, with a message.

The article [How didisc() computes the
estimate](https://dorleventer.github.io/didisc/articles/didisc-how-it-works.md)
rebuilds the type-continuity and composition-stability tests by hand and
runs all four on simulated designs in which one assumption fails at a
time, with the resulting bias in the estimate.

## Reference

Leventer, D. and D. Nevo (2024). Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data. arXiv:2408.05847.
<https://arxiv.org/abs/2408.05847>
