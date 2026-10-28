# Get started: estimate a diff-in-disc effect

## The setting

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
weighted average of the comparison-period jumps from the RD-period jump:
``` math
\mathrm{ATT}(t_{\mathrm{RD}}) = D_{t_{\mathrm{RD}}} - \sum_t w_t D_t ,
```
where $`t_{\mathrm{RD}}`$ is the RD period, $`D_t`$ is the jump in the
outcome at the cutoff in period $`t`$, the sum runs over the comparison
periods, and $`w_t`$ is the weight of comparison period $`t`$. How the
weights are set is the **confounding-trend assumption**: equal weights
if the confounding jump is the same in every period, a straight-line
extrapolation if it moves linearly in time.

## The data

``` r

library(didisc)
head(didisc_sim)
#>   id year          R V W           Y
#> 1  1    1 -0.7785626 0 0 -0.38918276
#> 2  2    1 -0.9019249 0 0 -0.75738418
#> 3  3    1  0.1621249 1 0  1.10233495
#> 4  4    1  0.1695525 1 0  1.50736287
#> 5  5    1 -0.7321141 0 0 -0.04729411
#> 6  6    1  0.8406841 1 0  1.59100391
```

`didisc_sim` is a simulated panel of 1,000 units observed in three
years, in long format: one row per unit and year.
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
uses four columns: the outcome `Y`, the running variable `R` (cutoff 0),
the period `year` and the unit `id`. The other two record the design:
the confounding policy `V` is on above the cutoff in every year and
raises `Y` by 0.5; the treatment of interest `W` is on above the cutoff
in year 3 only and raises `Y` by 1. So year 3 is the RD period, years 1
and 2 are comparison periods, and **the true effect is 1**.

With your own data,
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
needs the same shape: a data frame in long format, one row per unit and
period, with the outcome, the running variable, the period and, for a
panel, a unit identifier (the panel need not be balanced). The cutoff is
`c` (default 0), and a unit with `x >= c` is above it. The period column
is usually numeric, a year; a character column works with the default
constant trend, while `trend = "linear"` fits its line on the period
values and so needs a numeric column. By default the comparison periods
are all periods other than `t_rd`, in sorted order.

## Estimate

``` r

fit <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2694 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4136, 2 = 0.3901, 3 = 0.3984
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0922     0.1259    8.67   <0.001   [0.8454, 1.3391]
#>   Robust (bias-corrected)      1.1411     0.1488    7.67   <0.001   [0.8494, 1.4328]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

Line by line:

- **`Diff-in-disc estimate of the ATT in period 3`**: the estimand and
  the RD period (`t_rd = 3`). The ATT is the effect of the treatment on
  the units just above the cutoff, which are the ones treated in the RD
  period. If everybody is treated in the comparison periods,
  `estimand = "atu"` labels the same numbers as the ATU (see [the ATU
  section](#atu) below).
- **`Comparison periods: 1, 2 (constant confounding trend; weights 0.5, 0.5)`**:
  by default every period other than `t_rd` is a comparison period. The
  default `trend = "constant"` assumes the confounding jump is the same
  in every period, so each comparison period gets the same weight and
  their average stands in for the confounding jump in year 3.
- **`Sampling scheme`**: how the data were sampled, read off `id` and
  `R`: here a panel (the same units every year) in which no unit changes
  side of the cutoff, because the running variable does not move. The
  standard error takes into account that the same units enter every
  year’s fit. Under the default bandwidth rule the scheme can also move
  the bandwidth, and so the estimate; [Bandwidth rules and sampling
  schemes](https://dorleventer.github.io/didisc/articles/didisc-options.md)
  shows how.
- **`Bandwidth: common h`**: the default rule, `bwselect = "joint"`,
  uses one main bandwidth `h` in every period, chosen to minimize the
  asymptotic mean squared error of the diff-in-disc estimate.
  **`Pilot bandwidth b (period = value)`** lists the bandwidth of the
  bias correction in each period, in time order.
- **The two rows.** `Conventional` is the local-linear estimate with its
  conventional standard error. `Robust (bias-corrected)` is the
  bias-corrected estimate with its robust standard error (Calonico,
  Cattaneo and Titiunik, 2014). Each row has its z statistic, p-value
  and 95% confidence interval.

The true effect is 1; the two estimates are 1.09 and 1.14. A plain RD in
year 3 alone would estimate the effect plus the confounding jump, 1 +
0.5 = 1.5. The year-3 fit is stored in the object:

``` r

fit$fits[["3"]]
#> Single-period RD (p=1, h=0.2694, b=0.3984, kernel=triangular, n=1000)
#>   D (conventional)   = +1.664  (se 0.1635)
#>   D (bias-corrected) = +1.6856  (se 0.1924)
```

Its jump, 1.66, is that naive RD estimate, whose target is 1.5, not 1.
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
subtracts the average of the year-1 and year-2 jumps, (0.51 + 0.63) / 2
= 0.57, its estimate of the confounding jump of 0.5.

## Summary, coefficients and tables

``` r

summary(fit)
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2694 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4136, 2 = 0.3901, 3 = 0.3984
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0922     0.1259    8.67   <0.001   [0.8454, 1.3391]
#>   Robust (bias-corrected)      1.1411     0.1488    7.67   <0.001   [0.8494, 1.4328]
#> 
#>   Per-period local-linear fits (Conventional = sum of coef x jump; Robust = sum of coef x jump (bc)):
#>   period   role           coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison   -0.500   1000   0.2694   0.4136     0.5129    0.1683     0.4642    0.1942
#>   2        comparison   -0.500   1000   0.2694   0.3901     0.6306    0.1679     0.6249    0.2004
#>   3        RD            1.000   1000   0.2694   0.3984     1.6640    0.1635     1.6856    0.1924
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
#>     cross-section 0.2377   panel, no unit changes side 0.1488   panel, some change side 0.1488
```

[`summary()`](https://rdrr.io/r/base/summary.html) adds the per-period
fits. Each row gives the period’s role, its coefficient in the estimate
(+1 for the RD period, minus its weight for a comparison period), its
number of observations, its bandwidths, and its jump with standard
error, conventional (`jump`, `s.e.`) and bias-corrected (`jump (bc)`,
`s.e. (rb)`). The estimate is the sum of coefficient times jump: 1.664 -
0.5 × 0.513 - 0.5 × 0.631 = 1.092.

The last lines give the robust standard error under each of the three
sampling schemes. The printout uses the one for the scheme detected in
the data; the other two show what the formulas for the other schemes
give on the same data. Here the cross-section standard error, 0.238,
ignores that the same units appear in every year; the panel one, which
accounts for it, is 0.149.

``` r

coef(fit)
#> Conventional       Robust 
#>     1.092242     1.141077
confint(fit)
#>                  2.5 %   97.5 %
#> Conventional 0.8454070 1.339077
#> Robust       0.8493838 1.432769
tidy(fit)
#>           term estimate std.error statistic      p.value  conf.low conf.high
#> 1 Conventional 1.092242 0.1259386  8.672816 4.215617e-18 0.8454070  1.339077
#> 2       Robust 1.141077 0.1488256  7.667206 1.757833e-14 0.8493838  1.432769
```

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the two point
estimates and [`confint()`](https://rdrr.io/r/stats/confint.html) their
confidence intervals.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns one
row per estimate, so table makers such as modelsummary work with
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
fits; [A table for a paper](#a-table-for-a-paper) below has an example.

## A confounding jump that moves linearly

``` r

fit_lin <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                 trend = "linear")
fit_lin
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (linear confounding trend)
#>   Comparison weights: 1: -1.000, 2: 2.000
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.306 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4698, 2 = 0.4431, 3 = 0.4525
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9198     0.2423    3.80   <0.001   [0.4448, 1.3948]
#>   Robust (bias-corrected)      0.8801     0.2882    3.05    0.002   [0.3153, 1.4449]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

With `trend = "linear"` the confounding jump may move linearly in time.
The line through the year-1 and year-2 jumps is extrapolated to year 3,
which gives weights -1 and 2: the confounding jump in year 3 is
estimated by $`2 D_2 - D_1`$. This needs at least two comparison
periods. The standard error grows (0.242 against 0.126) because
extrapolating a line amplifies noise. In `didisc_sim` the confounding
jump is constant, so both estimates are valid; the linear one gives up
precision to allow for a linear trend. A numeric vector instead of
`"constant"` or `"linear"` sets the weights directly, one per comparison
period.

## Several RD periods

When the treatment of interest is in place in several periods, `t_rd`
can list them all.
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)
then estimates a weighted average of their effects, with weights
`rd_weights`: equal by default, or any finite weights, which define the
average being estimated. Each period’s jump is estimated once, and each
comparison period’s weight is the `rd_weights`-weighted sum of the
weights each RD period alone would give it. `didisc_sim` has one RD
period, so for illustration year 3 is copied as a second RD year, 4:

``` r

sim4 <- rbind(didisc_sim, transform(didisc_sim[didisc_sim$year == 3, ], year = 4))
fit_avg <- didisc(sim4, y = "Y", x = "R", time = "year", id = "id", t_rd = c(3, 4),
                 bwselect = "cct")
fit_avg
#> Diff-in-disc estimate of the ATT averaged over RD periods 3, 4
#>   RD-period weights: 3: 0.500, 4: 0.500
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights (summed over the RD periods): 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: per-period CCT MSE-optimal (rule "cct"): h (period = value) 1 = 0.2637, 2 = 0.3256, 3 = 0.4141, 4 = 0.4141
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0396     0.1170    8.88   <0.001   [0.8103, 1.2689]
#>   Robust (bias-corrected)      1.0964     0.1407    7.79   <0.001   [0.8205, 1.3722]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
fit_cct3 <- didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                  bwselect = "cct")
```

`bwselect = "cct"` gives each period its own MSE-optimal bandwidth (from
rdrobust’s `rdbwselect`; Calonico, Cattaneo and Titiunik, 2014), so the
two copies get the same bandwidth, the two RD-period jumps are
identical, and their average is the year-3 estimate (1.040 against
1.040). The standard error is the year-3 one too (0.117 against 0.117):
in a panel it includes the covariance between the RD-period jumps, which
here are perfectly correlated. The cross-section formula leaves that
covariance out, and so would report the average as more precise than
year 3 alone (0.148 against 0.176), crediting the copy with information
it does not have. Under the default bandwidth rule the two estimates
differ slightly: the rule penalizes noisy estimates of the bias (the
`regularize` option of
[`didisc()`](https://dorleventer.github.io/didisc/reference/didisc.md)),
and that penalty adds up the periods’ variances without their
covariances. When the RD periods hold different units or draws,
averaging them lowers the part of the variance that comes from the
RD-period jumps; under a constant trend the part that comes from the
comparison jumps, which every RD period shares, stays.

## A table for a paper

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives the
estimates and
[`glance()`](https://generics.r-lib.org/reference/glance.html) the
design (RD periods and their weights, comparison periods, trend,
bandwidth rule, `h`, scheme), which is what table makers read. With the
modelsummary package installed, the two fits side by side, with
confidence intervals:

``` r

modelsummary::modelsummary(
  list(Constant = fit, Linear = fit_lin), statistic = "conf.int",
  gof_map = list(list(raw = "h", clean = "h", fmt = 3),
                 list(raw = "nobs", clean = "Num.Obs.", fmt = 0)))
```

|              | Constant         | Linear           |
|--------------|------------------|------------------|
| Conventional | 1.092            | 0.920            |
|              | \[0.845, 1.339\] | \[0.445, 1.395\] |
| Robust       | 1.141            | 0.880            |
|              | \[0.849, 1.433\] | \[0.315, 1.445\] |
| h            | 0.269            | 0.306            |
| Num.Obs.     | 3000             | 3000             |

The same call with `output = "didisc_table.tex"` (or `.html`, `.md`,
`.docx`) writes the table to a file;
`write.csv(tidy(fit), "didisc_estimate.csv", row.names = FALSE)` writes
the estimate table to a plain CSV file.

`gof_map` sets which rows of
[`glance()`](https://generics.r-lib.org/reference/glance.html) appear
below the estimates and their digits (`fmt`): here `h` with three digits
and the number of observations. Without it, every column of
[`glance()`](https://generics.r-lib.org/reference/glance.html) is listed
and `h` prints with all its digits. The last line writes the estimate
table to a plain CSV file.

## When the running variable moves over time

In `didisc_sim` every unit stays on its side of the cutoff. In many
applications the running variable moves (a population count, a test
score, an income), and some units, the **switchers**, are above the
cutoff in one period and below it in another. The units at the cutoff in
the RD period can then be a different mix from the units at the cutoff
in a comparison period, and the difference of jumps can be biased. A
unit’s **type** records its side of the cutoff in periods other than the
one at hand. The tests use two versions:

| Test | A unit’s type in period *t* | Printed as |
|----|----|----|
| type continuity, composition stability | its sides of the cutoff in all the other periods, in time order | `++`, `+-`, `-+`, `--` (`+` above, `-` below) |
| homogeneous confounding, constant within-type confounding | its side of the cutoff in the RD period (the default, `type_by = "rd_side"`) | `+`, `-`; in plots “Above in 3”, “Below in 3” |

With two periods, both versions are the unit’s side of the cutoff in the
other period.

The diff-in-disc estimate identifies the ATT when the type distribution
is continuous at the cutoff and the confounding jump within each type is
constant over time, and in addition either the type composition is
stable across periods or the confounding jump is the same for every type
(Leventer and Nevo, 2024). Four tests check these assumptions.
`didisc_sim_pv` is the `didisc_sim` design with a running variable that
drifts between years. The estimate whose assumptions the four tests
below check:

``` r

fit_pv <- didisc(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_pv
#> Diff-in-disc estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Bandwidth: common h = 0.2436 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4001, 2 = 0.3958, 3 = 0.4047
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9698     0.1901    5.10   <0.001   [0.5972, 1.3424]
#>   Robust (bias-corrected)      0.9943     0.2256    4.41   <0.001   [0.5522, 1.4364]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

The sampling scheme is now a panel in which some units change side of
the cutoff between years.

**Type continuity.** H0: the share of each type jumps by zero at the
cutoff, in every period. A rejection means units sort around the cutoff
by their side in other periods, so the jump in the outcome partly
reflects who the units are; the estimate can then be biased, whatever
the other tests say.

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
#>     Period 1: sup-t (t critical values) max |t| = 1.232 over 3 jumps,  p = 0.500
#>     Period 2: sup-t (t critical values) max |t| = 0.582 over 3 jumps,  p = 0.902
#>     Period 3: sup-t (t critical values) max |t| = 1.048 over 3 jumps,  p = 0.633
#>   All joint statistics: x$joint_tests
```

Every test prints its null (H0) and the sampling scheme read off the
data. `Types: ++, +-, -+, --` are the types of the first row of the
table above. `Bandwidth: CCT MSE-optimal, chosen per cell` means that
each local-linear fit the test runs (one per period and type, a *cell*)
gets its own MSE-optimal bandwidth from
[`rdrobust::rdbwselect`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
(Calonico, Cattaneo and Titiunik, 2014; hence “CCT”).
`Joint sup-t (studentized) ... over 9 jumps` tests all the jumps in the
type shares at once: in each year the four shares sum to one, so three
jumps are free, and three years give 9. It compares the largest \|t\| of
the nine jumps with the sup-t critical value and names the year and type
with the largest \|t\| (`joint = "wald"` gives the chi-squared Wald
instead; see [Checking the identification
assumptions](https://dorleventer.github.io/didisc/articles/didisc-validation-tests.md)
for why sup-t is the default). Each `Period` line tests one year alone.

**Composition stability.** H0: the share of each type among the units
just above the cutoff is the same in the RD period and in each
comparison period. A rejection means the units at the cutoff in the RD
period are a different mix of types from those in a comparison period;
the estimate is then valid only if the next test’s null holds.

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

`Pair 3::1` compares year 3, the RD period, with year 1; its sup-t test
covers the three free jumps in the type shares. `n above the cutoff`
counts the units above the cutoff in the RD year, in the comparison
year, and in both. `Joint over pairs` tests the jumps of both pairs
together, with the covariance that comes from the year-3 units both
pairs share.

**Homogeneous confounding.** H0: in each comparison period the
confounding jump is the same for every type. Here, and in the next test,
the types are those of the second row of the table above: a unit’s side
of the cutoff in the RD period. A rejection means the confounding jump
differs by type; together with a rejection of composition stability, the
comparison-period jumps do not measure the confounding jump of the units
at the RD-period cutoff, and the estimate can be biased.

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

The table gives the confounding jump of each type in each comparison
year, with its standard error and number of observations, all at one
bandwidth (`Bandwidth`: the CCT bandwidth of the outcome with the
comparison years pooled). The Wald statistic has one contrast per year,
type `+` against the type marked `(reference)`, so two comparison years
give `chi-squared(2)`.

**Constant within-type confounding.** H0: within each type, the
confounding jump is the same in every comparison period. This is a
pre-trends check: the assumption concerns the RD period, where the
confounding jump cannot be seen apart from the effect, so the test asks
whether it is stable across the comparison periods. A rejection means
the constant-trend weights do not cancel the confounding jump.
`trend = "linear"` allows a confounding jump that moves linearly in
time; `rd_trendcell(trend = "linear")` tests that version (three or more
comparison periods).

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

The jumps are those of
[`rd_homog()`](https://dorleventer.github.io/didisc/reference/rd_homog.md),
grouped by type. The Wald statistic has one contrast per type: year 2
against the year marked `(reference)`, the type’s first comparison year.

One table for all four:

``` r

do.call(rbind, lapply(list(tc, cs, hg, tr), tidy))
#>            test                       assumption statistic df    p.value
#> 1   rd_typecont                  type continuity 1.2317753 NA 0.86933131
#> 2 rd_compstable            composition stability 3.0857080 NA 0.01639984
#> 3      rd_homog          homogeneous confounding 0.2474976  2 0.88360178
#> 4  rd_trendcell constant within-type confounding 1.4453616  2 0.48544912
#>                     method
#> 1 sup-t, t critical values
#> 2 sup-t, t critical values
#> 3      Wald chi-squared(2)
#> 4      Wald chi-squared(2)
```

The composition-stability row is the joint test over pairs; the per-pair
tests are in its printout above. For the two share tests the statistic
column holds the largest \|t\| (sup-t), for the two outcome tests the
chi-squared statistic.

In `didisc_sim_pv` three of the four nulls are true by construction:
types are continuous at the cutoff, and the confounding jump is 0.5 for
every unit in every year. Composition stability is false: units at the
year-3 cutoff were mostly on the same side in years 1 and 2, while units
at the year-1 cutoff are spread evenly over the four types. In this
sample the four p-values are 0.869, 0.016, 0.884 and 0.485 (in the order
of the table). A p-value below 0.05 for a true null is a false
rejection, which a test at the 5% level is designed to make in about one
sample in twenty. Because the other three assumptions hold, homogeneous
confounding among them, the failure of composition stability does not
bias `fit_pv`: its estimate is 0.97, against a true effect of 1. The
article on the identification assumptions shows how each statistic is
built and the options of the tests.

## Comparison periods where everybody is treated: the ATU

When everybody is treated in the comparison periods, rather than nobody,
the same difference of jumps estimates the ATU: the effect for the units
just below the cutoff, which are untreated in the RD period.
`estimand = "atu"` labels the output as the ATU: the estimate, standard
errors and bandwidths are the same as without it (among the tests, it
changes only
[`rd_compstable()`](https://dorleventer.github.io/didisc/reference/rd_compstable.md)).
For illustration, make everybody treated in years 1 and 2 of
`didisc_sim`, with the same effect of 1:

``` r

sim_atu <- transform(didisc_sim, W = ifelse(year < 3, 1, W), Y = Y + (year < 3))
didisc(sim_atu, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu")
#> Diff-in-disc estimate of the ATU in period 3
#>   Comparison periods: 1, 2   (constant confounding trend)
#>   Comparison weights: 1: 0.500, 2: 0.500
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2694 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4136, 2 = 0.3901, 3 = 0.3984
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0922     0.1259    8.67   <0.001   [0.8454, 1.3391]
#>   Robust (bias-corrected)      1.1411     0.1488    7.67   <0.001   [0.8494, 1.4328]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

The numbers are those of `fit`, now labeled the ATU; the true ATU here
is 1.

## What the paper reports

- **The estimate.** The paper reports both rows. It bases its confidence
  intervals and p-values on the `Robust (bias-corrected)` row, following
  Calonico, Cattaneo and Titiunik (2014): at a bandwidth chosen to
  minimize the mean squared error, the conventional interval is centered
  on a biased estimate and can under-cover.
- **The design.** The choices behind an estimate (the RD period, the
  comparison periods, the trend assumption with its weights, the
  bandwidth rule and `h`, the sampling scheme) are in the first lines of
  the printout; `glance(fit)` returns them as a one-row table.
- **The tests.** When the running variable moves over time, the table
  above collects the four tests of the assumptions the estimate relies
  on.

## Next

- [Checking the identification
  assumptions](https://dorleventer.github.io/didisc/articles/didisc-validation-tests.md):
  how each test statistic is built, the tests under the ATU design, and
  their options.
- [Plots of the estimate and the
  checks](https://dorleventer.github.io/didisc/articles/didisc-plots.md):
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit and
  on each test, the switchers picture, and saving a figure.
- [Bandwidth rules and sampling
  schemes](https://dorleventer.github.io/didisc/articles/didisc-options.md):
  `bwselect`, fixed bandwidths, and the standard error under each
  scheme.
- [How didisc() computes the
  estimate](https://dorleventer.github.io/didisc/articles/didisc-how-it-works.md):
  the estimate and the tests rebuilt by hand, for readers who want to
  check the computation.

## References

`citation("didisc")` returns the reference for the package, the paper
below, with a BibTeX entry.

Leventer, D. and D. Nevo (2024). Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data. arXiv:2408.05847.
<https://arxiv.org/abs/2408.05847>

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295–2326.
