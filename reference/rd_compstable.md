# Test composition stability across periods (assumption `ass:comp-stable`)

Tests the composition-stability assumption (`ass:comp-stable`, Section
4.4 of Leventer and Nevo): \$\$\pi\_{t\_{\mathrm{RD}},(+)}(\mathbf{u},
b) = \pi\_{t_0,(+)}(\mathbf{u}, b) \quad \forall\\(\mathbf{u}, b),\$\$
where \\\mathbf{u} = \mathbf{v}\_{-\\t_0, t\_{\mathrm{RD}}\\}\\ are the
sides of the OTHER periods (shared between the two confounding objects)
and \\b \in \\0,1\\\\ is the "partner" side. This is a cross-period
covariate-continuity statement: the above-cutoff \\(\mathbf{u},b)\\ type
mix must be the same at the RD period and at each comparison period.

## Usage

``` r
rd_compstable(
  data,
  x,
  time,
  id,
  t_rd,
  comparisons = NULL,
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  q = NULL,
  S = 499L,
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  bc = TRUE,
  ...
)
```

## Arguments

- data:

  A long data frame, one row per unit-period (balanced or unbalanced
  panel).

- x:

  Column name (string) for the running variable.

- time:

  Column name (string) for the period indicator.

- id:

  Column name (string) for the unit identifier.

- t_rd:

  Value of `time` identifying the RD period.

- comparisons:

  Values of `time` to use as comparison periods. If `NULL` (default),
  all periods except `t_rd` are used.

- c:

  Cutoff for the running variable (default 0).

- h:

  Bandwidth. If `NULL` (default), the bandwidth is determined by
  `bwselect`; an explicit numeric value overrides `bwselect` and is used
  directly.

- bwselect:

  Bandwidth selection rule when `h = NULL`: `"cct"` (default) computes a
  per-cell CCT MSE-optimal bandwidth via
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  for each type indicator RD in the reflected space; `"rot"` uses \\0.5
  \times \mathrm{IQR}(x)\\ as a rule-of-thumb applied to the full
  sample. Ignored when `h` is supplied explicitly.

- q:

  Number of observations nearest the artificial cutoff on each side for
  the Canay-Kamat permutation test. `NULL` (default) selects `q` per
  \\(t\_{\mathrm{RD}}, t_0)\\ pair by the Canay & Kamat (2018) rule of
  thumb (see
  [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md));
  this is the recommended choice, since a fixed `q` over-rejects in
  finite samples when the type distribution varies steeply in the
  running variable at the cutoff. Supply an integer to force a fixed `q`
  on every pair. The per-pair `q` actually used is returned in
  `meta$q_used`.

- S:

  Number of permutation replications (default 499).

- kernel:

  Kernel for the local-linear RD: `"triangular"` (default),
  `"epanechnikov"`, or `"uniform"`.

- scheme:

  Covariance scheme for the Wald test:

  `"auto"`

  :   Detects whether any unit appears on both sides of the artificial
      cutoff (i.e., above the true cutoff in both periods). If yes, uses
      `"pv"` (time-varying panel); otherwise `"cs"`.

  `"cs"`

  :   Treats the two sides as independent.

  `"pc"`

  :   Includes same-side cross-period covariance only.

  `"pv"`

  :   Full panel with time-varying running variable: includes same-side
      minus opposite-side cross-period covariance.

- bc:

  Use robust bias-corrected jumps and variances in the LL-Wald
  (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the
  test with the bias-corrected
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  estimator; `FALSE` uses the conventional local-linear jumps and
  variances. Does not affect the Canay-Kamat permutation test.

- ...:

  Currently unused.

## Value

An object of class `"rd_compstable"`, a named list with:

- `pairs`:

  A list, one element per \\(t\_{\mathrm{RD}}, t_0)\\ pair (named
  `"trd::t0"`), each containing:

  `ll_wald`

  :   list with `stat`, `df`, `p`.

  `ck_perm`

  :   list with `stat` (observed sum of \|mean diffs\|) and `p`.

  `type_values`

  :   Character vector of \\(\mathbf{u},b)\\ type labels present in this
      pair.

  `scheme`

  :   Scheme actually used.

  `q`

  :   Number of nearest observations per side actually used.

  `n_trd`

  :   Number of above-cutoff units from \\t_RD\\.

  `n_t0`

  :   Number of above-cutoff units from \\t_0\\.

  `n_both`

  :   Number of units above the cutoff in both periods.

- `joint`:

  Joint result over all pairs (stacked Wald + minimum-p permutation
  envelope):

  `ll_wald`

  :   list with `stat`, `df`, `p` (sum of the per-pair statistics and
      df, which assumes independent pairs; the pairs share the RD-period
      above-cutoff group, so treat the joint as approximate — the
      paper's test is per pair).

  `ck_perm`

  :   list with `stat` (sum of per-pair stats) and `p`.

- `meta`:

  list with `t_rd`, `comparisons`, `h` (NA when `bwselect = "cct"`),
  `bwselect`, `q` (`"rot"` when the rule of thumb is used, otherwise the
  integer supplied), `q_used` (per-pair `q` actually used), `S`, `c`,
  `bc`.

## Details

### Reflection construction

For each RD-period / comparison-period pair \\(t\_{\mathrm{RD}}, t_0)\\,
take the above-cutoff units of each period. For \\t_0\\-above units flip
the centred running variable: \\x' = -(R\_{i,t_0} - c)\\, placing them
just BELOW an artificial cutoff at 0. For \\t\_{\mathrm{RD}}\\-above
units set \\x' = R\_{i,t\_{\mathrm{RD}}} - c\\, keeping them just ABOVE
0. Stack the two groups into one artificial cross-section. At the
artificial cutoff the left/right limits of the \\(\mathbf{u},b)\\ type
share are then \\\pi\_{t_0,(+)}(\mathbf{u},b)\\ and
\\\pi\_{t\_{\mathrm{RD}},(+)}(\mathbf{u},b)\\, so a jump at 0 equals the
composition difference. The same two tests from the type-continuity
assessment (`ass:type-cont`,
[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md))
then apply:

1.  **LL-Wald** (necessary AND sufficient): local-linear RD of each
    \\(\mathbf{u},b)\\ indicator on the reflected running variable;
    joint Wald that all jumps are zero, dropping one reference type (the
    type shares sum to 1, so the full set of jumps is rank-deficient).
    With binary types this is the single share-jump test of the paper's
    Section 4.4. This is the paper's test.

2.  **Canay-Kamat permutation** (necessary AND sufficient): approximate
    sign randomisation test comparing the two sides of the artificial
    cutoff, permuted at the **unit** level (see "Unit-level wrinkle"
    below).

Both tests are **necessary AND sufficient** for `ass:comp-stable`.
Section 4.4 of the paper states the LL-Wald; the permutation test is
reported in the paper's Section 6 validation table (Canay and Kamat
2018, rule-of-thumb `q`). See `dev/tests_map.md`.

### Unit-level wrinkle

A unit that is above the cutoff in BOTH periods \\t\_{\mathrm{RD}}\\ and
\\t_0\\ appears on BOTH sides of the artificial cutoff (as a
\\t\_{\mathrm{RD}}\\-above observation above the artificial 0 and a
\\t_0\\-above observation below it). This has two consequences handled
by the function:

\(i\) **Wald test**: the covariance matrix between the left-side and
right-side intercept estimates must include the id-matched cross-side
covariance term (the two g-vectors can share unit ids). The function
computes \\(\text{cov}\_{++} + \text{cov}\_{--} - \text{cov}\_{+-} -
\text{cov}\_{-+})\\ — the same formula as the PV scheme in the main
estimator — rather than assuming the two sides are independent.

\(ii\) **Permutation test**: permutes at the unit level. Each unit
contributes its observations (possibly \>1) as a block; the side label
is permuted across units, not across individual rows. This follows Amro
and Pauly (2017) and Derrick et al. (2022) (see References).

## Note

**Necessary and sufficient status:** Both the LL-Wald and the
Canay-Kamat permutation test are necessary AND sufficient for
`ass:comp-stable` (composition stability). See Leventer and Nevo for the
proof.

## References

Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
Designs." Working paper.

Amro, L. and Pauly, M. (2017). Permuting longitudinal data in spite of
the dependencies. *Journal of Statistical Computation and Simulation*,
87(15), 3033-3044.

Derrick, B., Broad, A., Ruck, A., and White, P. (2022). The impact of
repeated measures on the permutation test. *Journal of Applied
Quantitative Methods*, 17(1).

Canay, I. A. and Kamat, V. (2018). Approximate permutation tests and
induced order statistics in the regression discontinuity design. *Review
of Economic Studies*, 85(3), 1577-1608.

## See also

[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-period panel with no composition shift (null DGP).
set.seed(1)
n <- 500
eta <- rnorm(n)
dat <- data.frame(
  id   = rep(seq_len(n), 2),
  time = rep(1:2, each = n),
  R    = c(eta + rnorm(n), eta + rnorm(n))
)
rd_compstable(dat, x = "R", time = "time", id = "id", t_rd = 2,
              comparisons = 1, h = 0.5, S = 99)
} # }
```
