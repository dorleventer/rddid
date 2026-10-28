# print and summary of a fit

    Code
      print(fit)
    Output
      Diff-in-disc estimate of the ATT in period 3
        Comparison periods: 1, 2   (constant confounding trend)
        Comparison weights: 1: 0.500, 2: 0.500
        Sampling scheme: panel, no unit changes side of the cutoff (detected)
        Bandwidth: common h = 0.2694 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
        Pilot bandwidth b (period = value): 1 = 0.4136, 2 = 0.3901, 3 = 0.3984
      
                                   Estimate  Std. err.       z  p-value   95% CI
        Conventional                 1.0922     0.1259    8.67   <0.001   [0.8454, 1.3391]
        Robust (bias-corrected)      1.1411     0.1488    7.67   <0.001   [0.8494, 1.4328]
      
        summary() shows the per-period fits and the s.e. under every sampling scheme.

---

    Code
      print(summary(fit))
    Output
      Diff-in-disc estimate of the ATT in period 3
        Comparison periods: 1, 2   (constant confounding trend)
        Comparison weights: 1: 0.500, 2: 0.500
        Sampling scheme: panel, no unit changes side of the cutoff (detected)
        Bandwidth: common h = 0.2694 (rule "joint": one bandwidth, chosen for the diff-in-disc estimate)
        Pilot bandwidth b (period = value): 1 = 0.4136, 2 = 0.3901, 3 = 0.3984
      
                                   Estimate  Std. err.       z  p-value   95% CI
        Conventional                 1.0922     0.1259    8.67   <0.001   [0.8454, 1.3391]
        Robust (bias-corrected)      1.1411     0.1488    7.67   <0.001   [0.8494, 1.4328]
      
        Per-period local-linear fits (Conventional = sum of coef x jump; Robust = sum of coef x jump (bc)):
        period   role           coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
        1        comparison   -0.500   1000   0.2694   0.4136     0.5129    0.1683     0.4642    0.1942
        2        comparison   -0.500   1000   0.2694   0.3901     0.6306    0.1679     0.6249    0.2004
        3        RD            1.000   1000   0.2694   0.3984     1.6640    0.1635     1.6856    0.1924
      
        Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
          cross-section 0.2377   panel, no unit changes side 0.1488   panel, some change side 0.1488

---

    Code
      print(coef(fit))
    Output
      Conventional       Robust 
          1.092242     1.141077 

---

    Code
      print(confint(fit))
    Output
                       2.5 %   97.5 %
      Conventional 0.8454070 1.339077
      Robust       0.8493838 1.432769

# print of the four tests

    Code
      print(rd_typecont(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3))
    Output
      Test of type continuity  [rd_typecont()]
        H0: the share of each type jumps by zero at the cutoff, in every period
        Sampling scheme: panel, some units change side of the cutoff (detected)
        Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
        Type in a period: a unit's sides of the cutoff in the other periods, in time order
          ("+" at or above, "-" below).
      
        Joint sup-t (t critical values) max |t| = 1.232 over 9 jumps,  p = #.###  (largest: period:type 1:+-)
          Period 1: sup-t (t critical values) max |t| = 1.232 over 3 jumps,  p = #.###
          Period 2: sup-t (t critical values) max |t| = 0.582 over 3 jumps,  p = #.###
          Period 3: sup-t (t critical values) max |t| = 1.048 over 3 jumps,  p = #.###
        All joint statistics: x$joint_tests

---

    Code
      print(rd_compstable(didisc_sim_pv, x = "R", time = "year", id = "id", t_rd = 3))
    Output
      Test of composition stability  [rd_compstable()]
        H0: the share of each type among the units just above the cutoff is the same in the RD period and in each comparison period
        Sampling scheme: panel, some units change side of the cutoff (detected)
        RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
      
        Pair 3::1: sup-t (t critical values) max |t| = 3.086 over 3 jumps,  p = #.###
          n above the cutoff: 490 (RD period), 508 (comparison), 430 in both
        Pair 3::2: sup-t (t critical values) max |t| = 2.059 over 3 jumps,  p = #.###
          n above the cutoff: 490 (RD period), 497 (comparison), 418 in both
      
        Joint over pairs: sup-t (t critical values) max |t| = 3.086 over 6 jumps,  p = #.###  (largest: pair:type 3::1:-+)
        All joint statistics: x$joint_tests

---

    Code
      print(rd_homog(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3))
    Output
      Test of homogeneous confounding  [rd_homog()]
        H0: in each comparison period the confounding jump is the same for every type
        Sampling scheme: panel, some units change side of the cutoff (detected)
        Comparison periods: 1, 2   Bandwidth: h = 0.2104 in every cell (CCT on the pooled comparison periods)
      
        Type: a unit's side of the cutoff in RD period 3 ("+" at or above, "-" below).
        Wald chi-squared(2) = 0.247,  p = 0.884
      
        Per-cell local-linear jumps (comparison periods):
          Period     Type               jump       s.e.       n
          1          -                0.5347     0.2640     510  (reference)
          1          +                0.7070     0.4047     490
          2          -                1.0987     0.4183     510  (reference)
          2          +                0.8699     0.5244     490

---

    Code
      print(rd_trendcell(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id",
        t_rd = 3))
    Output
      Test of constant within-type confounding  [rd_trendcell()]
        H0: within each type, the confounding jump is the same in every comparison period
        Sampling scheme: panel, some units change side of the cutoff (detected)
        Comparison periods: 1, 2   Trend: constant
        Bandwidth: h = 0.2104 in every cell (CCT on the pooled comparison periods)
      
        Type: a unit's side of the cutoff in RD period 3 ("+" at or above, "-" below).
        Wald chi-squared(2) = 1.445,  p = 0.485
      
        Per-cell local-linear jumps (comparison periods):
          Type       Period             jump       s.e.       n
          +          1                0.7070     0.4047     490  (reference)
          +          2                0.8699     0.5244     490
          -          1                0.5347     0.2640     510  (reference)
          -          2                1.0987     0.4183     510

# the main error and message texts

    Code
      didisc(didisc_sim, y = "Yy", x = "R", time = "year", id = "id", t_rd = 3)
    Condition
      Error in `didisc()`:
      ! not a column of `data`: `y` = "Yy". The columns are id, year, R, V, W, Y.

---

    Code
      didisc(didisc_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
        comparisons = c(1, 7))
    Condition
      Error in `didisc()`:
      ! comparison period(s) 7 not in `year` (its periods: 1, 2, 3).

---

    Code
      rd_typecont(didisc_sim, x = "R", time = "year", id = "id")
    Condition
      Error in `rd_typecont()`:
      ! no unit changes side of the cutoff between the periods, so every unit has the same type in every period and the test is not defined. The tests are for a running variable that varies over time (see ?rd_typecont).

---

    Code
      didisc(didisc_sim, y = "Y", x = "R", time = "year", t_rd = 3, bwselect = "cct")$
        scheme
    Message
      didisc(): no `id` given, so every row is treated as a different unit (repeated cross-section standard errors).
    Output
      [1] "cs"

# print of rd_checks()

    Code
      print(rd_checks(didisc_sim_pv, y = "Y", x = "R", time = "year", id = "id",
        t_rd = 3))
    Output
      Tests of the identification assumptions  [rd_checks()]
      
        Assumption                          Statistic        p   Statistic type             Bandwidth
        type continuity                         1.232    #.###   sup-t, t critical values   CCT per cell
        composition stability                   3.086    #.###   sup-t, t critical values   CCT per cell
        homogeneous confounding                 0.247    0.884   Wald chi-squared(2)        pooled CCT, h = 0.21
        constant within-type confounding        1.445    0.485   Wald chi-squared(2)        pooled CCT, h = 0.21
      
        Reading at level 0.95, rejecting when p < 0.05 (not rejecting does not establish an assumption):
          type continuity                    not rejected
          composition stability              rejected     } identification needs
          homogeneous confounding            not rejected } only one of these two
          constant within-type confounding   not rejected (checked in the comparison periods)
        Each test's details: x$typecont, x$compstable, x$homog, x$trendcell.

