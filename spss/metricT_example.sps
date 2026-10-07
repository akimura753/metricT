* Encoding: UTF-8.
* ======================================================================.
* Metric-T in SPSS: self-contained example. No other file is read.
* Open this file in SPSS and choose Run > All.
*   Part 1  the macro definition (a copy of metricT.sps)
*   Part 2  SIMULATED data: 12 subjects (6 per group), 2 bands x 6 features.
*           The numbers are artificial and describe no real person.
*   Part 3  the analysis
* Expected result (924 label assignments, exact test):
*   Band 1  DC_wPLI 0.0  DC_MSC 100.0  T -100.0  p 0.121
*   Band 2  DC_wPLI 33.3  DC_MSC 83.3  T -50.0  p 0.418
* ======================================================================.

* ---------------- Part 1: macro definition ----------------.
DEFINE !METRICT (GROUP = !TOKENS(1)
                /G1 = !TOKENS(1)
                /WPLI = !ENCLOSE('(',')')
                /MSC = !ENCLOSE('(',')')
                /NBANDS = !DEFAULT(1) !TOKENS(1)
                /FILE = !TOKENS(1)
                /NFEAT = !TOKENS(1)
                /MAXEXACT = !DEFAULT(200000) !TOKENS(1)
                /NPERM = !DEFAULT(10000) !TOKENS(1)
                /SEED = !DEFAULT(42) !TOKENS(1))
PRESERVE.
SET MXLOOPS = 100000000.
SET SEED = !SEED.
!IF (!NFEAT !NE !NULL) !THEN
!LET !wl = !NULL
!LET !ml = !NULL
!DO !b = 1 !TO !NBANDS
!LET !wl = !CONCAT(!wl, " b", !b, "w1 TO b", !b, "w", !NFEAT)
!LET !ml = !CONCAT(!ml, " b", !b, "m1 TO b", !b, "m", !NFEAT)
!DOEND
GET DATA /TYPE = TXT
  /FILE = !FILE
  /DELIMITERS = ","
  /QUALIFIER = '"'
  /ARRANGEMENT = DELIMITED
  /FIRSTCASE = 2
  /VARIABLES = subject A40 group A40
!DO !b = 1 !TO !NBANDS
!DO !j = 1 !TO !NFEAT
    !CONCAT("b", !b, "w", !j) F12.9
!DOEND
!DO !j = 1 !TO !NFEAT
    !CONCAT("b", !b, "m", !j) F12.9
!DOEND
!DOEND
  .
COMPUTE mt_grp = 2.
IF (group EQ !G1) mt_grp = 1.
EXECUTE.
FREQUENCIES VARIABLES = group.
!IFEND
MATRIX.
PRINT /TITLE = "Metric-T macro, version 1.2".
!IF (!NFEAT !NE !NULL) !THEN
GET d /VARIABLES = mt_grp !wl !ml /MISSING = OMIT.
COMPUTE g1val = 1.
!ELSE
GET d /VARIABLES = !GROUP !WPLI !MSC /MISSING = OMIT.
COMPUTE g1val = !G1.
!IFEND
COMPUTE nb = !NBANDS.
COMPUTE maxex = !MAXEXACT.
COMPUTE nperm = !NPERM.
COMPUTE n = NROW(d).
COMPUTE p2 = NCOL(d) - 1.
COMPUTE p = p2 / 2.
COMPUTE k = p / nb.
COMPUTE gobs = (d(:,1) EQ g1val).
COMPUTE n1 = CSUM(gobs).
COMPUTE n2 = n - n1.
COMPUTE ok = 1.
DO IF (p NE TRUNC(p) OR k NE TRUNC(k)).
  COMPUTE ok = 0.
  PRINT /TITLE = "ERROR: the wPLI and MSC lists must have the same length, divisible by NBANDS.".
END IF.
DO IF (n1 LT 1 OR n2 LT 1).
  COMPUTE ok = 0.
  PRINT /TITLE = "ERROR: both groups need at least one case. Check GROUP and G1.".
END IF.
DO IF (ok EQ 1).
  COMPUTE x = d(:, 2:(p2 + 1)).
  COMPUTE tot = CSUM(x).
  COMPUTE ncomb = 1.
  LOOP j = 1 TO n1.
    COMPUTE ncomb = ncomb * (n - n1 + j) / j.
  END LOOP.
  COMPUTE ncomb = RND(ncomb).
  COMPUTE exact = (ncomb LE maxex).
  DO IF (exact EQ 1).
    COMPUTE nlab = ncomb.
  ELSE.
    COMPUTE nlab = nperm.
  END IF.
  COMPUTE s1 = T(gobs) * x.
  COMPUTE pos = ((s1 / n1) GT ((tot - s1) / n2)).
  COMPUTE cw = MAKE(1, nb, 0).
  COMPUTE cm = MAKE(1, nb, 0).
  LOOP b = 1 TO nb.
    COMPUTE cw(1, b) = RSUM(pos(1, ((b - 1) * k + 1):(b * k))).
    COMPUTE cm(1, b) = RSUM(pos(1, (p + (b - 1) * k + 1):(p + b * k))).
  END LOOP.
  COMPUTE dobs = ABS(cw - cm).
  COMPUTE dall = MAKE(nlab, nb, 0).
  COMPUTE c = MAKE(1, n1, 0).
  LOOP j = 1 TO n1.
    COMPUTE c(1, j) = j.
  END LOOP.
  COMPUTE onesr = MAKE(1, n, 1).
  COMPUTE onesc = MAKE(n, 1, 1).
  LOOP it = 1 TO nlab.
    DO IF (exact EQ 1).
      COMPUTE g = MAKE(n, 1, 0).
      LOOP j = 1 TO n1.
        COMPUTE g(c(1, j), 1) = 1.
      END LOOP.
    ELSE.
      COMPUTE u = UNIFORM(n, 1).
      COMPUTE g = (RSUM((u * onesr) GE (onesc * T(u))) LE n1).
    END IF.
    COMPUTE s1 = T(g) * x.
    COMPUTE pos = ((s1 / n1) GT ((tot - s1) / n2)).
    LOOP b = 1 TO nb.
      COMPUTE dall(it, b) = ABS(RSUM(pos(1, ((b - 1) * k + 1):(b * k)))
                              - RSUM(pos(1, (p + (b - 1) * k + 1):(p + b * k)))).
    END LOOP.
    DO IF (exact EQ 1).
      COMPUTE cp = 0.
      LOOP j = 1 TO n1.
        DO IF (c(1, j) LT n - n1 + j).
          COMPUTE cp = j.
        END IF.
      END LOOP.
      DO IF (cp GT 0).
        COMPUTE c(1, cp) = c(1, cp) + 1.
        DO IF (cp LT n1).
          LOOP j = cp + 1 TO n1.
            COMPUTE c(1, j) = c(1, j - 1) + 1.
          END LOOP.
        END IF.
      END IF.
    END IF.
  END LOOP.
  COMPUTE addone = 1 - exact.
  COMPUTE praw = MAKE(1, nb, 0).
  COMPUTE pmin = MAKE(1, nb, 0).
  LOOP b = 1 TO nb.
    COMPUTE praw(1, b) = (addone + CSUM(dall(:, b) GE dobs(1, b))) / (nlab + addone).
    COMPUTE cnt = CSUM(dall(:, b) GE k).
    DO IF (exact EQ 1 AND cnt LT 1).
      COMPUTE cnt = 1.
    END IF.
    COMPUTE pmin(1, b) = (addone + cnt) / (nlab + addone).
  END LOOP.
  COMPUTE pholm = MAKE(1, nb, 0).
  COMPUTE used = MAKE(1, nb, 0).
  COMPUTE runmax = 0.
  LOOP i = 1 TO nb.
    COMPUTE best = 2.
    COMPUTE bi = 0.
    LOOP j = 1 TO nb.
      DO IF (used(1, j) EQ 0 AND praw(1, j) LT best).
        COMPUTE best = praw(1, j).
        COMPUTE bi = j.
      END IF.
    END LOOP.
    COMPUTE used(1, bi) = 1.
    COMPUTE val = best * (nb - i + 1).
    DO IF (val GT 1).
      COMPUTE val = 1.
    END IF.
    DO IF (val GT runmax).
      COMPUTE runmax = val.
    END IF.
    COMPUTE pholm(1, bi) = runmax.
  END LOOP.
  COMPUTE ord = MAKE(1, nb, 0).
  COMPUTE used = MAKE(1, nb, 0).
  LOOP i = 1 TO nb.
    COMPUTE best = -1.
    COMPUTE bi = 0.
    LOOP j = 1 TO nb.
      DO IF (used(1, j) EQ 0 AND dobs(1, j) GT best).
        COMPUTE best = dobs(1, j).
        COMPUTE bi = j.
      END IF.
    END LOOP.
    COMPUTE ord(1, i) = bi.
    COMPUTE used(1, bi) = 1.
  END LOOP.
  COMPUTE adj = MAKE(1, nb, 0).
  COMPUTE cur = dall(:, ord(1, nb)).
  LOOP ii = 1 TO nb.
    COMPUTE i = nb - ii + 1.
    DO IF (i LT nb).
      COMPUTE nxt = dall(:, ord(1, i)).
      COMPUTE cur = (nxt GT cur) &* nxt + (nxt LE cur) &* cur.
    END IF.
    COMPUTE adj(1, i) = (addone + CSUM(cur GE dobs(1, ord(1, i)))) / (nlab + addone).
  END LOOP.
  COMPUTE pwy = MAKE(1, nb, 0).
  COMPUTE runmax = 0.
  LOOP i = 1 TO nb.
    DO IF (adj(1, i) GT runmax).
      COMPUTE runmax = adj(1, i).
    END IF.
    COMPUTE pwy(1, ord(1, i)) = runmax.
  END LOOP.
  COMPUTE dcw = 100 * cw / k.
  COMPUTE dcm = 100 * cm / k.
  COMPUTE tobs = dcw - dcm.
  COMPUTE rev = (((dcw - 50) &* (dcm - 50)) LT 0).
  COMPUTE bandno = MAKE(nb, 1, 0).
  LOOP b = 1 TO nb.
    COMPUTE bandno(b, 1) = b.
  END LOOP.
  COMPUTE info = {n1, n2, k, nb, nlab, exact}.
  PRINT info /TITLE = "Metric-T: design"
    /CLABELS = "n1", "n2", "Features", "Bands", "Labels", "Exact" /FORMAT = F9.0.
  COMPUTE res1 = {bandno, T(dcw), T(dcm), T(tobs), T(rev)}.
  PRINT res1 /TITLE = "Metric-T: T = DC(wPLI) - DC(MSC), percentage points (group 1 vs group 2)"
    /CLABELS = "Band", "DC_wPLI", "DC_MSC", "T", "Reversal" /FORMAT = F9.1.
  COMPUTE res2 = {bandno, T(praw), T(pmin), T(pholm), T(pwy)}.
  PRINT res2 /TITLE = "Metric-T: two-sided p-values over the label assignments"
    /CLABELS = "Band", "p", "p_min", "p_Holm", "p_WY" /FORMAT = F9.3.
  DO IF (exact EQ 1).
    PRINT /TITLE = "All label assignments were enumerated (exact test).".
  ELSE.
    PRINT /TITLE = "Monte Carlo test: p = (1 + count) / (NPERM + 1).".
  END IF.
END IF.
END MATRIX.
RESTORE.
!ENDDEFINE.

* ---------------- Part 2: simulated data ----------------.
DATA LIST FREE
  / subject (A8) grp
  b1w1 b1w2 b1w3 b1w4 b1w5 b1w6
  b2w1 b2w2 b2w3 b2w4 b2w5 b2w6
  b1m1 b1m2 b1m3 b1m4 b1m5 b1m6
  b2m1 b2m2 b2m3 b2m4 b2m5 b2m6.
BEGIN DATA
s01 1 0.3162 0.2686 0.2973 0.2382 0.2144 0.1911
      0.2595 0.3948 0.3850 0.3545 0.2567 0.3629
      0.3351 0.4624 0.3827 0.4471 0.4236 0.4345
      0.3784 0.4447 0.2868 0.3441 0.3834 0.4882
s02 1 0.3091 0.2215 0.2070 0.2268 0.2970 0.2314
      0.2797 0.2468 0.2409 0.3355 0.2914 0.1800
      0.2809 0.4670 0.5225 0.5221 0.4140 0.4635
      0.3626 0.4324 0.3786 0.4484 0.3410 0.4150
s03 1 0.3004 0.2228 0.2470 0.2742 0.2719 0.2331
      0.2591 0.3023 0.2522 0.3134 0.2350 0.2802
      0.5128 0.4504 0.4659 0.5361 0.4793 0.4212
      0.3659 0.4546 0.3475 0.5351 0.3403 0.3810
s04 1 0.2456 0.1245 0.1894 0.2259 0.1898 0.1687
      0.1649 0.2738 0.1842 0.3002 0.3563 0.2321
      0.4264 0.4963 0.4996 0.4954 0.4266 0.4983
      0.4299 0.4801 0.4159 0.3628 0.3706 0.5244
s05 1 0.3317 0.2287 0.3252 0.2678 0.2898 0.3458
      0.3847 0.3663 0.3375 0.2811 0.3817 0.4070
      0.4358 0.2986 0.5167 0.4045 0.4165 0.4907
      0.3163 0.3636 0.2612 0.4396 0.3636 0.3641
s06 1 0.1401 0.1993 0.1007 0.2303 0.1700 0.1862
      0.3236 0.2039 0.1677 0.2481 0.2072 0.2465
      0.4732 0.5709 0.4924 0.5171 0.5086 0.6007
      0.4486 0.4745 0.3731 0.3665 0.5235 0.4185
s07 2 0.3096 0.3663 0.3163 0.3699 0.3113 0.3392
      0.3243 0.2811 0.2745 0.2951 0.3595 0.2846
      0.4805 0.4421 0.3465 0.4573 0.3784 0.3866
      0.4601 0.4081 0.3454 0.3634 0.3901 0.3718
s08 2 0.3251 0.2199 0.2317 0.3414 0.3306 0.2255
      0.1872 0.3632 0.2712 0.2944 0.2962 0.3068
      0.3811 0.3618 0.4133 0.4613 0.3700 0.3719
      0.3750 0.4698 0.4218 0.3987 0.4315 0.4467
s09 2 0.2550 0.3148 0.3087 0.2204 0.2171 0.4290
      0.2561 0.3537 0.2455 0.2946 0.3594 0.2786
      0.3506 0.3690 0.4199 0.4333 0.4086 0.4462
      0.3209 0.4377 0.4393 0.4296 0.3383 0.4554
s10 2 0.2694 0.3395 0.3230 0.3617 0.3597 0.2965
      0.3247 0.2896 0.3143 0.3192 0.3290 0.2395
      0.3475 0.4105 0.4330 0.3930 0.3612 0.2924
      0.3382 0.3264 0.4637 0.4118 0.3877 0.4260
s11 2 0.2542 0.3074 0.3089 0.2449 0.2871 0.2380
      0.2481 0.2853 0.3596 0.2918 0.2136 0.2021
      0.4172 0.4202 0.3691 0.4375 0.4075 0.3559
      0.4336 0.4067 0.4966 0.3562 0.3474 0.3019
s12 2 0.2962 0.3366 0.3328 0.4333 0.4081 0.3014
      0.3025 0.2725 0.3995 0.3658 0.3028 0.3400
      0.4173 0.3575 0.3144 0.3301 0.4127 0.4513
      0.3413 0.3235 0.3984 0.3708 0.4107 0.3433
END DATA.

* ---------------- Part 3: analysis ----------------.
* Both bands in one call; p_Holm and p_WY are adjusted over the two bands.
!METRICT GROUP = grp G1 = 1
  WPLI = (b1w1 TO b1w6 b2w1 TO b2w6)
  MSC = (b1m1 TO b1m6 b2m1 TO b2m6)
  NBANDS = 2.

* Band 1 only.
!METRICT GROUP = grp G1 = 1 WPLI = (b1w1 TO b1w6) MSC = (b1m1 TO b1m6).
