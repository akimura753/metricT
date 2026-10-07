* Encoding: UTF-8.
* ======================================================================.
* Metric-T for IBM SPSS Statistics (syntax only; no R or Python needed).
*   T = DC(wPLI) - DC(MSC), in percentage points.
*   DC = percentage of features in which the mean of group 1 exceeds
*        the mean of group 2 (direction consistency).
* Reference: Kimura A (2026). Neuroscience Informatics 6:100286.
*            doi:10.1016/j.neuri.2026.100286
* Companion to the R package metricT (https://github.com/akimura753/metricT);
* same definitions as metric_t(..., exact = TRUE) for two groups.
*
* Run this file once (Run > All). It only defines the macro !METRICT,
* which can then be called in two ways.
*
* ----------------------------------------------------------------------.
* A. Any data set that is open in SPSS
*    (.sav, Excel, CSV ... opened through File > Open > Data or by syntax).
*
*   !METRICT GROUP = grp  G1 = 1  WPLI = (w1 TO w10)  MSC = (m1 TO m10)  NBANDS = 1.
*
*   Data layout: one row per subject.
*   GROUP    numeric group variable
*   G1       value of GROUP that defines group 1 (all other cases: group 2)
*   WPLI     wPLI variables, in parentheses
*   MSC      MSC variables, in parentheses (same number, same order:
*            feature j of WPLI pairs with feature j of MSC)
*   NBANDS   number of bands contained in the lists (default 1). With several
*            bands, list band 1 first, then band 2, ...; every band must have
*            the same number of features.
*
* ----------------------------------------------------------------------.
* B. A feature table written by the R package (mt_write_features_csv).
*    The macro reads the file and runs the analysis.
*
*   !METRICT FILE = "C:/data/features.csv"  NBANDS = 4  NFEAT = 10  G1 = "CB".
*
*   File layout: a header line, then one row per subject with the columns
*   subject, group, and for each band NFEAT wPLI columns followed by NFEAT
*   MSC columns.
*   FILE     the CSV file, in quotes
*   NBANDS   number of bands in the file
*   NFEAT    number of features (channel pairs) per band
*   G1       label of group 1 in the group column, in quotes
*            (all other subjects: group 2)
*   The variables are named b1w1 ... (band 1, wPLI, feature 1) and b1m1 ...;
*   the numeric group variable is mt_grp (1 = group 1, 2 = group 2).
*
* ----------------------------------------------------------------------.
* Optional in both forms
*   MAXEXACT largest number of label assignments enumerated exactly
*            (default 200000); above it a Monte Carlo test is used
*   NPERM    number of random assignments for the Monte Carlo test (10000)
*   SEED     random number seed for the Monte Carlo test (42)
*
* OUTPUT, one row per band
*   DC_wPLI DC_MSC T   direction consistency (%) and Metric-T (points)
*   Reversal           1 if the two DC values lie on opposite sides of 50
*   p                  two-sided p-value over the label assignments
*   p_min              smallest p-value attainable with these group sizes
*   p_Holm p_WY        adjusted over bands (Holm; Westfall-Young max-T)
* Cases with a missing value in any listed variable are omitted.
* ======================================================================.

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
