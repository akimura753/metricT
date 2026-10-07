# Metric-T for IBM SPSS Statistics

SPSS syntax for **Metric-T**, T = DC(wPLI) − DC(MSC), with the exact two-group
permutation test. It uses the SPSS matrix language only: no R, Python or
extension is needed. It is a companion to the R package in this repository and
uses the same definitions as `metric_t(..., exact = TRUE)`.

| File | Content |
|---|---|
| `metricT.sps` | Defines the macro `!METRICT`. Run it once per session. |
| `metricT_example.sps` | Self-contained example with simulated data. Open it and choose *Run > All*. |
| `METRICT.spd` | Dialog box: choose the variables with the mouse, no syntax to write. |

## Dialog box (no syntax)

Install `METRICT.spd` once: *Extensions > Utilities > Install Custom Dialog
(Compatibility mode)*, or double-click the file. The dialog then appears under
*Analyze > Correlate > Metric-T (wPLI vs MSC)*. Open a data set with one row
per subject, move the numeric group variable, the wPLI variables and the MSC
variables into the three boxes, give the value of group 1 and the number of
bands, and press *OK*. The dialog contains the whole computation; `metricT.sps`
is not needed for it. Labels are available in English and Japanese.

## Use with syntax

Run `metricT.sps` (*Run > All*); it only defines the macro and prints nothing.
The file can be run again at any time. Then call the macro in one of two ways,
for example by adding the call as the last line of the file.

**From a feature table written by the R package** (`mt_write_features_csv()`).
The macro reads the file itself, so nothing needs to be open:

```
!METRICT FILE = "C:/data/features.csv"  NBANDS = 4  NFEAT = 10  G1 = "CB".
```

`NFEAT` is the number of features (channel pairs) per band, and `G1` is the
label of group 1 in the `group` column; all other subjects form group 2.

**From any data set that is open in SPSS** (`.sav`, Excel, CSV, ...), one row per
subject, with a numeric group variable and the wPLI and MSC features as variables:

```
!METRICT GROUP = grp  G1 = 1  WPLI = (w1 TO w10)  MSC = (m1 TO m10)  NBANDS = 1.
```

`G1` is the value of the group variable that defines group 1; all other cases
form group 2. The two variable lists must have the same length and order. With
several bands, list band 1 first, then band 2, and so on, and give `NBANDS`.
The variable names must be those of the open data set; the dialog box avoids
typing them.

## Output

One row per band: DC(wPLI), DC(MSC), T, a reversal flag, the two-sided p-value
over the label assignments, the smallest attainable p-value, and p-values
adjusted over bands (Holm and Westfall–Young max-T). All assignments are
enumerated when there are at most 200,000 (`MAXEXACT`); otherwise a Monte Carlo
test with `NPERM` random assignments is used.

## Scope and checks

The syntax covers the two-group test. The complete labelling space
(`mt_label_space`), the stratified test for a continuous covariate
(`metric_t_trend`) and the computation of wPLI and MSC from EEG signals are
available in the R package only.

Results were compared with the R package (metricT 0.2.0) on simulated data and
on the feature table of the reference paper: DC, T, p, the attainable minimum
and the adjusted p-values agreed. Both call forms of the macro and the dialog box were
run in IBM SPSS Statistics 27 (version 1.2); the syntax also runs in GNU PSPP 2.0.0.

Reference: Kimura A (2026). Metric-T: a permutation-based diagnostic for
directional fragility in EEG functional connectivity analysis. *Neuroscience
Informatics* 6:100286. doi:10.1016/j.neuri.2026.100286
