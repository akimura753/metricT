# metricT

R implementation of **Metric-T**, a permutation-based diagnostic that flags when the
reported direction of a group difference in EEG functional connectivity depends on the
choice of metric (wPLI vs. MSC).

Method: Kimura A (2026). *Metric-T: A permutation-based diagnostic for directional
fragility in EEG functional connectivity analysis.* Neuroscience Informatics 6:100286.
doi:10.1016/j.neuri.2026.100286

## Install

```r
install.packages(c("remotes", "shiny"))        # shiny is only needed for the GUI
remotes::install_github("akimura753/metricT")
```

or, from the release file: `install.packages("metricT_0.1.0.tar.gz", repos = NULL, type = "source")`.

## Graphical interface

```r
metricT::run_metricT_app()
```

Data sources: (1) built-in simulated demo; (2) a feature file, either `features.npz`
from the published Python pipeline (read natively, no Python required) or a CSV table;
(3) one CSV of preprocessed EEG per subject plus a participants table.

## Script

```r
library(metricT)
dat <- mt_read_npz("metricT_out/features.npz")
res <- metric_t(dat$features, dat$group,
                comparisons = list(c("AD","CN"), c("FTD","CN"), c("AD","FTD")),
                control = "CN", n_perm = 10000, seed = 42)
res          # table: DC(wPLI), DC(MSC), T, reversal flag, p, adjusted p
plot(res)    # DC bars with the 50% line
```

## Exact enumeration for small samples

With few subjects the reference set of a permutation test is small enough to be
enumerated completely, so no Monte Carlo error remains and the smallest attainable
p-value is known for the design at hand. Three functions share this idea.

```r
## 1. Exact two-group test: every assignment of the group labels
res <- metric_t(dat$features, dat$group, comparisons = list(c("A", "B")),
                exact = TRUE, wy_family = "all")

## 2. The complete labelling space: Metric-T for all 2^n - 2 two-group labellings
sp <- mt_label_space(dat$features)
sp                                   # SD of T, effective number of features, reversal share
mt_locate(sp, dat$group == "A")      # the observed labelling: T, exact p, p_min
mt_locate(sp, session == "day1")     # any other split: a batch, a cut-point, an exclusion

## 3. A continuous covariate, permuted only within strata (e.g. recording sessions)
tr <- metric_t_trend(dat$features, years, strata = session)
tr                                   # T, exact p, family-wise adjusted p
```

* `mt_label_space()` places the observed labelling, alternative group definitions and
  nuisance splits (such as recording sessions) on one distribution. Its summary reports
  the effective number of independent features, `2500 / Var(DC)`, and how often DC(wPLI)
  and DC(MSC) fall on opposite sides of 50% under arbitrary labelling.
* `metric_t_trend()` defines direction consistency as the percentage of features
  positively associated with the covariate. Stratum contributions are additive, so the
  complete reference set (the product of the stratum factorials) is obtained from the
  sum of the stratum factorials evaluations.
* A runnable example on simulated data is in `inst/examples/label_space.R`.

## Notes

* Core functions use base R only (FFT via `stats::mvfft`).
* Preprocessing (filtering, artifact removal, re-referencing) is done upstream;
  the package starts from epoched signals or from wPLI/MSC feature matrices.
* Ties between permuted and observed |T| are resolved exactly (integer counts).
* `p_min` reports the smallest p-value attainable in a condition. DC keeps only the
  sign of each mean difference, so when channel pairs vary in unison across subjects
  even |T| = 100 pp may not reach significance; a large p is then not evidence that
  the direction is stable.
