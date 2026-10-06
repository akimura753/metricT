# metricT 0.2.0

* `metric_t(exact = TRUE)`: complete enumeration of two-group label assignments;
  p-values, attainable minimum p and Westfall-Young adjustment are computed on the
  complete set.
* `mt_label_space()` and `mt_locate()`: Metric-T for every two-group labelling of
  the subjects (2^n - 2), with the effective number of independent features and the
  share of labellings showing a nominal reversal.
* `metric_t_trend()`: Metric-T for a continuous covariate, with exact permutation
  restricted to strata (for example recording sessions) and joint adjustment over
  several covariates.
* `mt_clean_epochs()`: linear detrending and peak-to-peak rejection of epochs.
* `mt_labelset()` / `mt_labelings()`: enumeration of labellings for ambiguous
  group membership.

# metricT 0.1.0

* First release: wPLI/MSC features, direction consistency, Metric-T with aligned
  permutations, Westfall-Young max-T, 'shiny' interface, readers for the feature
  files of the published Python pipeline.
