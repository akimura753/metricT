## Complete labelling space and stratified trend test on simulated data.
## Runs in a few seconds; no external files are needed.
library(metricT)

sim <- mt_simulate(n = c(A = 6, B = 7), n_ch = 5, n_epochs = 12, seed = 7)

## 1. Exact two-group test: all choose(13, 6) = 1716 label assignments
res <- metric_t(sim$features, sim$group, comparisons = list(c("A", "B")),
                exact = TRUE, wy_family = "all")
print(res)

## 2. The complete labelling space: 2^13 - 2 = 8190 labellings
sp <- mt_label_space(sim$features)
print(sp)

## Where does the observed labelling sit, and where would other splits sit?
mt_locate(sp, sim$group == "A")          # observed groups
session <- rep(c("s1", "s2", "s3"), c(4, 4, 5))
mt_locate(sp, session == "s1")           # a recording session treated as a "group"

## 3. A continuous covariate, permuted within recording sessions only
set.seed(11)
covariate <- round(runif(13, 0, 40))
tr <- metric_t_trend(sim$features, covariate, strata = session)
print(tr)                                # 4! * 4! * 5! = 69120 arrangements, exact
