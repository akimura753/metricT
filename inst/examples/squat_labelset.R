# Simulation: Metric-T over ALL labellings of the intermediate cases (2^m).
# Usage: Rscript squat_labelset.R <share of intermediate cases> <replicates> <N>
suppressMessages(library(metricT))
a <- commandArgs(trailingOnly = TRUE); q <- as.numeric(a[1]); R <- as.integer(a[2]); N <- as.integer(a[3])
set.seed(20261006 + round(1000 * q) + N)
one <- function() {
  n_int <- round(q * N); n_fi <- n_int %/% 2; n_sa <- n_int - n_fi
  n_con <- N - n_int; n_fa <- n_con %/% 2; n_si <- n_con - n_fa
  cls <- rep(c("FA", "FI", "SA", "SI"), c(n_fa, n_fi, n_sa, n_si))
  zs <- ifelse(cls %in% c("FA", "FI"), 1, -1) + rnorm(N, 0, 0.4)
  za <- ifelse(cls %in% c("FA", "SA"), 1, -1) + rnorm(N, 0, 0.4)
  s <- mt_simulate_motion(zs, za, seed = NULL)
  fast <- s$covariates$cadence_hz > stats::median(s$covariates$cadence_hz)
  acc <- s$covariates$period_cv < stats::median(s$covariates$period_cv)
  lab <- ifelse(fast & acc, "good", ifelse(!fast & !acc, "poor", NA))
  if (sum(lab == "good", na.rm = TRUE) < 2 || sum(lab == "poor", na.rm = TRUE) < 2) return(NULL)
  L <- mt_labelset(s$features["movement"], lab, seed = sample.int(1e6, 1))
  amb <- L$ambiguous
  inf <- if (length(amb)) L$influence[, 1] else numeric(0)
  cbind(L$summary, exact = L$exact, any_reversal = any(L$reversal, na.rm = TRUE),
        width = L$summary$T_max - L$summary$T_min,
        infl_fast_irregular = if (length(amb)) mean(inf[fast[amb]]) else NA,
        infl_slow_precise = if (length(amb)) mean(inf[!fast[amb]]) else NA)
}
res <- do.call(rbind, lapply(seq_len(R), function(i) one()))
saveRDS(res, sprintf("set_q%02d_N%d.rds", round(100 * q), N))
f <- function(x) round(mean(x, na.rm = TRUE), 1)
cat(sprintf("q=%.1f N=%d reps=%d | m=%.1f labellings(median)=%s exact=%.0f%% | Tmin %.1f  Tmed %.1f  Tmax %.1f  width %.1f  E[T] %.1f | reversal among labellings %.2f%% | replicates with any reversal %.0f%% | influence FI %.1f  SA %.1f\n",
  q, N, nrow(res), mean(res$n_ambiguous), format(stats::median(res$n_labellings), big.mark = ","), 100 * mean(res$exact),
  f(res$T_min), f(res$T_median), f(res$T_max), f(res$width), f(res$T_expected), mean(res$reversal_pct),
  100 * mean(res$any_reversal), f(res$infl_fast_irregular), f(res$infl_slow_precise)))
