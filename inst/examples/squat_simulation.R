# Simulation: Metric-T under eight labelling rules for a paced squat task.
# Usage: Rscript squat_simulation.R <share of intermediate cases> <replicates> <N> [n_perm]
suppressMessages(library(metricT))
a <- commandArgs(trailingOnly = TRUE)
q <- as.numeric(a[1]); R <- as.integer(a[2]); N <- as.integer(a[3]); n_perm <- if (length(a) > 3) as.integer(a[4]) else 0L
set.seed(20261005 + round(1000 * q) + N)
rules <- c("speed", "accuracy", "strict_both", "lenient_either", "exclude_intermediate",
           "rater_speed", "rater_form", "outcome_classifier")
gp <- function(x) ifelse(x, "good", "poor")
one <- function() {
  n_int <- round(q * N); n_fi <- n_int %/% 2; n_sa <- n_int - n_fi
  n_con <- N - n_int; n_fa <- n_con %/% 2; n_si <- n_con - n_fa
  cls <- rep(c("FA", "FI", "SA", "SI"), c(n_fa, n_fi, n_sa, n_si))
  zs <- ifelse(cls %in% c("FA", "FI"), 1, -1) + rnorm(N, 0, 0.4)
  za <- ifelse(cls %in% c("FA", "SA"), 1, -1) + rnorm(N, 0, 0.4)
  s <- mt_simulate_motion(zs, za, seed = NULL)
  cad <- s$covariates$cadence_hz; cv <- s$covariates$period_cv
  fast <- cad > stats::median(cad); acc <- cv < stats::median(cv)
  y <- 0.5 * zs + 0.5 * za + rnorm(N, 0, 0.7); Y <- as.integer(y > stats::median(y))
  fit <- suppressWarnings(stats::glm(Y ~ scale(cad) + scale(-log(cv)), family = stats::binomial))
  r1 <- 0.8 * zs + 0.2 * za + rnorm(N, 0, 0.5); r2 <- 0.2 * zs + 0.8 * za + rnorm(N, 0, 0.5)
  lab <- list(speed = gp(fast), accuracy = gp(acc), strict_both = gp(fast & acc),
              lenient_either = gp(fast | acc),
              exclude_intermediate = ifelse(fast & acc, "good", ifelse(!fast & !acc, "poor", NA)),
              rater_speed = gp(r1 > stats::median(r1)), rater_form = gp(r2 > stats::median(r2)),
              outcome_classifier = gp(stats::fitted(fit) > 0.5))
  t <- mt_labelings(s$features["movement"], lab, n_perm = n_perm, seed = sample.int(1e6, 1))
  t$disagree_speed_accuracy <- mean(fast != acc)
  t$disagree_raters <- mean(lab$rater_speed != lab$rater_form)
  t
}
res <- do.call(rbind, lapply(seq_len(R), function(i) cbind(rep = i, one())))
saveRDS(res, sprintf("res_q%02d_N%d_p%d.rds", round(100 * q), N, n_perm))
sm <- do.call(rbind, lapply(rules, function(r) { x <- res[res$rule == r, ]
  data.frame(rule = r, n_good = mean(x$n_good), DC_wPLI = mean(x$DC_wPLI, na.rm = TRUE), DC_MSC = mean(x$DC_MSC, na.rm = TRUE),
             T_mean = mean(x$T, na.rm = TRUE), T_sd = stats::sd(x$T, na.rm = TRUE), reversal_pct = 100 * mean(x$reversal, na.rm = TRUE),
             na = sum(is.na(x$T)),
             p05 = if (n_perm > 0) 100 * mean(x$p_raw <= 0.05, na.rm = TRUE) else NA,
             pmin_med = if (n_perm > 0) stats::median(x$p_min, na.rm = TRUE) else NA) }))
cat(sprintf("\nq = %.1f, N = %d, R = %d; speed/accuracy rules disagree on %.0f%% of subjects; raters on %.0f%%\n", q, N, R,
            100 * mean(res$disagree_speed_accuracy), 100 * mean(res$disagree_raters)))
print(format(sm, digits = 3, nsmall = 1), row.names = FALSE)
