#' Frequency bands for paced repetitive movement (Hz)
#'
#' \code{movement} covers the repetition frequency of a paced squat
#' (roughly 0.4 to 0.9 Hz); \code{harmonic} covers its second and third harmonics.
#' @return Named list of length-2 numeric vectors.
#' @examples
#' mt_motion_bands()
#' @export
mt_motion_bands <- function() list(movement = c(0.35, 1.0), harmonic = c(1.0, 2.2))

#' Simulate pose-estimated repetitive movement
#'
#' Generates joint-angle time series for a paced repetitive task such as the
#' squat, one recording per subject, from two latent qualities: speed and
#' accuracy (both standardised; 0 is average). The generative assumptions are
#' (i) faster performers move their joints more simultaneously (smaller
#' inter-joint phase lags) and with less joint-specific noise, and (ii) more
#' accurate performers repeat the same period and the same inter-joint timing.
#' Under these assumptions MSC rises mainly with speed, while wPLI rises with
#' accuracy and falls with speed. The assumptions are for studying the
#' behaviour of Metric-T, not established facts about squatting.
#'
#' @param z_speed,z_accuracy Latent speed and accuracy of each subject.
#' @param fs Frame rate (Hz).
#' @param duration Task duration (s).
#' @param f0 Repetition frequency of an average performer (Hz).
#' @param n_joints Number of joint-angle channels (2 to 6).
#' @param lag_unit Inter-joint phase lag step (rad) of an average performer.
#' @param lag_jitter SD (rad) of the slowly varying lag deviation of an average performer.
#' @param period_cv Repetition-period variability of an average performer.
#' @param noise Joint-specific noise SD of an average performer.
#' @param epoch_len,overlap Windowing for the cross-spectrum.
#' @param bands Named list of bands.
#' @param seed Random seed (NULL leaves the generator state untouched).
#' @param keep_signals Also return the simulated recordings.
#' @return List with \code{features} and \code{covariates}: observed cadence
#'   (Hz) and repetition-period CV per subject, as video analysis would give.
#' @examples
#' z_speed <- rep(c(1, -1), each = 4)
#' z_accuracy <- rep(c(1, -1), times = 4)
#' sim <- mt_simulate_motion(z_speed, z_accuracy, duration = 30, seed = 1)
#' sim$covariates
#' dim(sim$features$movement$wpli)
#' @export
mt_simulate_motion <- function(z_speed, z_accuracy, fs = 30, duration = 60, f0 = 0.63,
                               n_joints = 6L, lag_unit = 0.35, lag_jitter = 0.5, period_cv = 0.08,
                               noise = 0.8, epoch_len = 8, overlap = 0.5,
                               bands = mt_motion_bands(), seed = 1L, keep_signals = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(z_speed)
  if (length(z_accuracy) != n) stop("'z_speed' and 'z_accuracy' must have the same length.")
  pattern <- c(0, 1, 2, -1, 0.5, 1.5)[seq_len(n_joints)]
  jn <- c("hip", "knee", "ankle", "trunk", "com", "arm")[seq_len(n_joints)]
  tt <- seq(0, duration - 1 / fs, by = 1 / fs)
  cadence <- cv <- numeric(n)
  sig <- vector("list", n)
  for (s in seq_len(n)) {
    f <- f0 * exp(0.15 * z_speed[s])
    sp <- period_cv * exp(-0.6 * z_accuracy[s])
    per <- 1 / (f * exp(sp * rnorm(ceiling(duration * f * 2) + 5L)))
    ends <- cumsum(per)
    k <- findInterval(tt, c(0, ends))                       # repetition index (1-based)
    start <- c(0, ends)[k]
    phi <- 2 * pi * ((k - 1) + (tt - start) / per[k])
    done <- sum(ends <= duration)
    cadence[s] <- done / ends[done]
    cv[s] <- stats::sd(per[seq_len(done)]) / mean(per[seq_len(done)])
    lag <- lag_unit * exp(-0.5 * z_speed[s]) * pattern
    sl <- lag_jitter * exp(-0.7 * z_accuracy[s])
    sn <- noise * exp(-0.6 * z_speed[s])
    knots <- seq(0, duration + 2, by = 2)
    x <- sapply(seq_len(n_joints), function(j) {
      dev <- stats::approx(knots, rnorm(length(knots), 0, sl), xout = tt)$y
      th <- phi - lag[j] - dev
      cos(th) + 0.25 * cos(2 * th) + rnorm(length(tt), 0, sn)
    })
    colnames(x) <- jn
    sig[[s]] <- x
  }
  names(sig) <- sprintf("sub-%03d", seq_len(n))
  ep <- lapply(sig, mt_epoch, sfreq = fs, epoch_len = epoch_len, overlap = overlap)
  out <- list(features = mt_features_from_epochs(ep, fs, bands),
              covariates = data.frame(subject = names(sig), cadence_hz = cadence, period_cv = cv,
                                      stringsAsFactors = FALSE))
  if (keep_signals) out$signals <- sig
  out
}

#' Metric-T under alternative labelling rules
#'
#' Applies Metric-T to the same features under several ways of assigning
#' subjects to two classes, so that the dependence of the conclusion on the
#' labelling rule can be read from one table.
#'
#' @param features As for \code{\link{metric_t}}.
#' @param labelings Data frame or named list; each element gives one label
#'   per subject (\code{NA} excludes the subject under that rule).
#' @param good,poor The two labels; DC is the percentage of features with good > poor.
#' @param n_perm Permutations per rule (0 skips the test; p-values are then NA).
#' @param seed Random seed.
#' @return Data frame with one row per rule and band.
#' @examples
#' sim <- mt_simulate_motion(rep(c(1, -1), each = 6), rep(c(1, -1), times = 6),
#'                           duration = 30, seed = 1)
#' fast <- sim$covariates$cadence_hz > median(sim$covariates$cadence_hz)
#' steady <- sim$covariates$period_cv < median(sim$covariates$period_cv)
#' rules <- list(by_speed = ifelse(fast, "good", "poor"),
#'               by_accuracy = ifelse(steady, "good", "poor"))
#' mt_labelings(sim$features, rules, n_perm = 200)
#' @export
mt_labelings <- function(features, labelings, good = "good", poor = "poor", n_perm = 2000L, seed = 42L) {
  labelings <- as.list(labelings)
  if (is.null(names(labelings))) names(labelings) <- paste0("rule", seq_along(labelings))
  res <- lapply(names(labelings), function(r) {
    lab <- as.character(labelings[[r]])
    keep <- !is.na(lab) & lab %in% c(good, poor)
    ng <- sum(lab[keep] == good); np <- sum(lab[keep] == poor)
    if (ng < 2L || np < 2L)
      return(data.frame(rule = r, band = names(features), n_good = ng, n_poor = np, DC_wPLI = NA_real_,
                        DC_MSC = NA_real_, T = NA_real_, reversal = NA, p_raw = NA_real_, p_min = NA_real_,
                        stringsAsFactors = FALSE))
    sub <- lapply(features, function(f) list(wpli = as.matrix(f$wpli)[keep, , drop = FALSE],
                                             msc = as.matrix(f$msc)[keep, , drop = FALSE]))
    t <- metric_t(sub, lab[keep], comparisons = list(c(good, poor)), n_perm = n_perm, seed = seed,
                  keep_null = FALSE)$table
    if (n_perm < 1L) t$p_raw <- t$p_min <- NA_real_
    data.frame(rule = r, band = t$band, n_good = ng, n_poor = np, DC_wPLI = t$DC_wPLI, DC_MSC = t$DC_MSC,
               T = t$T, reversal = t$reversal, p_raw = t$p_raw, p_min = t$p_min, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

#' Metric-T over the complete set of admissible labellings
#'
#' Subjects whose class is clear keep their label; each of the \eqn{m}
#' ambiguous subjects is assigned to either class in every possible way,
#' giving \eqn{2^m} labellings. Metric-T is evaluated on each one. The range
#' of T and the proportion of labellings with a reversal describe how far the
#' conclusion depends on the treatment of ambiguous cases.
#'
#' @param features As for \code{\link{metric_t}}.
#' @param labels One label per subject; ambiguous subjects may be \code{NA}
#'   or any value other than \code{good} and \code{poor}.
#' @param ambiguous Logical or index vector of ambiguous subjects (default:
#'   every subject whose label is neither \code{good} nor \code{poor}).
#' @param good,poor The two class labels.
#' @param prob Probability that each ambiguous subject belongs to
#'   \code{good} (length 1 or \eqn{m}), e.g. from a sigmoid membership
#'   function. Used as weights for the expected T; 0.5 weights all
#'   labellings equally.
#' @param max_exact Largest number of labellings enumerated exactly; above
#'   this, \code{max_exact} labellings are drawn at random using \code{prob}.
#' @param seed Random seed for the sampled case.
#' @return Object of class \code{mt_labelset}: \code{summary} (one row per
#'   band), \code{T} (labellings by bands), \code{reversal}, \code{n_good},
#'   \code{influence} (mean change in T when a subject is moved from poor to
#'   good) and \code{exact}.
#' @examples
#' sim <- mt_simulate_motion(rep(c(1, -1), each = 6), rep(c(1, -1), times = 6),
#'                           duration = 30, seed = 1)
#' fast <- sim$covariates$cadence_hz > median(sim$covariates$cadence_hz)
#' steady <- sim$covariates$period_cv < median(sim$covariates$period_cv)
#' lab <- ifelse(fast & steady, "good", ifelse(!fast & !steady, "poor", NA))
#' table(lab, useNA = "ifany")
#' mt_labelset(sim$features, lab)
#' @export
mt_labelset <- function(features, labels, ambiguous = NULL, good = "good", poor = "poor",
                        prob = 0.5, max_exact = 65536L, seed = 42L) {
  labels <- as.character(labels); N <- length(labels)
  amb <- if (is.null(ambiguous)) is.na(labels) | !(labels %in% c(good, poor)) else {
    a <- logical(N); a[ambiguous] <- TRUE; a }
  if (any(!amb & !(labels %in% c(good, poor)))) stop("Non-ambiguous subjects need a 'good' or 'poor' label.")
  M <- which(amb); m <- length(M)
  prob <- rep_len(prob, m)
  exact <- 2^m <= max_exact
  if (exact) {
    G <- if (m) as.matrix(expand.grid(rep(list(0:1), m))) else matrix(0L, 1L, 0L)
    w <- if (m) exp(G %*% log(prob) + (1 - G) %*% log1p(-prob))[, 1] else 1
  } else {
    set.seed(seed)
    G <- matrix(as.integer(stats::runif(max_exact * m) < rep(prob, each = max_exact)), max_exact, m)
    w <- rep(1, max_exact)
  }
  storage.mode(G) <- "double"
  fg <- !amb & labels == good
  ng <- sum(fg) + rowSums(G); np <- N - ng
  ok <- ng >= 1 & np >= 1
  dc <- function(X) {
    X <- as.matrix(X)
    Sg <- sweep(G %*% X[M, , drop = FALSE], 2, colSums(X[fg, , drop = FALSE]), "+")
    St <- colSums(X)
    100 * rowMeans(Sg / ng > (rep(St, each = nrow(G)) - Sg) / np)
  }
  bands <- names(features)
  Tm <- Rv <- matrix(NA_real_, nrow(G), length(bands), dimnames = list(NULL, bands))
  for (b in bands) {
    dw <- dc(features[[b]]$wpli); dm <- dc(features[[b]]$msc)
    Tm[, b] <- dw - dm
    Rv[, b] <- (dw > 50 & dm < 50) | (dw < 50 & dm > 50)
  }
  Tm[!ok, ] <- NA; Rv[!ok, ] <- NA
  ww <- w * ok; ww <- ww / sum(ww)
  sm <- data.frame(band = bands, n_ambiguous = m, n_labellings = sum(ok),
                   T_min = apply(Tm, 2, min, na.rm = TRUE), T_median = apply(Tm, 2, stats::median, na.rm = TRUE),
                   T_max = apply(Tm, 2, max, na.rm = TRUE),
                   T_expected = colSums(Tm * ww, na.rm = TRUE),
                   reversal_pct = 100 * colSums(Rv * ww, na.rm = TRUE),
                   row.names = NULL, stringsAsFactors = FALSE)
  infl <- if (m) sapply(bands, function(b) sapply(seq_len(m), function(j)
    mean(Tm[G[, j] == 1, b], na.rm = TRUE) - mean(Tm[G[, j] == 0, b], na.rm = TRUE))) else NULL
  if (m) infl <- matrix(infl, m, length(bands), dimnames = list(names(features[[1]]$wpli[, 1])[M] %||% M, bands))
  structure(list(summary = sm, T = Tm, reversal = Rv == 1, n_good = ng, assignments = G, ambiguous = M,
                 influence = infl, exact = exact), class = "mt_labelset")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.mt_labelset <- function(x, ...) {
  cat(sprintf("Metric-T over %s labellings of %d ambiguous subjects (%s)\n",
              format(x$summary$n_labellings[1], big.mark = ","), x$summary$n_ambiguous[1],
              if (x$exact) "complete enumeration" else "random sample"))
  print(format(x$summary[, -(2:3)], digits = 3, nsmall = 1), row.names = FALSE)
  invisible(x)
}
