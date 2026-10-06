#' Simulate EEG with tunable coupling-mode composition
#'
#' Each channel mixes a phase-lagged oscillation (seen by wPLI) and a common
#' zero-lag oscillation (seen by MSC but not by wPLI) plus white noise. With
#' the defaults, group A has stronger but mostly zero-lag coupling and group B
#' weaker but mostly phase-lagged coupling, so MSC ranks A above B while wPLI
#' ranks B above A in the band containing \code{f0}: the situation Metric-T
#' is designed to flag. Intended for demonstration and software testing.
#'
#' @param n Named integer vector of group sizes.
#' @param zero_lag Zero-lag fraction per group (same order as \code{n}).
#' @param strength Amplitude of the oscillatory (coupled) part per group.
#' @param subject_sd Between-subject SD of the zero-lag fraction and (relative) strength.
#' @param gain_sd SD (log scale) of subject-by-channel gains of the coupled part; values above 0
#'   make channel pairs vary less in unison across subjects.
#' @param n_ch,n_epochs,sfreq,epoch_len Recording layout.
#' @param f0 Oscillation frequency (Hz).
#' @param noise Noise standard deviation.
#' @param bands Named list of bands.
#' @param seed Random seed.
#' @param keep_epochs Also return the simulated epochs.
#' @return List with \code{features}, \code{group} and optionally \code{epochs}.
#' @export
mt_simulate <- function(n = c(A = 20, B = 20), zero_lag = c(0.58, 0.42), strength = c(1, 0.9),
                        subject_sd = 0.05, gain_sd = 0.7, n_ch = 8L,
                        n_epochs = 30L, sfreq = 256, epoch_len = 2, f0 = 10, noise = 0.5,
                        bands = mt_bands(), seed = 1L, keep_epochs = FALSE) {
  set.seed(seed)
  if (is.null(names(n))) names(n) <- LETTERS[seq_along(n)]
  zero_lag <- rep_len(zero_lag, length(n))
  strength <- rep_len(strength, length(n))
  amp <- rep(strength, times = n) * exp(rnorm(sum(n), 0, subject_sd))
  n_t <- as.integer(round(sfreq * epoch_len))
  tt <- (0:(n_t - 1L)) / sfreq
  group <- rep(names(n), times = n)
  zl <- pmin(pmax(rep(zero_lag, times = n) + rnorm(length(group), 0, subject_sd), 0), 1)
  ep <- lapply(seq_along(group), function(s) {
    a <- array(0, dim = c(n_epochs, n_ch, n_t))
    gain <- exp(rnorm(n_ch, 0, gain_sd))
    for (e in seq_len(n_epochs)) {
      common <- sin(2 * pi * f0 * tt + runif(1, 0, 2 * pi))
      base <- runif(1, 0, 2 * pi)
      for (ch in seq_len(n_ch)) {
        lag <- sin(2 * pi * f0 * tt + base + 0.3 * ch)
        a[e, ch, ] <- amp[s] * gain[ch] * ((1 - zl[s]) * lag + zl[s] * common) + noise * rnorm(n_t)
      }
    }
    a
  })
  names(ep) <- sprintf("sub-%03d", seq_along(group))
  out <- list(features = mt_features_from_epochs(ep, sfreq, bands), group = group)
  if (keep_epochs) out$epochs <- ep
  out
}

#' Launch the Metric-T graphical interface
#'
#' @param ... Passed to \code{shiny::runApp}.
#' @export
run_metricT_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("Package 'shiny' is required for the GUI: install.packages('shiny')")
  shiny::runApp(system.file("shiny", "metricT", package = "metricT"), ...)
}
