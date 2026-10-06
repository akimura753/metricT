#' Default frequency bands (Hz), as in Kimura (2026)
#'
#' A band c(lo, hi) keeps Fourier bins with lo <= f < hi.
#' @return Named list of length-2 numeric vectors.
#' @examples
#' mt_bands()
#' @export
mt_bands <- function() {
  list(theta = c(4, 8), alpha = c(8, 13), beta = c(13, 30), broadband = c(1, 40))
}

#' Cut a continuous recording into fixed-length epochs
#'
#' @param x Numeric matrix, samples in rows and channels in columns.
#' @param sfreq Sampling frequency (Hz).
#' @param epoch_len Epoch length in seconds (default 2).
#' @param overlap Fractional overlap between successive epochs (default 0.5).
#' @return Array of dimension epochs x channels x samples.
#' @examples
#' x <- matrix(rnorm(128 * 10 * 3), ncol = 3)   # 10 s of 3 channels at 128 Hz
#' ep <- mt_epoch(x, sfreq = 128, epoch_len = 2, overlap = 0.5)
#' dim(ep)                                      # epochs x channels x samples
#' @export
mt_epoch <- function(x, sfreq, epoch_len = 2, overlap = 0.5) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric (samples x channels).")
  n_t <- as.integer(round(epoch_len * sfreq))
  step <- as.integer(round(n_t * (1 - overlap)))
  if (step < 1L) stop("'overlap' must be smaller than 1.")
  if (nrow(x) < n_t) stop("Recording is shorter than one epoch.")
  starts <- seq.int(1L, nrow(x) - n_t + 1L, by = step)
  out <- array(0, dim = c(length(starts), ncol(x), n_t))
  for (e in seq_along(starts)) {
    out[e, , ] <- t(x[starts[e]:(starts[e] + n_t - 1L), , drop = FALSE])
  }
  dimnames(out) <- list(NULL, colnames(x), NULL)
  out
}

#' wPLI and MSC from the same cross-spectrum (one subject)
#'
#' Each epoch is Hann-tapered and Fourier transformed. The cross-spectrum
#' S_ij(f) = X_i(f) Conj(X_j(f)) is averaged over the bins of a band within
#' each epoch; epochs are the observations:
#' wPLI_ij = |mean(Im S_ij)| / mean(|Im S_ij|) and
#' MSC_ij = |mean(S_ij)|^2 / (mean(S_ii) mean(S_jj)).
#' Both metrics therefore come from identical cross-spectral estimates.
#'
#' @param epochs Array, epochs x channels x samples.
#' @param sfreq Sampling frequency (Hz).
#' @param bands Named list of bands, see \code{\link{mt_bands}}.
#' @return Named list (one element per band) with numeric vectors
#'   \code{wpli} and \code{msc} over channel pairs in \code{combn} order.
#' @examples
#' ep <- mt_epoch(matrix(rnorm(128 * 20 * 3), ncol = 3), sfreq = 128)
#' con <- mt_connectivity(ep, sfreq = 128)
#' str(con$alpha)
#' @export
mt_connectivity <- function(epochs, sfreq, bands = mt_bands()) {
  d <- dim(epochs)
  if (length(d) != 3L) stop("'epochs' must be an array: epochs x channels x samples.")
  n_ep <- d[1]; n_ch <- d[2]; n_t <- d[3]
  if (n_ch < 2L) stop("At least two channels are required.")
  win <- 0.5 - 0.5 * cos(2 * pi * (0:(n_t - 1L)) / (n_t - 1L))   # symmetric Hann
  n_f <- n_t %/% 2L + 1L
  freqs <- (0:(n_f - 1L)) * sfreq / n_t
  sel <- lapply(bands, function(b) which(freqs >= b[1] & freqs < b[2]))
  empty <- names(sel)[lengths(sel) == 0L]
  if (length(empty)) stop("No Fourier bins in band(s): ", paste(empty, collapse = ", "))

  acc <- lapply(bands, function(b) list(S = matrix(0+0i, n_ch, n_ch),
                                        absIm = matrix(0, n_ch, n_ch),
                                        P = numeric(n_ch)))
  for (e in seq_len(n_ep)) {
    x <- matrix(epochs[e, , ], nrow = n_ch)              # channels x samples
    X <- mvfft(t(x) * win)[seq_len(n_f), , drop = FALSE] # bins x channels
    for (b in names(bands)) {
      Xb <- X[sel[[b]], , drop = FALSE]
      S <- (t(Xb) %*% Conj(Xb)) / nrow(Xb)               # S[i,j] = mean_f X_i Conj(X_j)
      acc[[b]]$S <- acc[[b]]$S + S
      acc[[b]]$absIm <- acc[[b]]$absIm + abs(Im(S))
      acc[[b]]$P <- acc[[b]]$P + Re(diag(S))
    }
  }
  pr <- combn(n_ch, 2L)
  idx <- cbind(pr[1, ], pr[2, ])
  chn <- dimnames(epochs)[[2]]
  if (is.null(chn)) chn <- paste0("ch", seq_len(n_ch))
  pair_names <- paste(chn[pr[1, ]], chn[pr[2, ]], sep = "-")
  lapply(acc, function(a) {
    S <- a$S[idx]; aim <- a$absIm[idx]
    den <- a$P[pr[1, ]] * a$P[pr[2, ]]
    wpli <- ifelse(aim > 0, abs(Im(S)) / aim, 0)
    msc <- ifelse(den > 0, Mod(S)^2 / den, 0)
    names(wpli) <- names(msc) <- pair_names
    list(wpli = wpli, msc = msc)
  })
}

#' Subject-by-feature matrices from per-subject epoch arrays
#'
#' @param epochs_list List of arrays (epochs x channels x samples), one per subject.
#' @param sfreq Sampling frequency (Hz); a single value or one per subject.
#' @param bands Named list of bands.
#' @param progress Optional function called with the subject index.
#' @return Named list by band of \code{list(wpli, msc)} matrices
#'   (subjects x channel pairs).
#' @examples
#' ep <- lapply(1:4, function(i) mt_epoch(matrix(rnorm(128 * 20 * 3), ncol = 3), sfreq = 128))
#' names(ep) <- paste0("s", 1:4)
#' feats <- mt_features_from_epochs(ep, sfreq = 128)
#' dim(feats$alpha$wpli)                        # subjects x channel pairs
#' @export
mt_features_from_epochs <- function(epochs_list, sfreq, bands = mt_bands(), progress = NULL) {
  n <- length(epochs_list)
  sfreq <- rep_len(sfreq, n)
  res <- vector("list", n)
  for (s in seq_len(n)) {
    res[[s]] <- mt_connectivity(epochs_list[[s]], sfreq[s], bands)
    if (is.function(progress)) progress(s)
  }
  np <- unique(vapply(res, function(r) length(r[[1]]$wpli), 1L))
  if (length(np) != 1L) stop("Subjects have different numbers of channels.")
  feats <- lapply(names(bands), function(b) {
    w <- do.call(rbind, lapply(res, function(r) r[[b]]$wpli))
    m <- do.call(rbind, lapply(res, function(r) r[[b]]$msc))
    rownames(w) <- rownames(m) <- names(epochs_list)
    list(wpli = w, msc = m)
  })
  names(feats) <- names(bands)
  feats
}

#' Features from one CSV file of preprocessed EEG per subject
#'
#' Each file holds one recording: samples in rows, channels in columns, with a
#' header row of channel names. Preprocessing (filtering, artifact removal,
#' re-referencing) is expected to have been done upstream.
#'
#' @param files Character vector of CSV paths.
#' @param ids Subject identifiers (default: file names without extension).
#' @param sfreq Sampling frequency (Hz).
#' @param epoch_len,overlap Passed to \code{\link{mt_epoch}}.
#' @param bands Named list of bands.
#' @param drop_cols Column names to ignore (e.g. a time column).
#' @param progress Optional function called with the subject index.
#' @return Feature list as in \code{\link{mt_features_from_epochs}}.
#' @examples
#' files <- file.path(tempdir(), c("s1.csv", "s2.csv"))
#' n <- 128 * 20                                # 20 s at 128 Hz
#' for (f in files) {
#'   d <- data.frame(time = seq_len(n) / 128, C3 = rnorm(n), C4 = rnorm(n), Pz = rnorm(n))
#'   write.csv(d, f, row.names = FALSE)
#' }
#' feats <- mt_features_from_eeg_csv(files, sfreq = 128)
#' dim(feats$alpha$wpli)
#' unlink(files)
#' @export
mt_features_from_eeg_csv <- function(files, ids = NULL, sfreq, epoch_len = 2, overlap = 0.5,
                                     bands = mt_bands(), drop_cols = c("time", "Time", "t"),
                                     progress = NULL) {
  if (is.null(ids)) ids <- sub("\\.[^.]*$", "", basename(files))
  ep <- lapply(files, function(f) {
    x <- read.csv(f, check.names = FALSE)
    x <- x[, !(names(x) %in% drop_cols), drop = FALSE]
    num <- vapply(x, is.numeric, TRUE)
    mt_epoch(as.matrix(x[, num, drop = FALSE]), sfreq, epoch_len, overlap)
  })
  names(ep) <- ids
  mt_features_from_epochs(ep, sfreq, bands, progress)
}

#' Detrend epochs and drop those with large excursions
#'
#' Minimal preprocessing for raw recordings: each channel of each epoch is
#' linearly detrended (removing offset and slow drift), and an epoch is dropped
#' when the peak-to-peak range of any channel exceeds \code{max_ptp}.
#'
#' @param epochs Array, epochs x channels x samples (see \code{\link{mt_epoch}}).
#' @param detrend Remove a least-squares line from each channel of each epoch.
#' @param max_ptp Peak-to-peak rejection threshold in the unit of the data
#'   (e.g. microvolts); \code{NULL} keeps all epochs.
#' @return The cleaned array, with attributes \code{n_total} and \code{n_kept}.
#' @examples
#' x <- matrix(rnorm(128 * 10 * 3), ncol = 3)
#' x[300, 1] <- 50                              # one artefact
#' ep <- mt_clean_epochs(mt_epoch(x, sfreq = 128), detrend = TRUE, max_ptp = 20)
#' c(kept = attr(ep, "n_kept"), total = attr(ep, "n_total"))
#' @export
mt_clean_epochs <- function(epochs, detrend = TRUE, max_ptp = NULL) {
  d <- dim(epochs)
  if (length(d) != 3L) stop("'epochs' must be an array: epochs x channels x samples.")
  n_t <- d[3]
  if (detrend) {
    tt <- seq_len(n_t) - (n_t + 1) / 2
    stt <- sum(tt^2)
    x <- matrix(epochs, nrow = d[1] * d[2])               # (epoch, channel) rows x samples
    x <- x - rowMeans(x)
    x <- x - tcrossprod(as.vector(x %*% tt) / stt, tt)
    epochs <- array(x, dim = d, dimnames = dimnames(epochs))
  }
  keep <- rep(TRUE, d[1])
  if (!is.null(max_ptp)) {
    ptp <- apply(epochs, c(1, 2), function(v) max(v) - min(v))
    keep <- apply(ptp, 1, max) <= max_ptp
  }
  out <- epochs[keep, , , drop = FALSE]
  attr(out, "n_total") <- d[1]; attr(out, "n_kept") <- sum(keep)
  out
}
