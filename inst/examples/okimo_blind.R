## Metric-T on raw EMOTIV Insight recordings (5 channels, 128 Hz):
## congenitally blind (CB) versus acquired blind (AB) during an electrotactile task.
## Usage: Rscript okimo_blind.R <folder with n*-AB-*.csv / n*-CB-*.csv> [output folder]
## All analysis settings are fixed here, before any result is inspected.
suppressMessages(library(metricT))
args <- commandArgs(trailingOnly = TRUE)
dir_in  <- if (length(args) >= 1) args[1] else "."
dir_out <- if (length(args) >= 2) args[2] else "."

CHANNELS <- c("EEG.AF3", "EEG.T7", "EEG.Pz", "EEG.T8", "EEG.AF4")
SFREQ    <- 128
SKIP_SEC <- 26            # initial rest period, not part of the task
EPOCH    <- 2; OVERLAP <- 0.5
PTP      <- c(primary = 250, strict = 150, lenient = 400)   # microvolt, peak-to-peak
BANDS    <- mt_bands()    # theta, alpha, beta, broadband (Kimura 2026)

files <- list.files(dir_in, pattern = "^n[0-9]+-(AB|CB)-.*\\.csv$", full.names = TRUE)
id    <- sub("^(n[0-9]+)-.*", "\\1", basename(files))
grp   <- sub("^n[0-9]+-(AB|CB)-.*", "\\1", basename(files))
keep  <- !duplicated(id); files <- files[keep]; id <- id[keep]; grp <- grp[keep]
ord   <- order(as.integer(sub("n", "", id))); files <- files[ord]; id <- id[ord]; grp <- grp[ord]

read_task <- function(f) {
  x <- read.csv(f, skip = 1, check.names = FALSE)[, CHANNELS]
  x <- as.matrix(x[stats::complete.cases(x), ])
  colnames(x) <- sub("EEG.", "", colnames(x), fixed = TRUE)
  x[-seq_len(SKIP_SEC * SFREQ), , drop = FALSE]
}
raw <- lapply(files, read_task); names(raw) <- id
ep0 <- lapply(raw, mt_epoch, sfreq = SFREQ, epoch_len = EPOCH, overlap = OVERLAP)

run <- function(ptp) {
  ep <- lapply(ep0, mt_clean_epochs, detrend = TRUE, max_ptp = ptp)
  qc <- data.frame(id = id, group = grp, task_sec = sapply(raw, nrow) / SFREQ,
                   epochs = sapply(ep, attr, "n_total"), kept = sapply(ep, attr, "n_kept"))
  feats <- mt_features_from_epochs(ep, SFREQ, BANDS)
  res <- metric_t(feats, grp, comparisons = list(c("CB", "AB")), exact = TRUE, wy_family = "all")
  list(qc = qc, features = feats, result = res)
}
out <- lapply(PTP, run)

cat("Participants:", length(id), " CB:", sum(grp == "CB"), " AB:", sum(grp == "AB"), "\n\n")
print(out$primary$qc, row.names = FALSE)
m <- unlist(lapply(out$primary$features, `[[`, "msc"))
cat(sprintf("\nMSC range %.3f-%.3f, wPLI range %.3f-%.3f\n\n", min(m), max(m),
            min(unlist(lapply(out$primary$features, `[[`, "wpli"))), max(unlist(lapply(out$primary$features, `[[`, "wpli")))))
for (nm in names(out)) { cat("== rejection threshold", PTP[nm], "uV (", nm, ") ==\n"); print(out[[nm]]$result); cat("\n") }

mt_write_features_csv(out$primary$features, grp, file.path(dir_out, "okimo_features.csv"), subject = id)
saveRDS(out, file.path(dir_out, "okimo_result.rds"))
