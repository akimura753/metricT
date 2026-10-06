# Dementia example (OpenNeuro ds004504: AD n=36, FTD n=23, CN n=29; 19 channels).
#
# Input: features.npz written by the published Python pipeline
#   (metricT_ds004504_pipeline.py --bids /path/to/ds004504 --out metricT_out),
#   i.e. the subject x channel-pair wPLI and MSC matrices behind Table I of
#   Kimura (2026), Neuroscience Informatics 6:100286.
#
# Usage:  Rscript ds004504_example.R /path/to/metricT_out/features.npz
library(metricT)
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "metricT_out/features.npz"
dat <- mt_read_npz(path)
print(table(dat$group))
res <- metric_t(dat$features, dat$group,
                comparisons = list(c("AD", "CN"), c("FTD", "CN"), c("AD", "FTD")),
                control = "CN", n_perm = 10000, seed = 42)
print(res)
write.csv(as.data.frame(res), "metricT_R_tableI.csv", row.names = FALSE)
grDevices::pdf("metricT_R_fig_DC.pdf", width = 7, height = 8.5); plot(res); grDevices::dev.off()
cat("\nWrote metricT_R_tableI.csv and metricT_R_fig_DC.pdf\n")
