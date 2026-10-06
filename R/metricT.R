#' Direction consistency (DC)
#'
#' Percentage of features for which the mean of group 1 exceeds the mean of
#' group 2. DC removes the native unit of a connectivity metric, so that
#' different metrics can be compared on one scale ("unit erasure").
#'
#' @param g1,g2 Numeric matrices (subjects x features) for the two groups.
#' @return A single number between 0 and 100.
#' @examples
#' g1 <- matrix(rnorm(5 * 10, mean = 0.5), nrow = 5)   # 5 subjects, 10 features
#' g2 <- matrix(rnorm(6 * 10), nrow = 6)               # 6 subjects
#' mt_dc(g1, g2)
#' @export
mt_dc <- function(g1, g2) {
  100 * mean(colMeans(as.matrix(g1)) > colMeans(as.matrix(g2)))
}

.mt_check_features <- function(features, n) {
  if (!is.list(features) || is.null(names(features)))
    stop("'features' must be a named list: features[[band]]$wpli / $msc.")
  for (b in names(features)) {
    f <- features[[b]]
    if (!all(c("wpli", "msc") %in% names(f)))
      stop("Band '", b, "' needs both 'wpli' and 'msc' matrices.")
    w <- as.matrix(f$wpli); m <- as.matrix(f$msc)
    if (!identical(dim(w), dim(m))) stop("wPLI and MSC dimensions differ in band '", b, "'.")
    if (nrow(w) != n) stop("Band '", b, "': rows (", nrow(w), ") != length(group) (", n, ").")
    if (anyNA(w) || anyNA(m)) stop("Missing values in band '", b, "'.")
  }
  invisible(TRUE)
}

#' Westfall-Young max-T step-down adjustment
#'
#' @param T_obs Observed statistics for the m conditions of one family.
#' @param T_perm m x B matrix of statistics under B aligned permutations
#'   (column b is one relabelling applied to every condition).
#' @param exact \code{TRUE} when the columns of \code{T_perm} are the complete
#'   set of label assignments (the observed one included).
#' @return Adjusted p-values, in the order of \code{T_obs}.
#' @examples
#' T_obs <- c(40, -10, 20)
#' T_perm <- matrix(sample(seq(-100, 100, by = 10), 3 * 500, replace = TRUE), nrow = 3)
#' mt_wy_maxt(T_obs, T_perm)
#' @export
mt_wy_maxt <- function(T_obs, T_perm, exact = FALSE) {
  so <- abs(T_obs); sp <- abs(T_perm)
  if (is.null(dim(sp))) sp <- matrix(sp, nrow = 1L)
  m <- nrow(sp); B <- ncol(sp)
  ord <- order(-so)
  so_o <- so[ord]; sp_o <- sp[ord, , drop = FALSE]
  u <- sp_o
  if (m > 1L) for (i in (m - 1L):1L) u[i, ] <- pmax(u[i + 1L, ], sp_o[i, ])
  adj <- if (exact) rowSums(u >= so_o) / B else (1 + rowSums(u >= so_o)) / (B + 1)
  adj <- cummax(adj)
  out <- numeric(m); out[ord] <- adj
  out
}

#' Metric-T: diagnostic for directional fragility
#'
#' For every comparison x band condition, computes DC(wPLI), DC(MSC) and
#' T = DC(wPLI) - DC(MSC) in percentage points, and tests T by permuting the
#' group labels of all subjects. One relabelling is shared by all conditions
#' (aligned permutations), which preserves their dependence and allows
#' Westfall-Young max-T adjustment.
#'
#' @param features Named list by band; each element is a list with matrices
#'   \code{wpli} and \code{msc} (subjects x channel pairs).
#' @param group Group label of each subject (row order of the matrices).
#' @param comparisons List of length-2 character vectors \code{c(g1, g2)}.
#'   Default: each group against \code{control} first, then remaining pairs.
#' @param control Optional reference group. If given, the Westfall-Young
#'   family is, within each band, the set of comparisons against it.
#' @param n_perm Number of permutations (default 10000).
#' @param seed Random seed.
#' @param wy_family \code{"control"}, \code{"band"} or \code{"all"}.
#' @param perm_index Optional n_perm x n matrix of permutation indices
#'   (1-based); used to replay permutations generated elsewhere.
#' @param exact If \code{TRUE} (two groups only), every distinct assignment of
#'   the group labels is enumerated instead of sampled, and p-values are exact
#'   proportions over that complete reference set. Feasible for small samples
#'   (e.g. 4 versus 9 subjects gives 715 assignments).
#' @param keep_null Keep the permutation matrix in the result.
#' @param progress Optional function called with the fraction completed.
#' @return Object of class \code{metricT}: \code{table} (one row per
#'   condition), \code{T_perm}, and settings. Column \code{p_min} is the
#'   smallest p-value attainable in that condition (the permutation
#'   probability of |T| = 100 pp); because DC keeps only the sign of each
#'   mean difference, strongly dependent features can make even a complete
#'   reversal non-significant, and \code{p_min} shows when that is the case.
#' @references Kimura A (2026). Metric-T: A permutation-based diagnostic for
#'   directional fragility in EEG functional connectivity analysis.
#'   Neuroscience Informatics 6:100286. doi:10.1016/j.neuri.2026.100286
#' @examples
#' sim <- mt_simulate(n = c(A = 12, B = 12), n_ch = 5, n_epochs = 10, seed = 1)
#' res <- metric_t(sim$features, sim$group, n_perm = 200)
#' res
#' @export
metric_t <- function(features, group, comparisons = NULL, control = NULL,
                     n_perm = 10000L, seed = 42L,
                     wy_family = c("control", "band", "all"),
                     perm_index = NULL, exact = FALSE, keep_null = TRUE, progress = NULL) {
  group <- as.character(group)
  n <- length(group)
  .mt_check_features(features, n)
  lev <- unique(group)
  if (length(lev) < 2L) stop("At least two groups are required.")
  if (!is.null(control) && !(control %in% lev)) stop("'control' is not a group label.")
  if (is.null(comparisons)) {
    if (!is.null(control)) {
      oth <- setdiff(lev, control)
      comparisons <- lapply(oth, function(g) c(g, control))
      if (length(oth) > 1L)
        comparisons <- c(comparisons, combn(oth, 2L, simplify = FALSE))
    } else {
      comparisons <- combn(lev, 2L, simplify = FALSE)
    }
  }
  if (is.matrix(comparisons)) comparisons <- lapply(seq_len(nrow(comparisons)), function(i) comparisons[i, ])
  cmp <- do.call(rbind, comparisons)
  if (!all(cmp %in% lev)) stop("Unknown group in 'comparisons'.")
  wy_family <- match.arg(wy_family)
  if (wy_family == "control" && is.null(control)) wy_family <- "band"

  bands <- names(features)
  K <- nrow(cmp); nb <- length(bands)
  W <- do.call(cbind, lapply(features, function(f) unname(as.matrix(f$wpli))))
  M <- do.call(cbind, lapply(features, function(f) unname(as.matrix(f$msc))))
  pb <- vapply(features, function(f) ncol(as.matrix(f$wpli)), 1L)
  Bm <- matrix(0, sum(pb), nb)
  Bm[cbind(seq_len(sum(pb)), rep(seq_len(nb), pb))] <- 1
  gi <- match(group, lev)
  ng <- tabulate(gi, length(lev))
  a <- match(cmp[, 1], lev); b <- match(cmp[, 2], lev)

  counts <- function(g) {
    mw <- rowsum(W, g, reorder = TRUE) / ng
    mm <- rowsum(M, g, reorder = TRUE) / ng
    cw <- (mw[a, , drop = FALSE] > mw[b, , drop = FALSE]) %*% Bm   # K x nb
    cm <- (mm[a, , drop = FALSE] > mm[b, , drop = FALSE]) %*% Bm
    list(cw = cw, cm = cm)
  }
  scale <- matrix(100 / pb, K, nb, byrow = TRUE)
  obs <- counts(gi)
  DCw <- obs$cw * scale; DCm <- obs$cm * scale
  T_obs <- as.vector((obs$cw - obs$cm) * scale)          # band-major, as in the reference

  lab_mat <- NULL
  if (isTRUE(exact)) {
    if (length(lev) != 2L) stop("'exact = TRUE' requires exactly two groups.")
    n_lab <- choose(n, ng[1])
    if (n_lab > 2e6) stop("Too many label assignments (", n_lab, ") for exact enumeration.")
    lab_mat <- utils::combn(n, ng[1])                     # members of group 1, one column per assignment
    n_perm <- ncol(lab_mat)
  } else if (!is.null(perm_index)) {
    perm_index <- as.matrix(perm_index)
    if (ncol(perm_index) != n) stop("'perm_index' must have one column per subject.")
    n_perm <- nrow(perm_index)
  } else {
    set.seed(seed)
  }
  n_perm <- as.integer(n_perm)
  T_perm <- matrix(0, K * nb, n_perm)
  tick <- max(1L, n_perm %/% 50L)
  for (i in seq_len(n_perm)) {
    g <- if (!is.null(lab_mat)) { z <- rep(2L, n); z[lab_mat[, i]] <- 1L; z }
         else if (is.null(perm_index)) gi[sample.int(n)] else gi[perm_index[i, ]]
    cc <- counts(g)
    T_perm[, i] <- as.vector((cc$cw - cc$cm) * scale)
    if (is.function(progress) && i %% tick == 0L) progress(i / n_perm)
  }
  if (isTRUE(exact)) {
    p_raw <- rowSums(abs(T_perm) >= abs(T_obs)) / n_perm
    p_min <- pmax(rowSums(abs(T_perm) >= 100 - 1e-9), 1) / n_perm
  } else {
    p_raw <- (1 + rowSums(abs(T_perm) >= abs(T_obs))) / (n_perm + 1)
    # smallest p-value any observed T could attain in this condition (|T| = 100 pp)
    p_min <- (1 + rowSums(abs(T_perm) >= 100 - 1e-9)) / (n_perm + 1)
  }

  cmp_lab <- paste(cmp[, 1], "vs", cmp[, 2])
  tab <- data.frame(
    comparison = rep(cmp_lab, times = nb),
    band = rep(bands, each = K),
    n1 = rep(ng[a], times = nb), n2 = rep(ng[b], times = nb),
    n_features = rep(pb, each = K),
    DC_wPLI = as.vector(DCw), DC_MSC = as.vector(DCm),
    T = T_obs,
    reversal = as.vector((DCw - 50) * (DCm - 50) < 0),
    p_raw = p_raw,
    p_min = p_min,
    p_bonferroni = p.adjust(p_raw, "bonferroni"),
    p_holm = p.adjust(p_raw, "holm"),
    q_BH = p.adjust(p_raw, "BH"),
    p_WY = NA_real_,
    stringsAsFactors = FALSE
  )
  in_ctrl <- rep(cmp[, 1] == control | cmp[, 2] == control, times = nb)
  fam <- switch(wy_family,
    control = ifelse(in_ctrl, tab$band, NA_character_),
    band = tab$band,
    all = rep("all", nrow(tab)))
  for (f in unique(fam[!is.na(fam)])) {
    ix <- which(fam == f)
    tab$p_WY[ix] <- mt_wy_maxt(T_obs[ix], T_perm[ix, , drop = FALSE], exact = isTRUE(exact))
  }
  rownames(T_perm) <- paste(tab$comparison, tab$band, sep = " | ")
  structure(list(table = tab, T_perm = if (keep_null) T_perm else NULL,
                 n_perm = n_perm, exact = isTRUE(exact), seed = if (is.null(perm_index) && !isTRUE(exact)) seed else NA,
                 control = control, wy_family = wy_family,
                 groups = stats::setNames(ng, lev), bands = bands),
            class = "metricT")
}

#' @export
as.data.frame.metricT <- function(x, ...) x$table

#' @export
print.metricT <- function(x, digits = 1, ...) {
  cat("Metric-T: T = DC(wPLI) - DC(MSC), percentage points\n")
  cat("Groups: ", paste0(names(x$groups), " (n=", x$groups, ")", collapse = ", "), "\n", sep = "")
  if (isTRUE(x$exact)) {
    cat("Label assignments: all ", x$n_perm, " enumerated (exact, two-sided)\n", sep = "")
  } else {
    cat("Permutations: ", x$n_perm, " (aligned label permutation, two-sided)",
        if (!is.na(x$seed)) paste0(", seed ", x$seed), "\n", sep = "")
  }
  cat("Westfall-Young family: ",
      switch(x$wy_family, control = paste0("per band, comparisons against ", x$control),
             band = "per band", all = "all conditions"), "\n\n", sep = "")
  t <- x$table
  out <- data.frame(Comparison = t$comparison, Band = t$band,
                    `DC(wPLI)%` = formatC(t$DC_wPLI, format = "f", digits = digits),
                    `DC(MSC)%` = formatC(t$DC_MSC, format = "f", digits = digits),
                    `T(pp)` = formatC(t$T, format = "f", digits = digits, flag = "+"),
                    Reversal = ifelse(t$reversal, "YES", ""),
                    p = formatC(t$p_raw, format = "f", digits = 3),
                    p_min = formatC(t$p_min, format = "f", digits = 3),
                    p_Holm = formatC(t$p_holm, format = "f", digits = 3),
                    q_BH = formatC(t$q_BH, format = "f", digits = 3),
                    p_WY = ifelse(is.na(t$p_WY), "", formatC(t$p_WY, format = "f", digits = 3)),
                    check.names = FALSE, stringsAsFactors = FALSE)
  print(out, row.names = FALSE)
  nrev <- sum(t$reversal)
  cat("\nDirectional reversals (DC values on opposite sides of 50%): ", nrev, " of ",
      nrow(t), " conditions\n", sep = "")
  nf <- sum(t$p_min > 0.05)
  if (nf) cat("Note: in ", nf, " condition(s) even |T| = 100 pp could not reach p <= 0.05 (p_min > 0.05):\n",
              "      channel pairs vary in unison across subjects, so a large p there is not evidence of stability.\n", sep = "")
  invisible(x)
}

#' Plot direction consistency for both metrics
#'
#' One panel per comparison; bars show DC(wPLI) and DC(MSC) by band, the
#' dashed line marks 50%, and shaded bands mark directional reversals.
#'
#' @param x A \code{metricT} object.
#' @param col Two colours for wPLI and MSC.
#' @param ... Unused.
#' @return The object \code{x}, invisibly. Called for its side effect of drawing a plot.
#' @examples
#' sim <- mt_simulate(n = c(A = 6, B = 6), n_ch = 4, n_epochs = 8, seed = 1)
#' res <- metric_t(sim$features, sim$group, exact = TRUE)
#' plot(res)
#' @export
plot.metricT <- function(x, col = c("#3b6fb0", "#cf5b4e"), ...) {
  t <- x$table
  cmps <- unique(t$comparison)
  op <- par(mfrow = c(length(cmps), 1L), mar = c(3, 4.5, 2.2, 1), oma = c(1.6, 0, 0, 0))
  on.exit(par(op))
  for (k in seq_along(cmps)) {
    s <- t[t$comparison == cmps[k], ]
    h <- rbind(s$DC_wPLI, s$DC_MSC)
    mids <- barplot(h, beside = TRUE, ylim = c(0, 112), col = col, border = NA,
                    names.arg = s$band, ylab = "Direction consistency (%)", las = 1)
    cx <- colMeans(mids)
    for (j in which(s$reversal)) {
      graphics::rect(cx[j] - 1.4, 0, cx[j] + 1.4, 112, col = adjustcolor("#cf5b4e", 0.10), border = NA)
      text(cx[j], max(h[, j]) + 6, "reversal", col = "#b03030", cex = 0.85, font = 2)
    }
    text(cx, pmax(h[1, ], h[2, ]) + ifelse(s$reversal, 13, 6),
         sprintf("T = %+.1f", s$T), cex = 0.8)
    abline(h = 50, lty = 2)
    mtext(sprintf("%s  (n = %d vs %d)", cmps[k], s$n1[1], s$n2[1]), side = 3, adj = 0, font = 2, line = 0.5)
    if (k == 1L) legend("topright", c("DC(wPLI)", "DC(MSC)"), fill = col, border = NA,
                        bty = "n", horiz = TRUE, cex = 0.9)
    box(bty = "l")
  }
  mtext(sprintf("Two-sided label permutation, %d iterations", x$n_perm),
        side = 1, outer = TRUE, cex = 0.75, line = 0.3)
  invisible(x)
}

#' Permutation null distribution of T for one condition
#'
#' @param x A \code{metricT} object computed with \code{keep_null = TRUE}.
#' @param condition Row number of \code{x$table}.
#' @return The permutation values of T for the chosen condition, invisibly.
#'   Called for its side effect of drawing a histogram.
#' @examples
#' sim <- mt_simulate(n = c(A = 6, B = 6), n_ch = 4, n_epochs = 8, seed = 1)
#' res <- metric_t(sim$features, sim$group, exact = TRUE)
#' mt_null_plot(res, condition = 1)
#' @export
mt_null_plot <- function(x, condition = 1L) {
  if (is.null(x$T_perm)) stop("Null distribution was not kept (keep_null = FALSE).")
  t <- x$table[condition, ]
  v <- x$T_perm[condition, ]
  lim <- max(abs(c(v, t$T))) * 1.05
  hist(v, breaks = 41, col = "grey80", border = "white", xlim = c(-lim, lim),
       main = sprintf("%s, %s", t$comparison, t$band), xlab = "T under label permutation (pp)", las = 1)
  abline(v = t$T, col = "#b03030", lwd = 2)
  abline(v = -t$T, col = "#b03030", lwd = 1, lty = 2)
  mtext(sprintf("observed T = %+.1f pp,  p = %.3f", t$T, t$p_raw), side = 3, line = 0.2, cex = 0.9)
  invisible(v)
}
