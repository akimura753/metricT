## Complete labelling space and stratified trend version of Metric-T.

.mt_feature_mats <- function(features) {
  bands <- names(features)
  list(bands = bands,
       W = lapply(features, function(f) as.matrix(f$wpli)),
       M = lapply(features, function(f) as.matrix(f$msc)))
}

## integer counts of features with mean(group 1) > mean(group 2), for every row of L
.mt_count_gt <- function(L, X, n1, n) {
  s1 <- L %*% X
  s0 <- matrix(colSums(X), nrow(L), ncol(X), byrow = TRUE) - s1
  rowSums(s1 / n1 > s0 / (n - n1))
}

#' Enumerate the complete two-group labelling space
#'
#' Computes Metric-T for every way of splitting the \eqn{n} subjects into two
#' non-empty groups (\eqn{2^n - 2} labellings). The result is the finite set
#' on which exact permutation tests (labellings of one size), attainable
#' significance levels, and alternative group definitions (cut-points,
#' exclusions, batches) can all be located.
#'
#' The spread of direction consistency (DC) over the labelling space measures
#' how many features carry independent information: with \eqn{m} independent
#' features the standard deviation of DC is \eqn{50 / \sqrt{m}}{50 / sqrt(m)}
#' percentage points, so \eqn{m_{eff} = 2500 / Var(DC)}{m_eff = 2500 / Var(DC)}.
#' Features whose group means tie (for example features that are constant
#' across subjects) never count towards DC and inflate this number.
#'
#' @param features Named list by band; each element is a list with matrices
#'   \code{wpli} and \code{msc} (subjects x features), as for \code{\link{metric_t}}.
#' @param max_n Largest number of subjects for which the space is enumerated.
#' @return An object of class \code{"mt_label_space"}: a list with
#'   \item{labelings}{Logical matrix, one row per labelling; \code{TRUE} marks group 1.}
#'   \item{n1}{Size of group 1 in each labelling.}
#'   \item{T, DC_wPLI, DC_MSC}{Matrices (labellings x bands), in percentage points.}
#'   \item{summary}{Data frame per band: SD of T, effective number of features for each
#'     metric, correlation between the two DC values, and the share of labellings in which
#'     the two DC values fall on opposite sides of 50 percent.}
#' @seealso \code{\link{mt_locate}}, \code{\link{metric_t}}
#' @examples
#' sim <- mt_simulate(n = c(A = 5, B = 5), n_ch = 4, n_epochs = 8, seed = 1)
#' sp <- mt_label_space(sim$features)
#' sp
#' mt_locate(sp, sim$group == "A")
#' @export
mt_label_space <- function(features, max_n = 16L) {
  fm <- .mt_feature_mats(features)
  n <- nrow(fm$W[[1L]])
  .mt_check_features(features, n)
  if (n < 2L) stop("At least two subjects are needed.")
  if (n > max_n) stop("n = ", n, " gives ", 2^n - 2, " labellings; raise 'max_n' to enumerate them.")
  L <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
  dimnames(L) <- list(NULL, rownames(fm$W[[1L]]))
  n1 <- rowSums(L)
  keep <- n1 >= 1L & n1 <= n - 1L
  L <- L[keep, , drop = FALSE]; n1 <- n1[keep]
  Ld <- L + 0
  dw <- dm <- matrix(0, nrow(L), length(fm$bands), dimnames = list(NULL, fm$bands))
  for (b in fm$bands) {
    k <- ncol(fm$W[[b]])
    dw[, b] <- 100 * .mt_count_gt(Ld, fm$W[[b]], n1, n) / k
    dm[, b] <- 100 * .mt_count_gt(Ld, fm$M[[b]], n1, n) / k
  }
  Tm <- round(dw - dm, 10)
  pv <- function(v) mean((v - mean(v))^2)
  summ <- data.frame(band = fm$bands,
                     n_features = vapply(fm$W, ncol, 1L),
                     SD_T = apply(Tm, 2L, function(v) sqrt(pv(v))),
                     m_eff_wPLI = 2500 / apply(dw, 2L, pv),
                     m_eff_MSC = 2500 / apply(dm, 2L, pv),
                     cor_DC = vapply(fm$bands, function(b) suppressWarnings(stats::cor(dw[, b], dm[, b])), 0),
                     reversal_share = vapply(fm$bands, function(b) mean((dw[, b] - 50) * (dm[, b] - 50) < 0), 0),
                     row.names = NULL, stringsAsFactors = FALSE)
  structure(list(labelings = L, n1 = n1, T = Tm, DC_wPLI = dw, DC_MSC = dm,
                 summary = summ, n = n, bands = fm$bands),
            class = "mt_label_space")
}

#' @export
print.mt_label_space <- function(x, digits = 2, ...) {
  cat("Complete labelling space: ", x$n, " subjects, ", nrow(x$labelings),
      " two-group labellings\n\n", sep = "")
  s <- x$summary
  num <- vapply(s, is.numeric, TRUE)
  s[num] <- lapply(s[num], round, digits)
  print(s, row.names = FALSE)
  cat("\nm_eff: effective number of independent features (2500 / Var(DC)).\n",
      "reversal_share: labellings with DC(wPLI) and DC(MSC) on opposite sides of 50%.\n", sep = "")
  invisible(x)
}

#' Locate a labelling in the labelling space
#'
#' Reports Metric-T for one labelling and the share of labellings whose
#' \eqn{|T|} is at least as large, both in the complete space and among the
#' labellings in which group 1 has the same size. The latter is the exact
#' two-sided permutation p-value of \code{metric_t(exact = TRUE)}; its smallest
#' attainable value is also returned.
#'
#' @param space Result of \code{\link{mt_label_space}}.
#' @param group1 Logical vector (one element per subject) marking group 1, or
#'   the integer positions of the subjects in group 1.
#' @return Data frame with one row per band: \code{T}, the two DC values,
#'   \code{share_all} (complete space), \code{p_exact} (same group sizes),
#'   \code{p_min} (smallest attainable \code{p_exact}: the share of those
#'   labellings with |T| = 100 pp, as in \code{\link{metric_t}}) and the number
#'   of labellings with the same group sizes.
#' @examples
#' sim <- mt_simulate(n = c(A = 5, B = 5), n_ch = 4, n_epochs = 8, seed = 1)
#' sp <- mt_label_space(sim$features)
#' mt_locate(sp, sim$group == "A")   # the observed labelling
#' mt_locate(sp, 1:3)                # any other split, e.g. a batch or a cut-point
#' @export
mt_locate <- function(space, group1) {
  if (!inherits(space, "mt_label_space")) stop("'space' must come from mt_label_space().")
  g <- rep(FALSE, space$n)
  if (is.logical(group1)) {
    if (length(group1) != space$n) stop("'group1' must have one element per subject.")
    g <- group1
  } else g[group1] <- TRUE
  if (!any(g) || all(g)) stop("Both groups must be non-empty.")
  row <- which(colSums(t(space$labelings) == g) == space$n)
  same <- space$n1 == sum(g)
  tol <- 1e-9
  out <- lapply(space$bands, function(b) {
    a <- abs(space$T[, b]); t0 <- a[row]
    data.frame(band = b, T = space$T[row, b], DC_wPLI = space$DC_wPLI[row, b],
               DC_MSC = space$DC_MSC[row, b],
               share_all = mean(a >= t0 - tol),
               p_exact = mean(a[same] >= t0 - tol),
               p_min = max(sum(a[same] >= 100 - tol), 1) / sum(same),
               n_same_size = sum(same),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.mt_perms <- function(m) {
  if (m <= 1L) return(matrix(1L, 1L, 1L))
  p <- .mt_perms(m - 1L)
  do.call(rbind, lapply(seq_len(m), function(j) cbind(j, p + (p >= j), deparse.level = 0)))
}

.mt_zrank <- function(v, strata) {
  stats::ave(v, strata, FUN = function(u) {
    r <- rank(u)
    s <- stats::sd(r)
    if (length(u) < 2L || is.na(s) || s == 0) r * 0 else (r - mean(r)) / s
  })
}

#' Metric-T for a continuous covariate, with stratified exact permutation
#'
#' Extends Metric-T from two groups to a numeric covariate (for example years
#' of visual experience). For each metric, direction consistency is the
#' percentage of features positively associated with the covariate, and
#' \eqn{T = DC(wPLI) - DC(MSC)}. The covariate is permuted among subjects
#' within strata only (for example recording sessions), so that differences
#' between strata are not treated as chance.
#'
#' With \code{method = "rank"} the covariate and every feature are replaced by
#' their standardised ranks within strata, which makes T invariant to any
#' monotone transformation. With \code{method = "linear"} raw values are
#' centred within strata; for a 0/1 covariate and a single stratum this
#' reproduces the two-group statistic of \code{\link{metric_t}}.
#'
#' Each feature statistic is a sum of stratum contributions, so the complete
#' reference set (\eqn{\prod_s n_s!}{prod(n_s!)} arrangements) is obtained from
#' \eqn{\sum_s n_s!}{sum(n_s!)} evaluations followed by additions.
#'
#' @param features Named list by band, as for \code{\link{metric_t}}.
#' @param x Numeric vector (one value per subject), or a numeric matrix or data
#'   frame with one column per covariate. Several covariates are permuted
#'   jointly, which allows a family-wise adjustment across them.
#' @param strata Optional factor defining the strata; \code{NULL} for none.
#' @param method \code{"rank"} (default) or \code{"linear"}.
#' @param exact \code{TRUE} to enumerate every within-stratum arrangement,
#'   \code{FALSE} to sample; \code{NULL} enumerates when the number of
#'   arrangements does not exceed \code{max_exact}.
#' @param n_perm Number of random arrangements when not exact.
#' @param max_exact Largest reference set enumerated when \code{exact = NULL}.
#' @param sets Optional named list of logical or integer vectors selecting
#'   features (columns), tested in addition to the bands as pooled sets across bands.
#' @param seed Random seed for sampled arrangements.
#' @return An object of class \code{"metricT_trend"}: a list with \code{table}
#'   (covariate, band, DC values, T, two-sided \code{p}, smallest attainable
#'   \code{p_min}, and \code{p_family}, the max-|T| adjusted p over all bands and
#'   covariates), \code{n_perm}, \code{exact} and the stratum sizes.
#' @seealso \code{\link{metric_t}}, \code{\link{mt_label_space}}
#' @examples
#' sim <- mt_simulate(n = c(A = 4, B = 4), n_ch = 4, n_epochs = 8, seed = 2)
#' set.seed(1)
#' years <- round(runif(8, 0, 40))
#' session <- rep(c("day1", "day2"), each = 4)
#' metric_t_trend(sim$features, years, strata = session)
#' @export
metric_t_trend <- function(features, x, strata = NULL, method = c("rank", "linear"),
                           exact = NULL, n_perm = 10000L, max_exact = 2e5,
                           sets = NULL, seed = 42L) {
  method <- match.arg(method)
  fm <- .mt_feature_mats(features)
  n <- nrow(fm$W[[1L]])
  .mt_check_features(features, n)
  X <- if (is.null(dim(x))) matrix(as.numeric(x), ncol = 1L, dimnames = list(NULL, "x")) else as.matrix(x)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (nrow(X) != n) stop("'x' must have one value (row) per subject.")
  if (anyNA(X)) stop("Missing values in 'x'.")
  strata <- if (is.null(strata)) rep("all", n) else as.character(strata)
  if (length(strata) != n) stop("'strata' must have one element per subject.")
  std <- if (method == "rank") .mt_zrank else function(v, s) stats::ave(v, s, FUN = function(u) u - mean(u))
  Zx <- apply(X, 2L, std, strata)
  dim(Zx) <- dim(X); colnames(Zx) <- colnames(X)
  Y <- do.call(cbind, c(fm$W, fm$M))
  Zy <- apply(Y, 2L, std, strata)
  kb <- vapply(fm$W, ncol, 1L)
  bandcol <- rep(fm$bands, kb)
  isW <- rep(c(TRUE, FALSE), each = length(bandcol))
  fsets <- lapply(fm$bands, function(b) bandcol == b); names(fsets) <- fm$bands
  if (!is.null(sets)) {
    if (is.null(names(sets))) stop("'sets' must be a named list.")
    if (length(unique(kb)) != 1L) stop("'sets' needs the same features in every band.")
    for (nm in names(sets)) {
      s <- rep(FALSE, kb[1L]); s[sets[[nm]]] <- TRUE
      fsets[[nm]] <- rep(s, length(fm$bands))
    }
  }
  lev <- unique(strata)
  idx <- lapply(lev, function(s) which(strata == s))
  sizes <- vapply(idx, length, 1L)
  n_arr <- prod(factorial(sizes))
  if (is.null(exact)) exact <- n_arr <= max_exact
  if (exact && n_arr > 5e6) stop("Too many arrangements (", n_arr, ") for exact enumeration.")
  tol <- 1e-9
  if (exact) {
    P <- lapply(sizes, .mt_perms)
    grid <- as.matrix(expand.grid(lapply(P, function(p) seq_len(nrow(p)))))
  } else {
    set.seed(seed)
    P <- lapply(sizes, function(m) rbind(seq_len(m), t(replicate(n_perm, sample.int(m)))))
    grid <- matrix(seq_len(n_perm + 1L), n_perm + 1L, length(lev))
  }
  B <- nrow(grid)
  Tall <- NULL; obs <- NULL
  for (j in seq_len(ncol(Zx))) {
    tot <- 0
    for (s in seq_along(lev)) {
      r <- idx[[s]]; Ps <- P[[s]]
      zmat <- matrix(Zx[r, j][Ps], nrow(Ps), length(r))
      tot <- tot + (zmat %*% Zy[r, , drop = FALSE])[grid[, s], , drop = FALSE]
    }
    S <- tot > tol
    for (nm in names(fsets)) {
      cw <- which(isW)[fsets[[nm]]]; cm <- which(!isW)[fsets[[nm]]]
      dw <- 100 * rowSums(S[, cw, drop = FALSE]) / length(cw)
      dm <- 100 * rowSums(S[, cm, drop = FALSE]) / length(cm)
      Tall <- cbind(Tall, round(dw - dm, 10))
      obs <- rbind(obs, data.frame(covariate = colnames(Zx)[j], band = nm, n_features = length(cw),
                                   DC_wPLI = dw[1L], DC_MSC = dm[1L], stringsAsFactors = FALSE))
    }
  }
  A <- abs(Tall)
  cnt <- function(v, t0) if (exact) mean(v >= t0 - tol) else (sum(v[-1L] >= t0 - tol) + 1) / B
  obs$T <- Tall[1L, ]
  obs$p <- vapply(seq_len(ncol(A)), function(i) cnt(A[, i], A[1L, i]), 0)
  obs$p_min <- vapply(seq_len(ncol(A)), function(i) max(cnt(A[, i], 100), 1 / B), 0)
  mx <- apply(A, 1L, max)
  obs$p_family <- vapply(A[1L, ], function(t0) cnt(mx, t0), 0)
  structure(list(table = obs, n_perm = if (exact) B else B - 1L, exact = exact, method = method,
                 strata = stats::setNames(sizes, lev), T_perm = Tall),
            class = "metricT_trend")
}

#' @export
print.metricT_trend <- function(x, digits = 3, ...) {
  cat("Metric-T for a continuous covariate (", x$method, " association)\n", sep = "")
  cat("Strata: ", paste(sprintf("%s (n=%d)", names(x$strata), x$strata), collapse = ", "), "\n", sep = "")
  if (x$exact) cat("Arrangements: all ", x$n_perm, " enumerated within strata (exact, two-sided)\n\n", sep = "")
  else cat("Arrangements: ", x$n_perm, " sampled within strata (two-sided)\n\n", sep = "")
  t <- x$table
  num <- vapply(t, is.numeric, TRUE)
  t[num] <- lapply(t[num], function(v) round(v, digits))
  print(t, row.names = FALSE)
  cat("\np_family: adjusted over all rows by the maximum |T| of each arrangement.\n")
  invisible(x)
}
