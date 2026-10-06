#' Read a NumPy .npy file (base R, no Python needed)
#'
#' Supports little-endian float64/float32/int32/int64/bool and fixed-width
#' unicode arrays, which covers the files written by the Python reference
#' pipeline. Pickled object arrays are skipped with a warning.
#'
#' @param path Path to a .npy file.
#' @return An R array (or vector) with the NumPy shape, or NULL if unsupported.
#' @export
mt_read_npy <- function(path) {
  raw <- readBin(path, "raw", n = file.info(path)$size)
  if (!identical(raw[1:6], as.raw(c(0x93, utf8ToInt("NUMPY"))))) stop("Not a .npy file: ", path)
  major <- as.integer(raw[7])
  if (major == 1L) {
    hlen <- readBin(raw[9:10], "integer", size = 2L, signed = FALSE, endian = "little"); off <- 10L
  } else {
    hlen <- readBin(raw[9:12], "integer", size = 4L, endian = "little"); off <- 12L
  }
  header <- rawToChar(raw[(off + 1L):(off + hlen)])
  descr <- sub(".*'descr': *'([^']+)'.*", "\\1", header)
  fortran <- grepl("'fortran_order': *True", header)
  shp <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", header)
  shape <- as.integer(Filter(nzchar, strsplit(gsub(" ", "", shp), ",")[[1]]))
  n <- if (length(shape)) prod(shape) else 1L
  dat <- raw[(off + hlen + 1L):length(raw)]
  v <- switch(descr,
    "<f8" = readBin(dat, "double", n = n, size = 8L, endian = "little"),
    "<f4" = readBin(dat, "double", n = n, size = 4L, endian = "little"),
    "<i4" = readBin(dat, "integer", n = n, size = 4L, endian = "little"),
    "<i8" = {
      x <- readBin(dat, "integer", n = 2L * n, size = 4L, endian = "little")
      lo <- x[c(TRUE, FALSE)]; hi <- x[c(FALSE, TRUE)]
      hi * 4294967296 + ifelse(lo < 0, lo + 4294967296, lo)
    },
    "|b1" = as.logical(as.integer(dat[seq_len(n)])),
    NULL)
  if (is.null(v) && grepl("^[<|=]U[0-9]+$", descr)) {
    k <- as.integer(sub("^.U", "", descr))
    cp <- matrix(readBin(dat, "integer", n = n * k, size = 4L, endian = "little"), nrow = k)
    v <- apply(cp, 2L, function(z) intToUtf8(z[z > 0L]))
  }
  if (is.null(v)) {
    warning("Unsupported dtype '", descr, "' in ", basename(path), "; skipped.")
    return(NULL)
  }
  if (length(shape) <= 1L) return(v)
  if (fortran) return(array(v, dim = shape))
  aperm(array(v, dim = rev(shape)), rev(seq_along(shape)))
}

#' Read the feature archive written by the Python reference pipeline
#'
#' The published pipeline saves \code{features.npz} with arrays named
#' \code{<band>__wpli} and \code{<band>__msc} (subjects x channel pairs) and
#' \code{subj_group}. This reads the archive directly in R.
#'
#' @param path Path to a .npz file.
#' @return List with \code{features} and \code{group}, ready for \code{\link{metric_t}}.
#' @export
mt_read_npz <- function(path) {
  tmp <- tempfile("npz"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  files <- unzip(path, exdir = tmp)
  nm <- sub("\\.npy$", "", basename(files))
  arr <- suppressWarnings(lapply(files, mt_read_npy)); names(arr) <- nm
  if (!("subj_group" %in% nm)) stop("Archive has no 'subj_group' array.")
  keys <- grep("__(wpli|msc)$", nm, value = TRUE)
  bands <- unique(sub("__(wpli|msc)$", "", keys))
  if (!length(bands)) stop("No '<band>__wpli' / '<band>__msc' arrays found.")
  pref <- c("delta", "theta", "alpha", "beta", "gamma", "broadband")
  bands <- c(intersect(pref, bands), setdiff(bands, pref))
  feats <- lapply(bands, function(b) list(wpli = arr[[paste0(b, "__wpli")]],
                                          msc = arr[[paste0(b, "__msc")]]))
  names(feats) <- bands
  list(features = feats, group = as.character(arr$subj_group))
}

#' Write and read features as a plain CSV table
#'
#' Wide layout: one row per subject; columns \code{subject}, \code{group},
#' then one column per feature named \code{<band>__<metric>__<pair>}. Any
#' software that can compute wPLI and MSC can produce this file.
#'
#' @param features,group As for \code{\link{metric_t}}.
#' @param path CSV path.
#' @param subject Optional subject identifiers.
#' @export
mt_write_features_csv <- function(features, group, path, subject = NULL) {
  n <- length(group)
  if (is.null(subject)) subject <- rownames(features[[1]]$wpli)
  if (is.null(subject)) subject <- sprintf("sub-%03d", seq_len(n))
  cols <- list(subject = subject, group = group)
  for (b in names(features)) for (m in c("wpli", "msc")) {
    x <- as.matrix(features[[b]][[m]])
    pn <- colnames(x); if (is.null(pn)) pn <- seq_len(ncol(x))
    colnames(x) <- paste(b, m, pn, sep = "__")
    cols[[paste(b, m)]] <- x
  }
  out <- do.call(cbind.data.frame, c(unname(cols[1:2]), unname(cols[-(1:2)])))
  names(out)[1:2] <- c("subject", "group")
  write.csv(out, path, row.names = FALSE)
  invisible(path)
}

#' @rdname mt_write_features_csv
#' @return \code{mt_read_features_csv}: list with \code{features}, \code{group}, \code{subject}.
#' @export
mt_read_features_csv <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("subject", "group") %in% names(d))) stop("CSV needs 'subject' and 'group' columns.")
  fn <- setdiff(names(d), c("subject", "group"))
  parts <- strsplit(fn, "__", fixed = TRUE)
  if (any(lengths(parts) < 3L)) stop("Feature columns must be named <band>__<metric>__<pair>.")
  band <- vapply(parts, `[`, "", 1L); metric <- tolower(vapply(parts, `[`, "", 2L))
  pair <- vapply(parts, function(p) paste(p[-(1:2)], collapse = "__"), "")
  bands <- unique(band)
  feats <- lapply(bands, function(b) {
    g <- function(m) {
      ix <- which(band == b & metric == m)
      x <- as.matrix(d[, fn[ix], drop = FALSE]); colnames(x) <- pair[ix]; x
    }
    w <- g("wpli"); m <- g("msc")
    if (!ncol(w) || ncol(w) != ncol(m)) stop("Band '", b, "' needs matching wpli and msc columns.")
    list(wpli = w, msc = m[, colnames(w), drop = FALSE])
  })
  names(feats) <- bands
  list(features = feats, group = as.character(d$group), subject = as.character(d$subject))
}
