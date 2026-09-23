# bed_cor(A, B = NULL, threads = 4L): signed correlations of REF (BIM A2)
# dosages of all variants in BED A (x B); missing calls take the variant's
# lower median. Files are prepared upstream (see README).
bed_cor <- function(A, B = NULL, threads = 4L) {
  A <- .bed_cor_path(A, "A")
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1 || threads != round(threads)) {
    stop("threads must be a single positive integer.", call. = FALSE)
  }
  threads <- as.integer(threads)
  self <- is.null(B)
  if (!self) B <- .bed_cor_path(B, "B")
  bed_cor_cpp(A, if (self) A else B, self, threads)
}

.bed_cor_path <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(arg, " must be one BED file path.", call. = FALSE)
  }
  if (!grepl("\\.bed$", x)) x <- paste0(x, ".bed")
  if (!file.exists(x)) stop("BED file not found: ", x, call. = FALSE)
  normalizePath(x, winslash = "/")
}
