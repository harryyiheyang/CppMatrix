# pgen_cor(A, B = NULL, threads = 4L): as bed_cor for hardcall-only PGEN
# files with a plain-text .pvar, read with the bundled pgenlib.
pgen_cor <- function(A, B = NULL, threads = 4L) {
  threads <- .pgen_cor_threads(threads)
  self <- is.null(B)
  files_A <- .pgen_cor_files(A)
  files_B <- if (self) files_A else .pgen_cor_files(B)
  samples <- .pgen_cor_samples(files_A$psam)
  if (!self && !identical(samples, .pgen_cor_samples(files_B$psam))) {
    stop("A and B must have the same samples in the same order.",
         call. = FALSE)
  }
  n <- length(samples)
  if (n < 2L) stop("PGEN correlation requires at least two samples.",
                   call. = FALSE)
  variants_A <- .pgen_cor_variants(files_A$pvar)
  variants_B <- if (self) variants_A else .pgen_cor_variants(files_B$pvar)
  m_A <- length(variants_A$id)
  m_B <- length(variants_B$id)
  if (m_A > 100000000 / m_B) {
    stop("Requested correlation matrix is too large; subset variants first.",
         call. = FALSE)
  }
  R <- pgen_cor_cpp(files_A$pgen, files_B$pgen, self, threads, n, m_A, m_B,
                    variants_A$allele_ct, variants_B$allele_ct)
  dimnames(R) <- list(variants_A$id, variants_B$id)
  R
}

.pgen_cor_files <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("A and B must be single PGEN file paths or prefixes.", call. = FALSE)
  }
  prefix <- sub("\\.pgen$", "", x, ignore.case = TRUE)
  pgen <- paste0(prefix, ".pgen")
  pvar <- paste0(prefix, ".pvar")
  psam <- paste0(prefix, ".psam")
  if (!file.exists(pvar) && file.exists(paste0(pvar, ".zst"))) {
    stop("Decompress the .pvar.zst first, e.g. plink2 --pfile ", prefix,
         " vzs --make-just-pvar --out ", prefix, call. = FALSE)
  }
  for (path in c(pgen, pvar, psam)) {
    if (!file.exists(path)) stop("Required PGEN companion file is missing: ",
                                 path, call. = FALSE)
  }
  list(pgen = normalizePath(pgen, mustWork = TRUE),
       pvar = normalizePath(pvar, mustWork = TRUE),
       psam = normalizePath(psam, mustWork = TRUE))
}

# First line that does not start with "##" and the number of lines up to and
# including it.
.pgen_cor_header <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return(list(line = NULL, skip = skip))
    skip <- skip + 1L
    if (!startsWith(line, "##")) return(list(line = line, skip = skip))
  }
}

# Only the named columns of a header-described text file, as character.
.pgen_cor_columns <- function(path, header, cols, keep, sep) {
  classes <- rep("NULL", length(cols))
  classes[cols %in% keep] <- "character"
  utils::read.table(path, sep = sep, skip = header$skip, header = FALSE,
                    colClasses = classes, col.names = cols,
                    comment.char = "", quote = "",
                    na.strings = character(0), check.names = FALSE)
}

.pgen_cor_samples <- function(path) {
  header <- .pgen_cor_header(path)
  if (is.null(header$line)) {
    stop("PSAM must contain a header and at least one sample: ", path,
         call. = FALSE)
  }
  cols <- strsplit(trimws(header$line), "[[:space:]]+")[[1L]]
  iid_col <- if ("#IID" %in% cols) "#IID" else "IID"
  if (!iid_col %in% cols) {
    stop("PSAM is missing its IID column: ", path, call. = FALSE)
  }
  fid_col <- if ("#FID" %in% cols) "#FID" else "FID"
  tab <- .pgen_cor_columns(path, header, cols, c(iid_col, fid_col), "")
  if (!nrow(tab)) {
    stop("PSAM must contain a header and at least one sample: ", path,
         call. = FALSE)
  }
  fid <- if (fid_col %in% cols) tab[[fid_col]] else rep("0", nrow(tab))
  paste(fid, tab[[iid_col]], sep = "\t")
}

# Variant IDs and allele counts (REF + ALTs) of a plain-text .pvar, in file
# order. "##" lines are skipped; the "#CHROM" header line names the columns.
.pgen_cor_variants <- function(path) {
  header <- .pgen_cor_header(path)
  if (is.null(header$line) || !startsWith(header$line, "#CHROM")) {
    stop("PVAR must have a #CHROM header line: ", path, call. = FALSE)
  }
  cols <- strsplit(sub("^#", "", header$line), "\t", fixed = TRUE)[[1L]]
  if (!all(c("ID", "ALT") %in% cols)) {
    stop("PVAR header is missing its ID or ALT column: ", path, call. = FALSE)
  }
  tab <- .pgen_cor_columns(path, header, cols, c("ID", "ALT"), "\t")
  if (!nrow(tab)) stop("PVAR file has no variants: ", path, call. = FALSE)
  alt <- tab[["ALT"]]
  n_alt <- nchar(alt) - nchar(gsub(",", "", alt, fixed = TRUE)) + 1L
  list(id = tab[["ID"]], allele_ct = as.integer(1L + n_alt))
}

.pgen_cor_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      !is.finite(threads) || threads < 1 || threads != round(threads) ||
      threads > .Machine$integer.max) {
    stop("threads must be a single positive integer.", call. = FALSE)
  }
  as.integer(threads)
}
