#' Correlation of variants in PGEN files
#'
#' Compute signed, unphased dosage correlations directly from PLINK 2 PGEN
#' files. As in PLINK 2 `--r-unphased`, each pair uses founders with both
#' dosages observed; means and variances are recomputed on that pair's common
#' samples. A pair with fewer than two common samples or zero variance returns
#' `NA`. At least 50 founders are required, as with PLINK 2 LD. Results are
#' ordered by chromosome, position, then variant ID. The
#' corresponding `.pvar` (or `.pvar.zst`) and `.psam` files must be present.
#' Chromosome X is currently unsupported because PLINK 2 uses special
#' sex-aware weighting; calls selecting X error rather than returning a
#' nonmatching value.
#'
#' @param A Path to a `.pgen` file or its prefix.
#' @param B Optional path to another `.pgen` file or its prefix.
#' @param alt_A Target ALT strings for `A`: a vector in full-file order, a
#'   vector aligned with `snp` before sorting, or a named vector keyed by SNP
#'   ID. `NA` selects the sole ALT at a biallelic site; multiallelic sites
#'   require an explicit ALT.
#' @param alt_B As `alt_A`, for `B`. Must be `NULL` when `B` is omitted.
#' @param keep Optional path to a sample-ID file. Its first column is IID, or
#'   its first two columns are FID and IID. A header is optional. Matching
#'   samples retain PSAM order; nonfounders are then excluded.
#' @param snp Optional SNP IDs or path to a one-column SNP-ID file. For a
#'   cross matrix, use `list(A = ..., B = ...)` to select each side separately;
#'   `NULL` on a side selects all variants in that file.
#' @return A numeric matrix with variants from `A` in rows and variants from
#'   `B` in columns. Both axes have SNP-ID names and chromosome-position order.
#' @export
pgen_cor <- function(A, B = NULL, alt_A = NULL, alt_B = NULL,
                     keep = NULL, snp = NULL) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("pgen_cor() requires the pgenlibr package. Install it from CRAN.",
         call. = FALSE)
  }
  if (is.null(B) && !is.null(alt_B)) {
    stop("alt_B must be NULL when B is omitted.", call. = FALSE)
  }
  if (is.null(B) && is.list(snp)) {
    stop("snp must be one list of IDs or a file path when B is omitted.",
         call. = FALSE)
  }
  if (!is.null(B) && !is.null(snp) &&
      (!is.list(snp) || !identical(sort(names(snp)), c("A", "B")))) {
    stop("For a cross matrix, snp must be list(A = ..., B = ...).",
         call. = FALSE)
  }
  snp_A <- if (is.list(snp)) snp$A else snp
  snp_B <- if (is.list(snp)) snp$B else NULL

  files_A <- .pgen_cor_files(A)
  files_B <- if (is.null(B)) files_A else .pgen_cor_files(B)
  samples_A <- .pgen_cor_samples(files_A$psam)
  samples_B <- if (is.null(B)) samples_A else .pgen_cor_samples(files_B$psam)
  keep_ids <- .pgen_cor_keep(keep)
  sample_A <- .pgen_cor_sample_index(samples_A, keep_ids)
  sample_B <- if (is.null(B)) sample_A else
    .pgen_cor_sample_index(samples_B, keep_ids)
  if (!identical(samples_A$key[sample_A], samples_B$key[sample_B])) {
    stop("A and B must have identical selected founders in the same order.",
         call. = FALSE)
  }
  n <- length(sample_A)
  if (n < 50L) {
    stop("At least 50 founder samples are required, as in PLINK 2 LD.",
         call. = FALSE)
  }

  pvar_A <- pgenlibr::NewPvar(files_A$pvar)
  on.exit(pgenlibr::ClosePvar(pvar_A), add = TRUE)
  pgen_A <- pgenlibr::NewPgen(files_A$pgen, pvar = pvar_A,
                              sample_subset = sample_A)
  on.exit(pgenlibr::ClosePgen(pgen_A), add = TRUE, after = FALSE)
  if (is.null(B)) {
    pvar_B <- pvar_A
    pgen_B <- pgen_A
  } else {
    pvar_B <- pgenlibr::NewPvar(files_B$pvar)
    on.exit(pgenlibr::ClosePvar(pvar_B), add = TRUE, after = FALSE)
    pgen_B <- pgenlibr::NewPgen(files_B$pgen, pvar = pvar_B,
                                sample_subset = sample_B)
    on.exit(pgenlibr::ClosePgen(pgen_B), add = TRUE, after = FALSE)
  }

  if (pgenlibr::GetRawSampleCt(pgen_A) != length(samples_A$key) ||
      pgenlibr::GetRawSampleCt(pgen_B) != length(samples_B$key)) {
    stop("PGEN sample count does not match its PSAM file.", call. = FALSE)
  }
  if (pgenlibr::GetVariantCt(pvar_A) != pgenlibr::GetVariantCt(pgen_A) ||
      pgenlibr::GetVariantCt(pvar_B) != pgenlibr::GetVariantCt(pgen_B)) {
    stop("PGEN variant count does not match its PVAR file.", call. = FALSE)
  }
  variants_A <- .pgen_cor_variants(pvar_A, snp_A, alt_A, "alt_A")
  variants_B <- if (is.null(B)) variants_A else
    .pgen_cor_variants(pvar_B, snp_B, alt_B, "alt_B")
  chr <- toupper(sub("^chr", "", c(variants_A$chr, variants_B$chr),
                       ignore.case = TRUE))
  sex_A <- samples_A$sex[sample_A]
  sex_B <- samples_B$sex[sample_B]
  if ("X" %in% chr) {
    stop("Chromosome X is not supported by pgen_cor; PLINK 2 uses special sex-aware weighting.",
         call. = FALSE)
  }
  if ("Y" %in% chr && any(sex_A != "1" | sex_B != "1")) {
    stop("Chromosome Y with nonmale founders is not yet supported.",
         call. = FALSE)
  }
  m_A <- length(variants_A$index)
  m_B <- length(variants_B$index)
  if (m_A > 100000000 / m_B) {
    stop("Requested correlation matrix is too large; subset SNPs first.",
         call. = FALSE)
  }

  block <- max(1L, min(256L, floor(268435456 / (16 * n))))
  R <- matrix(NA_real_, m_A, m_B,
              dimnames = list(variants_A$id, variants_B$id))
  self <- is.null(B)
  for (j0 in seq.int(1L, m_B, by = block)) {
    jb <- j0:min(m_B, j0 + block - 1L)
    XB <- .pgen_cor_block(pgen_B, variants_B$index[jb],
                          variants_B$allele[jb], n)
    for (i0 in seq.int(1L, if (self) j0 else m_A, by = block)) {
      ib <- i0:min(m_A, i0 + block - 1L)
      XA <- if (self && i0 == j0) XB else
        .pgen_cor_block(pgen_A, variants_A$index[ib],
                        variants_A$allele[ib], n)
      V <- pgen_pair_cor_cpp(XA, XB, self && i0 == j0)
      R[ib, jb] <- V
      if (self && i0 != j0) R[jb, ib] <- t(V)
    }
  }
  R
}

#' Extend a PGEN correlation matrix with new variants
#'
#' `R_old` must use the same sample selection, founder rule, allele directions,
#' and missing-value rule as [pgen_cor()]. SNPs already in `R_old` are removed
#' from `snp`. The old/new blocks are assembled first; a single coordinate
#' permutation is then applied to both axes.
#'
#' @param R_old Existing square matrix with matching SNP-ID row and column names.
#' @param A PGEN path or prefix containing the old SNPs.
#' @param snp New SNP IDs or a one-column SNP-ID file.
#' @param B PGEN path or prefix containing the new SNPs. Defaults to `A`.
#' @param alt_A,alt_B Target ALT vectors as in [pgen_cor()]. When `B` is
#'   omitted, `alt_B` defaults to `alt_A`.
#' @param keep Sample-ID file as in [pgen_cor()].
#' @return A numeric correlation matrix ordered by chromosome, position, ID.
#' @export
pgen_cor_update <- function(R_old, A, snp, B = NULL,
                            alt_A = NULL, alt_B = NULL, keep = NULL) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("pgen_cor_update() requires the pgenlibr package.", call. = FALSE)
  }
  if (missing(snp) || is.null(snp)) {
    stop("snp must list the SNPs to add.", call. = FALSE)
  }
  if (!is.matrix(R_old) || !is.numeric(R_old) ||
      nrow(R_old) != ncol(R_old) || is.null(rownames(R_old)) ||
      !identical(rownames(R_old), colnames(R_old)) ||
      anyDuplicated(rownames(R_old))) {
    stop("R_old must be a square numeric matrix with identical unique SNP-ID axes.",
         call. = FALSE)
  }
  old <- rownames(R_old)
  add <- setdiff(.pgen_cor_snp_ids(snp), old)
  source_B <- if (is.null(B)) A else B
  target_B <- if (is.null(B) && is.null(alt_B)) alt_A else alt_B

  files_A <- .pgen_cor_files(A)
  pvar_A <- pgenlibr::NewPvar(files_A$pvar)
  on.exit(pgenlibr::ClosePvar(pvar_A), add = TRUE)
  old_info <- .pgen_cor_variants(pvar_A, old, alt_A, "alt_A")
  if (length(add)) {
    R_cross <- pgen_cor(A, source_B, alt_A = alt_A, alt_B = target_B,
                        keep = keep, snp = list(A = old, B = add))
    R_new <- pgen_cor(source_B, alt_A = target_B, keep = keep, snp = add)
    add <- rownames(R_new)
    R_cross <- R_cross[old, add, drop = FALSE]
    out <- rbind(cbind(R_old, R_cross),
                 cbind(t(R_cross), R_new))
    pvar_B <- pgenlibr::NewPvar(.pgen_cor_files(source_B)$pvar)
    on.exit(pgenlibr::ClosePvar(pvar_B), add = TRUE, after = FALSE)
    new_info <- .pgen_cor_variants(pvar_B, add, target_B, "alt_B")
    chr <- c(old_info$chr[match(old, old_info$id)],
             new_info$chr[match(add, new_info$id)])
    bp <- c(old_info$bp[match(old, old_info$id)],
            new_info$bp[match(add, new_info$id)])
  } else {
    out <- R_old
    chr <- old_info$chr[match(old, old_info$id)]
    bp <- old_info$bp[match(old, old_info$id)]
  }
  id <- rownames(out)
  o <- .pgen_cor_order(chr, bp, id)
  out[o, o, drop = FALSE]
}

.pgen_cor_files <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("A and B must be single PGEN file paths or prefixes.", call. = FALSE)
  }
  prefix <- sub("\\.pgen$", "", x, ignore.case = TRUE)
  pgen <- paste0(prefix, ".pgen")
  pvar <- paste0(prefix, ".pvar")
  if (!file.exists(pvar)) pvar <- paste0(pvar, ".zst")
  psam <- paste0(prefix, ".psam")
  for (path in c(pgen, pvar, psam)) {
    if (!file.exists(path)) stop("Required PGEN companion file is missing: ",
                                 path, call. = FALSE)
  }
  list(pgen = normalizePath(pgen, mustWork = TRUE),
       pvar = normalizePath(pvar, mustWork = TRUE),
       psam = normalizePath(psam, mustWork = TRUE))
}

.pgen_cor_samples <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^##", lines) & nzchar(trimws(lines))]
  if (length(lines) < 2L) {
    stop("PSAM must contain a header and at least one sample: ", path,
         call. = FALSE)
  }
  tab <- utils::read.table(text = paste(lines, collapse = "\n"),
                           header = TRUE, check.names = FALSE,
                           colClasses = "character", comment.char = "",
                           quote = "")
  iid_col <- if ("#IID" %in% names(tab)) "#IID" else "IID"
  if (!iid_col %in% names(tab)) {
    stop("PSAM is missing its IID column: ", path, call. = FALSE)
  }
  fid_col <- if ("#FID" %in% names(tab)) "#FID" else "FID"
  fid <- if (fid_col %in% names(tab)) as.character(tab[[fid_col]]) else
    rep("0", nrow(tab))
  iid <- as.character(tab[[iid_col]])
  pat <- if ("PAT" %in% names(tab)) as.character(tab$PAT) else
    rep("0", nrow(tab))
  mat <- if ("MAT" %in% names(tab)) as.character(tab$MAT) else
    rep("0", nrow(tab))
  sex <- if ("SEX" %in% names(tab)) as.character(tab$SEX) else
    rep("0", nrow(tab))
  sex[is.na(sex)] <- "0"
  if (anyNA(iid) || anyDuplicated(paste(fid, iid, sep = "\t"))) {
    stop("PSAM sample IDs must be present and unique: ", path,
         call. = FALSE)
  }
  list(fid = fid, iid = iid, sex = sex,
       key = paste(fid, iid, sep = "\t"),
       founder = pat %in% c("0", "NA") & mat %in% c("0", "NA"))
}

.pgen_cor_keep <- function(path) {
  if (is.null(path)) return(NULL)
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("keep must be an existing sample-ID file.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^##", lines)]
  if (!length(lines)) stop("keep file is empty.", call. = FALSE)
  tab <- utils::read.table(text = paste(lines, collapse = "\n"),
                           header = FALSE, colClasses = "character",
                           comment.char = "", quote = "", fill = TRUE)
  first <- as.character(tab[[1L]])
  if (first[1L] %in% c("#IID", "IID", "#FID", "FID")) {
    tab <- tab[-1L, , drop = FALSE]
  }
  if (!nrow(tab)) stop("keep file has no sample IDs.", call. = FALSE)
  if (ncol(tab) == 1L) {
    list(iid = as.character(tab[[1L]]), key = NULL)
  } else {
    list(iid = NULL, key = paste(tab[[1L]], tab[[2L]], sep = "\t"))
  }
}

.pgen_cor_sample_index <- function(samples, keep) {
  idx <- which(samples$founder)
  if (!is.null(keep)) {
    selected <- if (is.null(keep$key)) samples$iid %in% keep$iid else
      samples$key %in% keep$key
    idx <- idx[selected[idx]]
  }
  idx
}

.pgen_cor_snp_ids <- function(snp) {
  if (is.null(snp)) return(NULL)
  if (!is.character(snp) || anyNA(snp)) {
    stop("snp must contain SNP IDs or name a one-column ID file.",
         call. = FALSE)
  }
  if (length(snp) == 1L && file.exists(snp)) {
    snp <- scan(snp, what = character(), quiet = TRUE)
  }
  if (!length(snp) || any(!nzchar(snp)) || anyDuplicated(snp)) {
    stop("snp must contain unique, nonempty SNP IDs.", call. = FALSE)
  }
  snp
}

.pgen_cor_variants <- function(pvar, snp, alt, arg) {
  m <- pgenlibr::GetVariantCt(pvar)
  ids <- .pgen_cor_snp_ids(snp)
  if (is.null(ids)) {
    index <- seq_len(m)
    ids <- vapply(index, function(i) pgenlibr::GetVariantId(pvar, i),
                  character(1))
  } else {
    hit <- lapply(ids, function(id) pgenlibr::GetVariantsById(pvar, id))
    if (any(lengths(hit) != 1L)) {
      stop("Every selected SNP ID must occur exactly once in PVAR.",
           call. = FALSE)
    }
    index <- as.integer(unlist(hit, use.names = FALSE))
  }
  if (anyDuplicated(ids)) {
    stop("Selected PVAR SNP IDs must be unique.", call. = FALSE)
  }
  if (is.null(alt)) {
    targets <- rep(NA_character_, length(index))
  } else if (!is.character(alt)) {
    stop(arg, " must be a character vector.", call. = FALSE)
  } else if (!is.null(names(alt))) {
    if (anyDuplicated(names(alt)) || !all(ids %in% names(alt))) {
      stop(arg, " names must uniquely cover every selected SNP ID.",
           call. = FALSE)
    }
    targets <- unname(alt[ids])
  } else if (length(alt) == m) {
    targets <- alt[index]
  } else if (length(alt) == length(index)) {
    targets <- alt
  } else {
    stop(arg, " must align with the full PVAR, selected SNPs, or SNP names.",
         call. = FALSE)
  }
  out <- integer(length(index))
  for (j in seq_along(index)) {
    i <- index[j]
    k <- pgenlibr::GetAlleleCt(pvar, i)
    target <- targets[j]
    if (is.na(target)) {
      if (k > 2L) {
        stop(arg, " must specify the target ALT at multiallelic variant ",
             ids[j], ".", call. = FALSE)
      }
      out[j] <- 2L
    } else {
      choices <- vapply(seq.int(2L, k), function(j)
        pgenlibr::GetAlleleCode(pvar, i, j), character(1))
      hit <- match(target, choices)
      if (is.na(hit)) {
        stop(arg, " target is not an ALT allele of variant ", ids[j],
             ".", call. = FALSE)
      }
      out[j] <- hit + 1L
    }
  }
  chr <- vapply(index, function(i) as.character(pgenlibr::GetVariantChrom(pvar, i)),
                character(1))
  bp <- vapply(index, function(i) as.numeric(pgenlibr::GetVariantPos(pvar, i)),
               numeric(1))
  o <- .pgen_cor_order(chr, bp, ids)
  list(index = index[o], id = ids[o], allele = out[o],
       chr = chr[o], bp = bp[o])
}

.pgen_cor_order <- function(chr, bp, id) {
  label <- sub("^chr", "", chr, ignore.case = TRUE)
  cnum <- suppressWarnings(as.numeric(label))
  special <- match(toupper(label), c("X", "Y", "XY", "MT", "M"))
  cnum[is.na(cnum) & !is.na(special)] <- 22 + special[is.na(cnum) & !is.na(special)]
  order(is.na(cnum), cnum, label, bp, id, method = "radix")
}

.pgen_cor_block <- function(pgen, idx, allele, n) {
  X <- matrix(NA_real_, n, length(idx))
  buf <- pgenlibr::Buf(pgen)
  for (j in seq_along(idx)) {
    pgenlibr::Read(pgen, buf, idx[j], allele[j])
    X[, j] <- buf
  }
  X
}
