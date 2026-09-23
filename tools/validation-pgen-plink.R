library(CppMatrix)
library(pgenlibr)

src <- file.path(system.file("extdata", package = "pgenlibr"),
                 "chr21_phase3_start")
plink <- Sys.getenv("PLINK2")
if (!nzchar(plink) || !file.exists(plink)) stop("Set PLINK2 to plink2.exe")
work <- Sys.getenv("CPP_MATRIX_VALIDATION_DIR",
                   file.path(tempdir(), "cppmatrix-pgen-plink"))
dir.create(work, showWarnings = FALSE)

v <- NewPvar(paste0(src, ".pvar.zst"))
p <- NewPgen(paste0(src, ".pgen"), pvar = v)
m <- GetVariantCt(v)
n <- GetRawSampleCt(p)
ids <- vapply(seq_len(m), function(i) GetVariantId(v, i), character(1))
alts <- rep(NA_character_, m)
target <- character(m)
missing <- numeric(m)
variance <- numeric(m)
buf <- Buf(p)
for (i in seq_len(m)) {
  k <- GetAlleleCt(v, i)
  if (k > 2L) alts[i] <- GetAlleleCode(v, i, min(3L, k))
  target[i] <- if (k > 2L) alts[i] else GetAlleleCode(v, i, 2L)
  Read(p, buf, i, if (k > 2L) min(3L, k) else 2L)
  missing[i] <- mean(is.na(buf))
  variance[i] <- stats::var(buf[seq_len(min(1000L, n))], na.rm = TRUE)
}
ClosePgen(p)
ClosePvar(v)

priority <- unique(c(which(!is.na(alts)), order(missing, decreasing = TRUE)[1:25]))
chosen <- sort(c(priority, head(setdiff(seq_len(m), priority), 200L - length(priority))))
id <- ids[chosen]
alt <- alts[chosen]
names(alt) <- id
writeLines(id, file.path(work, "ids.txt"))
write.table(data.frame(ALT = target[chosen], ID = id),
            file.path(work, "ref.txt"), quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")
psam <- read.table(paste0(src, ".psam"), header = TRUE,
                   check.names = FALSE, comment.char = "")
nonfounder <- which(as.character(psam$PAT) != "0" |
                    as.character(psam$MAT) != "0")
keep_id <- unique(c(psam[[1]][seq_len(min(1000L, n))],
                    psam[[1]][nonfounder]))
cat("selected samples:", length(keep_id), "including nonfounders:",
    length(intersect(keep_id, psam[[1]][nonfounder])), "\n")
cat("multiallelic target ALT2+ variants:", sum(!is.na(alt)), "\n")
writeLines(c("#IID", keep_id), file.path(work, "keep.txt"))

run <- function(args) {
  result <- system2(plink, args, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) stop(paste(result, collapse = "\n"))
}

# The AoU reference: select IDs, force each chosen ALT to REF, then use
# PLINK2's signed unphased REF-based correlation.
ref <- file.path(work, "ref")
ld <- file.path(work, "ld")
run(c("--pfile", shQuote(src), "vzs", "--keep", shQuote(file.path(work, "keep.txt")),
      "--extract", shQuote(file.path(work, "ids.txt")),
      "--ref-allele", "force", shQuote(file.path(work, "ref.txt")), "1", "2",
      "--make-pgen", "--out", shQuote(ref)))
run(c("--pfile", shQuote(ref), "--r-unphased", "ref-based", "square", "bin4",
      "--out", shQuote(ld)))
vars <- scan(paste0(ld, ".unphased.vcor1.bin.vars"), what = character(), quiet = TRUE)
con <- file(paste0(ld, ".unphased.vcor1.bin"), "rb")
plink_R <- matrix(readBin(con, "numeric", n = length(vars)^2L, size = 4L),
                  length(vars), length(vars), dimnames = list(vars, vars))
close(con)

R <- pgen_cor(src, alt_A = alt, keep = file.path(work, "keep.txt"), snp = rev(id))
P <- plink_R[rownames(R), colnames(R)]
diff <- abs(R - P)
sign_flip <- sum(is.finite(R) & is.finite(P) & R * P < 0)
cat("self max abs diff:", max(diff, na.rm = TRUE), "\n")
cat("self sign flips:", sign_flip, "\n")
cat("self NA pattern equal:", identical(is.na(R), is.na(P)), "\n")
stopifnot(max(diff, na.rm = TRUE) <= 1e-6,
          sign_flip == 0L, identical(is.na(R), is.na(P)),
          identical(rownames(R), id))

old <- id[seq(1L, length(id), by = 2L)]
new <- id[seq(2L, length(id), by = 2L)]
C <- pgen_cor(src, src, alt_A = alt, alt_B = alt,
              keep = file.path(work, "keep.txt"),
              snp = list(A = rev(old), B = rev(new)))
PC <- plink_R[rownames(C), colnames(C)]
cat("cross max abs diff:", max(abs(C - PC), na.rm = TRUE), "\n")
cat("cross sign flips:", sum(is.finite(C) & is.finite(PC) & C * PC < 0), "\n")
stopifnot(max(abs(C - PC), na.rm = TRUE) <= 1e-6)

R_old <- pgen_cor(src, alt_A = alt, keep = file.path(work, "keep.txt"), snp = old)
R_up <- pgen_cor_update(R_old, src, c(rev(new), old[1L]), alt_A = alt,
                        keep = file.path(work, "keep.txt"))
cat("incremental max abs diff:", max(abs(R_up - R), na.rm = TRUE), "\n")
cat("incremental identical axes:", identical(rownames(R_up), rownames(R)), "\n")
stopifnot(identical(rownames(R_up), rownames(R)),
          max(abs(R_up - R), na.rm = TRUE) <= 1e-12)

top <- head(order(missing[chosen], decreasing = TRUE), 5L)
report <- data.frame(SNP = id[top], missing_rate = missing[chosen][top],
                     max_abs_diff = vapply(id[top], function(x)
                       if (all(is.na(R[x, ] - P[x, ]))) NA_real_ else
                         max(abs(R[x, ] - P[x, ]), na.rm = TRUE), numeric(1)))
print(report, row.names = FALSE)

# This public fixture has no missing calls among selected SNPs. Introduce
# different missingness masks in five variants and repeat the PLINK comparison.
raw <- file.path(work, "raw")
run(c("--pfile", shQuote(src), "vzs", "--keep", shQuote(file.path(work, "keep.txt")),
      "--extract", shQuote(file.path(work, "ids.txt")),
      "--export", "vcf", "--out", shQuote(raw)))
vcf <- readLines(paste0(raw, ".vcf"))
record <- which(!grepl("^#", vcf))
rates <- c(0.50, 0.35, 0.20, 0.10, 0.05)
mutated <- head(which(variance[chosen] > 0), length(rates))
stopifnot(length(mutated) == length(rates))
for (j in seq_along(rates)) {
  fields <- strsplit(vcf[record[mutated[j]]], "\t", fixed = TRUE)[[1L]]
  count <- length(fields) - 9L
  set.seed(j)
  sample <- sample.int(count, floor(count * rates[j]))
  fields[9L + sample] <- "./."
  vcf[record[mutated[j]]] <- paste(fields, collapse = "\t")
}
writeLines(vcf, file.path(work, "missing.vcf"))
missing_src <- file.path(work, "missing")
missing_ref <- file.path(work, "missing-ref")
missing_ld <- file.path(work, "missing-ld")
run(c("--vcf", shQuote(file.path(work, "missing.vcf")), "--make-pgen",
      "--out", shQuote(missing_src)))
run(c("--pfile", shQuote(missing_src), "--ref-allele", "force",
      shQuote(file.path(work, "ref.txt")), "1", "2", "--make-pgen",
      "--out", shQuote(missing_ref)))
run(c("--pfile", shQuote(missing_ref), "--r-unphased", "ref-based", "square",
      "bin4", "--out", shQuote(missing_ld)))
mv <- scan(paste0(missing_ld, ".unphased.vcor1.bin.vars"),
           what = character(), quiet = TRUE)
con <- file(paste0(missing_ld, ".unphased.vcor1.bin"), "rb")
MP <- matrix(readBin(con, "numeric", n = length(mv)^2L, size = 4L),
             length(mv), length(mv), dimnames = list(mv, mv))
close(con)
MR <- pgen_cor(missing_src, alt_A = alt, snp = id)
MP <- MP[rownames(MR), colnames(MR)]
cat("missing self max abs diff:", max(abs(MR - MP), na.rm = TRUE), "\n")
cat("missing self sign flips:",
    sum(is.finite(MR) & is.finite(MP) & MR * MP < 0), "\n")
cat("missing self NA pattern equal:", identical(is.na(MR), is.na(MP)), "\n")
stopifnot(max(abs(MR - MP), na.rm = TRUE) <= 1e-6,
          identical(is.na(MR), is.na(MP)))
MC <- pgen_cor(missing_src, missing_src, alt_A = alt, alt_B = alt,
               snp = list(A = old, B = new))
cat("missing cross max abs diff:",
    max(abs(MC - MP[rownames(MC), colnames(MC)]), na.rm = TRUE), "\n")
cat("missing cross sign flips:",
    sum(is.finite(MC) & is.finite(MP[rownames(MC), colnames(MC)]) &
        MC * MP[rownames(MC), colnames(MC)] < 0), "\n")
stopifnot(max(abs(MC - MP[rownames(MC), colnames(MC)]),
              na.rm = TRUE) <= 1e-6)
missing_report <- data.frame(SNP = id[mutated],
                             missing_rate = rates,
                             max_abs_diff = vapply(id[mutated],
                               function(x) if (all(is.na(MR[x, ] - MP[x, ])))
                                 NA_real_ else max(abs(MR[x, ] - MP[x, ]),
                                                   na.rm = TRUE),
                               numeric(1)))
print(missing_report, row.names = FALSE)
MR_old <- pgen_cor(missing_src, alt_A = alt, snp = old)
MR_up <- pgen_cor_update(MR_old, missing_src, c(rev(new), old[1L]),
                         alt_A = alt)
cat("missing incremental max abs diff:",
    max(abs(MR_up - MR), na.rm = TRUE), "\n")
stopifnot(max(abs(MR_up - MR), na.rm = TRUE) <= 1e-12)
cat("validation directory:", work, "\n")
