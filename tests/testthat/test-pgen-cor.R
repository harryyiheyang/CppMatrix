test_that("pgen_cor selects the requested ALT at multiallelic sites", {
  skip_if_not_installed("pgenlibr")
  prefix <- file.path(system.file("extdata", package = "pgenlibr"),
                      "chr21_phase3_start")
  pvar <- pgenlibr::NewPvar(paste0(prefix, ".pvar.zst"))
  on.exit(pgenlibr::ClosePvar(pvar), add = TRUE)
  pgen <- pgenlibr::NewPgen(paste0(prefix, ".pgen"), pvar = pvar)
  on.exit(pgenlibr::ClosePgen(pgen), add = TRUE, after = FALSE)

  m <- pgenlibr::GetVariantCt(pvar)
  counts <- vapply(seq_len(m), function(i) pgenlibr::GetAlleleCt(pvar, i),
                   integer(1))
  multi <- which(counts > 2L)
  expect_true(length(multi) > 0L)
  expect_error(pgen_cor(prefix), "multiallelic variant")

  alt_A <- rep(NA_character_, m)
  alt_B <- alt_A
  alt_A[multi] <- vapply(multi, function(i)
    pgenlibr::GetAlleleCode(pvar, i, 3L), character(1))
  alt_B[multi] <- vapply(multi, function(i)
    pgenlibr::GetAlleleCode(pvar, i, 2L), character(1))

  R <- pgen_cor(prefix, alt_A = alt_A)
  expect_identical(dim(R), c(m, m))
  expect_equal(R, t(R))
  expect_equal(unname(diag(R)[!is.na(diag(R))]),
               rep(1, sum(!is.na(diag(R)))), tolerance = 1e-10)
  expect_true(anyNA(diag(R)))

  id <- pgenlibr::GetVariantId(pvar, multi[1L])
  C <- pgen_cor(prefix, prefix, alt_A = alt_A, alt_B = alt_B,
                snp = list(A = id, B = id))
  buf_A <- pgenlibr::Buf(pgen)
  buf_B <- pgenlibr::Buf(pgen)
  pgenlibr::Read(pgen, buf_A, multi[1L], 3L)
  pgenlibr::Read(pgen, buf_B, multi[1L], 2L)
  psam <- utils::read.table(paste0(prefix, ".psam"), header = TRUE,
                            check.names = FALSE, comment.char = "")
  founder <- psam$PAT == 0 & psam$MAT == 0
  expect_equal(C[id, id],
               unname(stats::cor(buf_A[founder], buf_B[founder],
                                 use = "pairwise.complete.obs")),
               tolerance = 1e-12)
  expect_gt(abs(C[id, id] - 1), 1e-3)
})

test_that("PGEN keep, SNP order, cross blocks, and incremental update agree", {
  skip_if_not_installed("pgenlibr")
  prefix <- file.path(system.file("extdata", package = "pgenlibr"),
                      "chr21_phase3_start")
  pvar <- pgenlibr::NewPvar(paste0(prefix, ".pvar.zst"))
  on.exit(pgenlibr::ClosePvar(pvar), add = TRUE)
  ids <- vapply(seq_len(5L), function(i) pgenlibr::GetVariantId(pvar, i),
                character(1))
  psam <- utils::read.table(paste0(prefix, ".psam"), header = TRUE,
                            check.names = FALSE, comment.char = "")
  keep <- tempfile()
  on.exit(unlink(keep), add = TRUE)
  writeLines(c("#IID", psam[["#IID"]][seq_len(60L)]), keep)

  R <- pgen_cor(prefix, snp = rev(ids), keep = keep)
  expect_identical(rownames(R), ids)
  expect_identical(colnames(R), ids)
  C <- pgen_cor(prefix, prefix,
                snp = list(A = rev(ids[1:2]), B = rev(ids[3:5])), keep = keep)
  expect_equal(C, R[ids[1:2], ids[3:5], drop = FALSE])
  old <- R[ids[1:2], ids[1:2], drop = FALSE]
  updated <- pgen_cor_update(old, prefix, c(rev(ids[3:5]), ids[1]),
                             keep = keep)
  expect_equal(updated, R)
})

test_that("PGEN sample IDs preserve leading zeros", {
  psam <- tempfile()
  keep <- tempfile()
  on.exit(unlink(c(psam, keep)), add = TRUE)
  writeLines(c("#IID\tPAT\tMAT\tSEX", "001\t0\t0\t1",
               "002\t0\t0\t2"), psam)
  writeLines(c("#IID", "001"), keep)
  samples <- CppMatrix:::.pgen_cor_samples(psam)
  selected <- CppMatrix:::.pgen_cor_sample_index(
    samples, CppMatrix:::.pgen_cor_keep(keep))
  expect_identical(samples$iid, c("001", "002"))
  expect_identical(selected, 1L)
})
