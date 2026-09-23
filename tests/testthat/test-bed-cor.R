write_test_bed <- function(path, X, ids = paste0("s", seq_len(nrow(X)))) {
  base <- sub("\\.bed$", "", path)
  write.table(cbind(ids, ids, 0, 0, 0, -9), paste0(base, ".fam"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(cbind(1, paste0("v", seq_len(ncol(X))), 0, seq_len(ncol(X)), "A", "G"),
              paste0(base, ".bim"), quote = FALSE, row.names = FALSE,
              col.names = FALSE)
  codes <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    bits <- ifelse(is.na(x), 1L, codes[as.character(x)])
    bits <- c(bits, rep(0L, (-length(bits)) %% 4L))
    packed <- colSums(matrix(bits * rep(c(1L, 4L, 16L, 64L),
                                    length.out = length(bits)), nrow = 4L))
    writeBin(as.raw(packed), con)
  }
}

test_that("BED correlations match mean-imputed Pearson correlation", {
  X <- cbind(c(0, 1, 2, NA, 1), c(2, 1, 0, 1, NA),
             c(1, 1, 1, 1, 1), c(NA, NA, 0, NA, NA))
  Y <- cbind(c(2, 2, 1, 0, NA), c(NA, 0, 1, 2, 0))
  a <- tempfile(fileext = ".bed")
  b <- tempfile(fileext = ".bed")
  write_test_bed(a, X)
  write_test_bed(b, Y)
  fill <- function(M) {
    for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
    M
  }
  XA <- fill(X)
  YB <- fill(Y)
  expected <- matrix(NA_real_, ncol(X), ncol(X))
  expected[1:2, 1:2] <- cor(XA[, 1:2])
  expect_equal(unname(bed_cor(a)), unname(expected), tolerance = 1e-12)
  cross <- matrix(NA_real_, ncol(X), ncol(Y))
  cross[1:2, ] <- cor(XA[, 1:2], YB)
  expect_equal(unname(bed_cor(a, b)), unname(cross), tolerance = 1e-12)
})

test_that("BED rejects mismatched sample order and damaged files", {
  X <- cbind(c(0, 1, 2, 0, 1), c(2, 1, 0, 2, 1))
  a <- tempfile(fileext = ".bed")
  b <- tempfile(fileext = ".bed")
  write_test_bed(a, X)
  write_test_bed(b, X, rev(paste0("s", seq_len(nrow(X)))))
  expect_error(bed_cor(a, b), "same FID/IID")
  con <- file(b, "r+b")
  writeBin(as.raw(c(0, 0, 1)), con)
  close(con)
  expect_error(bed_cor(b), "SNP-major")
})
