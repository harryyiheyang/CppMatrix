# CppMatrix

The `CppMatrix` package provides a collection of efficient matrix operation functions implemented using Rcpp and RcppArmadillo. These functions offer significant performance improvements compared to base R implementations, making it suitable for computationally intensive tasks involving large matrices.

## Installation

You can install the `CppMatrix` package from GitHub using the `devtools` package:

```R
devtools::install_github("harryyiheyang/CppMatrix")
```

## Functions

The package includes the following matrix operation functions:

- matrixInverse: Computes the inverse of a square matrix.
- matrixGeneralizedInverse: Computes the generalized inverse of a square matrix.
- matrixMultiply: Performs matrix multiplication, with optional transpose flags for either input.
- matrixListProduct: Performs matrix multiplication on a list of matrices sequentially.
- matrixVectorMultiply: Performs matrix-vector multiplication.
- matrixSolve: Solves a linear system AX = B and always treats B as a matrix.
- matrixScale: Centers and standardizes matrix columns, with an optional winsorized robust scale.
- matrixCor: Computes the sample correlation matrix from a data matrix.
- bed_cor: Computes variant correlations directly from PLINK BED files.
- pgen_cor: Computes variant correlations directly from PLINK 2 PGEN files, with explicit ALT selection for multiallelic variants (requires `pgenlibr`).
- pgen_cor_update: Extends an existing PGEN correlation matrix with new variants and restores chromosome-position order.
- matrixEigen: Computes the eigenvalue decomposition of a symmetric matrix.
- matrixSVD: Computes the singular values decomposition of a matrix.
- matrixKronecker: Computes the Kronecker product of two matrices.
- matrixSylvester: Solve the Sylvester equation: AX+XB-C=0
- matrixSylvesterEigen: Solve the Sylvester equation: AX+XB-C=0 but use the eigenvalue decomposition of A and B to reduce computational cost.

Most wrapper functions protect matrix arguments by converting non-matrix inputs with `as.matrix()` only when needed, so data frames and vectors can be passed to matrix-based routines more conveniently.

## Examples

```R
A <- matrix(1:6, nrow = 2)
B <- matrix(1:12, nrow = 4)

# A %*% t(B), without explicitly transposing B in R
matrixMultiply(A, B, transB = TRUE)

S <- matrix(c(3, 1, 1, 2), nrow = 2)
b <- c(1, 2)

# b is converted to an n x 1 matrix, and the result is returned as a matrix
matrixSolve(S, b)

X <- data.frame(x = c(1, 2, 100), y = c(4, 5, 6))

# Classical column centering and standardization
matrixScale(X)

# Robust centering and winsorized standardization
matrixScale(X, robust = TRUE)
```

For genotype files that already contain the desired variants and samples, `bed_cor(A)` returns a variant-by-variant correlation matrix. `bed_cor(A, B)` returns the correlations between variants in `A` and `B`. Each BED file needs matching `.bim` and `.fam` files, and cross-file calls require identical samples in the same order.

```R
Rkk <- Rold[keep, keep, drop = FALSE]
Rka <- bed_cor("keep.bed", "add.bed")
Raa <- bed_cor("add.bed")
Rnew <- rbind(cbind(Rkk, Rka), cbind(t(Rka), Raa))
```

`pgen_cor()` requires matching `.pvar` and `.psam` files. `keep` selects samples from a sample-ID file, and `snp` selects variants by ID (use `list(A = ..., B = ...)` for a cross matrix). Founders are used, and each pair is calculated over samples with both dosages observed, matching PLINK 2 `--r-unphased`. For a multiallelic variant, supply its target ALT string in `alt_A` or `alt_B`; these can be full-file vectors, vectors aligned with selected SNPs, or named by SNP ID. Biallelic variants default to ALT1. Both axes are sorted by chromosome, position, and SNP ID. Constant variants return `NA` correlations.

PLINK 2 applies special sex-aware weighting to chromosome X. For now, `pgen_cor()` reports an error for any selected X variant and for Y with nonmale founders, so it does not silently return a matrix that differs from PLINK 2.

```R
R_old <- pgen_cor("region.pgen", snp = old_ids, alt_A = target_alt, keep = "keep.txt")
R_new <- pgen_cor_update(R_old, "region.pgen", snp = new_ids,
                         alt_A = target_alt, keep = "keep.txt")
```

## Performance

The functions in this package are implemented using Rcpp and RcppArmadillo, which allow for efficient computation in C++. This can lead to significant performance improvements compared to equivalent implementations in base R, especially for large matrices.

## License

This package is licensed under the MIT License.

## Contact

Yihe Yang
Email: yxy1234@case.edu
