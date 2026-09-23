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
- pgen_cor: Computes variant correlations directly from PLINK 2 PGEN files.
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

## Input requirements for `bed_cor()` and `pgen_cor()`

Both functions compute signed correlations of **REF-allele dosages** over all samples and all variants in the given files. The package does not select samples or variants and never chooses or flips alleles, so prepare each file upstream: subset samples and variants, and make the allele that defines each variant's direction REF.

```sh
# BED (REF = BIM column 6, A2)
plink2 --bfile in --extract snps.txt --keep samples.txt \
  --ref-allele force target.txt 1 2 --make-bed --out out

# PGEN (multiallelic variants allowed, e.g. an ALT2 as the target)
plink2 --pfile in --extract snps.txt --keep samples.txt \
  --ref-allele force target.txt 1 2 --make-pgen --out out
```

`target.txt` has the allele to become REF in column 1 and the variant ID in column 2. For PGEN, the target becomes REF, the old REF becomes ALT1, and the remaining ALTs follow.

- `A` and `B` must contain the same samples in the same order.
- Missing calls are filled with each variant's median genotype; constant variants give `NA`.
- Rows and columns follow the variant order in the `.bim` / `.pvar` files.
- PGEN files must be hardcall-only with an uncompressed `.pvar` (`plink2 --make-pgen` output).

```R
R    <- bed_cor("gene.bed", threads = 16)           # or pgen_cor("gene.pgen")
Rkk  <- Rold[keep, keep, drop = FALSE]
Rka  <- bed_cor("keep.bed", "add.bed", threads = 16)
Raa  <- bed_cor("add.bed", threads = 16)
Rnew <- rbind(cbind(Rkk, Rka), cbind(t(Rka), Raa))
```

## Performance

The functions in this package are implemented using Rcpp and RcppArmadillo, which allow for efficient computation in C++. This can lead to significant performance improvements compared to equivalent implementations in base R, especially for large matrices.

## License

This package is licensed under the MIT License, except `src/pgenlib/` (unmodified pgenlib sources from [plink-ng](https://github.com/chrchang/plink-ng) by Christopher Chang), which is LGPL (>= 3); `src/simde/` holds MIT-licensed SIMDe headers used on non-x86 platforms.

## Contact

Yihe Yang
Email: yxy1234@case.edu
