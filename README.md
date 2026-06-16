# CppMatrix

The `CppMatrixFunction` package provides a collection of efficient matrix operation functions implemented using Rcpp and RcppArmadillo. These functions offer significant performance improvements compared to base R implementations, making it suitable for computationally intensive tasks involving large matrices.

## Installation

You can install the `CppMatrixFunction` package from GitHub using the `devtools` package:

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

## Performance

The functions in this package are implemented using Rcpp and RcppArmadillo, which allow for efficient computation in C++. This can lead to significant performance improvements compared to equivalent implementations in base R, especially for large matrices.

## License

This package is licensed under the MIT License.

## Contact

Yihe Yang
Email: yxy1234@case.edu
