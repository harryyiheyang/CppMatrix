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
- matrixMultiply: Performs matrix multiplication.
- matrixVectorMultiply: Performs matrix-vector multiplication.
- matrixCor: Computes the sample correlation matrix from a data matrix.
- matrixEigen: Computes the eigenvalue decomposition of a symmetric matrix.
- matrixKronecker: Computes the Kronecker product of two matrices.

## Performance

The functions in this package are implemented using Rcpp and RcppArmadillo, which allow for efficient computation in C++. This can lead to significant performance improvements compared to equivalent implementations in base R, especially for large matrices.

## License

This package is licensed under the MIT License.

## Contact

Yihe Yang
Email: yxy1234@case.edu
