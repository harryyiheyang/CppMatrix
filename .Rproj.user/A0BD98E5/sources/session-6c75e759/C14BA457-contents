#' Matrix-Vector Multiplication
#'
#' This function performs matrix-vector multiplication using RcppArmadillo.
#'
#' @param A The matrix.
#' @param b The vector.
#' @return The product of matrix A and vector b.
#' @export
matrixVectorMultiply <- function(A, b) {
  invisible(.Call(`_CppMatrix_matrixVectorMultiply`, A, b))
}
