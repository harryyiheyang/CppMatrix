#' Matrix-Vector Multiplication
#'
#' This function performs matrix-vector multiplication using RcppArmadillo.
#'
#' @param A The matrix.
#' @param b The vector.
#' @return The product of matrix A and vector b.
#' @export
matrixVectorMultiply <- function(A, b) {
  A <- .as_matrix_if_needed(A)
  a <- matrixVectorMultiply_cpp(A, b)
  return(as.vector(a))
}
