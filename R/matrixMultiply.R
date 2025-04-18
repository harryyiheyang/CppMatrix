#' Matrix Multiplication
#'
#' This function performs matrix multiplication using RcppArmadillo.
#'
#' @param A The first matrix.
#' @param B The second matrix.
#' @return The product of matrices A and B.
#' @export
matrixMultiply <- function(A, B) {
  return(invisible(.Call(`_CppMatrix_matrixMultiply`, A, B)))
}
