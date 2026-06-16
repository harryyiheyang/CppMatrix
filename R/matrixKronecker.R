#' Kronecker Product
#'
#' This function computes the Kronecker product of two matrices using RcppArmadillo.
#'
#' @param A The first matrix.
#' @param B The second matrix.
#' @return The Kronecker product of matrices A and B.
#' @export
matrixKronecker <- function(A, B) {
  A <- .as_matrix_if_needed(A)
  B <- .as_matrix_if_needed(B)
  return(invisible(matrixKronecker_cpp(A, B)))
}
