#' Generalized Inverse of a Matrix
#'
#' This function computes the generalized inverse (Moore-Penrose pseudoinverse) of a matrix using RcppArmadillo.
#'
#' @param A A numeric matrix.
#' @param tol Tolerance for singular values. Default is 1e-6.
#' @return The generalized inverse of the input matrix.
#' @export
matrixGeneralizedInverse <- function(A, tol = 5e-16) {
  A <- .as_matrix_if_needed(A)
  return(invisible(matrixGeneralizedInverse_cpp(A, tol)))
}
