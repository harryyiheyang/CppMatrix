#' Eigenvalue Decomposition
#'
#' This function computes the eigenvalue decomposition of a symmetric matrix using RcppArmadillo.
#'
#' @param A A symmetric matrix.
#' @return A list containing the eigenvalues and eigenvectors of matrix A.
#' @export
matrixEigen <- function(A) {
  A <- .as_matrix_if_needed(A)
  fit <- matrixEigen_cpp(A)
  fit$values <- c(fit$values)
  return(fit)
}
