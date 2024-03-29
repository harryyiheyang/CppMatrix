#' Matrix Inversion
#'
#' This function computes the inverse of a square matrix using RcppArmadillo.
#'
#' @param A A square matrix.
#' @return The inverse of matrix A.
#' @export
matrixInverse <- function(A) {
  invisible(.Call(`_CppMatrix_matrixInverse`, A))
}
