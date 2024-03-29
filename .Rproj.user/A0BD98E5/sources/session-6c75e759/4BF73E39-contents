#' Matrix Correlation
#'
#' This function computes the sample correlation matrix from a data matrix using RcppArmadillo.
#'
#' @param A The data matrix.
#' @return The sample correlation matrix.
#' @export
matrixCor <- function(A) {
  invisible(.Call(`_CppMatrix_matrixCor`, A))
}
