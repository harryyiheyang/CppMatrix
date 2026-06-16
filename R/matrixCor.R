#' Matrix Correlation
#'
#' This function computes the sample correlation matrix from a data matrix using RcppArmadillo.
#'
#' @param A The data matrix.
#' @return The sample correlation matrix.
#' @export
matrixCor <- function(A) {
  A <- .as_matrix_if_needed(A)
  return(invisible(matrixCor_cpp(A)))
}
