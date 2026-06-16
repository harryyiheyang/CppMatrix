#' Singular Value Decomposition
#'
#' This function performs singular value decomposition using RcppArmadillo.
#'
#' @param A The matrix to decompose.
#' @return A list containing d (singular values), u (left singular vectors), and v (right singular vectors).
#' @export
matrixSVD <- function(A) {
  A <- .as_matrix_if_needed(A)
  result <- matrixSVD_cpp(A)
  return(result)
}
