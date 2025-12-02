#' Singular Value Decomposition
#'
#' This function performs singular value decomposition using RcppArmadillo.
#'
#' @param A The matrix to decompose.
#' @return A list containing d (singular values), u (left singular vectors), and v (right singular vectors).
#' @export
matrixSVD <- function(A) {
  result <- .Call(`_CppMatrix_matrixSVD`, A)
  return(result)
}
