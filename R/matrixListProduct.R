#' Product of a List of Matrices
#'
#' This function computes the product of a list of matrices in the order they appear using RcppArmadillo.
#'
#' @param matrixList A list of numeric matrices.
#' @return The product of the input matrices.
#' @export
matrixListProduct <- function(matrixList) {
  invisible(.Call(`_CppMatrix_matrixListProduct`, matrixList))
}
