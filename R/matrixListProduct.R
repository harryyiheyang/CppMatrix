#' Product of a List of Matrices
#'
#' This function computes the product of a list of matrices in the order they appear using RcppArmadillo.
#'
#' @param matrixList A list of numeric matrices.
#' @return The product of the input matrices.
#' @export
matrixListProduct <- function(matrixList) {
  matrixList <- lapply(matrixList, .as_matrix_if_needed)
  return(invisible(matrixListProduct_cpp(matrixList)))
}
