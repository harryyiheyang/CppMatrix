#' Solve a Linear System
#'
#' This function solves AX = B using RcppArmadillo.
#'
#' @param A A square numeric matrix.
#' @param B A numeric vector or matrix.
#' @return A matrix.
#' @export
matrixSolve <- function(A, B) {
  A <- .as_matrix_if_needed(A)
  B <- .as_matrix_if_needed(B)

  result <- matrixSolveMat_cpp(A, B)
  dimnames(result) <- list(colnames(A), colnames(B))
  result
}
