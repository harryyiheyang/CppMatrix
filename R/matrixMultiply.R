#' Matrix Multiplication
#'
#' This function performs matrix multiplication using RcppArmadillo.
#'
#' @param A The first matrix.
#' @param B The second matrix.
#' @param transA Logical; if TRUE, use the transpose of A.
#' @param transB Logical; if TRUE, use the transpose of B.
#' @return The product of matrices A and B.
#' @export
matrixMultiply <- function(A, B, transA = FALSE, transB = FALSE) {
  A <- .as_matrix_if_needed(A)
  B <- .as_matrix_if_needed(B)
  transA <- .as_logical_flag(transA, "transA")
  transB <- .as_logical_flag(transB, "transB")
  return(invisible(matrixMultiply_cpp(A, B, transA, transB)))
}
