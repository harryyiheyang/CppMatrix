#' Solve Sylvester Equation using Armadillo
#'
#' This function solves the Sylvester equation AX + XB = C using Armadillo's built-in solver.
#'
#' @param A A square matrix.
#' @param B A square matrix.
#' @param C A square matrix.
#' @return The solution matrix X.
#' @export
matrixSylvester <- function(A, B, C) {
  A <- .as_matrix_if_needed(A)
  B <- .as_matrix_if_needed(B)
  C <- .as_matrix_if_needed(C)
  x <- matrixSylvester_cpp(A, B, C)
  return(x)
}

