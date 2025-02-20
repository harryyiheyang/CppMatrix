#' Solve Sylvester Equation using Armadillo
#'
#' This function solves the Sylvester equation AX + XB = C using Armadillo's built-in solver.
#'
#' @param A A square matrix.
#' @param B A square matrix.
#' @param C A square matrix.
#' @return The solution matrix X.
#' @export
matrixSylvesterEigen <- function(A, B, C) {
  x <- .Call(`_CppMatrix_matrixSylvester`, A, B, C)
  return(x)
}

