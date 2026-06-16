#' Solve Sylvester Equation using Eigen Decomposition
#'
#' This function solves the Sylvester equation AX + XB = C using precomputed eigen decompositions of A and B.
#'
#' @param UA Eigenvector matrix of A.
#' @param DA Eigenvalues of A.
#' @param UB Eigenvector matrix of B.
#' @param DB Eigenvalues of B.
#' @param C A square matrix.
#' @return The solution matrix X.
#' @export
matrixSylvesterEigen <- function(UA, DA, UB, DB, C) {
  UA <- .as_matrix_if_needed(UA)
  UB <- .as_matrix_if_needed(UB)
  C <- .as_matrix_if_needed(C)
  x <- matrixSylvesterEigen_cpp(UA, DA, UB, DB, C)
  return(x)
}
