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
  x <- .Call(`_CppMatrix_matrixSylvesterEigen`, UA, DA, UB, DB, C)
  return(x)
}
