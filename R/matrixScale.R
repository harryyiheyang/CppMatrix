#' Scale a Matrix
#'
#' This function centers and standardizes matrix columns using either classical or winsorized scale estimates.
#'
#' @param X A numeric matrix.
#' @param center Logical; if TRUE, center each column.
#' @param standardized Logical; if TRUE, standardize each column.
#' @param robust Logical; if TRUE, center by column medians and scale by winsorized standard deviations.
#' @return A scaled matrix with scaled:center and scaled:scale attributes when requested.
#' @export
matrixScale <- function(X, center = TRUE, standardized = TRUE, robust = FALSE) {
  X <- .as_matrix_if_needed(X)

  center <- .as_logical_flag(center, "center")
  standardized <- .as_logical_flag(standardized, "standardized")
  robust <- .as_logical_flag(robust, "robust")

  fit <- matrixScale_cpp(X, center, standardized, robust)
  result <- fit$x
  dimnames(result) <- dimnames(X)

  if (center) {
    attr(result, "scaled:center") <- .with_column_names(fit$center, X)
  }

  if (standardized) {
    attr(result, "scaled:scale") <- .with_column_names(fit$scale, X)
  }

  result
}
