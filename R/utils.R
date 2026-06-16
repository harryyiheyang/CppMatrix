.as_matrix_if_needed <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }

  as.matrix(x)
}

.as_logical_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }

  x
}

.with_column_names <- function(x, matrix) {
  x <- as.vector(x)
  names(x) <- colnames(matrix)
  x
}
