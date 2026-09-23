#' Correlation from PLINK BED files
#'
#' Calculates Pearson correlations between variants in SNP-major PLINK BED
#' files. Genotypes count copies of BIM allele 1. Missing genotypes are replaced
#' with their variant's observed mean, and variants with fewer than two observed
#' calls or zero variance produce `NA` correlations. `A` and `B` must have the
#' same FID/IID samples in the same order.
#'
#' @param A Path to a `.bed` file with matching `.bim` and `.fam` files.
#' @param B Optional path to another `.bed` file with matching sidecars.
#' @return A numeric matrix with A variants in rows and B variants in columns.
#'   When `B` is `NULL`, returns the symmetric A-by-A correlation matrix.
#' @export
bed_cor <- function(A, B = NULL) {
  if (!is.character(A) || length(A) != 1L || is.na(A) || !nzchar(A)) {
    stop("A must be one BED file path.")
  }
  self <- is.null(B)
  if (!self && (!is.character(B) || length(B) != 1L || is.na(B) || !nzchar(B))) {
    stop("B must be NULL or one BED file path.")
  }
  bed_cor_cpp(A, if (self) A else B, self)
}
