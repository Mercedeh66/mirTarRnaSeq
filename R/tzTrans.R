## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#'  tzTransTranspose and z-score transformation
#'
#' Transposes and z-score transforms a matrix or data.frame.
#' @param x matrix of miRNA or mRNA or the data frame to be transformed
#' @return transposed and transformed version of x as a matrix.
#' @export
#' @keywords zscore, scale, transform, transpose
#' @examples
#' x <- tzTrans(miRNA)
tzTrans <- function(x) {
  tx <- t(x)
  tnorm1 <- apply(tx, 1, scale)
  col1 <- colnames(tx)
  rownames(tnorm1) <- col1
  tnorm <- as.data.frame(tnorm1)
  return(tnorm)
}
