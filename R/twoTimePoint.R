## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' twoTimePoint miRNA and mRNA interrelation in two timepoints
#'
#' This function uses the output of one2OneRnaMiRNA and retruns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @return miRNA mRNA interelation dataframe
#' @export
#' @keywords mRNA miRNA interelation
#' @examples
#' \donttest{
#' x <- twoTimePoint(mRNA, miRNA)
#' }
twoTimePoint <- function(mRNA, miRNA) {
  assertthat::assert_that(ncol(mRNA) == 1 && ncol(miRNA) == 1)
  cc <- matrix(0, nrow = nrow(miRNA), ncol = nrow(mRNA))
  rownames(cc) <- rownames(miRNA)
  colnames(cc) <- rownames(mRNA)
  for (i in 1:nrow(miRNA)) {
    for (j in 1:nrow(mRNA)) {
      cc[i, j] <- abs(miRNA[i, 1] - mRNA[j, 1])
    }
  }
  mmycc <- reshape2::melt(cc)
  names(mmycc) <- c("V1", "V2", "value")
  mmycc$V1 <- as.character(mmycc$V1)
  mmycc$V2 <- as.character(mmycc$V2)
  return(as.data.frame(mmycc))
}
