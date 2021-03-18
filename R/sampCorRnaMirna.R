## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' sampCorRnaMirna sampling for correlation for miRNA and mRNA
#'
#' This function uses the output of one2OneRnaMiRNA and returns a sampled from original file
#' correlation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from fold changes (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from fold changes (FC) obj of the one2OneRnaMiRNA.
#' @param method Default is "pearson" else use "kendall" or "spearman".
#' @param Shrounds number of shuffling over the FC data, default is 100.
#' @param Srounds number of sampling from the shuffled data, default is 1000.
#' @return Correlation data frame
#' @export
#' @keywords sampling, sampling, correlation, shuffling
#' @examples
#' x <- sampCorRnaMirna(mRNA_fc, miRNA_fc, method = "pearson", Shrounds = 10, Srounds = 10)

sampCorRnaMirna <- function(mRNA, miRNA, method = "pearson", Shrounds = 100, Srounds = 1000) {
  outs <- c()
  for (i in 1:Shrounds) {
    print(i)
    shuffled_mrna <- mRNA
    shuffled_mirna <- miRNA
    for (col in 1:ncol(shuffled_mrna)) {
      shuffled_mrna[, col] <- sample(shuffled_mrna[, col]) # this shuffles all values in column _col_
      shuffled_mirna[, col] <- sample(shuffled_mirna[, col]) # this shuffles all values in column _col_
    }
    cc <- corMirnaRna(shuffled_mrna, shuffled_mirna, method = method) # run correlation on shuffled data
    outs <- c(outs, sample(cc$value, Srounds, replace = TRUE)) # take a sample of the correlations and add to _outs_
  }
  return(outs)
}
