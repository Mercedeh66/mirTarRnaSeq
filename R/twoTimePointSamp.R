## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' twoTimePointSamp miRNA and mRNA interrelation in two timepoints sampling
#'
#' This function uses the output of one2OneRnaMiRNA and retruns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param Shrounds number of shufflings over the FC data, default is 100.
#' @param Srounds number of sampling from the shuffled data, default is 1000.
#' @return miRNA mRNA interelation dataframe
#' @export
#' @keywords sampling, sampling, correlation, shuffling
#' @examples
#' \donttest{
#' x <- twoTimePointSamp(mRNA, miRNA, Shrounds = 100, Srounds = 1000)
#' }
twoTimePointSamp <- function(mRNA, miRNA, Shrounds = 100, Srounds = 1000) {
  outs <- c()
  for (i in 1:Shrounds) {
    print(i)
    shuffled_mrna <- mRNA
    shuffled_mirna <- miRNA
    for (col in 1:ncol(shuffled_mrna)) {
      shuffled_mrna[, col] <- sample(shuffled_mrna[, col]) # this shuffles all values in column _col_
      shuffled_mirna[, col] <- sample(shuffled_mirna[, col]) # this shuffles all values in column _col_
    }
    cc <- twoTimePoint(shuffled_mrna, shuffled_mirna)
    outs <- c(outs, sample(cc$value, Srounds, replace = T)) # take a sample of the corralations and add to _outs_
  }
  return(outs)
}