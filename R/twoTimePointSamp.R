## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' twoTimePointSamp miRNA and mRNA interrelation in two timepoints sampling
#'
#' This function uses the output of one2OneRnaMiRNA and returns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from fold changes (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from fold changes (FC) obj of the one2OneRnaMiRNA.
#' @param Shrounds number of shuffling over the FC data, default is 100.
#' @param Srounds number of sampling from the shuffled data, default is 1000.
#' @return miRNA mRNA interrelation dataframe
#' @export
#' @keywords sampling, sampling, correlation, shuffling
#' @examples
#' \donttest{
#' x <- twoTimePointSamp(mRNA, miRNA, Shrounds = 10, Srounds = 10)
#' }
twoTimePointSamp <- function(mRNA, miRNA, Shrounds = 100, Srounds = 1000) {
    outs <- c()
    for (i in seq_len(Shrounds)) {
        print(i)
        shuffled_mrna <- mRNA
        shuffled_mirna <- miRNA
        for (col in seq_len(ncol(shuffled_mrna))) {
            shuffled_mrna[, col] <- sample(shuffled_mrna[, col]) # this shuffles all values in column _col_
            shuffled_mirna[, col] <- sample(shuffled_mirna[, col]) # this shuffles all values in column _col_
        }
        cc <- twoTimePoint(shuffled_mrna, shuffled_mirna)
        outs <- c(outs, sample(cc$value, Srounds, replace = TRUE)) # take a sample of the correlations and add to _outs_
    }
    return(outs)
}
