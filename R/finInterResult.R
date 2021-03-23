## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' finInterResult miRNA and mRNA interrelation in two-time points  results in a dataframe.
#'
#' This function uses the output of one2OneRnaMiRNA and returns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param results Results from mirandaIntersectInter
#' @return miRNA mRNA interelation dataframe
#' @export
#' @keywords Results dataframe
#' @examples
#' x <- finInterResult(results)
finInterResult <- function(results) {
  final_results <- data.frame(
    mRNA = rownames(results$mrna),
    miRNA = rownames(results$mirna),
    FC_mRNA = results$mrna$FC1,
    FC_miRNA = results$mirna$FC1,
    pvalue = results$corrs$pvalue
  )
}
