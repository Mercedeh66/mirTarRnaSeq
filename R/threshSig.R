## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats quantile
NULL

#' threshSig Using shuffling threshold finds appropriate significant miRNA-mRNA correlation
#'
#' This function uses the sampCorRnaMirna shuffled output to determine an appropriate threshold
#' for significant mRNA and miRNA relationship of the dataset and shows all those with significant
#' relationships.
#' @param corr0 data.frame results of corMirnaRna function.
#' @param corrS vector of correlations, from the sampCorRnaMirna function.
#' @param pvalue The p value threshold to be used on the sampled data.
#' @return A dataframe of Significant mRNA and miRNA
#' @export
#' @keywords Significance, Threshold
#' @examples
#' \donttest{
#' x <- mirRnaHeatmap(outs, corr_0)
#' }
#'
threshSig <- function(corr0, corrS, pvalue = 0.05) {
  threshold <- quantile(corrS, pvalue) # 5% default
  sig_corrs <- corr0[corr0$value < threshold, ]
  return(sig_corrs)
}
