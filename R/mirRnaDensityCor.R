## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats density quantile
#' @importFrom graphics lines abline
NULL

#' mirRnaDensityCor for miRTarRNASeq miRNA and mRNA correlation real data versus sampled data
#'
#' This function draws density plots for miRNA and mRNA correlation while
#' comparing real data vs sampled data. It mainly illustrates the where the lower %5 (sig)
#' relationships lie.
#' @param corr0 data.frame results of corMirnaRna function.
#' @param outs data.frame results from the sampCorRnaMirna function.
#' @param pvalue The p value threshold to be used on the data density plot default is 0.05.
#' @return Density plot
#' @export
#' @keywords Density plot
#' @examples
#' x <- mirRnaDensityCor(corr_0, outs, pvalue = 0.05)

mirRnaDensityCor <- function(corr0, corrS, pvalue = 0.05) {
  dens_corr_0 <- density(corr0$value)
  dens_outs <- density(corrS) # this is the distribution for correlations for shuffled data
  mx_y <- max(dens_corr_0$y, dens_outs$y)

  plot(
    dens_outs,
    xlim = c(-1.5, 1.5),
    ylim = c(0, mx_y),
    col = "grey80",
    lwd = 2
  )
  lines(dens_corr_0, col = "red", lwd = 2) # this is what we got for our original data
  abline(
    v = c(-1, 1),
    col = "grey90",
    lty = 2
  ) # indicate -1 and 1 in the plot.
  # find a threshold based on the shuffled data
  threshold <- quantile(corrS, pvalue) # 5% if else user can define
  abline(v = threshold, col = "blue", lty = 2) # add theshold line to density plots
}
