## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats density quantile
#' @importFrom graphics lines abline
NULL

#' mirRnaDensityInter for mirTarRnaSeq miRNA and mRNA Interrelation real data versus sampled data
#'
#' This function draws density plots for miRNA and mRNA Interrelation while
#' comparing real data vs sampled data. It mainly illustrates the where the lower %5 (sig)
#' relationships lie.
#' @param Inter0 data.frame results of twoTimePoint function.
#' @param OUTS data.frame results from the twoTimePointSamp function.
#' @param pvalue The p value threshold to be used on the data density plot default is 0.05.
#' @return Density plot
#' @export
#' @keywords Density plot
#' @examples
#' \donttest{
#' x <- mirRnaDensityInter(Inter0, OUTS, pvalue = 0.05)
#' }
#'
mirRnaDensityInter <- function(Inter0, OUTS, pvalue = 0.05) {
  dens_Inter_0 <- density(Inter0$value)
  dens_outs <- density(OUTS) # this is the distribution for correlations for shuffled data
  mx_y <- max(dens_Inter_0$y, dens_outs$y)
  plot(
    dens_outs,
    ylim = c(0, mx_y),
    col = "grey80",
    lwd = 2
  )
  lines(dens_Inter_0, col = "red", lwd = 2) # this is what we got for our original data
  # find a threshold based on the shuffled data
  threshold <- quantile(OUTS, 1 - pvalue) # 5% if else user can define
  abline(v = threshold, col = "blue", lty = 2) # add theshold line to density plots
}
