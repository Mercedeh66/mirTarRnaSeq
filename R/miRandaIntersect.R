## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' miRandaIntersect Looks for Intersection of Significant output results with miRanda Results from getInputSpeciesDF
#' function
#'
#' Compares and looks for intersection if significant output results with miRanda Results from getInputSpeciesDF and outputs a final
#' filterd ourput for only those pairs of miRNA and mRNA which have actually been predicted to be targets
#' in miRanda file
#' function
#' @param sig_corrs correlation matrix, produced by threshSig
#' @param corrS vector of correlations, from the sampCorRnaMirna function.
#' @param getInputSpeciesDF miranda data, produced by getInputSpecies
#' @return An object containing data.frames of significant mRNA, miRNA
#'             and correlation matrix filtered by miranda input.
#' @export
#' @keywords Signficance, Threshold, intersect
#' @examples
#' \donttest{
#' x <- miRandaIntersect(sig_corrs, miranda)
#' }
miRandaIntersect <- function(sig_corrs, corrS, mRNA, miRNA, getInputSpeciesDF) {
  result_corrs <- dplyr::inner_join(sig_corrs, getInputSpeciesDF, by = c("V1", "V2"))
  result_mrna <- mRNA[result_corrs$V2, , drop = F]
  result_mirna <- miRNA[result_corrs$V1, , drop = F]
  # calculate "p-values" for correlations. (TODO check this.)
  # count how many values in corrS are > each correlation and calculate 1 - "percentage".
  pvalueiods <- 1 - sapply(result_corrs$value, function(x) {
    length(which(corrS > x))
  }) / length(corrS)

  # add p-value for result_corrs df
  result_corrs$pvalue <- pvalueiods
  return(list(mirna = result_mirna, mrna = result_mrna, corrs = result_corrs))
}
