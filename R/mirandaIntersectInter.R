## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' mirandaIntersectInter Looks for Intersection of Significant output results with miRanda Results from getInputSpeciesDF
#' function
#'
#' Compares and looks for intersection if significant output results with miRanda Results from getInputSpeciesDF and outputs a final
#' filterd ourput for only those pairs of miRNA and mRNA which have actually been predicted to be targets
#' in miRanda file
#' function
#' @param sig_corrs correlation matrix, produced by threshSig
#' @param corrS vector of Differences/Correlations, from the sampCorRnaMirna function.
#' @param mRNA mRNA FC matrix.
#' @param miRNA miRNA FC matrix.
#' @param getInputSpeciesDF miranda data, produced by getInputSpecies.
#' @return An object containing data.frames of significant mRNA, miRNA
#'             and correlation matrix filtered by miranda input.
#' @export
#' @keywords mirandaIntersectInter Threshold intersect miRanda
#' @examples
#' x <- mirandaIntersectInter(sig_InterR, outs2, mRNA_fc2, miRNA_fc2, miRandaM)
mirandaIntersectInter <- function(sig_corrs, corrS, mRNA, miRNA, getInputSpeciesDF) {
    result_corrs <- dplyr::inner_join(sig_corrs, getInputSpeciesDF, by = c("V1", "V2"))
    if (nrow(result_corrs) == 0) {
        stop("no common mRNA/miRNAs found.")
    }
    result_mrna <- mRNA[result_corrs$V2, , drop = FALSE]
    result_mirna <- miRNA[result_corrs$V1, , drop = FALSE]
    # count how many values in corrS are > each correlation and calculate 1 - "percentage".
    pvalueiods <- 1 - sapply(result_corrs$value, function(x) {
        length(which(corrS < x))
    }) / length(corrS)

    # add p-value for result_corrs df
    result_corrs$pvalue <- pvalueiods
    return(list(mirna = result_mirna, mrna = result_mrna, corrs = result_corrs))
}
