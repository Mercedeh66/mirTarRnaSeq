## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom grDevices colorRampPalette
#' @importFrom graphics par abline
#' @importFrom stats density
NULL

#' drawInterPlots for finInterResult miRNA and mRNA Interrelation real data
#'
#' This function draws miRNA, mRNA density plots for miRNA and mRNA Interrelation while
#' comparing in addition to overall FC_miRNA and FC_mRNA plots from the finInterResult dataframe function.
#' @param mrna mRNA results of twoTimePoint function.
#' @param mirna miRNA results of twoTimePoint function.
#' @param final_results finInterResult miRNA and mRNA interrelation in two timepoints  results in a dataframe.
#' @return par plots
#' @export
#' @keywords plot
#' @examples
#' x <- drawInterPlots(mRNA_fc2, miRNA_fc2, final_results)
drawInterPlots <- function(mrna, mirna, final_results) {
    par(mfrow = c(2, 2))

    # plot density of mRNA FCs
    plot(density(mrna$FC1))
    abline(v = 0, col = "grey80", lty = 2)

    # plot density of miRNA FCs
    plot(density(mirna$FC1))
    abline(v = 0, col = "grey80", lty = 2)

    # plot final results (i.e. p < 0.05) mRNA vs miRNA FCs
    # FC_mRNA and FC_miRNA are coloumns of final results
    with(final_results, plot(FC_mRNA, FC_miRNA))
}
