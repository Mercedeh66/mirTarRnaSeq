## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#'  geneVari Makes a list of gene names to be used in the runModels function
#'
#' This function defines the boudnaries of mRNA vs miRNAs of interest to be analysed by the runModels function
#' @param Combined the combined file for mRNA and selected miRNAs output of combiner function
#' @param miRNA_select The vector of selected miRNA/s
#' @return A vector of characters with defined mRNA dimensions
#' @export
#' @keywords dimensions
#' @examples
#' x <- geneVari(Combine, "ebv-mir-bart9-5p")
geneVari <- function(Combined, miRNA_select) {
    # geneVari1 <- colnames(Combined[, 1:(ncol(Combined) - length(miRNA_select))])
    geneVari1 <- setdiff(colnames(Combined), miRNA_select)
    return(geneVari1)
}
