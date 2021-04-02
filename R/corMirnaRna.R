## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats cor
NULL

## quiet concerns of R CMD check regarding unbound global variables (in dplyr::filter() calls)
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("V1", "V2"))
}

#' corMirnaRna correlation for miRNA and mRNA
#'
#' This function uses the output of one2OneRnaMiRNA and returns the correlation dataframe
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param method Default is "pearson" else use "kendall" or "spearman"
#' @return Correlation data.frame
#' @export
#' @keywords Correlation with miRanda, miRanda Threshold
#' @examples
#' x <- corMirnaRna(mRNA_fc, miRNA_fc, method = "spearman")
corMirnaRna <- function(mRNA, miRNA, method = "pearson") {
    tmRNA <- t(mRNA)
    tmiRNA <- t(miRNA)
    mycor <- cor(cbind(tmRNA, tmiRNA), method = method)
    mmycor <- reshape2::melt(mycor)
    names(mmycor) <- c("V1", "V2", "value")
    # Makes a square matrix the correlation should not be squared
    mycorf <- dplyr::filter(mmycor, !(V2 %in% rownames(miRNA)))
    mycorf <- dplyr::filter(mycorf, !(V1 %in% rownames(mRNA)))
    mycorf$V1 <- as.character(mycorf$V1)
    mycorf$V2 <- as.character(mycorf$V2)
    return(as.data.frame(mycorf))
}
