## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats cor
NULL

## quiet concerns of R CMD check regarding unbound global variables (in dplyr::filter() calls)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("V1", "V2", "value"))
}

#' corMirnaRnaMiranda correlation for miRNA and mRNA
#'
#' This function uses the output of one2OneRnaMiRNA and returns the correlation dataframe.
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param CorVal Correlation cut off.Example: If correlation -0.2 it would only return correlations with
#' smaller than this value correlation for miRNA and mRNA at various time points.
#' @param getInputSpeciesDF The dataframe generated from the getInputSpecies function.
#' @param method Default is "pearson" else use "kendall" or "spearman".
#' @return Correlation dataframe
#' @export
#' @keywords Correlation with miRanda, miRanda Threshold
#' @examples
#' x <- corMirnaRnaMiranda(mRNA_fc, miRNA_fc, Cor = -0.9, miRandaM)
corMirnaRnaMiranda <- function(mRNA, miRNA, CorVal, getInputSpeciesDF, method = "pearson") {
  tmRNA <- t(mRNA)
  tmiRNA <- t(miRNA)
  mycor <- cor(cbind(tmRNA, tmiRNA), method = method)
  mmycor <- reshape2::melt(mycor)
  names(mmycor) <- c("V1", "V2", "value")
  mmycor$V1 <- as.character(mmycor$V1)
  mmycor$V2 <- as.character(mmycor$V2)
  mycorf <- dplyr::filter(mmycor, !(V2 %in% rownames(miRNA)))
  mycorf <- dplyr::filter(mycorf, !(V1 %in% rownames(mRNA)))
  fmycorf <- dplyr::filter(mycorf, value < CorVal)
  Lfmycorf <- dplyr::left_join(getInputSpeciesDF, fmycorf, by = c("V1", "V2"))
  FLfmycorf <- Lfmycorf[is.finite(Lfmycorf$value), ]
  return(as.data.frame(FLfmycorf))
}
