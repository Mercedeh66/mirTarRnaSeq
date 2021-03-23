## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#'  miRanComp comparison of mRNAs present with miRanda file targets
#'
#' This function generates a dataframe consisting of mRNA or miRNAs present in
#'  miRanda generated file using the miRTarRNASeq:::getInputSpecies() function
#' @param miRNA Matrix or data.frame miRNA/RNA file or transformed diff expression file (generated
#' using TZtranz)
#' @param miRanda A dataframe of miRanda file with miRNA$V1 and miRNA targets miRNA$V2
#' @return An miRNA expression dataframe which includes only Genes/Targets present in miRanda file
#' @export
#' @keywords miRanda expression Intersection
#' @examples
#' x <- miRanComp(miRNA, miRanda)
miRanComp <- function(miRNA, miRanda) {
  subset <- intersect(rownames(miRNA), miRanda$V2)
  kept <- miRNA[subset, ]
  return(kept)
}
