## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#'  miRanComp comparison of mRNAs present with miRanda file targets
#'
#' This function genrates a dataframe consisting of mRNA or miRNAs present in
#'  miRanda genrated file using the miRTarRNASeq:::getInputSpecies() function
#' @param fl1 Matrix or data.frame mRNA/RNA from transforemed diff expression file (genrated
#' using TZtranz)
#' @param mirandadf A dataframe of miRanda file with miRNA$V1 and miRNA targets miRNA$V2
#' @return An mRNA expression dataframe which includes only Genes/Targets present in miRanda file
#' @export
#' @keywords miRanda expression Intersection
#' @examples
#' \donttest{
#' x <- miRanComp(TZmRNA, mirandadf)
#' }
miRanComp <- function(TZmRNA, mirandadf) {
  subset <- intersect(rownames(TZmRNA), mirandadf$V2)
  kept <- TZmRNA[subset, ]
  return(kept)
}
