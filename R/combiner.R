## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#'  combiner combines the miRNA and mRNA files
#'
#' This function makes and intersection dataframe for mRNA and miRNA/s of interest to be tested.
#' @param mRNA Matrix or data.frame mRNA/RNA from transforemed diff expression file (genrated
#' using TZtranz)
#' @param miRNA Matrix or data frame miRNA from transforemed diff file (genrated
#' using TZtranz)
#' @param miRNA_select A vector of character's for miRNAs which the user is interested in investigating if glm is use 1 miRNA
#' should be input. If multivariate several miRNAs should be imported, same goes for interaction determination for miRNAs. Note
#' we do not recommend more than 3-4 miRNAs at a time for the latter cases.
#' @return A dataframe which includes only mRNAs and miRNA intersection for the next estimation geneVari output.
#' @export
#' @keywords miRNA_mRNA_Intersection
#' @examples
#' miRNA_select <- c("ebv-mir-bart9-5p")
#' x <- combiner(mRNA, miRNA, miRNA_select)
combiner <- function(mRNA, miRNA, miRNA_select) {
  sub_miRNA <- miRNA[miRNA_select, ]
  common_samples <- intersect(names(mRNA), names(sub_miRNA))
  combination <- as.data.frame(t(rbind(mRNA[, common_samples], sub_miRNA[, common_samples])))
  combination <- combination[, !apply(is.na(combination), 2, any)]
  return(combination)
}
