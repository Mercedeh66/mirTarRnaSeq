## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' Sparse Partial Correlations On mRNA/miRNA Expression
#' We make mirTarRnaSeq compatible to SPONGE package in order to estimate
#' sparse matrix correlation (using elstic net) for prediction potential miRNA-mRNA
#' interaction. Note this function/method is suggested for miRNA/mRNA interactions
#' in many samples with a notable variance of mRNA/miRNA expression. This model also
#' only reports negative sparse partial correlation predictions.
#' @param mirna_exp miRna expression data.frame with miRNA for rows and samples for columns
#' @param diff_exp mRNA expression data.frame with mRNA for rows and samples for columns
#' @param miranda_sponge_predict miRanda sponge compatible matrix produced by miranda_sponge_predict function
#' @param non_null  The default for this parameter is TRUE, hence it returns only non-null estimated if
#' FALSE it would return all NULL and TRUE estimates.
#' @return matrix adjacency matrix with column names miRNA and row names mRNA
#' @keywords Sponge, mirTarRnaSeq, sparse partial correlation, ceRNA
#' @export
one2manySponge <- function(mirna_exp, diff_exp, miranda_sponge_predict,non_null=TRUE){
  # Attempt to load SPONGE
  if (!requireNamespace("SPONGE", quietly = TRUE)) {
    stop("SPONGE package unavailable")
  }
  mirna_exp_ad<-t(mirna_exp)
  mrna_exp_ad<-t(diff_exp)
  rowing_mRNA<-rownames(mrna_exp_ad)
  adequate_miRNA<- mirna_exp_ad[rownames(mirna_exp_ad) %in% rowing_mRNA, ]
  new_meth_mRNA <- mrna_exp_ad[ order(rownames(mrna_exp_ad)), ]
  new_meth_miRNA<-adequate_miRNA[ order(rownames(adequate_miRNA)), ]
  sponging<-SPONGE::sponge_gene_miRNA_interaction_filter(
    gene_expr = new_meth_mRNA,
    mir_expr = new_meth_miRNA,
    mir_predicted_targets = miranda_sponge_predict
  )
  if (non_null) {
    sponging <- sponging[!sapply(sponging, is.null)]
  }
  return(sponging)
}
