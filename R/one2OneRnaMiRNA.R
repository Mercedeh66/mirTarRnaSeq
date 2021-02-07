## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom dplyr filter
NULL

## quiet concerns of R CMD check regarding unbound global variables (in dplyr::filter() calls)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("pvalue"))
}

## Inputs mRNA or miRNA dif expression file

#' one2OneRnaMiRNA correlation for miRNA and mRNA  using differntial expression
#' fold change and if/when available p-value
#'
#' This function inputs accept a list of data.frames and returns an obj with
#' two dataframes calledd FC and p-value. FC with rownames == genes and columns are FC1, 2, 3, ...
#' (with fold-changes) - P-value with rownames == genes and columns are P1, 2, 3, ... (with p-values)
#' both data.frames have the same order dimensions.
#' @param files a list of dataframes either miRNAs or mRNAs from various time points.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param gene_colname Default is a vector character of length 1 "Gene" user can alter if they choose
#' This coloumn contains the gene names.
#' @param fc_colname Default "FC" is coloumn name for fold changes user can alter if they choose.
#' @param pval_colname Default is "pvalue" column name for p-values (in input).
#' @return Correlation dataframe
#' @export
#' @keywords correlation
#' @examples
#' \donttest{
#' x <- corAnmiRNAmRNA(mRNA, miRNA, Cor = -0.9, getInputSpeciesDF)
#' }
#'
one2OneRnaMiRNA <- function(files, gene_colname = "Gene", fc_colname = "FC",
                            pval_colname = "pvalue", pthreshold = NULL) {


  # asser file......
  # assert is.character(gene_colname) && length(gene_colname) == 1
  # assert is.character(fc_colname) && length(fc_colname) == 1
  # assert (is.numeric(pthreshold) && length(pthreshold) == 1) || is.null(pthreshold)
  # assert ... gene col is.character! not factor.

  # do some mininal cleaning here...
  files <- lapply(files, function(x) {
    x <- as.data.frame(x)
    x[, gene_colname] <- as.character(x[, gene_colname]) # account for users...
    return(x)
  })

  # grab fold changes and gene names from each data.frame supplied
  foldchanges <- lapply(files, function(x) {
    ret <- x[, c(gene_colname, fc_colname)]
    rownames(ret) <- ret[, gene_colname, drop = T]
    return(ret)
  })

  if (!is.null(pthreshold)) {
    # grab p-values and gene names from each data.frame supplied
    pvalues <- lapply(files, function(x) {
      ret <- x[, c(gene_colname, pval_colname)]
      rownames(ret) <- ret[, gene_colname, drop = T]
      return(ret)
    })

    # produce set of genes significant in any of the gene lists
    # all genes concatenated in one long unique list
    siggenes <- unique(unlist(sapply(pvalues, function(x) {
      ret <- filter(x, pvalue < pthreshold)[, gene_colname, drop = T]
      return(ret)
    }))) # for some reason sapply doesn't simplify here... hence the 'unlist'
  } else {
    # get set of genes regardless of p-value
    pvalues <- NULL
    siggenes <- unique(unlist(sapply(foldchanges, function(x) {
      x[, gene_colname, drop = T] # we're producing a long vector here...
    })))
  }

  # filter siggenes for gene names present in all input lists!
  for (x in foldchanges) {
    siggenes <- intersect(siggenes, x[, gene_colname, drop = T]) # we're intersecting _vectors_ here
  }

  # filter and reorder FCs by sig genes
  foldchanges <- lapply(foldchanges, function(x) {
    return(x[siggenes, fc_colname, drop = F])
  })

  # make one data frame with _just_ the FC columns.
  foldchanges <- do.call(cbind.data.frame, foldchanges)
  colnames(foldchanges) <- paste("FC", 1:ncol(foldchanges), sep = "") # cleanup column names

  if (!is.null(pthreshold)) {
    # filter and reorder pvalues by sig genes
    pvalues <- lapply(pvalues, function(x) {
      return(x[siggenes, pval_colname, drop = F])
    })

    # make one data frame with _just_ the pvalue columns.
    pvalues <- do.call(cbind.data.frame, pvalues)
    colnames(pvalues) <- paste("P", 1:ncol(pvalues), sep = "") # cleanup column names
  }

  # now 'foldchanges' and 'pvalues' are _data.frame_ with _just_ FC/P columns.
  return(list(foldchanges = foldchanges, pvalues = pvalues))
}
