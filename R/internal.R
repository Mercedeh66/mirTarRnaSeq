## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom assertthat assert_that
#' @importFrom dplyr filter select inner_join left_join `%>%`
#' @importFrom methods is
#' @importFrom viridis inferno
#' @importFrom reshape2 dcast melt
#' @importFrom pheatmap pheatmap
#' @importFrom corrplot corrplot
#' @importFrom graphics par plot abline lines points mtext
#' @importFrom stats cor quantile lm anova p.adjust density coefficients pf formula
#' @importFrom utils read.table
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach `%dopar%` foreach
#' @importFrom grDevices colorRampPalette
NULL

#'  tzTransTranspose and z-score transformation
#'
#' Transposes and z-score transforms a matrix or data.frame.
#' @param x matrix or data.frame to be transformed
#' @return transposed and transformed version of x as a matrix.
#' @export
#' @keywords zscore, scale, transform, transpose
#' @examples
#' \donttest{
#' x <- tzTrans(matrix(runif(100), nrow = 10))
#' }
tzTrans <- function(x) {
  tx <- t(x)
  tnorm1 <- apply(tx, 1, scale)
  col1 <- colnames(tx)
  rownames(tnorm1) <- col1
  tnorm <- as.data.frame(tnorm1)
  return(tnorm)
}


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
#' \donttest{
#' miRNA_select <- c("let7")
#' x <- miRanComp(mRNA, miRNA, miRNA_select)
#' }
combiner <- function(mRNA, miRNA, miRNA_select) {
  sub_miRNA <- miRNA[miRNA_select, ]
  common_samples <- intersect(names(mRNA), names(sub_miRNA))
  combination <- as.data.frame(t(rbind(mRNA[, common_samples], sub_miRNA[, common_samples])))
  combination <- combination[, !apply(is.na(combination), 2, any)]
  return(combination)
}


#'  geneVari Makes a list of gene names to be used in the runModels function
#'
#' This function defines the boudnaries of mRNA vs miRNAs of interest to be analysed by the runModels function
#' @param Combined the combined file for mRNA and selected miRNAs output of combiner function
#' @param miRNA_select The vector of selected miRNA/s
#' @return A vector of characters with defined mRNA dimensions
#' @export
#' @keywords dimensions
#' @examples
#' \donttest{
#' x <- geneVari(Combined, miRNA_select)
#' }
#'
geneVari <- function(Combined, miRNA_select) {
  geneVari1 <- colnames(Combined[, 1:(ncol(Combined) - length(miRNA_select))])
  return(geneVari1)
}


#' makeFormulaRightSide makes right hand side of formula for model variables: vector of indep. variables
#'
#' This function make right hand side of formula for model variables: vector of indep. variables (i.e. miRNAs)
#' mode: 'multi' for simple, 'inter' for model with interactions returns a string in the form "~ a + b", or "~ a + b + a * b"
#' @param variables The vector created by miRNA_select
#' @return data.frame containing Miranda data
#' @examples
#' \donttest{
#' x <- makeFormulaRightSide(variables, mode = "multi")
#' }
#'
makeFormulaRightSide <- function(variables, mode = "multi") {
  # be sure to properly quote variable names: `varname`.
  if (is.null(mode)) {
    mode <- "multi"
  }
  replace <- !grepl("^`.*`$", variables) # all that don't look like `somevar`
  variables[replace] <- paste("`", variables[replace], "`", sep = "")
  rightside <- paste(variables, collapse = " + ") # "~ a + b"
  if ("inter" == mode) {
    for (i in 1:(length(variables) - 1)) {
      for (j in (i + 1):length(variables)) {
        rightside <- paste(rightside, paste(variables[i], variables[j], sep = " * "), sep = " + ") # "~ a + b + a * b"
      }
    }
  }
  return(rightside)
}

#' runModels Makes a list of gene names to be used in the runModels function
#'
#' This function defines the boudnaries of mRNA vs miRNAs of interest to be analysed by the runModels function
#' @param combination the combined file for mRNA and selected miRNAs output of combiner function
#' @param geneVariable the output of
#' @param select_miRNA The vector of miRNA/s to be investigated
#' @param mode the mode of analysis if more than one miRNA is being investigated multivariate "multi"
#' or covariate/interaction analysis "inter" is being used
#' @return A list of p vlaues, annova, and significance for each gene and the miRNA/s of interest
#' @export
#' @keywords runModels multivariate interaction glm
#' @examples
#' \donttest{
#' x <- runModels(combination, geneVariable, miRNA_selectedVec, mode = "multi")
#' }
#'
runModels <- function(combination, geneVari, select_miRNA, mode = "multi") {
  is_significant <- c()
  all_models <- c()
  for (x in geneVari) {
    # construct formula by combining x with makeFormulaRightSide()'s output
    model_formula <- formula(paste(sprintf("`%s`", x), makeFormulaRightSide(select_miRNA, mode), sep = " ~ "))
    modelln <- lm(model_formula, data = combination)
    a <- anova(modelln)$`Pr(>F)`
    a <- a[!is.na(a)]
    a <- any(a < 0.05)
    is_significant <- c(is_significant, a)
    all_models <- c(all_models, list(modelln))
  }
  aaa <- lapply(all_models, anova) # retrieve the anovas for each model
  pvalues <- sapply(all_models, function(x) {
    k <- length(x$coefficients) - 1 # Number of mirnas in model!
    n <- nrow(combination)
    rs <- summary(x)$r.squared
    p <- pf((rs / k) / ((1 - rs) / (n - 1 - k)), k, n - 1 - k, lower.tail = F) # p value for R-squared
    return(p)
  })
  names(is_significant) <- geneVari
  names(all_models) <- geneVari
  names(aaa) <- geneVari
  names(pvalues) <- geneVari
  return(
    list(
      aaa = aaa,
      is_significant = is_significant,
      all_models = all_models,
      geneVari = geneVari,
      pvalues = pvalues
    )
  )
}


#' fdrSig Ruturns FDR significant miRNA/mRNA predictions
#'
#' This function performs FDR correction on the p_values generated by the runModels function list.
#' @param RMObj The output of runModels
#' @param value The FDR value default is 0.1
#' @param method The p-value adjustment method default is fdr. It could be either of the following "holm", "hochberg",
#' "hommel","bonferroni", "BH", "BY", or "fdr".
#' @param select_miRNA The vector of miRNA/s to be investigated.
#' @return A list of FDR corrected p vlaues, annova, and significance for each gene and the miRNA/s of interest
#' @export
#' @keywords p_adjust, correction
#' @examples
#' \donttest{
#' x <- fdrSig(RMObj, value = 0.1, method = "fdr")
#' }
#' # F9 fdrSig ruturns FDR significant genes
#' ## Desc: The User Can Here Define fdrSig value and get the FDR sig Samples, AllModelsGLMs annova value ,sig_genes names, and pvalues for all sig genes. In put is
#' ## the run model file obj and a value for FDR correction
fdrSig <- function(RMObj, value = 0.05, method = "fdr") {
  is_significant <- p.adjust(RMObj$pvalues, method = method) < value
  aaa <- RMObj$aaa[is_significant]
  all_models <- RMObj$all_models[is_significant]
  sig_genes <- RMObj$geneVari[is_significant]
  pvalues <- RMObj$pvalues[is_significant]
  return(list(
    sanova = aaa,
    FDR_significant = is_significant,
    # FDR_Value= p.adjust(RMObj$pvalues[is_significant], method=method),
    FDR_Value = p.adjust(RMObj$pvalues, method = method)[is_significant],
    all_models = all_models,
    sig_genes = sig_genes,
    pvalues = pvalues
  ))
}

#' importMirandaFile Read internal Miranda file
#'
#' Reads internal Miranda file from extdata and returns it as a data.frame
#' @param fn filename
#' @return data.frame containing Miranda data
#' @examples
#' \donttest{
#' x <- importMirandaFile("Mouse_miRanda.txt")
#' }
importMirandaFile <- function(fn) {
  ret1 <- read.table(system.file("extdata", "miRandaPrepFiles", fn,
    package = "mirTarRnaSeq", mustWork = T
  ),
  as.is = TRUE
  )
  return(ret1)
}


#' Return Miranda data for a given species.
#'
#' Reads Miranda file for a given speicies and returns it as
#' a data.frame, thresholded by percent identity.
#' Header options are Score (threshold), Energy-Kcal/Mol(energy),
#' Subject-IdentityPercent(targetIden), Query-IdentityPercent (mirnaIden)
#' @param selection Species
#' @return data.frame with Miranda data.
#' @keywords miranda, species
#' @export
#' @examples
#' \donttest{
#' x <- getInputSpecies("Epstein_Barr", threshold = 60, energy = -170, targetIden = 19, mirnaIden = 19) # Default is threshold 60
#' }
getInputSpecies <- function(selection, threshold = 60, energy = NULL, targetIden = NULL, mirnaIden = NULL) {
  if (selection == "Human") {
    ret <- importMirandaFile("Human_miRanda.txt.gz")
  } else if (selection == "Mouse") {
    ret <- importMirandaFile("Mouse_miRanda.txt.gz")
  } else if (selection == "Epstein_Barr") {
    ret <- importMirandaFile("Epstein_Barr_miRanda.txt.gz")
  } else if (selection == "C.elegans") {
    ret <- importMirandaFile("C.elegans_miRanda.txt.gz")
  } else if (selection == "Drosophila") {
    ret <- importMirandaFile("Drosophila_miRanda.txt.gz")
  } else if (selection == "Cytomegalovirus") {
    ret <- importMirandaFile("CMV_miRanda.txt.gz")
  } else if (selection == "Kaposi_Sarcoma") {
    ret <- importMirandaFile("Kaposi_miRanda.txt.gz")
  } else if (selection == "Epstein_Barr_Human") {
    ret <- importMirandaFile("EBV_Human_miRanda.txt.gz")
  } else if (selection == "CMV_Human") {
    ret <- importMirandaFile("CMV_Human_miRanda.txt.gz")
  } else if (selection == "KSHV_Human") {
    ret <- importMirandaFile("KSHV_Human_miRanda.txt.gz")
  }
  ret <- dplyr::filter(ret, V3 >= threshold)
  if (!is.null(energy)) {
    ret <- dplyr::filter(ret, V4 <= energy)
  }
  if (!is.null(targetIden)) {
    ret <- dplyr::filter(ret, V5 >= targetIden)
  }
  if (!is.null(mirnaIden)) {
    ret <- dplyr::filter(ret, V6 >= mirnaIden)
  }
  ret <- ret %>% dplyr::select(V1, V2, V3, V4, V5, V6)
  ret <- unique(ret)
  return(ret)
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

#' corMirnaRna correlation for miRNA and mRNA
#'
#' This function uses the output of one2OneRnaMiRNA and retruns the correlation dataframe
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param method Default is "pearson" else use "kendall" or "spearman"
#' @return Correlation data.frame
#' @export
#' @keywords Correlation with miRanda, miRanda Threshold
#' @examples
#' \donttest{
#' x <- corMirnaRna(mRNA, miRNA, method = "spearman")
#' }
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


#' corMirnaRnaMiranda correlation for miRNA and mRNA
#'
#' This function uses the output of one2OneRnaMiRNA and retruns the correlation dataframe
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA
#' @param CorVal Correlation cut off.Example: If correlation -0.2 it would only return correlations with
#' smaller than this value correlation for miRNA and mRNA at various time points
#' @param getInputSpeciesDF The dataframe generated from the getInputSpecies function
#' @param method Default is "pearson" else use "kendall" or "spearman"
#' @return Correlation dataframe
#' @export
#' @keywords Correlation with miRanda, miRanda Threshold
#' @examples
#' \donttest{
#' x <- corMirnaRnaMiranda(mRNA, miRNA, Cor = -0.9, getInputSpeciesDF)
#' }
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


#' mirRnaHeatmap pheatmap for miRTarRNASeq miRNA and mRNA correlation
#'
#' This function draws pheatmaps for miRNA and mRNA correlation while
#' using default and pheatmap for all other parameters
#' @param finalF data.frame results of corMirnaRnaMiranda or corMirnaRna function
#' @param ... arguments passed onto pheatmap
#' @param upper_bound is the upper_bound of the correlation pheatmap scale
#'  default is zero user can set to values based on output of correlation result (value)
#' @param main is the title of the pheatmap
#' @param color default inferno(50) from the library viridis R base,
#'  R colorbrewer and viridis compatible
#' @param fontsize default is 7 user adjustable
#' @return pheatmap Obj
#' @export
#' @keywords heatmap, pheatmap, color, correlation plot,correlation_plot
#' @examples
#' \donttest{
#' x <- mirRnaHeatmap(finalF, upper_bound = -0.1, color = rainbow(50), fontsize = 10)
#' }
mirRnaHeatmap <- function(finalF, ..., upper_bound = 0,
                          main = "Default mRNA miRNA heatmap",
                          color = c(viridis::inferno(50),"grey90"), fontsize = 7) {
  dfinalF <- dcast(finalF, V1 ~ V2, fun.aggregate = mean)
  dfinalF[is.na(dfinalF)] <- upper_bound
  rownames(dfinalF) <- dfinalF$V1
  dfinalF <- dfinalF %>% dplyr::select(-V1)
  p <- pheatmap::pheatmap(dfinalF,
    color = color, fontsize = fontsize,
    main = main, ...
  )
  return(p)
}



#' mirRnaHeatmapDiff pheatmap for miRTarRNASeq miRNA and mRNA correlation
#'
#' This function draws pheatmaps for miRNA and mRNA correlation while
#' using default and pheatmap for all other parameters
#' @param finalF data.frame results of corMirnaRnaMiranda or corMirnaRna function
#' @param ... arguments passed onto pheatmap
#' @param upper_bound is the upper_bound of the correlation pheatmap scale
#'  default is zero user can set to values based on output of correlation result (value)
#' @param main is the title of the pheatmap
#' @param color default inferno(50) from the library viridis R base,
#'  R colorbrewer and viridis compatible
#' @param fontsize default is 7 user adjustable
#' @return pheatmap Obj
#' @export
#' @keywords heatmap, pheatmap, color, correlation plot,correlation_plot
#' @examples
#' \donttest{
#' x <- mirRnaHeatmapDiff(finalF, upper_bound = -0.1, color = rainbow(50), fontsize = 10)
#' }
mirRnaHeatmapDiff <- function(finalF, ..., upper_bound = 0,
                          main = "Default mRNA miRNA heatmap",
                          color = c("grey90",viridis::inferno(50)), fontsize = 7) {
  dfinalF <- dcast(finalF, V1 ~ V2, fun.aggregate = mean)
  dfinalF[is.na(dfinalF)] <- upper_bound
  rownames(dfinalF) <- dfinalF$V1
  dfinalF <- dfinalF %>% dplyr::select(-V1)
  p <- pheatmap::pheatmap(dfinalF,
                          color = color, fontsize = fontsize,
                          main = main, ...
  )
  return(p)
}


#' sampCorRnaMirna sampling for correlation for miRNA and mRNA
#'
#' This function uses the output of one2OneRnaMiRNA and retruns a sampled from orig file
#' correlation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param method Default is "pearson" else use "kendall" or "spearman".
#' @param Shrounds number of shufflings over the FC data, default is 100.
#' @param Srounds number of sampling from the shuffled data, default is 1000.
#' @return Correlation dataframe
#' @export
#' @keywords sampling, sampling, correlation, shuffling
#' @examples
#' \donttest{
#' x <- sampCorRnaMirna(mRNA, miRNA, method = "pearson", Shrounds = 100, Srounds = 1000)
#' }
sampCorRnaMirna <- function(mRNA, miRNA, method = "pearson", Shrounds = 100, Srounds = 1000) {
  outs <- c()
  for (i in 1:Shrounds) {
    print(i)
    shuffled_mrna <- mRNA
    shuffled_mirna <- miRNA
    for (col in 1:ncol(shuffled_mrna)) {
      shuffled_mrna[, col] <- sample(shuffled_mrna[, col]) # this shuffles all values in column _col_
      shuffled_mirna[, col] <- sample(shuffled_mirna[, col]) # this shuffles all values in column _col_
    }
    cc <- corMirnaRna(shuffled_mrna, shuffled_mirna, method = method) # run correlation on shuffled data
    outs <- c(outs, sample(cc$value, Srounds, replace = T)) # take a sample of the corralations and add to _outs_
  }
  return(outs)
}


#' mirRnaDensityCor for miRTarRNASeq miRNA and mRNA correlation real data versus sampled data
#'
#' This function draws density plots for miRNA and mRNA correlation while
#' comparing real data vs sampled data. It mainly illustrates the where the lower %5 (sig)
#' relationships lie.
#' @param corr0 data.frame results of corMirnaRna function.
#' @param corrS data.frame results from the sampCorRnaMirna function.
#' @param pvalue The p value threshold to be used on the data density plot default is 0.05.
#' @return Density plot
#' @export
#' @keywords Density plot
#' @examples
#' \donttest{
#' x <- mirRnaDensityCor(corr0, corrS, pvalue = 0.05)
#' }
#'
mirRnaDensityCor <- function(corr0, corrS, pvalue = 0.05) {
  dens_corr_0 <- density(corr_0$value)
  dens_outs <- density(corrS) # this is the distribution for correlations for shuffled data
  mx_y <- max(dens_corr_0$y, dens_outs$y)

  plot(
    dens_outs,
    xlim = c(-1.5, 1.5),
    ylim = c(0, mx_y),
    col = "grey80",
    lwd = 2
  )
  lines(dens_corr_0, col = "red", lwd = 2) # this is what we got for our original data
  abline(
    v = c(-1, 1),
    col = "grey90",
    lty = 2
  ) # indicate -1 and 1 in the plot.
  # find a threshold based on the shuffled data
  threshold <- quantile(corrS, pvalue) # 5% if else user can define
  abline(v = threshold, col = "blue", lty = 2) # add theshold line to density plots
}

#' threshSig Using shuffling threshold finds appropriate significant miRNA-mRNA correlation
#'
#' This function uses the sampCorRnaMirna shuffled output to determine an appropriate thershold
#' for significant mRNA and miRNA relationship of the dataset and shows all those with significant
#' relationships.
#' @param corr0 data.frame results of corMirnaRna function.
#' @param corrS vector of correlations, from the sampCorRnaMirna function.
#' @param pvalue The p value threshold to be used on the sampled data.
#' @return A dataframe of Significant mRNA and miRNA
#' @export
#' @keywords Signficance, Threshold
#' @examples
#' \donttest{
#' x <- mirRnaHeatmap(corrS, corr0)
#' }
#'
threshSig <- function(corr0, corrS, pvalue = 0.05) {
  threshold <- quantile(corrS, pvalue) # 5% default
  sig_corrs <- corr0[corr0$value < threshold, ]
  return(sig_corrs)
}

#' threshSigInter Using shuffling threshold finds appropriate significant miRNA-mRNA correlation
#'
#' This function uses the sampCorRnaMirna shuffled output to determine an appropriate thershold
#' for significant mRNA and miRNA relationship of the dataset and shows all those with significant
#' relationships.
#' @param corr0 data.frame results of corMirnaRna function.
#' @param corrS vector of correlations, from the sampCorRnaMirna function.
#' @param pvalue The p value threshold to be used on the sampled data.
#' @return A dataframe of Significant mRNA and miRNA
#' @export
#' @keywords Signficance, Threshold
#' @examples
#' \donttest{
#' x <- mirRnaHeatmap(corrS, corr0)
#' }
#'
threshSigInter <- function(corr0, corrS, pvalue = 0.05) {
  threshold <- quantile(corrS, 1 - pvalue) # 5% default
  sig_corrs <- corr0[corr0$value > threshold, ]
  return(sig_corrs)
}

#' miRandaIntersect Looks for Intersection of Significant output results with miRanda Results from getInputSpeciesDF
#' function
#'
#' Compares and looks for intersection if significant output results with miRanda Results from getInputSpeciesDF and outputs a final
#' filterd ourput for only those pairs of miRNA and mRNA which have actually been predicted to be targets
#' in miRanda file
#' function
#' @param sig_corrs correlation matrix, produced by threshSig
#' @param corrS vector of correlations, from the sampCorRnaMirna function.
#' @param getInputSpeciesDF miranda data, produced by getInputSpecies
#' @return An object containing data.frames of significant mRNA, miRNA
#'             and correlation matrix filtered by miranda input.
#' @export
#' @keywords Signficance, Threshold, intersect
#' @examples
#' \donttest{
#' x <- miRandaIntersect(sig_corrs, miranda)
#' }
miRandaIntersect <- function(sig_corrs, corrS, mRNA, miRNA, getInputSpeciesDF) {
  result_corrs <- dplyr::inner_join(sig_corrs, getInputSpeciesDF, by = c("V1", "V2"))
  result_mrna <- mRNA[result_corrs$V2, , drop = F]
  result_mirna <- miRNA[result_corrs$V1, , drop = F]
  # calculate "p-values" for correlations. (TODO check this.)
  # count how many values in corrS are > each correlation and calculate 1 - "percentage".
  pvalueiods <- 1 - sapply(result_corrs$value, function(x) {
    length(which(corrS > x))
  }) / length(corrS)

  # add p-value for result_corrs df
  result_corrs$pvalue <- pvalueiods
  return(list(mirna = result_mirna, mrna = result_mrna, corrs = result_corrs))
}

#' miRandaIntersectInter Looks for Intersection of Significant output results with miRanda Results from getInputSpeciesDF
#' function
#'
#' Compares and looks for intersection if significant output results with miRanda Results from getInputSpeciesDF and outputs a final
#' filterd ourput for only those pairs of miRNA and mRNA which have actually been predicted to be targets
#' in miRanda file
#' function
#' @param sig_corrs correlation matrix, produced by threshSig
#' @param corrS vector of correlations, from the sampCorRnaMirna function.
#' @param getInputSpeciesDF miranda data, produced by getInputSpecies
#' @return An object containing data.frames of significant mRNA, miRNA
#'             and correlation matrix filtered by miranda input.
#' @export
#' @keywords Signficance, Threshold, intersect
#' @examples
#' \donttest{
#' x <- miRandaIntersect(sig_corrs, miranda)
#' }
miRandaIntersectInter <- function(sig_corrs, corrS, mRNA, miRNA, getInputSpeciesDF) {
  result_corrs <- dplyr::inner_join(sig_corrs, getInputSpeciesDF, by = c("V1", "V2"))
  if(nrow(result_corrs) == 0) {
    stop("no common mRNA/miRNAs found.")
  }
  result_mrna <- mRNA[result_corrs$V2, , drop = F]
  result_mirna <- miRNA[result_corrs$V1, , drop = F]
  # calculate "p-values" for correlations. (TODO check this.)
  # count how many values in corrS are > each correlation and calculate 1 - "percentage".
  pvalueiods <- 1 - sapply(result_corrs$value, function(x) {
    length(which(corrS < x))
  }) / length(corrS)

  # add p-value for result_corrs df
  result_corrs$pvalue <- pvalueiods
  return(list(mirna = result_mirna, mrna = result_mrna, corrs = result_corrs))
}



#' Plot model
#'
#' Plot 2D description
#' @param model linear model
#' @export
#' @keywords plot, 2d, model
#' @examples
#' \donttest{
#' plot2d(model)
#' }
#' # AddPValueSigAswell
plotFit <- function(model) {
  m <- model$model
  m$Fitted <- model$fitted.values
  last <- ncol(m)
  names(m)[last] <- "Fitted values\nlm(model_formula)"
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(m[, c(1, last)],
    type = "n", main = formula.tools:::as.character.formula(model$terms),
    cex.main = 1
  )
  abline(lm(m[, last] ~ m[, 1]), lwd = 2, col = "red")
  points(m[, c(1, last)])

  k <- length(model$coefficients) - 1 # Number of mirnas in model!
  n <- length(model$fitted.values)
  rs <- summary(model)$r.squared
  Pvalue <- pf((rs / k) / ((1 - rs) / (n - 1 - k)), k, n - 1 - k, lower.tail = F) # p value for R-squared

  mtext(sprintf("P value: %.2g", Pvalue), line = .74)
  mtext(sprintf("R-squared: %.2f", summary(model)$r.squared), line = .01)
}

#' R-squared from model
#'
#' Obtrain R-squared from a linear model
#' @param model linear model
#' @export
#' @keywords linear model, r-squared
#' @examples
#' \donttest{
#' x <- modelRsquared(model)
#' }
modelRsquared <- function(model) {
  return(summary(model)$r.squared)
}

#' Plot residuals
#'
#' Plot residuals description
#' @param model linear model
#' @export
#' @keywords plot, residuals, model
#' @examples
#' \donttest{
#' plotResiduals(model)
#' }
plotResiduals <- function(model) {
  plot(model, 1)
}

#' plotTerms
#'
#' Plot terms description
#' @param model linear model
#' @export
#' @keywords plot, residuals, model
#' @examples
#' \donttest{
#' plotTerms(model)
#' }
plotTerms <- function(model) {
  m <- model$model
  plot(m, main = formula.tools:::as.character.formula(model$terms))
  k <- length(model$coefficients) - 1 # Number of mirnas in model!
  n <- length(model$fitted.values)
  rs <- summary(model)$r.squared
  Pvalue <- pf((rs / k) / ((1 - rs) / (n - 1 - k)), k, n - 1 - k, lower.tail = F) # p value for R-squared

  mtext(sprintf("P value of model fit: %.2g", Pvalue), line = .74)
  mtext(sprintf("R-squared of model fit: %.2f", summary(model)$r.squared), line = .01)
}
#' twoTimePoint miRNA and mRNA interrelation in two timepoints
#'
#' This function uses the output of one2OneRnaMiRNA and retruns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @return miRNA mRNA interelation dataframe
#' @export
#' @keywords mRNA miRNA interelation
#' @examples
#' \donttest{
#' x <- twoTimePoint(mRNA, miRNA)
#' }
twoTimePoint <- function(mRNA, miRNA) {
  assertthat::assert_that(ncol(mRNA) == 1 && ncol(miRNA) == 1)
  cc <- matrix(0, nrow = nrow(miRNA), ncol = nrow(mRNA))
  rownames(cc) <- rownames(miRNA)
  colnames(cc) <- rownames(mRNA)
  for (i in 1:nrow(miRNA)) {
    for (j in 1:nrow(mRNA)) {
      cc[i, j] <- abs(miRNA[i, 1] - mRNA[j, 1])
    }
  }
  mmycc <- reshape2::melt(cc)
  names(mmycc) <- c("V1", "V2", "value")
  mmycc$V1 <- as.character(mmycc$V1)
  mmycc$V2 <- as.character(mmycc$V2)
  return(as.data.frame(mmycc))
}

#' twoTimePointSamp miRNA and mRNA interrelation in two timepoints sampling
#'
#' This function uses the output of one2OneRnaMiRNA and retruns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param mRNA mRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param miRNA miRNA file generated from foldchanges (FC) obj of the one2OneRnaMiRNA.
#' @param Shrounds number of shufflings over the FC data, default is 100.
#' @param Srounds number of sampling from the shuffled data, default is 1000.
#' @return miRNA mRNA interelation dataframe
#' @export
#' @keywords sampling, sampling, correlation, shuffling
#' @examples
#' \donttest{
#' x <- twoTimePointSamp(mRNA, miRNA, Shrounds = 100, Srounds = 1000)
#' }
twoTimePointSamp <- function(mRNA, miRNA, Shrounds = 100, Srounds = 1000) {
  outs <- c()
  for (i in 1:Shrounds) {
    print(i)
    shuffled_mrna <- mRNA
    shuffled_mirna <- miRNA
    for (col in 1:ncol(shuffled_mrna)) {
      shuffled_mrna[, col] <- sample(shuffled_mrna[, col]) # this shuffles all values in column _col_
      shuffled_mirna[, col] <- sample(shuffled_mirna[, col]) # this shuffles all values in column _col_
    }
    cc <- twoTimePoint(shuffled_mrna, shuffled_mirna)
    outs <- c(outs, sample(cc$value, Srounds, replace = T)) # take a sample of the corralations and add to _outs_
  }
  return(outs)
}


#' mirRnaDensityInter for miRTarRNASeq miRNA and mRNA Interrelation real data versus sampled data
#'
#' This function draws density plots for miRNA and mRNA Interrelation while
#' comparing real data vs sampled data. It mainly illustrates the where the lower %5 (sig)
#' relationships lie.
#' @param Inter0 data.frame results of twoTimePoint function.
#' @param OUTS data.frame results from the twoTimePointSamp function.
#' @param pvalue The p value threshold to be used on the data density plot default is 0.05.
#' @return Density plot
#' @export
#' @keywords Density plot
#' @examples
#' \donttest{
#' x <- mirRnaDensityInter(Inter0, OUTS, pvalue = 0.05)
#' }
#'
mirRnaDensityInter <- function(Inter0, OUTS, pvalue = 0.05) {
  dens_Inter_0 <- density(Inter0$value)
  dens_outs <- density(OUTS) # this is the distribution for correlations for shuffled data
  mx_y <- max(dens_Inter_0$y, dens_outs$y)
  plot(
    dens_outs,
    ylim = c(0, mx_y),
    col = "grey80",
    lwd = 2
  )
  lines(dens_Inter_0, col = "red", lwd = 2) # this is what we got for our original data
  # find a threshold based on the shuffled data
  threshold <- quantile(OUTS, 1 - pvalue) # 5% if else user can define
  abline(v = threshold, col = "blue", lty = 2) # add theshold line to density plots
}

#' finInterResult miRNA and mRNA interrelation in two timepoints  results in a dataframe.
#'
#' This function uses the output of one2OneRnaMiRNA and retruns a sampled from orig file
#' interrelation dataframe depending on user sampling selection.
#' @param results Results from miRandaIntersectInter
#' @return miRNA mRNA interelation dataframe
#' @export
#' @keywords Results dataframe
#' @examples
#' \donttest{
#' x <- finInterResult(results)
#' }
#'
finInterResult <- function(results) {
  final_results <- data.frame(
    mRNA = rownames(results$mrna),
    miRNA = rownames(results$mirna),
    FC_mRNA = results$mrna$FC1,
    FC_miRNA = results$mirna$FC1,
    pvalue = results$corrs$pvalue
  )
}

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
#' \donttest{
#' x <- drawInterPlots(mrna, mirna, final_results)
#' }
#'
drawInterPlots <- function(mrna, mirna, final_results) {
  par(mfrow = c(2, 2))

  # plot density of mRNA FCs
  plot(density(mrna$FC1))
  abline(v = 0, col = "grey80", lty = 2)

  # plot density of miRNA FCs
  plot(density(mirna$FC1))
  abline(v = 0, col = "grey80", lty = 2)

  # plot final results (i.e. p < 0.05) mRNA vs miRNA FCs
  with(final_results, plot(FC_mRNA, FC_miRNA))
}

#' runAllMirnaModels run_model for all miRNAs
#'
#' This function runs the "run_model" function for all miRNAs and mRNA combinations of two and returns a
#' list with significant genes and FDR models
#' @param mirnas vector of unique miRNAs under investigation.
#' @param DiffExpmRNA differentially/expressed mRNAs expression file.
#' @param DiffExpmiRNA differentially/expressed miRNAs expression file.
#' @param miranda_data getInputSpecies output file ( use low filters).
#' @param nPar number of threads for parallel processing.
#' @param prob user defined ratio for miranda distibution for miranda score selection default is 0.90.
#' @param method finInterResult miRNA and mRNA interrelation in two timepoints  results in a dataframe.
#' @param fdr_cutoff cutoff for FDR selection default is 0.1.
#' @param all_coeff if true only models with negative coefficient will be selected if false at least one
#' negative coefficient should be in the model; default is TRUE .
#' @param mode model mode, default is Null, can be changed to "multi" and "inter".
#' @return List of run models
#' @export
#' @keywords runAllMirnaModels
#' @examples
#' \donttest{
#' x <- runAllMirnaModels(mirnas, DiffExpmRNA, DiffExpmiRNA, miranda_data,
#'   nPar = 4, prob = 0.90, fdr_cutoff = 0.1,
#'   method = "fdr", all_coeff = TRUE, mode = "multi"
#' )
#' }
#'
runAllMirnaModels <- function(mirnas, DiffExpmRNA, DiffExpmiRNA, miranda_data,
                                 nPar = 4, prob = 0.75, fdr_cutoff = 0.1, method = "fdr",
                                 all_coeff = FALSE, mode = NULL) {
  # LoadLibrary(dopar)
  registerDoParallel(nPar)

  # make unique
  unique_mirnas <- unique(mirnas)

  # generate combinations of 2 miRNAs
  if (!is.null(mode)) {
    unique_mirnas <- caTools::combs(unique_mirnas, 2)
    unique_mirnas <- lapply(1:nrow(unique_mirnas), function(x) unique_mirnas[x, ])
  }


  # run models; parallellize by miRNA
  ret <- foreach(mirna = unique_mirnas) %dopar% {
    # ForDetermining Prob it needs to percent
    valmMedia <- quantile(miranda_data$V3, probs = prob)[1]
    # MirandaParsingOfmyRNA
    DiffExpmRNASub <- miRanComp(DiffExpmRNA, miranda_data)

    # everything <- lapply(uni_miRanda, function(mirna) {
    # combiner
    Combine <- combiner(DiffExpmRNASub, DiffExpmiRNA, mirna)

    # GeneVari
    geneVariant <- geneVari(Combine, mirna)

    # ActualModelRan
    MRun <- runModels(Combine, geneVariant, mirna, mode = mode)
    # #Bonferroni sig Genes
    FDRModel <- fdrSig(MRun, value = fdr_cutoff, method = method)

    if (length(FDRModel$all_models) > 0) {
      # look for negative coeffs; "-1" to exclude intercept while checking...#If true all if false any
      if (all_coeff) {
        valid_models <- sapply(FDRModel$all_models, function(m) all(coefficients(m)[-1] < 0))
      } else {
        valid_models <- sapply(FDRModel$all_models, function(m) any(coefficients(m)[-1] < 0))
      }

      # subset fdrmodel's things
      FDRModel$all_models <- FDRModel$all_models[valid_models]
      FDRModel$pvalues <- FDRModel$pvalues[valid_models]
      FDRModel$FDR_Value <- FDRModel$FDR_Value[valid_models]
      FDRModel$sig_genes <- FDRModel$sig_genes[valid_models]
      FDRModel$sanova <- FDRModel$sanova[valid_models]
    }

    # make siggenes data.frame
    SigFDRGenes <- data.frame(cbind(FDRModel$pvalues, FDRModel$FDR_Value))
    names(SigFDRGenes) <- c("P Value", "FDR")
    return(list(SigFDRGenes = SigFDRGenes, FDRModel = FDRModel))
  }
  names(ret) <- sapply(unique_mirnas, paste, collapse = " and ")
  return(ret)
}

#' rsquRes returns r matrix of correlation
#'
#' This function parses R correlation value of the from the list of significant FDR adjusted models with
#' negative/positive correlations.
#' @param FDRSigList list of FDR adjusted models with negative/positive correlations.
#' @return matrix of correlation
#' @export
#' @keywords R correlation
#' @examples
#' \donttest{
#' x <- rsquRes(FDRSigList)
#' }
#'
rsquRes <- function(FDRSigList) {
  # for all mirnas, for all the models get the r-squares
  rsquares <- lapply(FDRSigList, function(x) sapply(x[["FDRModel"]][["all_models"]], modelRsquared))

  # make rownames and colnames for the matrix we want to make
  colnames <- sort(unique(unlist(sapply(rsquares, names))))
  rownames <- names(rsquares)

  # make matrix with zeros
  m <- matrix(0.0, nrow = length(rownames), ncol = length(colnames))

  # set colnames and rownames
  rownames(m) <- rownames
  colnames(m) <- colnames

  # for all mirnas, put the rsqaures in the right places (by name)
  for (row in rownames) {
    m[row, names(rsquares[[row]])] <- rsquares[[row]]
  }

  m1 <- -sqrt(m)
  return(m1)
}


#' DrawCorPlot correlation plots for mRNA and miRNA regression results
#'
#' This function plots correlations for mRNA and miRNAs regression results (negative correlation for multi and
#'  individual interactions and positive and negative for interactions)
#' @param corMatrix correlation matrix obtained from rsquRes function
#' @param ...  parameters form the corrplot package
#' @return miRNA mRNA target correlation plot
#' @export
#' @keywords R correlation plot
#' @examples
#' \donttest{
#' x <- DrawCorPlot(corMatrix)
#' }
#'
DrawCorPlot <- function(corMatrix, ...) {
  col2 <- colorRampPalette(
    c(
      colorRampPalette(c(
        "black",
        "#543005",
        "#8c510a",
        "#f46d43",
        "#762a83",
        "#9970ab",
        "#c2a5cf",
        "#e7d4e8",
        "#f7f7f7",
        "#d9f0d3",
        "#a6dba0",
        "#5aae61",
        "#1b7837",
        "#a6d96a",
        "#9ecae1",
        "#3182bd",
        "#08306b"
      ))(50),
      "grey50",
      colorRampPalette(c("red", "#fde0dd", "#fcc5c0", "#fa9fb5", "pink"))(50)
    )
  )
  corrplot(corMatrix, ..., col = col2(211))
}
