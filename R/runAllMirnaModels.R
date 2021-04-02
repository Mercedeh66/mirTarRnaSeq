## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats quantile
#' @importFrom caTools combs
NULL

#' runAllMirnaModels runModel for all miRNAs
#'
#' This function runs the "runModel" function for all miRNAs and mRNA combinations of two and returns a
#' list with significant genes and FDR models
#' @param mirnas vector of unique miRNAs under investigation.
#' @param DiffExpmRNA deferentially/expressed mRNAs expression file.
#' @param DiffExpmiRNA deferentially/expressed miRNAs expression file.
#' @param miranda_data getInputSpecies output file ( use low filters).
#' @param prob user defined ratio for miRanda distribution for miRanda score selection default is 0.75.
#' @param fdr_cutoff cutoff for FDR selection default is 0.1.
#' @param method finInterResult miRNA and mRNA interrelation in two time points  results in a dataframe.
#' @param cutoff P-value cutoff of the model.
#' @param all_coeff if true only models with all negative coefficients will be selected if false at least one
#' negative coefficient should be in the model; default is TRUE.
#' @param family Default is glm_poisson(), for zero inflated negative binomial NB option use glm_zeroinfl(dist="negbin").
#' @param mode model mode, default is Null, can be changed to "multi" and "inter".
#' @param scale  if normalized data (FPKM,RPKM,TPM,CPM), scale to 10 etc., however the higher you go on
#' #scale the less accuracy your p-value estimate will be.
#' @return List of run models
#' @export
#' @keywords runAllMirnaModels all_miRNAs
#' @examples
#' mirnas <- c("ebv-mir-bart9-5p", "ebv-mir-bart6-3p")
#' x <- runAllMirnaModels(mirnas, mRNA, miRNA, miRanda,
#'     prob = 0.90, fdr_cutoff = 0.1, method = "fdr",
#'     all_coeff = TRUE, mode = "multi",
#'     family = glm_poisson(), scale = 100
#' )
runAllMirnaModels <- function(mirnas, DiffExpmRNA, DiffExpmiRNA, miranda_data,
                              prob = 0.75, fdr_cutoff = 0.1, method = "fdr", cutoff = 0.05,
                              all_coeff = FALSE, mode = NULL, family = glm_poisson(), scale = 1) {
    assertthat::assert_that(is.character(mirnas))

    # make unique
    unique_mirnas <- unique(mirnas)

    # generate combinations of 2 miRNAs
    if (!is.null(mode)) {
        # produces list of pair of mirnas, e.g.:
        #
        # ...
        # [[120]]
        # [1] "ebv-mir-bart11-3p" "ebv-mir-bart18-3p"
        #
        # [[121]]
        # [1] "ebv-mir-bart11-3p" "ebv-mir-bart18-5p"
        #
        # [[122]]
        # [1] "ebv-mir-bart11-3p" "ebv-mir-bart19-3p"
        # ...
        unique_mirnas <- caTools::combs(unique_mirnas, 2)
        unique_mirnas <- lapply(seq_len(nrow(unique_mirnas)), function(x) unique_mirnas[x, ])
    }

    # run models
    count <- 0
    ret <- lapply(unique_mirnas, function(mirna) {
        count <<- count + 1
        cat(sprintf("%d: %s\n", count, paste(mirna, collapse = ", ")))

        # For Determining Prob it needs to percent
        valmMedia <- quantile(miranda_data$V3, probs = prob)[1]
        # MirandaParsingOfmyRNA
        DiffExpmRNASub <- miRanComp(DiffExpmRNA, miranda_data)

        # everything <- lapply(uni_miRanda, function(mirna) {
        # combiner
        Combine <- combiner(DiffExpmRNASub, DiffExpmiRNA, mirna)

        # GeneVari
        geneVariant <- geneVari(Combine, mirna)

        MRun <- runModels(Combine, geneVariant, mirna, mode = mode, family = family, scale = scale, cutoff = cutoff)

        # #Bonferroni sig Genes
        FDRModel <- fdrSig(MRun, value = fdr_cutoff, method = method)

        if (length(FDRModel$sig_models) > 0) {
            # look for negative coeffs; "-1" to exclude intercept while checking...#If true all if false any
            if (all_coeff) {
                valid_models <- sapply(FDRModel$sig_models, function(m) all(modelCoefficients(m) < 0))
            } else {
                valid_models <- sapply(FDRModel$sig_models, function(m) any(modelCoefficients(m) < 0))
            }
            # subset fdrmodel's things
            FDRModel$sig_models <- FDRModel$sig_models[valid_models]
            FDRModel$sig_pvalues <- FDRModel$sig_pvalues[valid_models, , drop = FALSE]
            FDRModel$sig_p_adjusted <- FDRModel$sig_p_adjusted[valid_models, , drop = FALSE]
            FDRModel$sig_genes <- FDRModel$sig_genes[valid_models]
        }
        # make siggenes data.frame
        SigFDRGenes <- data.frame(cbind(FDRModel$sig_pvalues, FDRModel$sig_p_adjusted))

        names(SigFDRGenes) <- c("P Value", "FDR")

        return(list(SigFDRGenes = SigFDRGenes, FDRModel = FDRModel))
    })
    names(ret) <- sapply(unique_mirnas, paste, collapse = " and ")
    return(ret)
}
