## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats var formula
NULL

#' runModels runs miRNA mrna model model for various miRNA-mRNA data distributions
#'
#' This function defines the boudnaries of mRNA vs miRNAs of interest to be analysed by the runModels function
#' @param combination the combined file for mRNA and selected miRNAs output of combiner function
#' @param select_mRNA the output of
#' @param select_miRNA The vector of miRNA/s to be investigated
#' @param mode the mode of analysis if more than one miRNA is being investigated multivariate "multi"
#' or covariate/interaction analysis "inter" is being used
#' @param family gaussian or poisson
#' @param scale factor to scale input data (for genes) by, prior to rounding and model
#'              fitting. (\code{scale} must be greater than zero).
#' @param cutoff p-value cut off to call significance
#' @param all_coeff if true only models with all negative coefficients will be selected if false at least one
#' @return A list of p vlaues, annova, and significance for each gene and the miRNA/s of interest
#' @export
#' @keywords runModels univariate multivariate interaction glm
#' @examples
#' \donttest{
#' x <- runModels(combination, geneVariable, miRNA_selectedVec, mode = NULL, family = "poisson")
#' }
#'
runModels <- function(combination, select_mRNA, select_miRNA, mode = NULL,
                      family = glm_poisson(), scale = 1, cutoff = 0.05,
                      all_coeff=NULL) {
  assertthat::assert_that(is.numeric(scale) && 1 == length(scale) && scale > 0)
  assertthat::assert_that(is.numeric(cutoff) && 1 == length(cutoff) && cutoff > 0)
  assertthat::assert_that(is.null(all_coeff) || (is.logical(all_coeff) && 1 == length(all_coeff)))
  is_significant <- c()
  all_models <- list()
  genes <- c()
  pvalues <- c()
  AICvalues <- c()
  # remove zero-variance columns.
  combination <- combination[, apply(combination, 2, var) != 0, drop = F]
  select_mRNA <- intersect(select_mRNA, colnames(combination))

  for (gene in select_mRNA) {
    model_formula <- formula(paste(sprintf("`%s`", gene),
      makeFormulaRightSide(select_miRNA, mode = mode),
      sep = " ~ "
    ))

    # run the model for a mRNA
    x <- runModel(model_formula, data = combination, model = family, scale = scale)

    # check for negative coeffs (if 'all_coeff' is used)
    keep <- TRUE
    if (!is.null(all_coeff)) {
      if (identical(all_coeff, TRUE)) {
        if (!(all(modelCoefficients(x) < 0))) {
          keep <- FALSE
        }
      } else {
        if (!(any(modelCoefficients(x) < 0))) {
          keep <- FALSE
        }
      }
    }

    # if we were unable to make a model (i.e. it is NULL) or we didn't get
    # the negative coeffs we wanted, we don't add it to the lists here.
    if (!is.null(x) && keep) {
      ret <- modelTermPvalues(x)
      is_significant <- c(is_significant, any(is.finite(ret) & ret < cutoff))
      pvalues <- c(pvalues, list(ret))
      all_models <- c(all_models, list(x))
      genes <- c(genes, gene)
      AICvalues <- c(AICvalues, list(modelAIC(x)))
    }
  }
  names(is_significant) <- genes
  names(pvalues) <- genes
  names(all_models) <- genes
  names(AICvalues) <- genes

  # turn p-values into data.frame (with columns for every mirna)
  pvalues <- t(as.data.frame(pvalues))
  colnames(pvalues) <- gsub("^`?(.*)[^`]`?$", "\\1", colnames(pvalues))  # remove ` ` from colnames...
  rownames(pvalues) <- genes
  return(
    list(
      is_significant = is_significant, # logical vector TRUE if any term sign.; FALSE otherwise.
      genes = genes,
      pvalues = pvalues,
      all_models = all_models,
      AICvalues = AICvalues,
      select_mRNA = select_mRNA,
      select_miRNA = select_miRNA,
      combination = combination
    )
  )
}
