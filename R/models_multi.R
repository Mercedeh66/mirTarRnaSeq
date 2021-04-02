## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' Model functions for GLM with negative binomial family.
#'
#' Runs models 'glm_gaussian', 'glm_nb', 'glm_poisson',
#' 'glm_zeroinfl(poisson)', 'glm_zeroinfl(negbin)' and returns
#' mode with lowest AIC.
#' @param models Model type, one or more of glm_gaussian, glm_nb, glm_poisson, glm_zeroinfl_poisson or glm_zeroinfl_negbin
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_multi" in \code{model}.
#' @export
#' @examples
#' x <- glm_multi()
glm_multi <- function(models = c(
                          glm_gaussian, glm_nb, glm_poisson,
                          glm_zeroinfl_poisson, glm_zeroinfl_negbin
                      )) {

    # fit
    fit <- function(x, data, ...) {
        bestaic <- NULL
        bestmodel <- NULL
        for (m in models) {
            model <- runModel(x, data, model = m, ...)
            if (!is.null(model)) {
                currentaic <- modelAIC(model)
                if (is.null(bestaic) || currentaic < bestaic) {
                    bestaic <- currentaic
                    bestmodel <- model
                }
            }
        }
        return(bestmodel)
    }

    structure(list(
        fit = fit, coefficients = modelCoefficients, aic = modelAIC, data = modelData,
        pterm = modelTermPvalues, pmodel = modelModelPvalue, model = "glm_multi"
    ))
}
