## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats sd glm gaussian
NULL

#' Model functions for GLM with Gaussian model.
#'
#' Implements standardaized functions to fit the glm with
#' Gaussian family and to obtain coefficients, pvalues, etc.
#'
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_gaussian" in \code{model}.
#' @export
glm_gaussian <- function() {
  model_name <- "glm_gaussian"

  transform_ <- function(x, data, zscore) {
    mf <- modelframe_(x, data) # extract model data.frame
    if (identical(zscore, TRUE)) {
      mf[, 1] <- (mf[, 1] - mean(mf[, 1])) / sd(mf[, 1])
    }
    attr(mf, "zscore") <- zscore # add 'zscore' as attribute to data.frame
    return(mf)
  }

  fit_ <- function(x, data) {
    tryCatch(
      {
        g <- glm(x, data = data, family = gaussian())
        if (!is.null(g)) {
          attr(g, "model") <- model_name
        }
        return(g)
      },
      warning = function(e, ...) {
#        warning(e)
        return(NULL)
      },
      error = function(e, ...) {
#        warning(e)
        return(NULL)
      }
    )
  }

  # fit
  fit <- function(x, data, zscore = FALSE, ...) {
    data <- transform_(x, data, zscore)
    return(fit_(x, data))
  }

  pterm <- function(model) {
    nm <- termlabels_(model$formula)
    ret <- coefficients_(model)[, "Pr(>|t|)", drop = TRUE]
    names(ret) <- nm
    return(ret)
  }

  structure(list(
    fit = fit, coefficients = coefficients_, aic = aic_, data = data_,
    pterm = pterm, pmodel = pmodel_, model = model_name
  ))
}
