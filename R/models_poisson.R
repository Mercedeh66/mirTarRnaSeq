## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats poisson glm
NULL

#' Model functions for GLM with Poisson model.
#'
#' Implements standardaized functions to fit the glm with
#' Poisson family and to obtain coefficients, pvalues, etc.
#'
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_poisson" in \code{model}.
#' @export
#' @examples 
#' x <- glm_poisson()
glm_poisson <- function() {
  model_name <- "glm_poisson"

  transform_ <- function(x, data, scale) {
    mf <- modelframe_(x, data) # extract model data.frame
    mf[, 1] <- as.integer(round(mf[, 1] * scale)) # scale and round response variable
    attr(mf, "scale") <- scale # add 'scale' as attribute to data.frame
    return(mf)
  }

  fit_ <- function(x, data) {
    tryCatch(
      {
        g <- glm(x, data = data, family = poisson())
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
  fit <- function(x, data, scale = 1, ...) {
    data <- transform_(x, data, scale)
    return(fit_(x, data))
  }

  pterm <- function(model) {
    nm <- termlabels_(model$formula)
    ret <- coefficients_(model)[, "Pr(>|z|)", drop = TRUE]
    names(ret) <- nm
    return(ret)
  }

  structure(list(
    fit = fit, coefficients = coefficients_, aic = aic_, data = data_,
    pterm = pterm, pmodel = pmodel_, model = model_name
  ))
}
