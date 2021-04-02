## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom MASS glm.nb
#' @importFrom stats terms formula glm.control
NULL

#' Model functions for GLM with negative binomial family.
#'
#' Implements standardized functions to fit the negative
#' binomial GLM and to obtain coefficients, pvalues, etc.
#'
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_nb" in \code{model}.
#' @export
#' @examples
#' x <- glm_nb()
glm_nb <- function() {
    model_name <- "glm_nb"

    coefficients_ <- function(model) {
        nm <- termlabels_(model$formula)
        return(as.data.frame(summary(model)$coefficients)[nm, ])
    }

    transform_ <- function(x, data, scale) {
        mf <- modelframe_(x, data) # extract model data.frame
        mf[, 1] <- as.integer(round(mf[, 1] * scale)) # scale and round response variable
        attr(mf, "scale") <- scale # add 'scale' as attribute to data.frame
        return(mf)
    }

    fit_ <- function(x, data) {
        tryCatch(
            {
                g <- glm.nb(x, data = data, control = glm.control(maxit = 1000))
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
        ret <- fit_(x, data)
        if (!is.null(ret)) {
            ret$data <- data # attach data here; glm.nb doesn't do that for some reason...
            ret$formula <- formula(x) # attach formula here; glm.nb doesn't do that for some reason...
        }
        return(ret)
    }

    pterm <- function(model) {
        nm <- attr(terms(model), "term.labels")
        ret <- coefficients_(model)[, "Pr(>|z|)", drop = TRUE]
        names(ret) <- nm
        return(ret)
    }

    structure(list(
        fit = fit, coefficients = coefficients_, aic = aic_, data = data_,
        pterm = pterm, pmodel = pmodel_, model = model_name
    ))
}
