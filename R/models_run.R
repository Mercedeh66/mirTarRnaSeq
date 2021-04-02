## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom assertthat assert_that
#' @importFrom purrr is_formula
#' @importFrom stats terms
NULL

#' Run a model of a specific kind
#'
#' @param x model formula
#' @param data data.frame to run the model on
#' @param ... passed on to \code{fit()}
#' @param model model type
#' @return fitted model
#' @export
runModel <- function(x, data, ..., model = glm_gaussian()) {
    assert_that(is.character(x) || is_formula(x))
    assert_that(is.data.frame(data))
    # obtain model and extract member functions
    model <- canonicalModel_(model)
    fit_ <- model$fit
    ret <- fit_(x, data, ...)
    return(ret)
}

#' Obtain coefficients
#'
#' @param x fitted model
#' @return fitted model coefficients
#' @export
#' @examples
#' modelCoefficients(some_model)
modelCoefficients <- function(x) {
    assert_that(!is.vector(x), msg = "not a valid model object")
    assert_that(!is.null(x$model), msg = "not a valid model object")
    # obtain model and extract member functions
    model <- canonicalModel_(attr(x, "model"))
    coefficients_ <- model$coefficients
    x <- coefficients_(x)
    ret <- x[, "Estimate"]
    names(ret) <- rownames(x)
    return(ret)
}

#' Obtain model AIC
#'
#' @param x fitted model
#' @return AIC for model
#' @export
#' @examples
#' modelAIC(some_model)
modelAIC <- function(x) {
    assert_that(!is.vector(x), msg = "not a valid model object")
    assert_that(!is.null(x$model), msg = "not a valid model object")
    # obtain model and extract member functions
    model <- canonicalModel_(attr(x, "model"))
    return(model$aic(x))
}

#' Obtain model input data
#'
#' @param x fitted model
#' @return  Input data for the fitted model
#' @export
#' @examples
#' x <- modelData(some_model)
modelData <- function(x) {
    assert_that(!is.vector(x), msg = "not a valid model object")
    assert_that(!is.null(x$model), msg = "not a valid model object")
    # obtain model and extract member functions
    model <- canonicalModel_(attr(x, "model"))
    return(model$data(x))
}

#' Obtain p-values for terms in model formula
#'
#' @param x fitted model
#' @return  Pvalue for the terms in the fitted model
#' @export
#' @examples
#' modelTermPvalues(some_model)
modelTermPvalues <- function(x) {
    assert_that(!is.vector(x), msg = "not a valid model object")
    assert_that(!is.null(x$model), msg = "not a valid model object")
    model <- canonicalModel_(attr(x, "model"))
    return(model$pterm(x))
}

#' Obtain model p-value
#'
#' @param x fitted model
#' @return Pvalue for the model
#' @export
#' @examples
#' modelModelPvalue(some_model)
modelModelPvalue <- function(x) {
    assert_that(!is.vector(x), msg = "not a valid model object")
    assert_that(!is.null(x$model), msg = "not a valid model object")
    model <- canonicalModel_(attr(x, "model"))
    return(model$pmodel(x))
}

#' Obtain model name
#'
#' @param x fitted model
#' @return model name
#' @export
#' @examples
#' modelModelName(some_model)
modelModelName <- function(x) {
    assert_that(!is.vector(x), msg = "not a valid model object")
    assert_that(!is.null(x$model), msg = "not a valid model object")
    return(attr(x, "model"))
}


#' Decifer a 'model parameter' and run appropriate \code{glm_...} function.
#'
#' Return canonical model from model type string, function of object. Returns
#' a model as returned by \code{glm_gaussian()} and others, based on a string,
#' function or model type object (i.e. \code{"glm_gaussian"}, \code{glm_gaussian} or
#' \code{glm_gaussian()}).
#'
#' @param model string, function or object representing a model type.
#' @return model type object
canonicalModel_ <- function(model) {
    # figure out the type of model.
    assert_that(!is.null(model))
    if (is.function(model)) {
        # function was supplied: run the function
        x <- model()
    } else if (is.character(model)) {
        # character vector was supplied: try to obtain function of specified name
        assert_that(1 == length(model))
        model <- get(model, mode = "function", envir = parent.frame())
        x <- model()
    } else {
        # assume that what we got was the result of model function call
        x <- model
    }
    # by now x should be a structure with a "model" entry
    assert_that(!is.null(x$model))
    return(x)
}
