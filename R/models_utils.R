## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom purrr is_formula
#' @importFrom stats terms formula model.frame deviance pchisq df.residual
NULL

coefficients_ <- function(model) {
  nm <- attr(terms(model$formula), "term.labels")
  return(as.data.frame(summary(model)$coefficients)[nm, ])
}

pmodel_ <- function(model) {
  ret <- pchisq(model$null.deviance - deviance(model), model$df.null - df.residual(model), lower.tail = FALSE)
  # ret <- pchisq(model$deviance, model$df.residual, lower.tail = FALSE)
  return(ret)
}

aic_ <- function(model) {
  return(model$aic)
}

data_ <- function(model) {
  return(model$data)
}

modelframe_ <- function(x, data) {
  if (is_formula(x)) {
    x <- as.character.formula(x)
  }
  x <- formula(gsub("\\|", "+", x))
  model.frame(x, data) # extract model data.frame
}

termlabels_ <- function(x) {
  f <- formula(gsub("\\|.*$", "", as.character.formula(x)))
  nm <- attr(terms(f), "term.labels")
}

# notes on model evaluation
# https://www2.stat.duke.edu/courses/Fall17/sta521/knitr/Lec-9-Poisson-Checks/predictive-checking.pdf
# https://stats.stackexchange.com/questions/129958/glm-in-r-which-pvalue-represents-the-goodness-of-fit-of-entire-model
# http://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
# http://users.stat.umn.edu/~helwig/notes/msd-Notes.pdf
# https://www.scribbr.com/statistics/akaike-information-criterion/
