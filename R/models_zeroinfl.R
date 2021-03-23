## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom pscl zeroinfl zeroinfl.control
#' @importFrom assertthat assert_that
#' @importFrom purrr is_formula
#' @importFrom stats formula AIC
NULL


#' Model functions for zero inflated model using either Poisson or Negative Binomial distributions.
#'
#' Implements standardaized functions to fit the zero inflated model with
#' Poisson or Negative Binomial distribution, and to obtain coefficients,
#' pvalues, etc.
#'
#' @param dist either 'poisson' or 'negbin'
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_zeroinfl" in \code{model}.
#' @export
#' @examples 
#' x <- glm_zeroinfl("negbin")
glm_zeroinfl <- function(dist = "poisson") { # poisson, negbin
  assert_that(dist %in% c("poisson", "negbin"))
  model_name <- sprintf("glm_zeroinfl_%s", dist)

  transform_ <- function(x, data, scale) {
    mf <- modelframe_(x, data) # extract model data.frame
    mf[, 1] <- as.integer(round(mf[, 1] * scale)) # scale and round response variable
    attr(mf, "scale") <- scale # add 'scale' as attribute to data.frame
    return(mf)
  }

  fit_ <- function(x, data) {
    if (is_formula(x)) {
      x <- as.character.formula(x)
    }
    if (!grepl("\\|", x)) {
      x <- paste(x, "|", "1")
    }
    tryCatch(
      {
        g <- zeroinfl(formula(x), data = data, dist = dist,
                      control=zeroinfl.control(maxit=100000))
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
      ret$data <- data
    }
    return(ret)
  }

  pterm <- function(model) {
    nm <- termlabels_(model$formula)
    ret <- coefficients_(model)[, "Pr(>|z|)", drop = TRUE]
    names(ret) <- nm
    return(ret)
  }

  coefficients_ <- function(model) {
    nm <- termlabels_(model$formula)
    return(as.data.frame(summary(model)$coefficients$count)[nm, ])
  }

  aic_ <- function(model) {
    return(AIC(model))
  }

  structure(list(
    fit = fit, coefficients = coefficients_, aic = aic_, data = data_,
    pterm = pterm, pmodel = pmodel_, model = model_name
  ))
}

#' alias for glm_zeroinfl("poisson")
#' @param ... passed to \code{glm_zeroinfl}
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_zeroinfl" in \code{model}.
#' @export
#' @examples
#' x <- glm_zeroinfl_poisson()
glm_zeroinfl_poisson <- function(...) glm_zeroinfl("poisson")

#' alias for glm_zeroinfl("negbin")
#' @param ... passed to \code{glm_zeroinfl}
#' @return structure containing functions \code{fit}, \code{coefficients},
#'         \code{aic}, \code{data}, \code{pterm}, \code{pmodel}, and a
#'         character string "glm_zeroinfl" in \code{model}.
#' @export
#' @examples
#' x <- glm_zeroinfl_negbin()
glm_zeroinfl_negbin <- function(...) glm_zeroinfl("negbin")
