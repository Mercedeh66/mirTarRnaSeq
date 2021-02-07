#' ## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020
#'
#' #' @importFrom rsq rsq
#' NULL
#'
#' #
#' #' R-squared from model
#' #'
#' #' Obtrain R-squared from a linear model
#' #' @param model linear model
#' #' @export
#' #' @keywords linear model, r-squared
#' #' @examples
#' #' \donttest{
#' #' x <- modelRsquared(model)
#' #' }
#' modelRsquared <- function(model) {
#'   return(rsq(model))
#'  }
