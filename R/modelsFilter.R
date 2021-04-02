## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' modelsFilter Filter a list of models based on logical expression
#'
#' This function can be used to filter a list of models (such as returned by \code{runModelsZInf()})
#' based on a logical expression.
#' @param models list of models and related elemenets, such as returned by \code{runModelsZInf()}
#' @param expr expresion that yields a logical vector (evaluated in the environmnet of \code{model})
#' @param quiet suppress warnings
#' @return \code{models} but with all elements filtered by logical expression \code{expr}. Elements
#'         for which filter could not be applied (e.g. length mismatch between element and condition)
#'         are set to \code{NA}.
#' @export
#' @keywords R models filter
#' @examples
#' \donttest{
#' x <- modelsFilter(models, pvalues < 0.05)
#' x <- modelsFilter(models, is_significant)
#' x <- modelsFilter(models, is_significant == FALSE)
#' }
modelsFilter <- function(models, expr, quiet = FALSE) {
    assertthat::assert_that(is.list(models)) # `models` must be a list
    cond <- eval(substitute(expr), envir = models, enclos = parent.frame(n = 1))
    if (!is.logical(cond)) {
        stop(paste0("expression '", deparse(substitute(expr)), "' does not yield logical value"))
    }
    for (i in seq_len(length(models))) {
        if (length(cond) != length(models[[i]])) {
            if (identical(quiet, FALSE)) {
                name <- names(models)[i]
                if ("" != name) {
                    warning(sprintf("condition length mismatch ('%s')", name))
                } else {
                    warning(sprintf("condition length mismatch (i=%d)", i))
                }
            }
            models[[i]] <- NA
        } else {
            models[[i]] <- models[[i]][which(cond)]
        }
    }
    return(models)
}
