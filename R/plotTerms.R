## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom graphics mtext
#' @importFrom stats pf
NULL

#' plotTerms
#'
#' Plot terms description
#' @param model linear model
#' @return does not return value
#' @export
#' @keywords plot, residuals, model
#' @examples
#' plotTerms(lm(x ~ y, data = data.frame(x = runif(10), y = runif(10))))
plotTerms <- function(model) {
    m <- model$model
    plot(m, main = as.character.formula(model$terms))
    k <- length(model$coefficients) - 1 # Number of mirnas in model!
    n <- length(model$fitted.values)
    rs <- summary(model)$r.squared
    Pvalue <- pf((rs / k) / ((1 - rs) / (n - 1 - k)), k, n - 1 - k, lower.tail = FALSE) # p value for R-squared

    mtext(sprintf("P value of model fit: %.2g", Pvalue), line = .74)
    mtext(sprintf("R-squared of model fit: %.2f", summary(model)$r.squared), line = .01)
    return(invisible(NULL))
}
