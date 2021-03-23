## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom graphics par abline points mtext
#' @importFrom stats pf lm
NULL

#' Plot model
#'
#' Plot 2D description
#' @param model linear model
#' @return does not return value
#' @export
#' @keywords plot, 2d, model
#' @examples
#' plotFit(lm(x ~ y, data = data.frame(x = runif(10), y = runif(10))))
plotFit <- function(model) {
  m <- model$model
  m$Fitted <- model$fitted.values
  last <- ncol(m)
  names(m)[last] <- "Fitted values\nlm(model_formula)"
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(m[, c(1, last)],
    type = "n", main = as.character.formula(model$terms),
    cex.main = 1
  )
  abline(lm(m[, last] ~ m[, 1]), lwd = 2, col = "red")
  points(m[, c(1, last)])

  k <- length(model$coefficients) - 1 # Number of mirnas in model!
  n <- length(model$fitted.values)
  rs <- summary(model)$r.squared
  Pvalue <- pf((rs / k) / ((1 - rs) / (n - 1 - k)), k, n - 1 - k, lower.tail = FALSE) # p value for R-squared

  mtext(sprintf("P value: %.2g", Pvalue), line = .74)
  mtext(sprintf("R-squared: %.2f", summary(model)$r.squared), line = .01)
  return(invisible(NULL))
}
