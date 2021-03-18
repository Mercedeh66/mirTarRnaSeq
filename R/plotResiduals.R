## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' Plot residuals
#'
#' Plot residuals description
#' @param model linear model
#' @return does not return value
#' @export
#' @keywords plot, residuals, model
#' @examples
#' plotResiduals(lm(x ~ y, data=data.frame(x=runif(10), y=runif(10))))
plotResiduals <- function(model) {
  plot(model, 1)
  return(invisible(NULL))
}
