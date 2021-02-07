## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' Plot residuals
#'
#' Plot residuals description
#' @param model linear model
#' @export
#' @keywords plot, residuals, model
#' @examples
#' \donttest{
#' plotResiduals(model)
#' }
plotResiduals <- function(model) {
  plot(model, 1)
}
