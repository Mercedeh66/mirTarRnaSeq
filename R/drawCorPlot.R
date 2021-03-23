## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom grDevices colorRampPalette
#' @importFrom corrplot corrplot
NULL

#' drawCorPlot correlation plots for mRNA and miRNA regression results
#'
#' This function plots correlations for mRNA and miRNAs regression results (negative correlation for multi and
#'  individual interactions and positive and negative for interactions)
#' @param corMatrix Significant correlation matrix 
#' @param ...  parameters form the corrplot package
#' @return miRNA mRNA target correlation plot
#' @export
#' @keywords R correlation plot
#' @examples
#' \donttest{
#' x <- drawCorPlot(corMatrix)
#' }

drawCorPlot <- function(corMatrix, ...) {
  col2 <- colorRampPalette(
    c(
      colorRampPalette(c(
        "black",
        "#543005",
        "#8c510a",
        "#f46d43",
        "#762a83",
        "#9970ab",
        "#c2a5cf",
        "#e7d4e8",
        "#f7f7f7",
        "#d9f0d3",
        "#a6dba0",
        "#5aae61",
        "#1b7837",
        "#a6d96a",
        "#9ecae1",
        "#3182bd",
        "#08306b"
      ))(50),
      "grey50",
      colorRampPalette(c("red", "#fde0dd", "#fcc5c0", "#fa9fb5", "pink"))(50)
    )
  )
  corrplot(corMatrix, ..., col = col2(211))
}
