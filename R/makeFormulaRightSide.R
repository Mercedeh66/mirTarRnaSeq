## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' makeFormulaRightSide makes right hand side of formula for model variables: vector of indep. variables
#'
#' This function make right hand side of formula for model variables: vector of indep. variables (i.e. miRNAs)
#' mode: 'multi' for simple, 'inter' for model with interactions returns a string in the form "~ a + b", or "~ a + b + a * b"
#' @param variables The vector created by miRNA_select
#' @param mode One of "multi", "inter" or NULL
#' @return data.frame containing Miranda data
#' @examples
#' \donttest{
#' x <- makeFormulaRightSide(variables, mode = "multi")
#' }
#'
makeFormulaRightSide <- function(variables, mode = "multi") {
  # be sure to properly quote variable names: `varname`.
  if (is.null(mode)) {
    mode <- "multi"
  }
  replace <- !grepl("^`.*`$", variables) # all that don't look like `somevar`
  variables[replace] <- paste("`", variables[replace], "`", sep = "")
  rightside <- paste(variables, collapse = " + ") # "~ a + b"
  if ("inter" == mode) {
    for (i in 1:(length(variables) - 1)) {
      for (j in (i + 1):length(variables)) {
        rightside <- paste(rightside, paste(variables[i], variables[j], sep = " * "), sep = " + ") # "~ a + b + a * b"
      }
    }
  }
  return(rightside)
}
