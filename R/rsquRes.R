#' ## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020
#'
#'  #' rsquRes returns r matrix of correlation
#'  #'
#'  #' This function parses R correlation value of the from the list of significant FDR adjusted models with
#'  #' negative/positive correlations.
#'  #' @param FDRSigList list of FDR adjusted models with negative/positive correlations.
#'  #' @return matrix of correlation
#'  #' @export
#'  #' @keywords R correlation
#'  #' @examples
#'  #' \donttest{
#'  #' x <- rsquRes(FDRSigList)
#'  #' }
#'  #'
#' rsquRes <- function(FDRSigList) {
#'  # for all mirnas, for all the models get the r-squares
#'  rsquares <- lapply(FDRSigList, function(x) sapply(x[["all_models"]], modelRsquared))
#'
#' # make rownames and colnames for the matrix we want to make
#' colnames <- sort(unique(unlist(sapply(rsquares, names))))
#' rownames <- names(rsquares)
#' # make matrix with zeros
#' m <- matrix(0.0, nrow = length(rownames), ncol = length(colnames))
#'  # set colnames and rownames
#' rownames(m) <- rownames
#' colnames(m) <- colnames
#' # for all mirnas, put the rsqaures in the right places (by name)
#'  for (row in rownames) {
#'  m[row, names(rsquares[[row]])] <- rsquares[[row]]
#'  }
#'
#'  m1 <- -sqrt(m)
#'  return(m1)
#' }
