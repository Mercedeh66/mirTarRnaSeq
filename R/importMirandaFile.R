## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom utils read.table
NULL

#' importMirandaFile Read internal Miranda file
#'
#' Reads internal Miranda file from extdata and returns it as a data.frame
#' @param fn filename
#' @return data.frame containing Miranda data
#' @examples
#' \donttest{
#' x <- importMirandaFile("Mouse_miRanda.txt")
#' }
importMirandaFile <- function(fn) {
  fn <- system.file("extdata", "miRandaPrepFiles", fn, package = "mirTarRnaSeq", mustWork = T)
  ret1 <- read.table(fn, as.is = TRUE, sep = "\t")
  return(ret1)
}
