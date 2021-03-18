## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom utils read.table
#' @importFrom data.table fread
#' @import R.utils
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
  fn <- system.file("extdata", "miRandaPrepFiles", fn, package = "mirTarRnaSeq", mustWork = TRUE)
  ret1 <- read.table(fn, as.is = TRUE, sep = "\t")
  return(ret1)
}


#' downloadMirandaFile Read internal Miranda file
#'
#' Reads internal Miranda file from extdata and returns it as a data.frame
#' @param urlf URL of the specific chosen file
#' @return data.frame containing downloaded miRanda file
#' @examples
#' x <- downloadMirandaFile("https://zenodo.org/record/4615670/files/Mouse_miRanda.txt.gz?download=1")

downloadMirandaFile <- function(urlf) {
  #make effective statement error if it returns true or false for http with or without s
  assert_that(grepl("^https?://.*$", urlf))
  tmp <- tempfile(fileext=".gz") 
  download.file(urlf, tmp) 
  ret1 <- fread(tmp, sep="\t", header=FALSE) 
  unlink(tmp)
  return(as.data.frame(ret1))
}
#No need to download a second time if downloaded once we use memoization
.onLoad <- function(libname, pkgname) {
  # setup up memoization for expensive functions (if R.cache available)
  if (requireNamespace("R.cache", quietly = TRUE)) {
    assign("downloadMirandaFile", R.cache::addMemoization(downloadMirandaFile), envir = parent.env(environment()))
  }
}
