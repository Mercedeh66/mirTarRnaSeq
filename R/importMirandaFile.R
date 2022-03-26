## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom utils read.table download.file
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
#' \donttest{
#' x <- downloadMirandaFile("https://zenodo.org/record/4615670/files/Mouse_miRanda.txt.gz")
#' }
downloadMirandaFile <- function(urlf) {
    return(downloadMirandaFile_(urlf))
}

# use internal function here so we can R.cache it without changing how it looks like for the user.
downloadMirandaFile_ <- function(urlf) {
    # make effective statement error if it returns true or false for http with or without s
    assert_that(grepl("^https?://.*$", urlf))
    tmp <- tempfile(fileext = ".gz")
    op <- options(timeout = 99999)
    on.exit({
        unlink(tmp)
        options(op)
    })
    download.file(urlf, tmp)
    ret1 <- fread(tmp, sep = "\t", header = FALSE)
    return(as.data.frame(ret1))
}
