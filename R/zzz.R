# No need to download a second time if downloaded once we use memoization
.onLoad <- function(libname, pkgname) {
    # setup up memoization for expensive functions (if R.cache available)
    if (requireNamespace("R.cache", quietly = TRUE)) {
        assign("getInputSpecies_", R.cache::addMemoization(getInputSpecies_), envir = parent.env(environment()))
        assign("downloadMirandaFile_", R.cache::addMemoization(downloadMirandaFile_), envir = parent.env(environment()))
        assert_that(!is.null(getInputSpecies_) && is.function(getInputSpecies_)) # if NULL, memoization failed
        assert_that(!is.null(downloadMirandaFile_) && is.function(downloadMirandaFile_)) # if NULL, memoization failed
    }
}
