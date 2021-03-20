## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

#' @importFrom stats cor
#' @importFrom dplyr `%>%`
NULL

## quiet concerns of R CMD check regarding unbound global variables (in dplyr::filter() calls)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("V1", "V2", "V3", "V4", "V5", "V6"))
}

#' Return Miranda data for a given species.
#'
#' Reads Miranda file for a given speicies and returns it as
#' a data.frame, thresholded by percent identity.
#' Header options are Score (threshold), Energy-Kcal/Mol(energy),
#' Subject-IdentityPercent(targetIden), Query-IdentityPercent (mirnaIden)
#' @param selection Species (species selection are either for mature miRNA species "Human1","Mouse",
#'  "C.elegans", "Epstein_Barr", "Epstein_Barr_Human",
#'  "Drosophila","Kaposi_Sarcoma", "KSHV_Human",
#'  "Cytomegalovirus","CMV_Human") Note only EBV files are provided in the package for other files to be determined lol
#' @return data.frame with Miranda data.
#' @keywords miranda, species
#' @export
#' @examples
#' x <- getInputSpecies("Epstein_Barr", threshold = 60) # Default is threshold 60
getInputSpecies <- function(selection, threshold = 60, energy = NULL, targetIden = NULL, mirnaIden = NULL) {
    return(getInputSpecies_(selection=selection, threshold=threshold,
                            energy=energy, targetIden=targetIden,
                            mirnaIden=mirnaIden))
}

# use internal function here so we can R.cache it without changing how it looks like for the user.
getInputSpecies_ <- function(selection, threshold = 60, energy = NULL, targetIden = NULL, mirnaIden = NULL) {
  if (selection == "Human1") {
    ret <- downloadMirandaFile("https://zenodo.org/record/4615670/files/Human_miRanda.txt.gz?download=1")
  } else if (selection == "Mouse") {
    ret <- downloadMirandaFile("https://zenodo.org/record/4615670/files/Mouse_miRanda.txt.gz?download=1")
  } else if (selection == "Epstein_Barr") {
    ret <- importMirandaFile("Epstein_Barr_miRanda.txt.gz")
  } else if (selection == "C.elegans") {
    ret <- downloadMirandaFile("https://zenodo.org/record/4615670/files/C_elegans_miRanda.txt.gz?download=1")
  } else if (selection == "Drosophila") {
    ret <- downloadMirandaFile("https://zenodo.org/record/4615670/files/Drosophila_miRanda.txt.gz?download=1")
  } else if (selection == "Cytomegalovirus") {
    ret <- importMirandaFile("CMV_miRanda.txt.gz")
  } else if (selection == "Kaposi_Sarcoma") {
    ret <- importMirandaFile("Kaposi_miRanda.txt.gz")
  } else if (selection == "Epstein_Barr_Human") {
    ret <- importMirandaFile("EBV_Human_miRanda.txt.gz")
  } else if (selection == "CMV_Human") {
    ret <- importMirandaFile("CMV_Human_miRanda.txt.gz")
  } else if (selection == "KSHV_Human") {
    ret <- importMirandaFile("KSHV_Human_miRanda.txt.gz")
  }
  ret <- dplyr::filter(ret, V3 >= threshold)
  if (!is.null(energy)) {
    ret <- dplyr::filter(ret, V4 <= energy)
  }
  if (!is.null(targetIden)) {
    ret <- dplyr::filter(ret, V5 >= targetIden)
  }
  if (!is.null(mirnaIden)) {
    ret <- dplyr::filter(ret, V6 >= mirnaIden)
  }
  ret <- ret %>% dplyr::select(V1, V2, V3, V4, V5, V6)
  # ret <- unique(ret)
  return(ret)
}
