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
#' @param selection Species (species selection are either for mature miRNA species "Human1"
#'  for full length human miRNAs use "Human(1-10)", "Mouse",
#'  "C.elegans", "Epstein_Barr", "Epstein_Barr_Human",
#'  "Drosophila","Kaposi_Sarcoma", "KSHV_Human",
#'  "Cytomegalovirus","CMV_Human") Note only EBV files are provided in the package for other files to be determined lol
#' @return data.frame with Miranda data.
#' @keywords miranda, species
#' @export
#' @examples
#' x <- getInputSpecies("Epstein_Barr", threshold = 60) # Default is threshold 60

getInputSpecies <- function(selection, threshold = 60, energy = NULL, targetIden = NULL, mirnaIden = NULL) {
  if (selection == "Human1") {
    ret <- importMirandaFile("Human_miRanda1.txt.gz")
  } else if (selection == "Human2") {
    ret <- importMirandaFile("Human_miRanda2.txt.gz")
  } else if (selection == "Human3") {
    ret <- importMirandaFile("Human_miRanda3.txt.gz")
  } else if (selection == "Human4") {
    ret <- importMirandaFile("Human_miRanda4.txt.gz")
  } else if (selection == "Human5") {
    ret <- importMirandaFile("Human_miRanda5.txt.gz")
  } else if (selection == "Human6") {
    ret <- importMirandaFile("Human_miRanda6.txt.gz")
  }  else if (selection == "Human7") {
    ret <- importMirandaFile("Human_miRanda7.txt.gz")
  }  else if (selection == "Human8") {
    ret <- importMirandaFile("Human_miRanda8.txt.gz")
  }  else if (selection == "Human9") {
    ret <- importMirandaFile("Human_miRanda9.txt.gz")
  }  else if (selection == "Human10") {
    ret <- importMirandaFile("Human_miRanda10.txt.gz")
  } else if (selection == "Mouse") {
    ret <- importMirandaFile("Mouse_miRanda.txt.gz")
  } else if (selection == "Epstein_Barr") {
    ret <- importMirandaFile("Epstein_Barr_miRanda.txt.gz")
  } else if (selection == "C.elegans") {
    ret <- importMirandaFile("C.elegans_miRanda.txt.gz")
  } else if (selection == "Drosophila") {
    ret <- importMirandaFile("Drosophila_miRanda.txt.gz")
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
