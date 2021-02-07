## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Aug 2020

# convert formula to character (same as private formula.tools:::as.character.formula)
as.character.formula <- function(x, ...) {
  # see formula.tools:::as.character.formula
  form <- paste(deparse(x), collapse = " ")
  form <- gsub("\\s+", " ", form, perl = FALSE)
  return(form)
}
