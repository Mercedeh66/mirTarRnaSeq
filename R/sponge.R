#' @importFrom dplyr `%>%` filter select mutate
#' @importFrom reshape2 dcast
NULL

#' Transform miRanda data for relevant mRNA and miRNA to matrix form compatible with sponge
#'
#' Transforms miRanda data into adjacency matrix, with 1 indicating
#' presence of a relationship between a mRNA and miRNA, and 0
#' otherwise. miRanda input is filtered by miRNA and mRNA present in
#' `mirna_exp` and `diff_expr` row names, respectively.
#'
#' @param mirna_exp miRna expression data.frame with miRNA for rows and samples for columns
#' @param diff_exp mRNA expression data.frame with mRNA for rows and samples for columns
#' @param miranda_data miRanda data.frame with the first two columns having miRNA and mRNA names
#' @return matrix adjacency matrix with column names miRNA and row names mRNA
#' @keywords miRanda, sparse partial correlation, ceRNA, SPONGE
#' @export
miranda_sponge_predict <- function(mirna_exp, diff_exp, miranda_data) {
  # ensure all inputs are data.frames
  stopifnot("mirna_exp not a data.frame"=is.data.frame(mirna_exp))
  stopifnot("diff_exp not a data.frame"=is.data.frame(diff_exp))
  stopifnot("miranda_data not a data.frame"=is.data.frame(miranda_data))

  # Ensure column names for miRanda are what we expect
  miranda_data <- miranda_data[, c(1, 2)] %>% mutate(V3=1L)
  colnames(miranda_data) <- paste0("V", seq.int(1, ncol(miranda_data)))

  # Add rows for both mrna and mirna to miranda, to ensure all combinations
  # exists in final "dcasted" output
  miranda_data <- miranda_data %>%
                  rbind(data.frame(V1=rownames(mirna_exp),
                                   V2=rownames(diff_exp)[1],
                                   V3=0L)) %>%
                  rbind(data.frame(V1=rownames(mirna_exp)[1],
                                   V2=rownames(diff_exp),
                                   V3=0L))

  # dcast into adjacency matrix, only for combinations of miRNA and mRNA
  # that exist in expressio data
  sponge_mtx <- miranda_data %>%
    filter(V2 %in% rownames(diff_exp) & V1 %in% rownames(mirna_exp))

  if (nrow(sponge_mtx) == 0 || ncol(sponge_mtx) == 0) {
    stop("No interactions present for given input matrices")
  }

  # dcast into adjacency matrix
  sponge_mtx <- dcast(sponge_mtx, V2 ~ V1, max, fill=0L, value.var="V3")

  # Use V1 for rownames
  rownames(sponge_mtx) <- sponge_mtx$V2
  sponge_mtx <- sponge_mtx %>% select(-V2)

  if (nrow(sponge_mtx) == 0 || ncol(sponge_mtx) == 0) {
    stop("No interactions present for given input matrices")
  }

  return(as.matrix(sponge_mtx))
}



