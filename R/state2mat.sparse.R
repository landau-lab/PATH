#' Convert vector of cell states into matrix to compute phylogenetic correlations and perform transition inferences.
#'
#' @param states A vector of cell states from a phylogeny.
#' @param n Number of possible cell states. 
#' @export
state2mat.sparse <- function(states, n=NULL) {
  if(is.null(n)==T) {
    sparseMatrix(1:length(states), states)
  } else {
    sparseMatrix(1:length(states), states, dims = c(length(states), n))
  }
}
