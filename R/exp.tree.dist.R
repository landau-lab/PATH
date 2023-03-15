#' Get the (inverse) exponential tree distances between all cells
#'
#' This function gets the (inverse) exponential tree or phylogenetic distances, in terms of node distances or branch lengths, between all cell pairs. Exponential tree distances are normalized to sum to one.  
#' @param tree A phylogeny.
#' @param node Logical. Default TRUE. Tree distances are measured in terms of nodes. If set to FALSE distances in terms of branch lengths.
#' @param norm. Logical. Default TRUE. If TRUE, iverse tree distances are row normalized before sum normalizing.
#' @return A matrix of inverse phylogenetic distances between cells. 
#' @export exp.tree.dist
#' @usage exp.tree.dist(tree, node=TRUE, norm=TRUE)
exp.tree.dist <- function(tree, node=TRUE, norm=TRUE) {
  w <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = node)
  w <- exp(-w)
  diag(w) <- 0
  if(norm==TRUE) {
    w <- rowNorm(w)
  }
  return(w/sum(w))
}
